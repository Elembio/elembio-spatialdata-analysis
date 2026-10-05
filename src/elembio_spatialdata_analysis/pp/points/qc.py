"""Schema-aware resolution of polony-quality columns in points parquets.

Different upstream pipeline versions used different column names for the
per-polony mean base-quality score.  Rather than spread that knowledge
across every consumer, this module owns the lookup and exposes a single
``resolve_quality_column`` helper that returns the column name to use
for a given points ``LazyFrame``.

Recognized names (preference order):

- ``"MeanBaseQ"`` — canonical, matches sibling Q-columns (``MapQ``,
  ``MeanAlnBaseQ``) and the convention used in 3' Untargeted batches.
- ``"MeanBaseQuality"`` — legacy long-form, used in older OPS points
  parquets.  Accepted as a fallback so existing zarrs continue to work
  without rebuild; new writers should emit ``"MeanBaseQ"``.

Returns ``None`` when neither column is present (e.g. Protein points
in the current schema carry no per-polony quality at all).  Callers
should treat ``None`` as "polony quality unavailable for this
modality" and skip polony-quality metrics rather than fall back to a
synthetic value.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import TypedDict

import polars as pl

# Ordered tuple of recognised polony-quality column names.  First match
# wins in ``resolve_quality_column``.  Order encodes the preferred /
# canonical name first, legacy names after.
POLONY_QUALITY_COLUMNS: tuple[str, ...] = (
    "MeanBaseQ",
    "MeanBaseQuality",
)


def resolve_quality_column(
    lf_or_schema: pl.LazyFrame | pl.DataFrame | pl.Schema | Iterable[str],
) -> str | None:
    """Return the polony-quality column name in *lf_or_schema*, or None.

    Accepts a polars ``LazyFrame`` / ``DataFrame`` (schema is read
    lazily without materializing data), a ``Schema``, or any iterable
    of column names — whichever the caller has on hand.  Returns the
    first column name from :data:`POLONY_QUALITY_COLUMNS` that is
    present in the input; ``None`` if neither is.

    The schema check is cheap on a ``LazyFrame`` (no scan, no row
    materialization), so it is safe to call before constructing an
    aggregation expression — e.g. ``pl.col(resolve_quality_column(lf))``.
    """
    if isinstance(lf_or_schema, pl.LazyFrame):
        names: Iterable[str] = lf_or_schema.collect_schema().names()
    elif isinstance(lf_or_schema, pl.DataFrame):
        names = lf_or_schema.schema.names()
    elif isinstance(lf_or_schema, pl.Schema):
        names = lf_or_schema.names()
    else:
        names = lf_or_schema  # already an iterable of strings
    present = set(names)
    for candidate in POLONY_QUALITY_COLUMNS:
        if candidate in present:
            return candidate
    return None


class PolonyQcStats(TypedDict):
    """Aggregate polony-quality stats over a points scan."""

    n_polonies: int
    mean_q: float
    pct_q30: float


def polony_qc_stats(
    lf: pl.LazyFrame,
    *,
    cell_filter: Iterable[int] | None = None,
    use_in_cell_column: bool = False,
    q30_threshold: float = 30.0,
) -> PolonyQcStats | None:
    """Aggregate polony-level QC stats from a points ``LazyFrame``.

    Computes — in a single collect on a column-pruned scan — the
    number of polonies, mean base quality, and % of polonies whose
    mean base quality is ≥ ``q30_threshold``.  Returns ``None`` when
    the input has no recognised quality column (e.g. Protein points,
    which currently carry no per-polony quality field).

    Cell-scope filtering happens via one of two predicates, in order
    of preference:

    - ``use_in_cell_column=True`` — apply the boolean ``InCell``
      column from the points parquet (``pl.col("InCell")``).  This
      predicate pushes down to the parquet scan and avoids
      materialising the per-polony ``cell`` column.  Use this on the
      Raw leaf where the semantic is "all polonies that landed in any
      segmented cell" — equivalent to passing every cell id but
      O(1) memory instead of O(n).
    - ``cell_filter=<iterable of cell ids>`` — true membership filter.
      Required for the Filtered leaf where the caller has already
      decided which cells survive (a strict subset).  Pays the
      ``is_in()`` materialisation cost but typically run on smaller
      data.

    Pass neither for unfiltered stats (every polony in the parquet,
    intra- or extra-cellular).  Pass ``cell_filter=[]`` (empty list)
    to short-circuit to a zero-polonies result without scanning.

    ``q30_threshold`` defaults to 30 to match the canonical "%Q30"
    definition used by bases2fastq, with the important semantic
    difference that this metric counts *polonies whose mean is ≥ Q30*
    rather than individual bases at Q30+.  The polony-level metric is
    stricter (a polony with a few low-Q bases but high mean still
    counts).
    """
    qcol = resolve_quality_column(lf)
    # ``qcol is None`` means the parquet lacks per-polony quality scores —
    # typical for pre-alignment instrument zarrs whose points carry only
    # spatial coordinates + cell IDs.  Quality-derived fields will return
    # NaN and callers (e.g. the polony-QC strip renderer) skip those badges,
    # but the polony **count** is still useful and is reported.

    if use_in_cell_column:
        # Pushdown-friendly "polonies in any cell" predicate, picked
        # per the schema the build pipeline ships:
        # - 3'Untargeted points carry a boolean ``InCell`` column —
        #   use it directly.
        # - OPS points carry no ``InCell`` but use ``cell = 0`` as
        #   the extracellular sentinel (``cell`` is ``Int32``).
        # Either form pushes into the parquet scan, avoiding the
        # ``is_in(<677k cells>)`` materialisation that the previous
        # code used (which forced polars to materialise the per-polony
        # ``cell`` column to do the join).
        schema = lf.collect_schema()
        if "InCell" in schema:
            lf = lf.filter(pl.col("InCell"))
        elif "cell" in schema and schema["cell"].is_integer():
            lf = lf.filter(pl.col("cell") > 0)
        else:
            # String ``cell`` column without ``InCell`` — bail out
            # rather than scan everything and over-count extracellular
            # polonies.  Should not happen on supported modalities.
            return None
    elif cell_filter is not None:
        cells = list(cell_filter)
        if not cells:
            # Empty filter — no polonies pass.  Return a zero-stat
            # struct rather than None so callers can distinguish "0
            # polonies" from "modality has no quality column at all".
            return PolonyQcStats(n_polonies=0, mean_q=float("nan"), pct_q30=float("nan"))
        # Normalize ``cell`` dtype before the ``is_in`` join: OPS points
        # store ``cell`` as ``Int32``; 3'Untargeted points store it as
        # a zero-padded ``String`` (``"00000123"``).  Casting to Int64
        # parses both forms — polars will fail loudly if the string
        # column ever contains a non-numeric value, which would itself
        # be a real upstream regression worth surfacing.
        lf = lf.with_columns(pl.col("cell").cast(pl.Int64))
        lf = lf.filter(pl.col("cell").is_in(cells))

    # ``engine="streaming"`` keeps peak memory bounded when the scan
    # fans out over hundreds of partitioned parquet parts.  The default
    # in-memory engine parallelizes per-file, which on a 12-well × 6-batch
    # transcript polony scan (~1700 parts, ~10 GB) blew the resident set
    # past 20 GB; streaming holds it under 2 GB.
    if qcol is not None:
        result = lf.select(
            n_polonies=pl.len(),
            mean_q=pl.col(qcol).mean(),
            pct_q30=pl.col(qcol).ge(q30_threshold).mean() * 100.0,
        ).collect(engine="streaming")
    else:
        # Count-only path for parquets without quality scores (pre-alignment).
        result = lf.select(
            n_polonies=pl.len(),
            mean_q=pl.lit(None, dtype=pl.Float64),
            pct_q30=pl.lit(None, dtype=pl.Float64),
        ).collect(engine="streaming")
    if result.height == 0:
        return None
    n, mean_q, pct_q30 = result.row(0)
    return PolonyQcStats(
        n_polonies=int(n),
        mean_q=float(mean_q) if mean_q is not None else float("nan"),
        pct_q30=float(pct_q30) if pct_q30 is not None else float("nan"),
    )
