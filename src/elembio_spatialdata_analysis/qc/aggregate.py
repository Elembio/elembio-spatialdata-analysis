"""Per-batch_type and per-well QC rollups.

These helpers consolidate the various source-prefixed runstats DataFrames
in ``global-table-metadata.uns`` and the parsed ``metadata/inputs`` JSON
into small flat summary tables.
"""

from __future__ import annotations

import logging
import math
from typing import Any

import pandas as pd

_logger = logging.getLogger(__name__)


def join_batch_metadata(batches_df: pd.DataFrame | None, panel_assays: dict[str, str] | None) -> pd.DataFrame:
    """Return the run's batch table with a ``batch_type`` column guaranteed.

    ``batches_df`` is ``meta_uns["batches"]`` from the build.  ``panel_assays``
    is the mapping returned by ``PanelInfo.batch_assay_types()``.  When the
    build's ``batch_type`` column is missing or partial, panel assays back-fill.
    """
    if batches_df is None or batches_df.empty:
        return pd.DataFrame(columns=["batch_name", "batch_type", "batch_class"])

    df = batches_df.copy()
    if "batch_name" not in df.columns and df.index.name == "batch_name":
        df = df.reset_index()

    if "batch_type" not in df.columns:
        df["batch_type"] = None

    if panel_assays:
        backfill = df["batch_type"].isna() | (df["batch_type"] == "")
        if backfill.any() and "batch_name" in df.columns:
            df.loc[backfill, "batch_type"] = df.loc[backfill, "batch_name"].map(panel_assays)

    return df


def per_batch_type_summary(batches_df: pd.DataFrame) -> pd.DataFrame:
    """Group ``batches_df`` by ``batch_type``: count of batches, list of names."""
    if batches_df is None or batches_df.empty or "batch_type" not in batches_df.columns:
        return pd.DataFrame(columns=["batch_type", "n_batches", "batch_names"])

    out = (
        batches_df.dropna(subset=["batch_type"])
        .groupby("batch_type", as_index=False)
        .agg(
            n_batches=("batch_name", "count"),
            batch_names=("batch_name", lambda s: ", ".join(sorted(map(str, s)))),
        )
    )
    return out


def b2f_for_batches(
    b2f_per_batch: list[dict[str, Any]] | None, batch_names: list[str]
) -> dict[str, float | int | None]:
    """Aggregate bases2fastq per-batch rows for a given subset of batches.

    Returns ``{"NumPolonies": ..., "PercentQ30": ..., "QualityScoreMean": ...}``.
    Polonies and yield are summed, percentages and quality scores are
    polony-weighted means when polony counts are available.
    """
    if not b2f_per_batch or not batch_names:
        return {"NumPolonies": None, "PercentQ30": None, "QualityScoreMean": None, "TotalYield": None}

    name_set = set(batch_names)
    rows = [r for r in b2f_per_batch if r.get("BatchName") in name_set]
    if not rows:
        return {"NumPolonies": None, "PercentQ30": None, "QualityScoreMean": None, "TotalYield": None}

    total_polonies = _sum_or_none(r.get("NumPolonies") for r in rows)
    total_yield = _sum_or_none(r.get("TotalYield") for r in rows)

    weighted_q30 = _weighted_mean([(r.get("PercentQ30"), r.get("NumPolonies")) for r in rows])
    weighted_qmean = _weighted_mean([(r.get("QualityScoreMean"), r.get("NumPolonies")) for r in rows])

    return {
        "NumPolonies": total_polonies,
        "TotalYield": total_yield,
        "PercentQ30": weighted_q30,
        "QualityScoreMean": weighted_qmean,
    }


def per_tile_cell_count(adata: Any) -> dict[str, int]:
    """Return ``{tile_id: cell_count}`` from ``adata.obs.tile_id``.

    Modality-agnostic: works for any AnnData whose ``obs`` carries a
    ``tile_id`` column (the convention enforced by the loader for all
    cell-shaped tables).  Returns an empty dict when the column is
    absent so callers can branch cleanly.
    """
    if adata is None or "tile_id" not in adata.obs.columns:
        return {}
    counts = adata.obs.groupby("tile_id", observed=True).size()
    return {str(k): int(v) for k, v in counts.items()}


def per_tile_total_counts(adata: Any) -> dict[str, float]:
    """Return ``{tile_id: sum(Counts)}`` when both columns exist."""
    if adata is None or "tile_id" not in adata.obs.columns or "Counts" not in adata.obs.columns:
        return {}
    totals = adata.obs.groupby("tile_id", observed=True)["Counts"].sum()
    return {str(k): float(v) for k, v in totals.items()}


def per_well_overview(
    cyto_wells: pd.DataFrame | None,
    raw_reads_per_well: list[dict[str, Any]] | None,
) -> pd.DataFrame:
    """Build a one-row-per-well table combining cytometry and raw read counts.

    Returns columns: ``WellLocation, WellLabel, CellCount, MedianCellDiameter,
    PercentNucleatedCells, PercentConfluency, NumPoloniesBeforeTrimming,
    ReadsPerCell, ExtracellularRatio`` (subset of these per data availability).
    The ``ReadsPerCell`` / ``ExtracellularRatio`` read-disposition rates come
    from :func:`_well_raw_read_rates`.
    """
    parts: list[pd.DataFrame] = []

    if cyto_wells is not None and not cyto_wells.empty:
        cyto = cyto_wells.copy()
        if "WellLocation" not in cyto.columns and cyto.index.name == "WellLocation":
            cyto = cyto.reset_index()
        keep = [
            c
            for c in (
                "WellLocation",
                "WellLabel",
                "CellCount",
                "MedianCellDiameter",
                "PercentNucleatedCells",
                "PercentConfluency",
            )
            if c in cyto.columns
        ]
        if keep:
            parts.append(cyto[keep].copy())

    if raw_reads_per_well:
        rr_rows = [
            {
                "WellLocation": w.get("WellLocation"),
                "NumPoloniesBeforeTrimming": w.get("NumPoloniesBeforeTrimming"),
                **_well_raw_read_rates(w),
            }
            for w in raw_reads_per_well
            if w.get("WellLocation")
        ]
        if rr_rows:
            parts.append(pd.DataFrame(rr_rows))

    if not parts:
        return pd.DataFrame(columns=["WellLocation", "WellLabel"])

    out = parts[0]
    for extra in parts[1:]:
        out = out.merge(extra, on="WellLocation", how="outer")
    return out


def _well_raw_read_rates(well: dict[str, Any]) -> dict[str, float | None]:
    """Per-well pre-trim polony disposition from synthesized ``Wells[].Batches[]``.

    Surfaces *where a well's polonies land* — polonies/cell and the
    extracellular fraction — so an imaging-healthy well whose counts never
    convert into in-cell signal stands out against its pre-trim polony count
    (the "alive but empty" case).  Mirrors the run-level
    ``Cells2StatsStats.raw_read_aggregates_for_batches`` math per well:
    polonies/cell sums each batch's ``NumPoloniesPerCellBeforeTrimming``
    (batches deposit onto the same cells); the extracellular ratio is a
    count-weighted mean.  Returns ``None`` for a rate absent from the
    payload so the column drops cleanly when no well carries it.
    """
    reads_per_cell = 0.0
    have_rpc = False
    ec_num = 0.0
    ec_den = 0.0
    for batch in well.get("Batches") or []:
        mpc = _safe_float(batch.get("NumPoloniesPerCellBeforeTrimming"))
        if mpc is not None:
            reads_per_cell += mpc
            have_rpc = True
        rc = _safe_float(batch.get("NumPoloniesBeforeTrimming"))
        er = _safe_float(batch.get("ExtraCellularRatioBeforeTrimming"))
        if rc is not None and er is not None:
            ec_num += rc * er
            ec_den += rc
    return {
        "ReadsPerCell": reads_per_cell if have_rpc else None,
        "ExtracellularRatio": (ec_num / ec_den) if ec_den > 0 else None,
    }


def _sum_or_none(values: Any) -> float | None:
    nums = [_safe_float(v) for v in values]
    nums = [n for n in nums if n is not None]
    if not nums:
        return None
    return sum(nums)


def _weighted_mean(pairs: list[tuple[Any, Any]]) -> float | None:
    cleaned = [(_safe_float(v), _safe_float(w)) for v, w in pairs]
    cleaned = [(v, w) for v, w in cleaned if v is not None and w is not None and w > 0]
    if not cleaned:
        # Fall back to unweighted mean over present values
        plain = [_safe_float(v) for v, _ in pairs]
        plain = [n for n in plain if n is not None]
        return sum(plain) / len(plain) if plain else None
    total_w = sum(w for _, w in cleaned)
    if total_w == 0:
        return None
    return sum(v * w for v, w in cleaned) / total_w


def _safe_float(v: Any) -> float | None:
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    if math.isnan(f):
        return None
    return f


__all__ = [
    "b2f_for_batches",
    "join_batch_metadata",
    "per_batch_type_summary",
    "per_tile_cell_count",
    "per_tile_total_counts",
    "per_well_overview",
]
