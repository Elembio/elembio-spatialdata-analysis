"""Identify Cell Painting morphology channels for the imaging QC gate.

The imaging intensity gate must run on the structural Cell Painting stains
only, never on the protein / antibody panel, CellProfiler-derived features,
compartment variants, or background channels.

The reliable discriminator is the CellPaint *batch* (``batch_name`` starting
with ``"CP"``), not the per-feature ``target_type``.  On real IF panels every
imaging target is declared ``target_type == "CellPaint"`` in the panel (the
antibody markers, the phospho markers, even the IgG isotype controls), so
``target_type`` alone cannot separate structural stains from markers.  A
product build always ships the stains in a dedicated ``CP`` batch, so we anchor
on the batch and adapt to whatever ``var`` metadata the store format carries.
"""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

# CellProfiler-derived feature-name prefixes — used only when the store carries
# neither ``feature_class`` nor ``target_type`` to classify features.
_CELLPROFILER_NAME_RE = (
    r"^(?:Texture|Granularity|RadialDistribution|Location|AreaShape|Correlation|Intensity)_"
)
# Compartment / threshold variants of a base channel (drop; keep whole-cell only).
_VARIANT_SUFFIX_RE = r"_(?:Nuclear|Cytoplasmic|FractionAboveThreshold|MeanIntensityOfFractionAboveThreshold)$"
# Background-channel names (e.g. ``Bkg_RED``).
_BACKGROUND_NAME_RE = r"(?i)bkg|background"


def select_morphology_channels(
    var: pd.DataFrame,
    *,
    cp_batch_prefixes: Sequence[str] = ("CP",),
) -> list[str]:
    """Return the whole-cell Cell Painting morphology channels named in ``var``.

    These are the structural stains the imaging intensity gate should run on.
    The protein / antibody batch, CellProfiler-derived features, compartment
    variants, and background channels are all excluded.

    Selection (adaptive to store format):

    1. **CellPaint morphology batch(es) only** — ``batch_name`` starts with one
       of *cp_batch_prefixes* (default ``"CP"``).  This is what excludes the
       protein / antibody batch.
    2. **Drop CellProfiler-derived features** — ``feature_class == "base"`` when
       that column is present (current ``global-table-imaging`` format), else
       ``target_type == "CellPaint"`` (legacy ``global-table-CellPaint`` format),
       else a name-prefix fallback.
    3. **Whole-cell only** — ``suffix == ""`` when present, else drop the
       ``_Nuclear`` / ``_Cytoplasmic`` / threshold name variants.
    4. **Drop background channels** — ``is_background == False`` when present,
       else drop ``Bkg_*`` by name.

    Parameters
    ----------
    var : pandas.DataFrame
        An ``AnnData.var`` for an imaging table.
    cp_batch_prefixes : sequence of str
        Batch-name prefixes identifying the Cell Painting morphology batch(es).

    Returns
    -------
    list[str]
        The morphology channel names.  Empty (a deliberate, fail-safe no-op)
        when ``var`` has no ``batch_name`` column or no CellPaint batch — the
        gate must never fall back to the protein batch.
    """
    # Fail safe: without batch identity we cannot isolate the morphology batch.
    if "batch_name" not in var.columns:
        return []

    names = var.index.to_series()
    batch = var["batch_name"].astype("string").fillna("")
    prefixes = tuple(p.upper() for p in cp_batch_prefixes)
    mask = batch.str.upper().str.startswith(prefixes)

    # 2. Drop CellProfiler-derived features.
    if "feature_class" in var.columns:
        mask &= var["feature_class"] == "base"
    elif "target_type" in var.columns:
        mask &= var["target_type"] == "CellPaint"
    else:
        mask &= ~names.str.match(_CELLPROFILER_NAME_RE)

    # 3. Whole-cell only (drop compartment / threshold variants).
    if "suffix" in var.columns:
        mask &= var["suffix"].fillna("") == ""
    else:
        mask &= ~names.str.contains(_VARIANT_SUFFIX_RE, regex=True)

    # 4. Drop background channels.
    if "is_background" in var.columns:
        mask &= ~var["is_background"].fillna(False)
    else:
        mask &= ~names.str.contains(_BACKGROUND_NAME_RE, regex=True)

    return var.index[mask.to_numpy()].tolist()
