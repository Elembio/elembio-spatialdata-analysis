"""Core filtering functions: filter_cells, _base_filter_cells, filter_genes, compute_filter_diagnostics, get_filter_flags."""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Any

import anndata as ad
import pandas as pd

from elembio_spatialdata_analysis.pp.filter._constants import (
    _UNSET,
    FILTER_FLAGS_FULL_UNS_KEY,
    FILTER_FLAGS_OBSM_KEY,
    FILTER_PARAMS_UNS_KEY,
    QC_METRICS_OBSM_KEY,
)
from elembio_spatialdata_analysis.pp.filter._flags import (
    _flag_assignment_rate,
    _flag_cell_area,
    _flag_empty_cells,
    _flag_failed_batches,
    _flag_gapdh_outliers,
    _flag_genes_per_well,
    _flag_imaging_intensity,
    _flag_mt_outliers,
    _flag_nsb_outliers,
    _flag_spatial_artifacts,
    _flag_total_counts_per_batch,
    _flag_total_counts_summed,
    _get_primary_intensity_features,
)
from elembio_spatialdata_analysis.pp.filter._qc_metrics import _add_obsm_qc_matrix
from elembio_spatialdata_analysis.pp.filter._quantile import (
    QuantileSpec,
    _resolve_threshold_defaults,
)
from elembio_spatialdata_analysis.pp.filter._spatialdata import _update_spatialdata_attrs
from elembio_spatialdata_analysis.tl.constants import (
    SUPPORTED_TARGET_TYPES,
    normalize_target_type,
)

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _subset_obs_only(adata: ad.AnnData, mask: pd.Series) -> ad.AnnData:
    """Lightweight subset that shares X/layers rather than copying them.

    Flag functions only need ``obs``, ``obsm``, and ``var`` — never the
    expression matrix itself.  Using a view avoids materialising a full copy
    of X just to compute per-well quantiles or assignment rates.
    """
    idx = adata.obs_names[mask.reindex(adata.obs_names).fillna(False).astype(bool)]
    view = adata[idx]
    # AnnData views are backed by the parent matrix — no copy of X occurs.
    return view


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class FilterDiagnostics:
    """Container for filter results and the configuration used to compute them."""

    flags: pd.DataFrame
    thresholds: dict[str, Any]
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dataframe(self) -> pd.DataFrame:
        """Return a copy of the filter flags dataframe."""
        return self.flags.copy()


@dataclass(slots=True)
class ModalityFilterResult:
    """Pack filtered AnnData with diagnostics for a given modality."""

    adata: ad.AnnData
    diagnostics: FilterDiagnostics


# ---------------------------------------------------------------------------
# Modality-specific flag dispatchers
# ---------------------------------------------------------------------------


def _add_spatial_artifact_flag(
    filter_flags: dict[str, pd.Series],
    adata: ad.AnnData,
    *,
    target_type: str,
    enabled: bool,
    original_index: pd.Index,
) -> None:
    """Write ``passes_spatial_artifacts`` from the original (unsubset) adata."""
    try:
        artifact_flag = _flag_spatial_artifacts(
            adata,
            target_type=target_type,
            enabled=enabled,
        )
        filter_flags["passes_spatial_artifacts"] = artifact_flag.reindex(
            original_index, fill_value=not enabled
        )
    except Exception as e:
        _warn_flag_skipped(e)


def _warn_flag_skipped(exc: Exception) -> None:
    """Surface a swallowed QC-flag computation error.

    The modality flag builders compute each flag best-effort so one failure does
    not abort the whole filter. When a flag cannot be computed it is dropped from
    the result, meaning affected cells are *not* gated on that criterion, so the
    failure must be loud (a warning), never a silent ``print``.
    """
    warnings.warn(
        f"A QC filter flag could not be computed and was skipped ({exc}); "
        "affected cells are not gated on that criterion.",
        stacklevel=2,
    )


def _get_transcript_filter_flags(
    adata: ad.AnnData,
    *,
    batch_colname: str = "batch_name",
    well_colname: str = "Well",
    area_colname: str = "Area",
    min_total_count: float = 1.0,
    min_quant_count: QuantileSpec = 0.025,
    max_quant_count: QuantileSpec = None,
    min_quant_area: float = 0.01,
    max_quant_area: float = 0.99,
    max_quant_nsb: float = 0.99,
    max_quant_gapdh: QuantileSpec = None,
    min_assignment_rate: float | None = None,
    filter_spatial_artifacts: bool = False,
) -> pd.DataFrame:
    assert adata.obs_names.is_unique, "Cell names are not unique"
    assert adata.var_names.is_unique, "Gene names are not unique"

    original_index = adata.obs_names
    filter_flags: dict[str, pd.Series] = {}

    try:
        nan_flag = _flag_failed_batches(adata, batch_colname=batch_colname)
    except Exception as e:
        _warn_flag_skipped(e)
        nan_flag = pd.Series(False, index=original_index, name="passes_nan_filter")
    filter_flags["passes_nan_filter"] = nan_flag

    try:
        empty_flag = _flag_empty_cells(adata, min_total_counts=min_total_count)
        filter_flags["passes_empty_filter"] = empty_flag
    except Exception as e:
        _warn_flag_skipped(e)
        empty_flag = pd.Series(False, index=original_index, name="passes_empty_filter")

    filtered_adata = _subset_obs_only(adata, empty_flag & nan_flag)

    try:
        counts_flag = _flag_total_counts_per_batch(
            filtered_adata, batch_colname, well_colname, min_quant_count, max_quant_count
        )
        filter_flags["passes_counts"] = counts_flag.reindex(original_index, fill_value=False)
    except Exception as e:
        _warn_flag_skipped(e)

    try:
        assignment_flag = _flag_assignment_rate(filtered_adata, min_assignment_rate=min_assignment_rate)
        filter_flags["passes_assignment_rate"] = assignment_flag.reindex(original_index, fill_value=False)
    except Exception as e:
        _warn_flag_skipped(e)

    try:
        area_flag = _flag_cell_area(
            filtered_adata,
            well_colname=well_colname,
            area_colname=area_colname,
            min_quant_area=min_quant_area,
            max_quant_area=max_quant_area,
        )
        filter_flags["passes_area"] = area_flag.reindex(original_index, fill_value=False)
    except Exception as e:
        _warn_flag_skipped(e)

    try:
        nsb_flag = _flag_nsb_outliers(
            filtered_adata,
            batch_colname=batch_colname,
            well_colname=well_colname,
            max_quantile=max_quant_nsb,
        )
        filter_flags["passes_nsb_outlier"] = nsb_flag.reindex(original_index, fill_value=False)
    except Exception as e:
        _warn_flag_skipped(e)

    try:
        gapdh_flag = _flag_gapdh_outliers(
            filtered_adata,
            batch_colname=batch_colname,
            well_colname=well_colname,
            max_quantile=max_quant_gapdh,
        )
        filter_flags["passes_gapdh_outlier"] = gapdh_flag.reindex(original_index, fill_value=False)
    except Exception as e:
        _warn_flag_skipped(e)

    _add_spatial_artifact_flag(
        filter_flags,
        adata,
        target_type="transcript",
        enabled=filter_spatial_artifacts,
        original_index=original_index,
    )

    flag_df = pd.DataFrame(filter_flags).reindex(original_index, fill_value=False)
    flag_df["passes_all_filters"] = flag_df.all(axis=1)
    return flag_df


def _get_protein_filter_flags(
    adata: ad.AnnData,
    *,
    batch_colname: str = "batch_name",
    well_colname: str = "Well",
    area_colname: str = "Area",
    min_total_count: float = 1.0,
    min_quant_count: QuantileSpec = 0.025,
    max_quant_count: QuantileSpec = None,
    min_quant_area: float = 0.01,
    max_quant_area: float = 0.99,
    min_assignment_rate: float | None = None,
    filter_spatial_artifacts: bool = False,
    spatial_artifact_target_type: str = "protein",
) -> pd.DataFrame:
    assert adata.obs_names.is_unique, "Cell names are not unique"
    assert adata.var_names.is_unique, "Feature names are not unique"

    original_index = adata.obs_names
    filter_flags: dict[str, pd.Series] = {}

    try:
        nan_flag = _flag_failed_batches(adata, batch_colname=batch_colname)
    except Exception as e:
        _warn_flag_skipped(e)
        nan_flag = pd.Series(False, index=original_index, name="passes_nan_filter")
    filter_flags["passes_nan_filter"] = nan_flag

    try:
        empty_flag = _flag_empty_cells(adata, min_total_counts=min_total_count)
        filter_flags["passes_empty_filter"] = empty_flag
    except Exception as e:
        _warn_flag_skipped(e)
        empty_flag = pd.Series(False, index=original_index, name="passes_empty_filter")

    filtered_adata = _subset_obs_only(adata, empty_flag & nan_flag)

    try:
        counts_flag = _flag_total_counts_per_batch(
            filtered_adata, batch_colname, well_colname, min_quant_count, max_quant_count
        )
        filter_flags["passes_counts"] = counts_flag.reindex(original_index, fill_value=False)
    except Exception as e:
        _warn_flag_skipped(e)

    try:
        assignment_flag = _flag_assignment_rate(filtered_adata, min_assignment_rate=min_assignment_rate)
        filter_flags["passes_assignment_rate"] = assignment_flag.reindex(original_index, fill_value=False)
    except Exception as e:
        _warn_flag_skipped(e)

    try:
        area_flag = _flag_cell_area(
            filtered_adata,
            well_colname=well_colname,
            area_colname=area_colname,
            min_quant_area=min_quant_area,
            max_quant_area=max_quant_area,
        )
        filter_flags["passes_area"] = area_flag.reindex(original_index, fill_value=False)
    except Exception as e:
        _warn_flag_skipped(e)

    _add_spatial_artifact_flag(
        filter_flags,
        adata,
        target_type=spatial_artifact_target_type,
        enabled=filter_spatial_artifacts,
        original_index=original_index,
    )

    flag_df = pd.DataFrame(filter_flags).reindex(original_index, fill_value=False)
    flag_df["passes_all_filters"] = flag_df.all(axis=1)
    return flag_df


def _get_imaging_filter_flags(
    adata: ad.AnnData,
    *,
    batch_colname: str = "batch_name",
    well_colname: str = "Well",
    area_colname: str = "Area",
    min_total_count: float = 1.0,
    min_quant_count: QuantileSpec = 0.025,
    max_quant_count: QuantileSpec = None,
    min_quant_area: float = 0.01,
    max_quant_area: float = 0.99,
    filter_spatial_artifacts: bool = False,
) -> pd.DataFrame:
    assert adata.obs_names.is_unique, "Cell names are not unique"
    assert adata.var_names.is_unique, "Feature names are not unique"

    original_index = adata.obs_names
    filter_flags: dict[str, pd.Series] = {}

    # Resolve morphology channels once for NaN + empty + intensity (same CP* set).
    primary_feature_names = _get_primary_intensity_features(adata)

    try:
        # Scope NaNs to morphology only — IF / CellProfiler / compartment NaNs
        # must not drop cells when primary morphology is intact.
        nan_flag = _flag_failed_batches(
            adata,
            batch_colname=batch_colname,
            feature_names=primary_feature_names,
        )
    except Exception as e:
        warnings.warn(
            f"Imaging NaN/failed-batch filter could not be computed ({e}); "
            "'passes_nan_filter' set to False (those cells treated as failing).",
            stacklevel=2,
        )
        nan_flag = pd.Series(False, index=original_index, name="passes_nan_filter")
    filter_flags["passes_nan_filter"] = nan_flag

    try:
        # Same morphology scope as NaN — a NaN in Bkg_* / IF must not make the
        # full-matrix row sum NaN and wipe the well via the empty gate.
        empty_flag = _flag_empty_cells(
            adata,
            min_total_counts=min_total_count,
            feature_names=primary_feature_names,
        )
        filter_flags["passes_empty_filter"] = empty_flag
    except Exception as e:
        warnings.warn(
            f"Imaging empty-cell filter could not be computed ({e}); "
            "'passes_empty_filter' set to False (those cells treated as failing).",
            stacklevel=2,
        )
        empty_flag = pd.Series(False, index=original_index, name="passes_empty_filter")

    filtered_adata = _subset_obs_only(adata, empty_flag & nan_flag)

    try:
        imaging_flag = _flag_imaging_intensity(
            filtered_adata,
            well_colname=well_colname,
            min_quant_count=min_quant_count,
            max_quant_count=max_quant_count,
            primary_feature_names=primary_feature_names,
        )
        filter_flags["passes_imaging_intensity"] = imaging_flag.reindex(original_index, fill_value=False)
    except Exception as e:
        warnings.warn(
            f"Imaging intensity gate could not be computed ({e}); "
            "'passes_imaging_intensity' was NOT applied (those cells are not gated on intensity).",
            stacklevel=2,
        )

    try:
        area_flag = _flag_cell_area(
            filtered_adata,
            well_colname=well_colname,
            area_colname=area_colname,
            min_quant_area=min_quant_area,
            max_quant_area=max_quant_area,
        )
        filter_flags["passes_area"] = area_flag.reindex(original_index, fill_value=False)
    except Exception as e:
        warnings.warn(
            f"Imaging cell-area filter could not be computed ({e}); "
            "'passes_area' was NOT applied (those cells are not gated on area).",
            stacklevel=2,
        )

    _add_spatial_artifact_flag(
        filter_flags,
        adata,
        target_type="imaging",
        enabled=filter_spatial_artifacts,
        original_index=original_index,
    )

    flag_df = pd.DataFrame(filter_flags).reindex(original_index, fill_value=False)
    flag_df["passes_all_filters"] = flag_df.all(axis=1)
    return flag_df


# Backward-compat alias
_get_cell_paint_filter_flags = _get_imaging_filter_flags


def _get_threeprimeuntargeted_filter_flags(
    adata: ad.AnnData,
    *,
    batch_colname: str = "batch_name",
    well_colname: str = "Well",
    area_colname: str = "Area",
    min_total_count: float = 1.0,
    min_quant_count: QuantileSpec = 0.025,
    max_quant_count: QuantileSpec = None,
    min_quant_area: float = 0.01,
    max_quant_area: float = 0.99,
    max_pct_mt: float | None = None,
    min_assignment_rate: float | None = None,
    filter_spatial_artifacts: bool = False,
) -> pd.DataFrame:
    """Filter flags for ThreePrimeUntargeted data (no NSB/GAPDH filters)."""
    assert adata.obs_names.is_unique, "Cell names are not unique"
    assert adata.var_names.is_unique, "Gene names are not unique"

    original_index = adata.obs_names
    filter_flags: dict[str, pd.Series] = {}

    try:
        nan_flag = _flag_failed_batches(adata, batch_colname=batch_colname)
    except Exception as e:
        _warn_flag_skipped(e)
        nan_flag = pd.Series(False, index=original_index, name="passes_nan_filter")
    filter_flags["passes_nan_filter"] = nan_flag

    try:
        empty_flag = _flag_empty_cells(adata, min_total_counts=min_total_count)
        filter_flags["passes_empty_filter"] = empty_flag
    except Exception as e:
        _warn_flag_skipped(e)
        empty_flag = pd.Series(False, index=original_index, name="passes_empty_filter")

    filtered_adata = _subset_obs_only(adata, empty_flag & nan_flag)

    try:
        # ThreePrimeUntargeted batches are summed into a single expression matrix,
        # so filter on summed total counts per well rather than per-batch counts.
        # Per-batch filtering incorrectly drops cells with zero expression in some batches.
        counts_flag = _flag_total_counts_summed(
            filtered_adata,
            well_colname=well_colname,
            min_quant_count=min_quant_count,
            max_quant_count=max_quant_count,
        )
        filter_flags["passes_counts"] = counts_flag.reindex(original_index, fill_value=False)
    except (ValueError, TypeError):
        # Invalid threshold specs are caller errors; don't downgrade them to a skipped flag.
        raise
    except Exception as e:
        _warn_flag_skipped(e)

    try:
        assignment_flag = _flag_assignment_rate(filtered_adata, min_assignment_rate=min_assignment_rate)
        filter_flags["passes_assignment_rate"] = assignment_flag.reindex(original_index, fill_value=False)
    except Exception:
        filter_flags["passes_assignment_rate"] = pd.Series(True, index=original_index)

    try:
        area_flag = _flag_cell_area(
            filtered_adata,
            well_colname=well_colname,
            area_colname=area_colname,
            min_quant_area=min_quant_area,
            max_quant_area=max_quant_area,
        )
        filter_flags["passes_area"] = area_flag.reindex(original_index, fill_value=False)
    except Exception as e:
        _warn_flag_skipped(e)

    try:
        mt_flag = _flag_mt_outliers(filtered_adata, well_colname=well_colname, max_pct_mt=max_pct_mt)
        filter_flags["passes_mt_outlier"] = mt_flag.reindex(original_index, fill_value=False)
    except Exception as e:
        _warn_flag_skipped(e)

    _add_spatial_artifact_flag(
        filter_flags,
        adata,
        target_type="threeprimeuntargeted",
        enabled=filter_spatial_artifacts,
        original_index=original_index,
    )

    flag_df = pd.DataFrame(filter_flags).reindex(original_index, fill_value=False)
    flag_df["passes_all_filters"] = flag_df.all(axis=1)
    return flag_df


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def get_filter_flags(
    adata: ad.AnnData,
    *,
    target_type: str,
    batch_colname: str = "batch_name",
    well_colname: str = "Well",
    area_colname: str = "Area",
    min_total_count: float = 1.0,
    min_quant_count: QuantileSpec = 0.025,
    max_quant_count: QuantileSpec = None,
    min_quant_area: QuantileSpec = 0.01,
    max_quant_area: QuantileSpec = 0.99,
    max_quant_nsb: QuantileSpec = 0.99,
    max_quant_gapdh: QuantileSpec = None,
    max_pct_mt: float | None = None,
    min_assignment_rate: float | None = None,
    filter_spatial_artifacts: bool = False,
) -> pd.DataFrame:
    """Dispatch to the correct modality-specific flag function."""
    normalized_target_type = normalize_target_type(target_type)

    if normalized_target_type == "transcript":
        return _get_transcript_filter_flags(
            adata=adata,
            batch_colname=batch_colname,
            well_colname=well_colname,
            area_colname=area_colname,
            min_total_count=min_total_count,
            min_quant_count=min_quant_count,
            max_quant_count=max_quant_count,
            min_quant_area=min_quant_area,
            max_quant_area=max_quant_area,
            max_quant_nsb=max_quant_nsb,
            max_quant_gapdh=max_quant_gapdh,
            min_assignment_rate=min_assignment_rate,
            filter_spatial_artifacts=filter_spatial_artifacts,
        )

    if normalized_target_type == "threeprimeuntargeted":
        return _get_threeprimeuntargeted_filter_flags(
            adata=adata,
            batch_colname=batch_colname,
            well_colname=well_colname,
            area_colname=area_colname,
            min_total_count=min_total_count,
            min_quant_count=min_quant_count,
            max_quant_count=max_quant_count,
            min_quant_area=min_quant_area,
            max_quant_area=max_quant_area,
            max_pct_mt=max_pct_mt,
            min_assignment_rate=min_assignment_rate,
            filter_spatial_artifacts=filter_spatial_artifacts,
        )

    if normalized_target_type == "protein":
        return _get_protein_filter_flags(
            adata=adata,
            batch_colname=batch_colname,
            well_colname=well_colname,
            area_colname=area_colname,
            min_total_count=min_total_count,
            min_quant_count=min_quant_count,
            max_quant_count=max_quant_count,
            min_quant_area=min_quant_area,
            max_quant_area=max_quant_area,
            min_assignment_rate=min_assignment_rate,
            filter_spatial_artifacts=filter_spatial_artifacts,
        )

    if normalized_target_type == "imaging":
        return _get_imaging_filter_flags(
            adata=adata,
            batch_colname=batch_colname,
            well_colname=well_colname,
            area_colname=area_colname,
            min_total_count=min_total_count,
            min_quant_count=min_quant_count,
            max_quant_count=max_quant_count,
            min_quant_area=min_quant_area,
            max_quant_area=max_quant_area,
            filter_spatial_artifacts=filter_spatial_artifacts,
        )

    if normalized_target_type == "ops":
        # OPS shares Protein's obsm shape but lacks ``unassigned_counts``
        # by design — force ``min_assignment_rate=None`` so the
        # protein dispatcher doesn't try to gate on a slot that will be all
        # NaN.  See ``_get_ops_filter_flags`` for the formal entry point.
        return _get_ops_filter_flags(
            adata=adata,
            batch_colname=batch_colname,
            well_colname=well_colname,
            area_colname=area_colname,
            min_total_count=min_total_count,
            min_quant_count=min_quant_count,
            max_quant_count=max_quant_count,
            min_quant_area=min_quant_area,
            max_quant_area=max_quant_area,
            filter_spatial_artifacts=filter_spatial_artifacts,
        )

    raise ValueError(
        f"Unrecognized target_type '{target_type}'. Supported types: "
        "Transcript, Protein, Imaging (or CellPaint), ThreePrimeUntargeted, and OPS."
    )


def _get_ops_filter_flags(
    *,
    adata: ad.AnnData,
    batch_colname: str,
    well_colname: str | None,
    area_colname: str,
    min_total_count: int | None,
    min_quant_count: QuantileSpec,
    max_quant_count: QuantileSpec,
    min_quant_area: QuantileSpec,
    max_quant_area: QuantileSpec,
    filter_spatial_artifacts: bool = False,
) -> pd.DataFrame:
    """OPS filter flags — count + area quantile bounds, no assignment gate.

    OPS adata structurally mirrors Protein (``total_counts`` per cell,
    ``cell_geometry`` with ``AreaUm``) but carries no
    ``unassigned_counts`` / ``total_assignment_rate`` by design.  This wrapper
    delegates to the Protein flag dispatcher with the assignment-rate
    gate force-disabled so we don't generate a column of NaN flags.

    Defined as a thin wrapper rather than inlined in the dispatch so
    that future OPS-specific filters (e.g. ``MostAbundantTargetCount /
    TotalCount`` quality gates from ``obsm['target_assignment']``) have
    a clear home.
    """
    return _get_protein_filter_flags(
        adata=adata,
        batch_colname=batch_colname,
        well_colname=well_colname,
        area_colname=area_colname,
        min_total_count=min_total_count,
        min_quant_count=min_quant_count,
        max_quant_count=max_quant_count,
        min_quant_area=min_quant_area,
        max_quant_area=max_quant_area,
        min_assignment_rate=None,
        filter_spatial_artifacts=filter_spatial_artifacts,
        spatial_artifact_target_type="ops",
    )


def _require_target_type(target_type: str | None) -> None:
    if target_type is None:
        supported = ", ".join(sorted(SUPPORTED_TARGET_TYPES))
        raise ValueError(
            f"`target_type` is required. Choose one of: {supported}. "
            "Aliases such as 'RNA' or 'CellPaint' are also accepted."
        )


def compute_filter_diagnostics(
    adata: ad.AnnData,
    target_type: str | None = None,
    *,
    batch_colname: str = "batch_name",
    well_colname: str = "Well",
    area_colname: str = "Area",
    min_quant_count: QuantileSpec = _UNSET,
    max_quant_count: QuantileSpec = _UNSET,
    min_quant_area: QuantileSpec = _UNSET,
    max_quant_area: QuantileSpec = _UNSET,
    max_quant_nsb: QuantileSpec = _UNSET,
    max_quant_gapdh: QuantileSpec = _UNSET,
    max_pct_mt: float | None = _UNSET,
    min_assignment_rate: float | None = _UNSET,
    filter_spatial_artifacts: bool = _UNSET,
    min_cells_per_gene: int | None = None,
    min_counts_per_gene: float | None = None,
    inplace: bool = True,
    return_diagnostics: bool = False,
) -> FilterDiagnostics | None:
    """Compute QC metrics and filter flags without subsetting ``adata``.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object.
    target_type : str
        Modality identifier (e.g. ``"Transcript"``, ``"Protein"``, ``"Imaging"``).
    batch_colname : str
        Column in adata.obs identifying feature batches.
    well_colname : str
        Column in adata.obs identifying well assignments.
    area_colname : str
        Column in adata.obs containing cell area measurements.
    min_quant_count, max_quant_count, min_quant_area, max_quant_area,
    max_quant_nsb, max_quant_gapdh, max_pct_mt, min_assignment_rate,
    filter_spatial_artifacts :
        See ``filter_cells`` for full documentation.
    min_cells_per_gene : int, optional
        Minimum cells per gene per well.
    min_counts_per_gene : float, optional
        Minimum total counts per gene per well.
    inplace : bool
        If True, store flags in ``adata.obsm`` in place.
    return_diagnostics : bool
        If True, return a ``FilterDiagnostics`` object; otherwise return None.

    Returns
    -------
    FilterDiagnostics or None

    Raises
    ------
    ValueError
        If ``target_type`` is not provided.
    """
    _require_target_type(target_type)

    normalized_target_type = normalize_target_type(target_type)
    threshold_overrides = {
        "min_quant_count": min_quant_count,
        "max_quant_count": max_quant_count,
        "min_quant_area": min_quant_area,
        "max_quant_area": max_quant_area,
        "max_quant_nsb": max_quant_nsb,
        "max_quant_gapdh": max_quant_gapdh,
        "max_pct_mt": max_pct_mt,
        "min_assignment_rate": min_assignment_rate,
        "filter_spatial_artifacts": filter_spatial_artifacts,
        "min_cells_per_gene": min_cells_per_gene,
        "min_counts_per_gene": min_counts_per_gene,
    }
    resolved_thresholds = _resolve_threshold_defaults(normalized_target_type, threshold_overrides)

    # Skip recomputing if already present (user called compute_qc_metrics explicitly)
    if QC_METRICS_OBSM_KEY not in adata.obsm:
        _add_obsm_qc_matrix(adata, normalized_target_type, batch_colname)

    filter_flags = get_filter_flags(
        adata=adata,
        target_type=normalized_target_type,
        batch_colname=batch_colname,
        well_colname=well_colname,
        area_colname=area_colname,
        min_quant_count=resolved_thresholds["min_quant_count"],
        max_quant_count=resolved_thresholds["max_quant_count"],
        min_quant_area=resolved_thresholds["min_quant_area"],
        max_quant_area=resolved_thresholds["max_quant_area"],
        max_quant_nsb=resolved_thresholds["max_quant_nsb"],
        max_quant_gapdh=resolved_thresholds["max_quant_gapdh"],
        max_pct_mt=resolved_thresholds["max_pct_mt"],
        min_assignment_rate=resolved_thresholds["min_assignment_rate"],
        filter_spatial_artifacts=resolved_thresholds["filter_spatial_artifacts"],
    )

    if "passes_all_filters" not in filter_flags.columns:
        raise KeyError("Missing 'passes_all_filters' column in filter flags.")

    if not filter_flags.index.equals(adata.obs_names):
        filter_flags = filter_flags.reindex(adata.obs_names)

    keep = filter_flags["passes_all_filters"]
    if keep.isna().any():
        raise ValueError("`filter_flags` index does not match `adata.obs`.")

    thresholds: dict[str, Any] = {
        "target_type": target_type,
        "normalized_target_type": normalized_target_type,
        "batch_colname": batch_colname,
        "well_colname": well_colname,
        "area_colname": area_colname,
        "min_quant_count": resolved_thresholds["min_quant_count"],
        "max_quant_count": resolved_thresholds["max_quant_count"],
        "min_quant_area": resolved_thresholds["min_quant_area"],
        "max_quant_area": resolved_thresholds["max_quant_area"],
        "max_quant_nsb": resolved_thresholds["max_quant_nsb"],
        "max_quant_gapdh": resolved_thresholds["max_quant_gapdh"],
        "max_pct_mt": resolved_thresholds["max_pct_mt"],
        "min_assignment_rate": resolved_thresholds["min_assignment_rate"],
        "filter_spatial_artifacts": resolved_thresholds["filter_spatial_artifacts"],
        "min_cells_per_gene": resolved_thresholds.get("min_cells_per_gene"),
        "min_counts_per_gene": resolved_thresholds.get("min_counts_per_gene"),
        "inplace": inplace,
    }

    # Always write flags and params to the AnnData so downstream functions
    # (generate_filter_report, plotting) can read them without needing the
    # FilterDiagnostics object.
    adata.obsm[FILTER_FLAGS_OBSM_KEY] = filter_flags
    adata.uns[FILTER_PARAMS_UNS_KEY] = thresholds

    diagnostics = FilterDiagnostics(flags=filter_flags, thresholds=thresholds)

    if return_diagnostics:
        return diagnostics
    return None


def filter_cells(
    adata: ad.AnnData,
    target_type: str | None = None,
    *,
    batch_colname: str = "batch_name",
    well_colname: str = "Well",
    area_colname: str = "Area",
    min_quant_count: QuantileSpec = _UNSET,
    max_quant_count: QuantileSpec = _UNSET,
    min_quant_area: QuantileSpec = _UNSET,
    max_quant_area: QuantileSpec = _UNSET,
    max_quant_nsb: QuantileSpec = _UNSET,
    max_quant_gapdh: QuantileSpec = _UNSET,
    max_pct_mt: float | None = _UNSET,
    min_assignment_rate: float | None = _UNSET,
    filter_spatial_artifacts: bool = _UNSET,
    inplace: bool = True,
    diagnostics: bool = False,
) -> ad.AnnData | tuple[ad.AnnData, FilterDiagnostics]:
    """Filter cells in an AnnData object based on QC metrics.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object.
    target_type : str
        Modality identifier (e.g. ``"Transcript"``, ``"Protein"``, ``"Imaging"``).
    batch_colname : str
        Column in ``adata.obs`` identifying feature batches.
    well_colname : str
        Column in ``adata.obs`` identifying well assignments.
    area_colname : str
        Column in ``adata.obs`` containing cell area measurements.
    min_quant_count, max_quant_count : float or Mapping, optional
        Lower/upper quantile thresholds for per-cell counts within each batch.
    min_quant_area, max_quant_area : float or Mapping, optional
        Lower/upper quantile thresholds for cell area within each well.
    max_quant_nsb : float or Mapping, optional
        Upper quantile threshold for NSB control counts.
    max_quant_gapdh : float or Mapping, optional
        Upper quantile threshold for GAPDH positive-control counts (opt-in, default None).
    max_pct_mt : float, optional
        Maximum mitochondrial gene percentage (0–1).
    min_assignment_rate : float, optional
        Minimum fraction of assigned counts per cell.
    filter_spatial_artifacts : bool
        If True, drop cells flagged in the modality rollup column of
        ``obsm['spatial_artifacts']`` (``any-imaging``,
        ``any-ThreePrimeUntargeted``, or ``any-OPS``). Default False
        (flag written as all-True).
    inplace : bool
        If True, subset ``adata`` in place (default). If False, return a filtered copy
        leaving the original ``adata`` unchanged.
    diagnostics : bool
        If True, return ``(adata, FilterDiagnostics)``; otherwise return only ``adata``.

    Returns
    -------
    AnnData or (AnnData, FilterDiagnostics)

    Raises
    ------
    ValueError
        If ``target_type`` is not provided.
    """
    _require_target_type(target_type)

    normalized_target_type = normalize_target_type(target_type)
    threshold_overrides = {
        "min_quant_count": min_quant_count,
        "max_quant_count": max_quant_count,
        "min_quant_area": min_quant_area,
        "max_quant_area": max_quant_area,
        "max_quant_nsb": max_quant_nsb,
        "max_quant_gapdh": max_quant_gapdh,
        "max_pct_mt": max_pct_mt,
        "min_assignment_rate": min_assignment_rate,
        "filter_spatial_artifacts": filter_spatial_artifacts,
    }
    _resolve_threshold_defaults(normalized_target_type, threshold_overrides)

    # For inplace=False, make a single copy up front so that compute_filter_diagnostics
    # (which writes qc_metrics and filter_flags to obsm) does not mutate the caller's
    # adata. We then subset that copy in-place — no further full copies are needed.
    adata_work = adata if inplace else adata.copy()

    diagnostics_out = compute_filter_diagnostics(
        adata_work,
        target_type,
        batch_colname=batch_colname,
        well_colname=well_colname,
        area_colname=area_colname,
        min_quant_count=min_quant_count,
        max_quant_count=max_quant_count,
        min_quant_area=min_quant_area,
        max_quant_area=max_quant_area,
        max_quant_nsb=max_quant_nsb,
        max_quant_gapdh=max_quant_gapdh,
        max_pct_mt=max_pct_mt,
        min_assignment_rate=min_assignment_rate,
        filter_spatial_artifacts=filter_spatial_artifacts,
        inplace=True,
        return_diagnostics=True,
    )

    filter_flags = diagnostics_out.flags
    keep_mask = filter_flags["passes_all_filters"].to_numpy(dtype=bool)

    # Persist the full pre-subset flags in uns so that generate_filter_report(adata)
    # works correctly after inplace subsetting (obsm is trimmed, uns survives intact).
    adata_work.uns[FILTER_FLAGS_FULL_UNS_KEY] = filter_flags

    # Subset cells in-place on adata_work (which is either adata itself or its copy).
    adata_work._inplace_subset_obs(keep_mask)
    # Trim obsm DataFrames that were written by compute_filter_diagnostics to match new shape.
    adata_work.obsm[FILTER_FLAGS_OBSM_KEY] = filter_flags.loc[adata_work.obs_names]
    if QC_METRICS_OBSM_KEY in adata_work.obsm:
        adata_work.obsm[QC_METRICS_OBSM_KEY] = adata_work.obsm[QC_METRICS_OBSM_KEY].loc[adata_work.obs_names]

    _update_spatialdata_attrs(adata_work)

    if diagnostics:
        return adata_work, diagnostics_out
    return adata_work


def filter_genes(
    adata: ad.AnnData,
    target_type: str | None = None,
    *,
    well_colname: str = "Well",
    min_cells_per_gene: int | None = _UNSET,
    min_counts_per_gene: float | None = _UNSET,
    inplace: bool = True,
) -> ad.AnnData:
    """Filter genes from an AnnData object based on per-well expression thresholds.

    For ``ThreePrimeUntargeted`` data the thresholds default to the values in
    ``DEFAULT_BOUNDS_PER_MODALITY``; for all other modalities the defaults are
    ``None`` (no gene filtering unless explicitly provided).

    Parameters
    ----------
    adata : AnnData
        Input AnnData object (typically already cell-filtered via ``filter_cells``).
    target_type : str, optional
        Modality identifier used to look up default thresholds. Pass ``None`` to
        disable default lookup and rely entirely on the explicit parameters.
    well_colname : str
        Column in ``adata.obs`` identifying well assignments.
    min_cells_per_gene : int, optional
        Minimum number of cells expressing each gene within at least one well.
        Genes that do not meet this threshold in any well are removed.
    min_counts_per_gene : float, optional
        Minimum total counts per gene within at least one well.
        Genes that do not meet this threshold in any well are removed.
    inplace : bool
        If True, subset ``adata`` in place (default). If False, return a filtered copy.

    Returns
    -------
    AnnData
        AnnData with low-quality genes removed.
    """
    if target_type is not None:
        normalized_target_type = normalize_target_type(target_type)
        threshold_overrides = {
            "min_cells_per_gene": min_cells_per_gene,
            "min_counts_per_gene": min_counts_per_gene,
        }
        resolved = _resolve_threshold_defaults(normalized_target_type, threshold_overrides)
        resolved_min_cells = resolved.get("min_cells_per_gene")
        resolved_min_counts = resolved.get("min_counts_per_gene")
    else:
        resolved_min_cells = None if min_cells_per_gene is _UNSET else min_cells_per_gene
        resolved_min_counts = None if min_counts_per_gene is _UNSET else min_counts_per_gene

    adata_work = adata if inplace else adata.copy()

    gene_keep = _flag_genes_per_well(
        adata_work,
        well_colname=well_colname,
        min_cells_per_gene=resolved_min_cells,
        min_counts_per_gene=resolved_min_counts,
    )
    adata_work._inplace_subset_var(gene_keep.to_numpy(dtype=bool))
    _update_spatialdata_attrs(adata_work)

    return adata_work


def _base_filter_cells(
    adata: ad.AnnData,
    *,
    min_total_count: float = 1.0,
    inplace: bool = True,
    diagnostics: bool = False,
):
    """Minimal filtering: drop cells with NaN values or total counts below ``min_total_count``.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object.
    min_total_count : float
        Minimum total counts required to keep a cell.
    inplace : bool
        If True, mutates ``adata`` in place.
    diagnostics : bool
        If True, return ``(adata, filter_flags)``; otherwise return only ``adata``.

    Returns
    -------
    AnnData or (AnnData, pandas.DataFrame)
    """
    # Flags are computed on the original adata (read-only operations, no obsm writes).
    # The copy (if needed) is made once and subset in-place — avoids a second full copy.
    nan_flag = _flag_failed_batches(adata)
    empty_flag = _flag_empty_cells(adata, min_total_counts=min_total_count)

    flags = pd.DataFrame(
        {
            "passes_nan_filter": nan_flag,
            "passes_empty_filter": empty_flag,
        }
    )
    flags["passes_all_filters"] = flags.all(axis=1)

    keep_mask = flags["passes_all_filters"].to_numpy(dtype=bool)

    adata_work = adata if inplace else adata.copy()
    adata_work._inplace_subset_obs(keep_mask)
    _update_spatialdata_attrs(adata_work)

    if diagnostics:
        return adata_work, flags
    return adata_work
