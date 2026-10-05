"""Individual filter flag computation functions.

Each ``_flag_*`` function returns a boolean ``pd.Series`` indexed by cell names:
  - ``True``  → cell passes the filter (keep)
  - ``False`` → cell fails the filter (drop)

Data access conventions:
- Total counts per batch      → adata.obsm['total_counts']      (DataFrame)
- Unassigned counts per batch → adata.obsm['unassigned_counts'] (DataFrame)
- NSB / control variables     → adata.obsm['control_variables'] (DataFrame)
- GAPDH positive controls     → adata.X  (var['control_type'] == 'Positive Control 1')
- Cell geometry (Area, etc.)  → adata.obsm['cell_geometry']     (DataFrame, fallback: adata.obs)
- Data missing flags          → adata.obsm['data_missing']      (DataFrame)
- Spatial artifacts           → adata.obsm['spatial_artifacts'] (boolean DataFrame)
"""

from __future__ import annotations

import logging
import warnings
from collections.abc import Mapping, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from elembio_spatialdata_analysis.pp.filter._constants import (
    DEFAULT_BOUNDS_PER_MODALITY,
    OBSM_CELL_GEOMETRY,
    OBSM_CONTROL_VARIABLES,
    OBSM_DATA_MISSING,
    OBSM_SPATIAL_ARTIFACTS,
    OBSM_TOTAL_COUNTS,
    OBSM_UNASSIGNED_COUNTS,
    QC_METRICS_OBSM_KEY,
    SPATIAL_ARTIFACT_ROLLUP,
    VAR_CONTROL_TYPE_COL,
    VAR_GAPDH_CONTROL_TYPE,
)
from elembio_spatialdata_analysis.pp.filter._qc_metrics import (
    _col_sum,
    _obs_col_to_series,
    _obsm_to_df,
    _row_sum,
)
from elembio_spatialdata_analysis.pp.filter._quantile import (
    QuantileSpec,
    _resolve_well_batch_quantiles,
)

_logger = logging.getLogger(__name__)


def _groupby_apply_to_cell_series(result: pd.Series | pd.DataFrame) -> pd.Series:
    """Normalize ``groupby(...).apply`` results (pandas single-group quirk).

    Multiple groups → Series with MultiIndex ``(well, cell_id)`` — drop the well level.
    Single group → a one-row DataFrame ``(1, n_obs)`` — take ``iloc[0]``.
    """
    if isinstance(result, pd.DataFrame):
        if result.shape[0] == 1:
            out = result.iloc[0]
            if not isinstance(out, pd.Series):
                raise TypeError("Expected one row DataFrame to yield a Series.")
            return out
        raise ValueError(f"Unexpected groupby-apply shape: {result.shape}")

    if isinstance(result.index, pd.MultiIndex) and result.index.nlevels > 1:
        return result.droplevel(0)
    return result


# ---------------------------------------------------------------------------
# Structural / batch quality flags
# ---------------------------------------------------------------------------


def _flag_failed_batches(
    adata: ad.AnnData,
    batch_colname: str = "batch_name",
    feature_names: Sequence[str] | None = None,
) -> pd.Series:
    """Flag cells that have NaN values in ``adata.X`` (failed batches).

    Parameters
    ----------
    feature_names
        Optional column names to check. ``None`` (default) checks the full
        matrix. An empty sequence is a no-op (all cells pass) — used by the
        imaging path when no CellPaint morphology channels are present.
    """
    if feature_names is not None:
        if len(feature_names) == 0:
            return pd.Series(True, index=adata.obs_names, name="passes_nan_filter")
        missing = [n for n in feature_names if n not in adata.var_names]
        if missing:
            raise KeyError(f"feature_names not in adata.var_names: {missing}")
        col_indices = [adata.var_names.get_loc(name) for name in feature_names]
        X = adata.X[:, col_indices]
    else:
        X = adata.X
    if sp.issparse(X):
        X = sp.csr_matrix(X)
        nan_in_data = np.isnan(X.data)
        if not np.any(nan_in_data):
            return pd.Series(True, index=adata.obs_names, name="passes_nan_filter")
        has_nan = np.zeros(X.shape[0], dtype=bool)
        for i in range(X.shape[0]):
            start, end = X.indptr[i], X.indptr[i + 1]
            if np.any(nan_in_data[start:end]):
                has_nan[i] = True
        return pd.Series(~has_nan, index=adata.obs_names, name="passes_nan_filter")
    has_nan = np.isnan(np.asarray(X)).any(axis=1)
    return pd.Series(~has_nan, index=adata.obs_names, name="passes_nan_filter")


def _flag_data_missing(adata: ad.AnnData) -> pd.Series:
    """Flag cells that have any data-missing indicators in ``adata.obsm['data_missing']``.

    Returns ``True`` (pass) for cells where no missing-data flag is set.
    Returns all-True when the obsm key is absent.
    """
    if OBSM_DATA_MISSING not in adata.obsm:
        return pd.Series(True, index=adata.obs_names, name="passes_data_missing")

    dm = _obsm_to_df(adata.obsm[OBSM_DATA_MISSING], adata.obs_names)

    has_missing = dm.any(axis=1)
    return (~has_missing).rename("passes_data_missing")


def _flag_empty_cells(
    adata: ad.AnnData,
    min_total_counts: float = 1.0,
    feature_names: Sequence[str] | None = None,
) -> pd.Series:
    """Flag cells whose total count in ``adata.X`` is below ``min_total_counts``.

    Parameters
    ----------
    feature_names
        Optional column names to sum. ``None`` (default) sums the full matrix.
        An empty sequence is a no-op (all cells pass) — used by the imaging
        path when no CellPaint morphology channels are present.
    """
    if feature_names is not None:
        if len(feature_names) == 0:
            return pd.Series(True, index=adata.obs_names, name="passes_empty_filter")
        missing = [n for n in feature_names if n not in adata.var_names]
        if missing:
            raise KeyError(f"feature_names not in adata.var_names: {missing}")
        col_indices = [adata.var_names.get_loc(name) for name in feature_names]
        X = adata.X[:, col_indices]
    else:
        X = adata.X
    totals = _row_sum(X)
    return pd.Series(
        totals >= min_total_counts,
        index=adata.obs_names,
        name="passes_empty_filter",
    )


# ---------------------------------------------------------------------------
# Count / batch flags
# ---------------------------------------------------------------------------


def _check_per_well_spec(spec: QuantileSpec, *, spec_name: str, wells: pd.Index) -> None:
    """Reject per-batch forms of a quantile spec and warn on keys that match no well.

    Summed counts have no batch dimension, so only a scalar or a flat
    ``{well: quantile, "__default__": quantile}`` mapping is meaningful.
    """
    if not isinstance(spec, Mapping):
        return
    unknown = []
    for key, val in spec.items():
        key_str = str(key)
        if key_str == "__batch__" or isinstance(val, Mapping):
            raise ValueError(
                f"{spec_name}: per-batch thresholds are not supported for ThreePrimeUntargeted, "
                "whose batches are summed into one count matrix. Pass a number or a "
                "{well: quantile, '__default__': quantile} mapping."
            )
        if key_str not in {"__default__", "default"} and key_str not in wells:
            unknown.append(key_str)
    if unknown:
        warnings.warn(
            f"{spec_name}: keys {unknown} match no well in the data and are ignored.",
            stacklevel=3,
        )


def _flag_total_counts_summed(
    adata: ad.AnnData,
    well_colname: str = "Well",
    min_quant_count: QuantileSpec = 0.025,
    max_quant_count: QuantileSpec = None,
) -> pd.Series:
    """Flag cells based on their *summed* total counts (across all batches), per well.

    This is the correct filter for ThreePrimeUntargeted data, where counts are
    measured in multiple batches that are summed into a single expression matrix.
    Filtering per-batch would incorrectly remove cells that legitimately express
    zero genes from some batches.

    Quantiles are computed per well. Each bound is a number, ``None`` (disabled),
    or a per-well mapping ``{well: quantile, "__default__": quantile}``.
    """
    totals = pd.Series(_row_sum(adata.X), index=adata.obs_names, name="total_counts")
    well_series = _obs_col_to_series(adata, well_colname)
    wells = pd.Index(well_series.astype(str)).unique()
    modality_defaults = DEFAULT_BOUNDS_PER_MODALITY["threeprimeuntargeted"]

    levels = {}
    for name, spec in (("min_quant_count", min_quant_count), ("max_quant_count", max_quant_count)):
        _check_per_well_spec(spec, spec_name=name, wells=wells)
        # Wells a mapping doesn't name (and no "__default__") get the modality default.
        fallback = modality_defaults[name] if isinstance(spec, Mapping) else spec
        levels[name] = _resolve_well_batch_quantiles(
            spec,
            count_cols=["total_counts"],
            well_labels=well_series,
            default_value=None if fallback is None else float(fallback),
            spec_name=name,
        )

    flag = pd.Series(False, index=adata.obs_names, name="passes_counts")
    for well, group in totals.groupby(well_series, observed=True):
        min_q = levels["min_quant_count"][str(well)]["total_counts"]
        max_q = levels["max_quant_count"][str(well)]["total_counts"]
        lower = -np.inf if min_q is None else float(group.quantile(min_q))
        upper = np.inf if max_q is None else float(group.quantile(max_q))
        flag.loc[group.index] = (group > lower) & (group <= upper)
    return flag


def _flag_total_counts_per_batch(
    adata: ad.AnnData,
    batch_colname: str = "batch_name",
    well_colname: str = "Well",
    min_quant_count: QuantileSpec = 0.025,
    max_quant_count: QuantileSpec = None,
    suppress_extra_cols_warning: bool = False,
) -> pd.Series:
    """Flag cells whose per-batch total counts fall within the quantile bounds.

    Reads batch counts from ``adata.obsm['total_counts']``, falling back to
    ``adata.obs['total_counts_*']`` columns in older stores.
    """
    # --- Resolve count columns ---
    if OBSM_TOTAL_COUNTS in adata.obsm:
        tc = _obsm_to_df(adata.obsm[OBSM_TOTAL_COUNTS], adata.obs_names)
        # Exclude nuclear columns — nuclear counts are a subset of the per-batch
        # signal already captured in non-nuclear columns; including them double-counts.
        non_nuclear = [c for c in tc.columns if "nuclear" not in str(c).lower()]
        tc = tc[non_nuclear]
        # Rename columns to total_counts_<batch> to match quantile resolver expectations
        count_cols = [f"total_counts_{c}" for c in tc.columns]
        counts_per_batch = tc.copy()
        counts_per_batch.columns = count_cols
    else:
        # Fallback: older stores keep per-batch totals in adata.obs
        batch_ids = (
            pd.Index(adata.var[batch_colname].dropna().unique())
            if batch_colname in adata.var.columns
            else pd.Index([])
        )
        count_cols = [f"total_counts_{batch}" for batch in batch_ids if f"total_counts_{batch}" in adata.obs]
        extra_cols = sorted(c for c in adata.obs if c.startswith("total_counts_") and c not in count_cols)
        if extra_cols and not suppress_extra_cols_warning:
            import warnings

            warnings.warn(
                f"Ignoring total-count columns with no matching batch metadata: {extra_cols}", stacklevel=2
            )
        counts_per_batch = adata.obs[count_cols]

    if not count_cols:
        return pd.Series(True, index=adata.obs_names, name="passes_counts")

    if isinstance(min_quant_count, (int, float)):
        default_min = float(min_quant_count)
    elif min_quant_count is None:
        default_min = None
    else:
        default_min = 0.025

    if isinstance(max_quant_count, (int, float)):
        default_max = float(max_quant_count)
    elif max_quant_count is None:
        default_max = None
    else:
        default_max = None

    well_series = _obs_col_to_series(adata, well_colname)
    min_levels = _resolve_well_batch_quantiles(
        min_quant_count,
        count_cols=count_cols,
        well_labels=well_series,
        default_value=default_min,
        spec_name="min_quant_count",
    )
    max_levels = _resolve_well_batch_quantiles(
        max_quant_count,
        count_cols=count_cols,
        well_labels=well_series,
        default_value=default_max,
        spec_name="max_quant_count",
    )

    def _passes_quantiles(group: pd.DataFrame) -> pd.Series:
        well = str(group.name)
        min_targets = min_levels[well]
        max_targets = max_levels[well]
        lower = pd.Series(
            {
                col: (-np.inf if min_targets[col] is None else group[col].quantile(min_targets[col]))
                for col in count_cols
            },
            index=count_cols,
        )
        upper = pd.Series(
            {
                col: (
                    np.inf
                    if max_targets[col] is None or max_targets[col] >= 1.0
                    else group[col].quantile(max_targets[col])
                )
                for col in count_cols
            },
            index=count_cols,
        )
        return ((group > lower) & (group <= upper)).all(axis=1)

    flag = _groupby_apply_to_cell_series(
        counts_per_batch.groupby(well_series, observed=False).apply(_passes_quantiles),
    )
    flag.name = "passes_counts"
    return flag


# ---------------------------------------------------------------------------
# Assignment rate flag
# ---------------------------------------------------------------------------


def _flag_assignment_rate(
    adata: ad.AnnData,
    *,
    min_assignment_rate: float | None = 0.5,
) -> pd.Series:
    """Flag cells whose assignment rate meets or exceeds ``min_assignment_rate``.

    Reads from ``adata.obsm['total_counts']`` and ``adata.obsm['unassigned_counts']``.
    """
    if min_assignment_rate is None:
        return pd.Series(True, index=adata.obs_names, name="passes_assignment_rate")

    if not 0.0 <= min_assignment_rate <= 1.0:
        raise ValueError(f"`min_assignment_rate` must lie within [0, 1]; received {min_assignment_rate}.")

    if OBSM_TOTAL_COUNTS not in adata.obsm or OBSM_UNASSIGNED_COUNTS not in adata.obsm:
        return pd.Series(True, index=adata.obs_names, name="passes_assignment_rate")

    assigned = _obsm_to_df(adata.obsm[OBSM_TOTAL_COUNTS], adata.obs_names)
    unassigned = _obsm_to_df(adata.obsm[OBSM_UNASSIGNED_COUNTS], adata.obs_names)

    # Try direct column intersection first.
    common_cols = assigned.columns.intersection(unassigned.columns)
    if common_cols.empty:
        assigned_batch_cols = [c for c in assigned.columns if "nuclear" not in c.lower()]

        # Pattern A: unassigned col = batch + "_Unassigned"
        unassigned_map_a: dict[str, str] = {
            col[: -len("_Unassigned")]: col
            for col in unassigned.columns
            if col.endswith("_Unassigned") and "nuclear" not in col.lower()
        }
        common_batches_a = [b for b in assigned_batch_cols if b in unassigned_map_a]

        # Pattern B: single bare "Unassigned" column (no batch prefix, e.g. Protein)
        bare_unassigned = [c for c in unassigned.columns if c.lower() == "unassigned"]

        if common_batches_a:
            assigned = assigned[common_batches_a].copy()
            unassigned = unassigned[[unassigned_map_a[b] for b in common_batches_a]].copy()
            unassigned.columns = common_batches_a
        elif bare_unassigned and assigned_batch_cols:
            # Single unassigned column — compute rate globally across all batches
            assigned = assigned[assigned_batch_cols].copy()
            total_assigned = assigned.sum(axis=1)
            total_unassigned = unassigned[bare_unassigned[0]]
            denom = total_assigned.add(total_unassigned, fill_value=0)
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = total_unassigned.divide(denom)
            mask = denom == 0
            rate = 1.0 - ratio.where(~mask, 1.0)
            flag = (rate >= min_assignment_rate).rename("passes_assignment_rate")
            return flag
        else:
            return pd.Series(True, index=adata.obs_names, name="passes_assignment_rate")
    else:
        assigned = assigned[common_cols]
        unassigned = unassigned[common_cols]

    denom = assigned.add(unassigned, fill_value=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = unassigned.divide(denom)
    mask = denom == 0
    unassigned_rate = ratio.where(~mask, 1.0)
    assignment_rate = 1.0 - unassigned_rate

    flag = assignment_rate.min(axis=1) >= min_assignment_rate
    flag.name = "passes_assignment_rate"
    return flag


# ---------------------------------------------------------------------------
# Cell area flag
# ---------------------------------------------------------------------------


def _flag_cell_area(
    adata: ad.AnnData,
    well_colname: str = "Well",
    area_colname: str = "Area",
    min_quant_area: QuantileSpec = 0.01,
    max_quant_area: QuantileSpec = 0.99,
) -> pd.Series:
    """Flag cells whose area falls within the per-well quantile bounds.

    Reads area from ``adata.obsm['cell_geometry']['Area']`` when available,
    falling back to ``adata.obs[area_colname]``.
    """
    # Resolve area values
    if OBSM_CELL_GEOMETRY in adata.obsm:
        cg = _obsm_to_df(adata.obsm[OBSM_CELL_GEOMETRY], adata.obs_names)
        if area_colname in cg.columns:
            area_values = cg[area_colname]
        elif "Area" in cg.columns:
            area_values = cg["Area"]
        elif area_colname in adata.obs.columns:
            area_values = adata.obs[area_colname]
        else:
            return pd.Series(True, index=adata.obs_names, name="passes_area")
    elif area_colname in adata.obs.columns:
        area_values = adata.obs[area_colname]
    else:
        return pd.Series(True, index=adata.obs_names, name="passes_area")

    well_series = _obs_col_to_series(adata, well_colname)

    default_min = (
        float(min_quant_area)
        if isinstance(min_quant_area, (int, float))
        else (None if min_quant_area is None else 0.01)
    )
    default_max = (
        float(max_quant_area)
        if isinstance(max_quant_area, (int, float))
        else (None if max_quant_area is None else 0.99)
    )

    min_levels = _resolve_well_batch_quantiles(
        min_quant_area,
        count_cols=[area_colname],
        well_labels=well_series,
        default_value=default_min,
        spec_name="min_quant_area",
    )
    max_levels = _resolve_well_batch_quantiles(
        max_quant_area,
        count_cols=[area_colname],
        well_labels=well_series,
        default_value=default_max,
        spec_name="max_quant_area",
    )

    low_bounds = np.empty(len(area_values), dtype=np.float64)
    high_bounds = np.empty(len(area_values), dtype=np.float64)

    for well in pd.Index(well_series).unique():
        mask = (well_series == well).to_numpy()
        values_well = area_values[mask]
        min_q = min_levels[well][area_colname]
        max_q = max_levels[well][area_colname]
        if not mask.any():
            low_val = np.nan
            high_val = np.nan
        else:
            low_val = -np.inf if min_q is None else values_well.quantile(min_q)
            high_val = np.inf if max_q is None or max_q >= 1.0 else values_well.quantile(max_q)
        low_bounds[mask] = low_val
        high_bounds[mask] = high_val

    values_np = area_values.to_numpy(dtype=np.float64)
    flag_area = (values_np > low_bounds) & (values_np < high_bounds)
    return pd.Series(flag_area, index=adata.obs_names, name="passes_area")


# ---------------------------------------------------------------------------
# NSB / GAPDH outlier flags
# ---------------------------------------------------------------------------


def _flag_nsb_outliers(
    adata: ad.AnnData,
    batch_colname: str = "batch_name",
    well_colname: str = "Well",
    max_quantile: QuantileSpec = 0.99,
) -> pd.Series:
    """Flag cells with unusually high NSB counts within each well.

    Reads NSB columns from ``adata.obsm['control_variables']``.
    """
    if OBSM_CONTROL_VARIABLES not in adata.obsm:
        return pd.Series(True, index=adata.obs_names, name="passes_nsb_outlier")

    cv = _obsm_to_df(adata.obsm[OBSM_CONTROL_VARIABLES], adata.obs_names)

    nsb_cols = [col for col in cv.columns if "NSB" in col and "nuclear" not in col.lower()]
    if not nsb_cols:
        return pd.Series(True, index=adata.obs_names, name="passes_nsb_outlier")

    nsb_df = cv[nsb_cols].fillna(0.0)
    wells = _obs_col_to_series(adata, well_colname)

    default_max = float(max_quantile) if isinstance(max_quantile, (int, float)) else None
    quantile_levels = _resolve_well_batch_quantiles(
        max_quantile,
        count_cols=list(nsb_df.columns),
        well_labels=wells,
        default_value=default_max,
        spec_name="max_quant_nsb",
    )

    def _passes(group: pd.DataFrame) -> pd.Series:
        well = str(group.name)
        targets = quantile_levels[well]
        thresholds = pd.Series(
            {
                col: (
                    np.inf
                    if targets[col] is None or targets[col] >= 1.0
                    else group[col].quantile(targets[col])
                )
                for col in group.columns
            },
            index=group.columns,
        )
        return (group <= thresholds).all(axis=1)

    flag = _groupby_apply_to_cell_series(
        nsb_df.groupby(wells, observed=False).apply(_passes),
    ).reindex(adata.obs_names, fill_value=True)
    flag.name = "passes_nsb_outlier"
    return flag


def _flag_gapdh_outliers(
    adata: ad.AnnData,
    batch_colname: str = "batch_name",
    well_colname: str = "Well",
    max_quantile: QuantileSpec = 0.99,
) -> pd.Series:
    """Flag cells with unusually high GAPDH counts within each well.

    Reads GAPDH counts from ``adata.X`` for rows where
    ``adata.var['control_type'] == 'Positive Control 1'``.

    Returns all-True (no filtering) when:
    - ``max_quantile`` is ``None`` (opt-in disabled)
    - ``control_type`` column is absent from ``adata.var``
    - No GAPDH probes are found
    """
    if max_quantile is None:
        return pd.Series(True, index=adata.obs_names, name="passes_gapdh_outlier")

    if VAR_CONTROL_TYPE_COL not in adata.var.columns:
        return pd.Series(True, index=adata.obs_names, name="passes_gapdh_outlier")

    gapdh_mask = adata.var[VAR_CONTROL_TYPE_COL] == VAR_GAPDH_CONTROL_TYPE
    if not gapdh_mask.any():
        return pd.Series(True, index=adata.obs_names, name="passes_gapdh_outlier")

    X_gapdh = adata.X[:, gapdh_mask.to_numpy()]
    gapdh_probe_names = adata.var_names[gapdh_mask].tolist()

    gapdh_arr = X_gapdh.toarray() if sp.issparse(X_gapdh) else np.asarray(X_gapdh)

    gapdh_df = pd.DataFrame(gapdh_arr, index=adata.obs_names, columns=gapdh_probe_names)
    wells = _obs_col_to_series(adata, well_colname)

    default_max = float(max_quantile) if isinstance(max_quantile, (int, float)) else None
    quantile_levels = _resolve_well_batch_quantiles(
        max_quantile,
        count_cols=gapdh_probe_names,
        well_labels=wells,
        default_value=default_max,
        spec_name="max_quant_gapdh",
    )

    def _passes(group: pd.DataFrame) -> pd.Series:
        well = str(group.name)
        targets = quantile_levels[well]
        thresholds = pd.Series(
            {
                col: (
                    np.inf
                    if targets[col] is None or targets[col] >= 1.0
                    else group[col].quantile(targets[col])
                )
                for col in group.columns
            },
            index=group.columns,
        )
        return (group <= thresholds).all(axis=1)

    flag = _groupby_apply_to_cell_series(
        gapdh_df.groupby(wells, observed=False).apply(_passes),
    ).reindex(adata.obs_names, fill_value=True)
    flag.name = "passes_gapdh_outlier"
    return flag


# ---------------------------------------------------------------------------
# MT / nuclear fraction flags
# ---------------------------------------------------------------------------


def _flag_spatial_artifacts(
    adata: ad.AnnData,
    *,
    target_type: str,
    enabled: bool = False,
) -> pd.Series:
    """Flag cells that are not members of a spatial artifact region.

    Reads the modality rollup column from ``adata.obsm['spatial_artifacts']``:
    ``any-imaging``, ``any-ThreePrimeUntargeted``, or ``any-OPS``. ``True`` in
    that frame means the cell sits in a detected artifact; this flag inverts
    that so ``True`` means pass (keep).

    When ``enabled`` is False (the default), every cell passes. When the
    slot or rollup column is missing — stores built without spatial-artifact
    detection, older stores, or a modality with no detector (Transcript /
    Protein) — the gate is skipped with a warning and every cell passes.
    """
    name = "passes_spatial_artifacts"
    if not enabled:
        return pd.Series(True, index=adata.obs_names, name=name)

    rollup = SPATIAL_ARTIFACT_ROLLUP.get(target_type)
    if rollup is None:
        warnings.warn(
            f"Spatial-artifact filtering is not defined for target_type {target_type!r} "
            f"(no {list(SPATIAL_ARTIFACT_ROLLUP.values())} rollup); the gate was not applied.",
            stacklevel=2,
        )
        return pd.Series(True, index=adata.obs_names, name=name)

    if OBSM_SPATIAL_ARTIFACTS not in adata.obsm:
        warnings.warn(
            "obsm['spatial_artifacts'] is absent (store built without spatial-artifact "
            "detection, or an older store); the spatial-artifact gate was not applied.",
            stacklevel=2,
        )
        return pd.Series(True, index=adata.obs_names, name=name)

    artifacts = _obsm_to_df(adata.obsm[OBSM_SPATIAL_ARTIFACTS], adata.obs_names)
    if rollup not in artifacts.columns:
        warnings.warn(
            f"obsm['spatial_artifacts'] has no {rollup!r} column; the spatial-artifact gate was not applied.",
            stacklevel=2,
        )
        return pd.Series(True, index=adata.obs_names, name=name)

    flagged = artifacts[rollup].fillna(False).astype(bool)
    return (~flagged).rename(name)


def _flag_mt_outliers(
    adata: ad.AnnData,
    well_colname: str = "Well",
    max_pct_mt: float | None = None,
) -> pd.Series:
    """Flag cells with mitochondrial gene percentage above ``max_pct_mt``.

    Reads ``pct_counts_mt`` from ``adata.obsm['qc_metrics']``.
    """
    if max_pct_mt is None:
        return pd.Series(True, index=adata.obs_names, name="passes_mt_outlier")

    if not 0.0 <= max_pct_mt <= 1.0:
        raise ValueError(f"`max_pct_mt` must be between 0 and 1; received {max_pct_mt}.")

    if QC_METRICS_OBSM_KEY not in adata.obsm:
        return pd.Series(True, index=adata.obs_names, name="passes_mt_outlier")

    qc_metrics = adata.obsm[QC_METRICS_OBSM_KEY]
    if "pct_counts_mt" not in qc_metrics.columns:
        return pd.Series(True, index=adata.obs_names, name="passes_mt_outlier")

    mt_pct = qc_metrics["pct_counts_mt"].reindex(adata.obs_names)
    if mt_pct.isna().all():
        return pd.Series(True, index=adata.obs_names, name="passes_mt_outlier")

    # pct_counts_mt is stored as 0–100 (from scanpy); max_pct_mt is a 0–1 fraction.
    flag = (mt_pct <= max_pct_mt * 100).fillna(True)
    flag.name = "passes_mt_outlier"
    return flag


def _flag_nuclear_fraction(
    adata: ad.AnnData,
    *,
    min_nuclear_fraction: float | None = None,
    max_nuclear_fraction: float | None = None,
    nuclear_layer: str = "nuclear",
    x_layer: str | None = None,
) -> pd.Series:
    """Flag cells whose nuclear fraction lies within [min, max]."""
    for name, v in (
        ("min_nuclear_fraction", min_nuclear_fraction),
        ("max_nuclear_fraction", max_nuclear_fraction),
    ):
        if v is not None and not (0.0 <= v <= 1.0):
            raise ValueError(f"`{name}` must lie within [0, 1] or be None; received {v}.")

    if nuclear_layer not in adata.layers:
        raise KeyError(f"Expected nuclear counts in `adata.layers['{nuclear_layer}']`.")

    total_mat = adata.layers[x_layer] if x_layer is not None else adata.X
    nuc_mat = adata.layers[nuclear_layer]

    if total_mat.shape != nuc_mat.shape:
        raise ValueError(f"Shape mismatch: total {total_mat.shape} vs nuclear {nuc_mat.shape}.")

    total_sum = _row_sum(total_mat)
    nuclear_sum = _row_sum(nuc_mat)

    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.divide(
            nuclear_sum,
            total_sum,
            out=np.full_like(nuclear_sum, np.nan, dtype=float),
            where=total_sum > 0,
        )

    flag = pd.Series(True, index=adata.obs_names, name="passes_nuclear_fraction")
    flag &= ~pd.isna(frac)

    if min_nuclear_fraction is not None:
        flag &= frac >= float(min_nuclear_fraction)
    if max_nuclear_fraction is not None:
        flag &= frac <= float(max_nuclear_fraction)

    return flag


# ---------------------------------------------------------------------------
# Gene-level filtering
# ---------------------------------------------------------------------------


def _flag_genes_per_well(
    adata: ad.AnnData,
    *,
    well_colname: str = "Well",
    min_cells_per_gene: int | None = None,
    min_counts_per_gene: float | None = None,
) -> pd.Series:
    """Flag genes that meet expression thresholds in at least one well."""
    if min_cells_per_gene is None and min_counts_per_gene is None:
        return pd.Series(True, index=adata.var_names, name="passes_gene_filter")

    if min_cells_per_gene is not None and min_cells_per_gene < 0:
        raise ValueError(f"`min_cells_per_gene` must be >= 0; received {min_cells_per_gene}.")
    if min_counts_per_gene is not None and min_counts_per_gene < 0:
        raise ValueError(f"`min_counts_per_gene` must be >= 0; received {min_counts_per_gene}.")

    well_series = _obs_col_to_series(adata, well_colname)
    wells = well_series.unique()
    gene_passes = pd.Series(False, index=adata.var_names, name="passes_gene_filter")
    X = adata.X

    for well in wells:
        well_mask = (well_series == well).to_numpy()
        if not well_mask.any():
            continue
        X_well = X[well_mask, :]

        if min_cells_per_gene is not None:
            if sp.issparse(X_well):
                cells_expressing = np.asarray((X_well > 0).sum(axis=0)).ravel()
            else:
                cells_expressing = (X_well > 0).sum(axis=0)
            passes_cells = cells_expressing >= min_cells_per_gene
        else:
            passes_cells = np.ones(adata.n_vars, dtype=bool)

        if min_counts_per_gene is not None:
            counts_per_gene = _col_sum(X_well)
            passes_counts = counts_per_gene >= min_counts_per_gene
        else:
            passes_counts = np.ones(adata.n_vars, dtype=bool)

        gene_passes |= passes_cells & passes_counts

    return gene_passes


# ---------------------------------------------------------------------------
# Imaging-specific flags (CellPaint + IF Protein)
# ---------------------------------------------------------------------------


def _get_primary_intensity_features(adata: ad.AnnData) -> list[str]:
    """Identify the structural morphology channels for the imaging intensity gate.

    The gate answers "is this a real, well-segmented, well-stained cell?" and
    must run on the Cell Painting structural stains only, never on the protein /
    antibody panel, CellProfiler descriptors, or background channels.

    Delegates to
    :func:`elembio_spatialdata_analysis.tl.morphology.select_morphology_channels`,
    which anchors on the CellPaint batch (``batch_name`` starting with ``"CP"``)
    and adapts to the store format.  Returns ``[]`` (a fail-safe no-op) when no
    CellPaint batch is present; :func:`_flag_imaging_intensity` then passes all
    cells rather than gating on the protein batch.
    """
    from elembio_spatialdata_analysis.tl.morphology import select_morphology_channels

    names = select_morphology_channels(adata.var)
    if names:
        _logger.info(
            "Imaging intensity gate: gating on %d morphology channel(s): %s",
            len(names),
            names,
        )
    else:
        warnings.warn(
            "No CellPaint morphology channel found; skipping the imaging intensity gate "
            "(all cells pass). Imaging QC requires a CellPaint (CP) batch of structural stains.",
            stacklevel=2,
        )
    return names


def _flag_imaging_intensity(
    adata: ad.AnnData,
    *,
    well_colname: str = "Well",
    min_quant_count: QuantileSpec = 0.025,
    max_quant_count: QuantileSpec = None,
    primary_feature_names: Sequence[str] | None = None,
) -> pd.Series:
    """Flag imaging cells based on morphology-channel intensity quantiles per well.

    Runs only on the structural Cell Painting channels selected by
    :func:`_get_primary_intensity_features` (the CellPaint batch), never on the
    protein / antibody panel. Passes all cells when no morphology channel is
    found.

    Parameters
    ----------
    primary_feature_names
        Optional precomputed morphology channel names. When ``None``, channels
        are resolved via :func:`_get_primary_intensity_features`.
    """
    if adata.n_vars == 0:
        raise ValueError("No imaging features available for filtering.")

    well_series = _obs_col_to_series(adata, well_colname)

    default_min = (
        float(min_quant_count)
        if isinstance(min_quant_count, (int, float))
        else (None if min_quant_count is None else 0.025)
    )
    default_max = float(max_quant_count) if isinstance(max_quant_count, (int, float)) else None

    if primary_feature_names is None:
        primary_feature_names = _get_primary_intensity_features(adata)

    if not primary_feature_names:
        return pd.Series(True, index=adata.obs_names, name="passes_imaging_intensity")

    col_indices = [adata.var_names.get_loc(name) for name in primary_feature_names]
    X_primary = adata.X[:, col_indices]
    X_primary = X_primary.toarray() if sp.issparse(X_primary) else np.asarray(X_primary)
    primary_expression = pd.DataFrame(X_primary, index=adata.obs_names, columns=primary_feature_names)
    count_cols = primary_feature_names

    min_levels = _resolve_well_batch_quantiles(
        min_quant_count,
        count_cols=count_cols,
        well_labels=well_series,
        default_value=default_min,
        spec_name="min_quant_count",
    )
    max_levels = _resolve_well_batch_quantiles(
        max_quant_count,
        count_cols=count_cols,
        well_labels=well_series,
        default_value=default_max,
        spec_name="max_quant_count",
    )

    def _passes_primary_quantiles(group: pd.DataFrame) -> pd.Series:
        well = str(group.name)
        min_targets = min_levels[well]
        max_targets = max_levels[well]
        out = pd.Series(True, index=group.index)
        for col in count_cols:
            low = -np.inf if min_targets[col] is None else group[col].quantile(min_targets[col])
            high = (
                np.inf
                if max_targets[col] is None or max_targets[col] >= 1.0
                else group[col].quantile(max_targets[col])
            )
            out &= (group[col] > low) & (group[col] <= high)
        return out

    flag = _groupby_apply_to_cell_series(
        primary_expression.groupby(well_series, observed=False).apply(_passes_primary_quantiles),
    )
    flag.name = "passes_imaging_intensity"
    return flag.reindex(adata.obs_names, fill_value=False)


# Priority order for first-failing-reason attribution: each removed cell
# is attributed to the first ``passes_*`` flag it failed in this order;
# later flags it also failed are not double-counted.
FIRST_FAILING_ORDER: tuple[tuple[str, str], ...] = (
    ("passes_nan_filter", "NaN values"),
    ("passes_data_missing", "Missing data"),
    ("passes_empty_filter", "Empty cells"),
    ("passes_counts", "Total counts out of range"),
    ("passes_counts_per_batch", "Total counts out of range"),  # legacy compat (0.11.0)
    ("passes_imaging_intensity", "Channel intensity out of range"),
    ("passes_imaging_counts", "Channel intensity out of range"),  # legacy compat (0.10.0)
    ("passes_cell_paint_counts", "Channel intensity out of range"),  # legacy compat
    ("passes_assignment_rate", "Assignment rate below threshold"),
    ("passes_area", "Cell area out of range"),
    ("passes_nsb_outlier", "High NSB counts"),
    ("passes_gapdh_outlier", "High GAPDH counts"),
    ("passes_mt_outlier", "High mitochondrial %"),
    ("passes_nuclear_fraction", "Nuclear fraction out of range"),
    ("passes_spatial_artifacts", "In a spatial artifact region"),
)


def first_failing_reason(
    flags: pd.DataFrame,
    *,
    order: Sequence[tuple[str, str]] = FIRST_FAILING_ORDER,
) -> dict[str, int]:
    """Per-reason count of cells removed by the filter pipeline.

    Walks ``order`` in priority order; for each ``(column, label)`` pair
    present in ``flags``, counts cells where ``flags[column] == False``
    and which have not already been attributed to an earlier flag.
    Skips columns missing from ``flags`` and skips reasons with zero
    cells in the result, so the returned dict contains only non-empty
    rows in ``order`` priority.
    """
    if flags is None or flags.empty:
        return {}
    attributed = np.zeros(len(flags), dtype=bool)
    counts: dict[str, int] = {}
    for col, label in order:
        if col not in flags.columns:
            continue
        failed = (~flags[col].to_numpy(dtype=bool)) & (~attributed)
        n = int(failed.sum())
        if n:
            counts[label] = n
            attributed |= failed
    return counts
