"""QC metric computation functions for the filter module.

Data access conventions:
- Total assigned counts per batch  → adata.obsm['total_counts']   (DataFrame, cols = batch names)
- Unassigned counts per batch      → adata.obsm['unassigned_counts'] (DataFrame, cols = batch names)
- NSB / control variables          → adata.obsm['control_variables'] (DataFrame, cols = probe names)
- GAPDH positive controls          → adata.X  (rows where adata.var['control_type'] == 'Positive Control 1')
- Nuclear counts                   → adata.layers['nuclear']
- Cell geometry (area, etc.)       → adata.obsm['cell_geometry']
"""

from __future__ import annotations

import contextlib
import warnings
from collections.abc import Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from elembio_spatialdata_analysis.pp.filter._constants import (
    OBSM_CONTROL_VARIABLES,
    OBSM_TOTAL_COUNTS,
    OBSM_UNASSIGNED_COUNTS,
    QC_METRICS_OBSM_KEY,
    VAR_CONTROL_TYPE_COL,
    VAR_GAPDH_CONTROL_TYPE,
)
from elembio_spatialdata_analysis.tl.constants import normalize_target_type

# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------


def _row_sum(mat) -> np.ndarray:
    if sp.issparse(mat):
        return np.asarray(mat.sum(axis=1)).ravel().astype(float)
    return np.asarray(mat.sum(axis=1)).ravel().astype(float)


def _col_sum(mat) -> np.ndarray:
    if sp.issparse(mat):
        return np.asarray(mat.sum(axis=0)).ravel().astype(float)
    return np.asarray(mat.sum(axis=0)).ravel().astype(float)


def _obsm_to_df(value, obs_names: pd.Index) -> pd.DataFrame:
    """Normalise any obsm value to a ``pd.DataFrame`` indexed by ``obs_names``.

    Handles the three types that appear in practice:

    * ``pd.DataFrame`` — reindexed to ``obs_names`` and returned as-is.
    * ``anndata.Dataset2D`` (xarray-backed obsm from zarr) — materialised via
      ``.to_memory()``, which triggers the zarr read and returns a properly
      indexed ``pd.DataFrame``.  The index from zarr is then replaced with
      ``obs_names`` for consistent downstream alignment.
    * ``np.ndarray`` (1-D or 2-D) — reshaped to 2-D if needed, wrapped in
      a DataFrame indexed by ``obs_names``.

    ``obs_names`` must always match the row count of ``value``; callers
    should pass ``adata.obs_names`` (not the parent's obs_names when working
    on a view/subset).
    """
    if isinstance(value, pd.DataFrame):
        return value.reindex(obs_names)

    # Dataset2D: use the official .to_memory() API to trigger zarr reads and
    # get a properly-formed DataFrame.  Do NOT pass obs_names to pd.DataFrame()
    # here — the Dataset2D shape is (n_obs_in_view, n_cols) and injecting a
    # foreign-length index causes a shape mismatch ValueError.
    if hasattr(value, "to_memory"):
        df = value.to_memory()
        df.index = obs_names
        return df

    # Plain numpy array fallback
    arr = np.asarray(value)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    return pd.DataFrame(arr, index=obs_names)


def _obs_col_to_series(adata: ad.AnnData, col: str) -> pd.Series:
    """Materialize an obs column to a pandas Series.

    With ``read_lazy``, ``adata.obs`` is xarray-backed; accessing a column
    returns a DataArray that lacks ``.unique()`` and triggers per-element
    loads when iterated. This helper forces a single materialization.
    """
    obs_col = adata.obs[col]
    if isinstance(obs_col, pd.Series):
        return obs_col.astype(str)
    if hasattr(obs_col, "to_numpy"):
        arr = obs_col.to_numpy()
    elif hasattr(obs_col, "values"):
        arr = obs_col.values
    else:
        arr = np.asarray(obs_col, dtype=object)
    out = pd.Series(arr, index=adata.obs_names)
    return out.astype(str)


# ---------------------------------------------------------------------------
# Per-metric computation
# ---------------------------------------------------------------------------


def _compute_total_counts(
    adata: ad.AnnData,
    feature_names: Sequence[str] | None = None,
) -> pd.Series:
    """Sum assigned counts across all non-nuclear batches from ``adata.obsm['total_counts']``.

    Nuclear counts are a subset of the cytoplasmic+nuclear signal already captured
    in the non-nuclear batch columns, so nuclear columns (containing "nuclear" in
    their name) are excluded to avoid double-counting.

    Falls back to summing ``adata.X`` directly when the obsm key is absent (e.g.
    for imaging modalities that have no batch breakdown).

    Parameters
    ----------
    feature_names
        Optional ``adata.X`` column names to sum when falling back (imaging:
        primary morphology only). Ignored when ``obsm['total_counts']`` is
        present. An empty sequence yields all-NaN totals (no morphology
        channels to sum).
    """
    if feature_names is None and OBSM_TOTAL_COUNTS in adata.obsm:
        df = _obsm_to_df(adata.obsm[OBSM_TOTAL_COUNTS], adata.obs_names)
        non_nuclear_cols = [c for c in df.columns if "nuclear" not in str(c).lower()]
        totals = df[non_nuclear_cols].sum(axis=1).rename("total_counts")
        return totals
    if feature_names is not None:
        if len(feature_names) == 0:
            return pd.Series(np.nan, index=adata.obs_names, name="total_counts")
        missing = [n for n in feature_names if n not in adata.var_names]
        if missing:
            raise KeyError(f"feature_names not in adata.var_names: {missing}")
        col_indices = [adata.var_names.get_loc(name) for name in feature_names]
        totals = _row_sum(adata.X[:, col_indices])
        return pd.Series(totals, index=adata.obs_names, name="total_counts")
    # Fallback: sum X directly
    totals = _row_sum(adata.X)
    return pd.Series(totals, index=adata.obs_names, name="total_counts")


def _compute_n_genes_by_counts(adata: ad.AnnData) -> pd.Series:
    """Number of genes with non-zero counts per cell."""
    X = adata.X
    n_genes = np.asarray((X > 0).sum(axis=1)).ravel() if sp.issparse(X) else (X > 0).sum(axis=1)
    return pd.Series(n_genes.astype(np.float64), index=adata.obs_names, name="n_genes_by_counts")


def _compute_assignment_rate(adata: ad.AnnData) -> pd.Series:
    """Compute per-cell assignment rate from ``obsm['total_counts']`` and ``obsm['unassigned_counts']``.

    assignment_rate = total_assigned / (total_assigned + total_unassigned)

    The two DataFrames use different column naming conventions in the store:
    - ``total_counts`` columns: batch names like ``B02``, ``B03``, …
    - ``unassigned_counts`` columns: ``B02_Unassigned``, ``B03_Unassigned``, …

    We strip the ``_Unassigned`` suffix (and skip ``_Unassigned_nuclear`` columns)
    to align the two DataFrames, then compute per-batch assignment rates and
    return the *minimum* across all batches per cell.
    """
    if OBSM_TOTAL_COUNTS not in adata.obsm:
        return pd.Series(np.nan, index=adata.obs_names, name="total_assignment_rate")
    if OBSM_UNASSIGNED_COUNTS not in adata.obsm:
        return pd.Series(np.nan, index=adata.obs_names, name="total_assignment_rate")

    assigned = _obsm_to_df(adata.obsm[OBSM_TOTAL_COUNTS], adata.obs_names)
    unassigned = _obsm_to_df(adata.obsm[OBSM_UNASSIGNED_COUNTS], adata.obs_names)

    # Try direct column intersection first (works if naming already matches).
    common_cols = assigned.columns.intersection(unassigned.columns)
    if common_cols.empty:
        # Naming mismatch — two observed patterns:
        #   Pattern A (Transcript): total_counts=[B02,B03,...], unassigned=[B02_Unassigned,...]
        #   Pattern B (Protein):    total_counts=[B01,...],     unassigned=[Unassigned,...]
        # Non-nuclear assigned columns (source of truth for batch names):
        assigned_batch_cols = [c for c in assigned.columns if "nuclear" not in c.lower()]

        # Pattern A: unassigned col = batch + "_Unassigned"
        unassigned_map_a: dict[str, str] = {
            col[: -len("_Unassigned")]: col
            for col in unassigned.columns
            if col.endswith("_Unassigned") and "nuclear" not in col.lower()
        }
        common_batches_a = [b for b in assigned_batch_cols if b in unassigned_map_a]

        # Pattern B: single bare "Unassigned" column (no batch prefix)
        bare_unassigned = [c for c in unassigned.columns if c.lower() == "unassigned"]

        if common_batches_a:
            # Pattern A match
            assigned = assigned[common_batches_a].copy()
            unassigned = unassigned[[unassigned_map_a[b] for b in common_batches_a]].copy()
            unassigned.columns = common_batches_a
        elif bare_unassigned and assigned_batch_cols:
            # Pattern B: broadcast single unassigned column across all assigned batch columns
            assigned = assigned[assigned_batch_cols].copy()
            total_assigned = assigned.sum(axis=1)
            total_unassigned = unassigned[bare_unassigned[0]]
            denom = total_assigned.add(total_unassigned, fill_value=0)
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = total_unassigned.divide(denom)
            mask = denom == 0
            result = (1.0 - ratio.where(~mask, 1.0)).rename("total_assignment_rate")
            return result
        else:
            return pd.Series(np.nan, index=adata.obs_names, name="total_assignment_rate")
    else:
        assigned = assigned[common_cols]
        unassigned = unassigned[common_cols]

    denom = assigned.add(unassigned, fill_value=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = unassigned.divide(denom)
    mask = denom == 0
    unassigned_rate = ratio.where(~mask, 1.0)
    assignment_rate = 1.0 - unassigned_rate
    result = assignment_rate.min(axis=1)
    result.name = "total_assignment_rate"
    return result


def _compute_total_nsb_counts(adata: ad.AnnData, batch_colname: str) -> pd.Series:
    """Sum NSB probe counts from ``adata.obsm['control_variables']``.

    Columns are matched by the substring ``NSB`` appearing anywhere in the name
    (e.g. ``B02_NSB_1``, ``NSB_1``) while excluding nuclear variants (containing
    ``Nuclear`` or ``nuclear``).  Returns zeros when the obsm key is absent or
    no NSB columns are found.
    """
    if OBSM_CONTROL_VARIABLES not in adata.obsm:
        return pd.Series(0.0, index=adata.obs_names, name="total_nsb_counts")

    cv = _obsm_to_df(adata.obsm[OBSM_CONTROL_VARIABLES], adata.obs_names)

    nsb_cols = [col for col in cv.columns if "NSB" in col and "nuclear" not in col.lower()]
    if not nsb_cols:
        return pd.Series(0.0, index=adata.obs_names, name="total_nsb_counts")

    totals = cv[nsb_cols].sum(axis=1)
    totals.name = "total_nsb_counts"
    return totals


def _compute_total_gapdh_counts(adata: ad.AnnData, batch_colname: str) -> pd.Series:
    """Sum GAPDH positive-control counts from ``adata.X``.

    Rows in ``adata.var`` where ``control_type == 'Positive Control 1'`` are
    treated as GAPDH probes.  Returns zeros when no such rows exist.
    """
    if VAR_CONTROL_TYPE_COL not in adata.var.columns:
        return pd.Series(0.0, index=adata.obs_names, name="total_gapdh_counts")

    gapdh_mask = adata.var[VAR_CONTROL_TYPE_COL] == VAR_GAPDH_CONTROL_TYPE
    if not gapdh_mask.any():
        return pd.Series(0.0, index=adata.obs_names, name="total_gapdh_counts")

    X_gapdh = adata.X[:, gapdh_mask.to_numpy()]
    totals = _row_sum(X_gapdh)
    return pd.Series(totals, index=adata.obs_names, name="total_gapdh_counts")


def _compute_nuclear_fraction(adata: ad.AnnData) -> pd.Series:
    """Compute per-cell nuclear fraction from ``adata.layers['nuclear']``."""
    nuclear_layer = "nuclear"
    if nuclear_layer not in adata.layers:
        raise KeyError(f"Expected nuclear counts in `adata.layers['{nuclear_layer}']`.")

    total_mat = adata.X
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
    return pd.Series(frac, index=adata.obs_names, name="nuclear_frac")


def _extract_rnaseq_qc_metrics_from_obs(adata: ad.AnnData) -> dict[str, pd.Series]:
    """Extract RNA-seq QC metrics from ``adata.obs`` or compute via scanpy."""
    metrics: dict[str, pd.Series] = {}

    qc_metric_names = [
        "n_genes_by_counts",
        "pct_counts_in_top_10_genes",
        "pct_counts_in_top_20_genes",
        "total_counts_mt",
        "pct_counts_mt",
        "total_counts_chrm",
        "pct_counts_chrm",
        "total_counts_ribo",
        "pct_counts_ribo",
        "total_counts_hb",
        "pct_counts_hb",
    ]

    for metric_name in qc_metric_names:
        if metric_name in adata.obs.columns:
            metrics[metric_name] = adata.obs[metric_name]

    missing_metrics = [m for m in qc_metric_names if m not in metrics]

    if missing_metrics and adata.var_names.size > 0:
        try:
            import scanpy as sc

            if "mt" not in adata.var.columns:
                # Literature / scanpy idiom: MT-* gene symbols only.
                adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
            if "chrm" not in adata.var.columns:
                # All assigned chrM in the 3' matrix: MT-* genes plus the
                # synthetic catch-all row (matrix GeneID stays chrM_umbrella).
                adata.var["chrm"] = adata.var["mt"] | (adata.var_names.astype(str) == "chrM_umbrella")
            if "ribo" not in adata.var.columns:
                adata.var["ribo"] = adata.var_names.str.match(r"^RP[SL]\d+", case=False)
            if "hb" not in adata.var.columns:
                # Matches the scanpy-tutorial / single-cell-best-practices idiom:
                # HB-prefixed genes excluding HBP*. Note this admits two known
                # non-hemoglobin false positives (HBEGF, HBS1L) for parity with
                # the convention; an explicit canonical list is the alternative.
                adata.var["hb"] = adata.var_names.str.contains(r"^HB[^(P)]", case=False)

            # scanpy requires float32/float64; cast integer sparse matrices first
            # to avoid "Output dtype not compatible with inputs" from scipy.
            X = adata.X
            if sp.issparse(X) and not np.issubdtype(X.dtype, np.floating):
                adata.X = X.astype(np.float32)
            elif not sp.issparse(X) and not np.issubdtype(np.asarray(X).dtype, np.floating):
                adata.X = np.asarray(X, dtype=np.float32)

            sc.pp.calculate_qc_metrics(
                adata,
                percent_top=None,
                log1p=False,
                inplace=True,
                qc_vars=["mt", "chrm", "ribo", "hb"],
            )

            for metric_name in qc_metric_names:
                if metric_name in adata.obs.columns:
                    metrics[metric_name] = adata.obs[metric_name]
        except ImportError:
            if not metrics:
                warnings.warn(
                    "scanpy is not available. QC metrics were not found in adata.obs and "
                    "cannot be computed automatically. Install scanpy and call "
                    "sc.pp.calculate_qc_metrics() before filtering.",
                    UserWarning,
                    stacklevel=3,
                )
        except Exception as e:
            if not metrics:
                warnings.warn(
                    f"Failed to compute QC metrics automatically: {e}.",
                    UserWarning,
                    stacklevel=3,
                )

    return metrics


def compute_qc_metrics(
    adata: ad.AnnData,
    target_type: str,
    *,
    batch_colname: str = "batch_name",
) -> None:
    """Compute QC metrics and store them in ``adata.obsm['qc_metrics']``.

    This is an explicit, always-inplace step that should be called before
    plotting or filtering — analogous to ``sc.pp.calculate_qc_metrics``.
    Calling this before ``filter_cells`` allows pre-filter distribution plots
    to work correctly; ``filter_cells`` will skip recomputing if metrics are
    already present.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object. Modified in place.
    target_type : str
        Modality identifier (e.g. ``"Transcript"``, ``"Protein"``, ``"Imaging"``).
        Controls which metrics are computed.
    batch_colname : str
        Column in ``adata.obs`` identifying feature batches. Default ``"batch_name"``.

    Writes
    ------
    ``adata.obsm['qc_metrics']`` : DataFrame
        Per-cell QC metrics. Columns vary by modality:

        - All modalities: ``total_counts`` (Imaging: sum of primary morphology
          channels only, matching the NaN/empty/intensity gates)
        - Transcript: + ``n_genes_by_counts``, ``total_assignment_rate``,
          ``total_nsb_counts``, ``total_gapdh_counts``, ``nuclear_frac``
        - ThreePrimeUntargeted: + ``n_genes_by_counts`` (always from ``.X``), RNA-seq
          metrics from Scanpy when available (``pct_counts_mt`` = MT-* genes;
          ``pct_counts_chrm`` = MT-* ∪ ``chrM_umbrella``; ``pct_counts_ribo`` /
          ``pct_counts_hb``), ``total_assignment_rate``, ``nuclear_frac``
        - Protein: + ``total_assignment_rate``, ``nuclear_frac``
    """
    _add_obsm_qc_matrix(adata, normalize_target_type(target_type), batch_colname)


def _add_obsm_qc_matrix(adata: ad.AnnData, target_type: str, batch_colname: str) -> None:
    """Compute QC metrics once and persist them under ``adata.obsm['qc_metrics']``."""
    normalized_type = normalize_target_type(target_type)

    metrics: dict[str, pd.Series] = {}

    if normalized_type == "imaging":
        # Align with imaging NaN/empty/intensity gates: sum primary morphology
        # only so Bkg_*/IF NaNs do not poison qc_metrics['total_counts'].
        from elembio_spatialdata_analysis.tl.morphology import select_morphology_channels

        primary = select_morphology_channels(adata.var)
        metrics["total_counts"] = _compute_total_counts(adata, feature_names=primary)
    else:
        metrics["total_counts"] = _compute_total_counts(adata)

    if normalized_type == "transcript":
        metrics["n_genes_by_counts"] = _compute_n_genes_by_counts(adata)
        with contextlib.suppress(Exception):
            metrics["total_assignment_rate"] = _compute_assignment_rate(adata)
        metrics["total_nsb_counts"] = _compute_total_nsb_counts(adata, batch_colname)
        metrics["total_gapdh_counts"] = _compute_total_gapdh_counts(adata, batch_colname)
        with contextlib.suppress(KeyError):
            metrics["nuclear_frac"] = _compute_nuclear_fraction(adata)

    elif normalized_type == "threeprimeuntargeted":
        rnaseq_metrics = _extract_rnaseq_qc_metrics_from_obs(adata)
        metrics.update(rnaseq_metrics)
        # Match Transcript: always compute genes detected from .X. For 3′ we otherwise
        # depended on Scanpy inside _extract_rnaseq_qc_metrics_from_obs, which can fail
        # (e.g. dtype compatibility on large sparse matrices) and then omit this column.
        if "n_genes_by_counts" not in metrics:
            metrics["n_genes_by_counts"] = _compute_n_genes_by_counts(adata)
        with contextlib.suppress(Exception):
            metrics["total_assignment_rate"] = _compute_assignment_rate(adata)
        with contextlib.suppress(KeyError):
            metrics["nuclear_frac"] = _compute_nuclear_fraction(adata)

    elif normalized_type == "protein":
        with contextlib.suppress(Exception):
            metrics["total_assignment_rate"] = _compute_assignment_rate(adata)
        with contextlib.suppress(KeyError):
            metrics["nuclear_frac"] = _compute_nuclear_fraction(adata)

    elif normalized_type == "ops":
        # OPS adata only carries ``total_counts`` (already populated above)
        # and ``cell_geometry``.  OPS has no ``extracellular_counts`` or
        # ``unassigned_counts`` by design — the counting model treats every
        # assigned read as belonging to exactly one target, so there is
        # no "unassigned" bucket.  ``_compute_assignment_rate`` will
        # short-circuit to NaN when ``unassigned_counts`` is absent.
        # Nuclear fraction is also not modelled for OPS.
        with contextlib.suppress(Exception):
            metrics["total_assignment_rate"] = _compute_assignment_rate(adata)

    qc_df = pd.DataFrame(metrics).reindex(adata.obs_names)
    adata.obsm[QC_METRICS_OBSM_KEY] = qc_df

    if normalized_type == "transcript":
        adata.obs["total_counts"] = qc_df["total_counts"].values
        if "n_genes_by_counts" in qc_df.columns:
            adata.obs["n_genes_by_counts"] = qc_df["n_genes_by_counts"].values
