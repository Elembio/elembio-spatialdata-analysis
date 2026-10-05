"""Private normalization primitives.

All functions are prefixed with ``_`` and are not part of the public API.
Public-facing wrappers live in ``modalities.py`` and ``normalize.py``.

Four building blocks:

* ``_normalize_agg_counts``  — library-size scaling per batch
* ``_normalize_median_of_ratios`` — well-level median-of-ratios correction
* ``_log1p``  — element-wise log(1 + x), sparse-safe
* ``_clr``    — Centered Log-Ratio transformation, sparse-safe
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from elembio_spatialdata_analysis.pp.normalize._resolve import _resolve_well_groups

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _var_col_to_series(adata, col: str) -> pd.Series:
    """Materialize a var column to a pandas Series.

    With lazy loading, ``adata.var`` is xarray-backed; column access returns
    a DataArray that lacks ``.unique()`` and other pandas methods.
    """
    var_col = adata.var[col]
    if isinstance(var_col, pd.Series):
        return var_col
    arr = var_col.to_numpy() if hasattr(var_col, "to_numpy") else np.asarray(var_col)
    return pd.Series(arr, index=adata.var_names)


def _get_matrix(adata, layer: str | None):
    """Return the matrix referenced by *layer* (or X if None)."""
    return adata.layers[layer] if layer is not None else adata.X


def _set_matrix(adata, layer: str | None, matrix) -> None:
    if layer is None:
        adata.X = matrix
    else:
        adata.layers[layer] = matrix


def _snapshot(adata, matrix, snapshot_layer: str | None) -> None:
    """Write *matrix* to ``adata.layers[snapshot_layer]`` when the key is not None."""
    if snapshot_layer is None:
        return
    if sp.issparse(matrix):
        adata.layers[snapshot_layer] = matrix.copy()
    else:
        adata.layers[snapshot_layer] = np.array(matrix, copy=True, dtype=np.float64)


def _to_csr_float64(matrix):
    """Return a float64 CSR copy, converting from other sparse formats if needed."""
    if not isinstance(matrix, sp.csr_matrix):
        matrix = matrix.tocsr()
    if matrix.dtype != np.float64:
        matrix = matrix.astype(np.float64)
    return matrix.copy()


def _replace_batch_cols(matrix, batch_indices: np.ndarray, scaled_batch):
    """Return *matrix* with columns *batch_indices* replaced by *scaled_batch*.

    Works for both sparse and dense matrices.
    """
    if not sp.issparse(matrix):
        matrix[:, batch_indices] = scaled_batch
        return matrix

    batch_indices_sorted = np.sort(batch_indices)
    is_contiguous = len(batch_indices_sorted) == 1 or np.all(np.diff(batch_indices_sorted) == 1)

    if is_contiguous and len(batch_indices) == matrix.shape[1]:
        return scaled_batch

    if is_contiguous:
        start_col = int(batch_indices_sorted[0])
        end_col = int(batch_indices_sorted[-1]) + 1
        parts = []
        if start_col > 0:
            parts.append(matrix[:, :start_col])
        parts.append(scaled_batch)
        if end_col < matrix.shape[1]:
            parts.append(matrix[:, end_col:])
        return sp.hstack(parts, format="csr")

    # Non-contiguous: modify data in COO representation
    coo = matrix.tocoo()
    batch_col_set = set(batch_indices)
    # Build a row → scale lookup from scaled_batch
    # scaled_batch is (n_cells × len(batch_indices)); we need per-cell per-col values.
    # For the non-contiguous case we reconstruct via COO of scaled_batch.
    scaled_coo = scaled_batch.tocoo() if sp.issparse(scaled_batch) else None

    if scaled_coo is not None:
        # Map local column index → original column index
        local_to_orig = dict(enumerate(batch_indices_sorted))
        # Build lookup: (row, orig_col) → new_value
        new_val_lookup: dict[tuple[int, int], float] = {}
        for r, c, v in zip(scaled_coo.row, scaled_coo.col, scaled_coo.data, strict=False):
            orig_col = local_to_orig.get(c)
            if orig_col is not None:
                new_val_lookup[(int(r), int(orig_col))] = float(v)

        new_data = coo.data.copy().astype(np.float64)
        for i, (r, c) in enumerate(zip(coo.row, coo.col, strict=False)):
            if c in batch_col_set:
                new_data[i] = new_val_lookup.get((int(r), int(c)), 0.0)
        return sp.csr_matrix((new_data, (coo.row, coo.col)), shape=matrix.shape, dtype=np.float64)

    # scaled_batch is dense
    new_data = coo.data.copy().astype(np.float64)
    col_local = {orig: local for local, orig in enumerate(batch_indices_sorted)}
    for i, (r, c) in enumerate(zip(coo.row, coo.col, strict=False)):
        if c in batch_col_set:
            new_data[i] = float(scaled_batch[int(r), col_local[int(c)]])
    return sp.csr_matrix((new_data, (coo.row, coo.col)), shape=matrix.shape, dtype=np.float64)


# ---------------------------------------------------------------------------
# _normalize_agg_counts
# ---------------------------------------------------------------------------


def _normalize_agg_counts(
    adata,
    *,
    batch_colname: str = "batch_name",
    well_colname: str = "Well",
    run_colname: str | None = None,
    layer: str | None = None,
    snapshot_layer: str | None = None,
    target_counts: float | None = None,
    inplace: bool = True,
):
    """Scale each cell's counts to a target total per batch.

    Parameters
    ----------
    adata:
        AnnData object to normalize.
    batch_colname:
        ``adata.var`` column identifying feature batches.  If the column is
        absent (e.g. 3' untargeted data where all genes form a single panel),
        all genes are treated as one batch — equivalent to a global
        total-count normalization.
    well_colname:
        ``adata.obs`` column identifying wells.
    run_colname:
        ``adata.obs`` column identifying the run/experiment.  Used to
        disambiguate wells that share the same label across runs.
        ``None`` (default) → auto-detect: tries ``"run_name"`` (written by
        ``Loader``) then ``"run_id"`` (legacy); silently skips run
        disambiguation when neither is present (single-run data).
    layer:
        Layer to read from and write to.  ``None`` → operate on ``adata.X``.
    snapshot_layer:
        When not ``None``, the input matrix is copied to this layer before
        scaling is applied.
    target_counts:
        Fixed normalization target (Seurat/scanpy ``target_sum`` style).
        ``None`` → use the per-well median total counts per batch (data-driven).
    inplace:
        Modify ``adata`` in place when ``True``; return a copy when ``False``.
    """
    target = adata if inplace else adata.copy()
    matrix_in = _get_matrix(target, layer)
    is_sparse = sp.issparse(matrix_in)

    _snapshot(target, matrix_in, snapshot_layer)

    matrix = _to_csr_float64(matrix_in) if is_sparse else np.array(matrix_in, copy=True, dtype=np.float64)

    _, well_row_indices = _resolve_well_groups(target.obs, well_colname, run_colname)

    # If no batch column exists (e.g. 3' untargeted, single unified panel),
    # treat all genes as one batch — equivalent to a global total-count norm.
    if batch_colname not in target.var.columns:
        batch_series = pd.Series("__all__", index=target.var_names)
    else:
        batch_series = _var_col_to_series(target, batch_colname)

    for batch in batch_series.unique():
        batch_mask = batch_series.isna() if pd.isna(batch) else batch_series == batch
        batch_cols = np.where(batch_mask)[0]
        if batch_cols.size == 0:
            continue

        if target_counts is not None:
            target_value = float(target_counts)
        else:
            medians = []
            for rows in well_row_indices:
                if rows.size == 0:
                    medians.append(np.nan)
                    continue
                sub = matrix[rows][:, batch_cols]
                counts = np.asarray(sub.sum(axis=1)).ravel() if is_sparse else np.nansum(sub, axis=1)
                medians.append(np.nanmedian(counts))
            target_value = float(np.nanmedian(medians))

        sub = matrix[:, batch_cols]
        norm_values = np.asarray(sub.sum(axis=1)).ravel() if is_sparse else np.nansum(sub, axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            scale = np.divide(
                target_value,
                norm_values,
                out=np.zeros_like(norm_values, dtype=np.float64),
                where=norm_values > 0,
            )

        if is_sparse:
            scale_diag = sp.diags(scale, format="csr", dtype=np.float64)
            scaled_batch = scale_diag @ matrix[:, batch_cols]
            matrix = _replace_batch_cols(matrix, batch_cols, scaled_batch)
        else:
            matrix[:, batch_cols] *= scale[:, None]

    _set_matrix(target, layer, matrix)
    return None if inplace else target


# ---------------------------------------------------------------------------
# _normalize_median_of_ratios
# ---------------------------------------------------------------------------


def _normalize_median_of_ratios(
    adata,
    *,
    batch_colname: str = "batch_name",
    well_colname: str = "Well",
    run_colname: str | None = None,
    layer: str | None = None,
    snapshot_layer: str | None = None,
    inplace: bool = True,
):
    """Apply well-level median-of-ratios correction per batch.

    Computes one scaling factor per well (per batch) and divides every cell
    in that well by its well's factor.  Sparse-safe.

    If ``batch_colname`` is absent from ``adata.var`` (e.g. 3' untargeted
    data), all genes are treated as a single batch and one global scaling
    factor is computed per well.

    The algorithm:
    1. Compute mean expression per well for each batch (one vector per well).
    2. Derive a pseudo-reference as the geometric mean across wells.
    3. For each well: ``scale = exp(median(log(well_mean) - log(pseudo_ref)))``.
    4. Divide every cell in the well by its ``scale``.
    """
    target = adata if inplace else adata.copy()
    matrix_in = _get_matrix(target, layer)
    is_sparse = sp.issparse(matrix_in)

    _snapshot(target, matrix_in, snapshot_layer)

    matrix = _to_csr_float64(matrix_in) if is_sparse else np.array(matrix_in, copy=True, dtype=np.float64)

    _, well_row_indices = _resolve_well_groups(target.obs, well_colname, run_colname)

    # If no batch column exists (e.g. 3' untargeted, single unified panel),
    # treat all genes as one batch — MOR correction becomes a single
    # global well-level scaling factor.
    if batch_colname not in target.var.columns:
        batch_series = pd.Series("__all__", index=target.var_names)
    else:
        batch_series = _var_col_to_series(target, batch_colname)

    for batch in batch_series.unique():
        batch_mask = batch_series.isna() if pd.isna(batch) else batch_series == batch
        batch_cols = np.where(batch_mask)[0]
        if batch_cols.size == 0:
            continue

        well_means = []
        for rows in well_row_indices:
            if rows.size == 0:
                well_means.append(np.full(len(batch_cols), np.nan))
                continue
            sub = matrix[rows][:, batch_cols]
            if is_sparse:
                col_sums = np.asarray(sub.sum(axis=0)).ravel()
                well_means.append(col_sums / len(rows))
            else:
                well_means.append(np.nanmean(sub, axis=0))

        log_means = np.log(np.vstack(well_means) + 1e-8)
        pseudo_ref = np.nanmean(log_means, axis=0)
        log_ratios = log_means - pseudo_ref
        scaling_factors = np.exp(np.nanmedian(log_ratios, axis=1))

        if is_sparse:
            cell_scale = np.ones(matrix.shape[0], dtype=np.float64)
            for rows, sf in zip(well_row_indices, scaling_factors, strict=False):
                if rows.size > 0:
                    cell_scale[rows] = 1.0 / sf
            scale_diag = sp.diags(cell_scale, format="csr", dtype=np.float64)
            scaled_batch = scale_diag @ matrix[:, batch_cols]
            matrix = _replace_batch_cols(matrix, batch_cols, scaled_batch)
        else:
            for rows, sf in zip(well_row_indices, scaling_factors, strict=False):
                if rows.size > 0:
                    matrix[np.ix_(rows, batch_cols)] /= sf

    _set_matrix(target, layer, matrix)
    return None if inplace else target


# ---------------------------------------------------------------------------
# _log1p
# ---------------------------------------------------------------------------


def _log1p(
    adata: ad.AnnData,
    *,
    layer: str | None = None,
    snapshot_layer: str | None = None,
    base: float | None = None,
    inplace: bool = True,
) -> ad.AnnData | None:
    """Element-wise log(1 + x) transformation, sparse-safe.

    For sparse matrices only the stored non-zero values are transformed,
    which preserves the sparsity structure (log1p(0) = 0).
    """
    if base is not None and base <= 0:
        raise ValueError("Logarithm base must be positive.")

    target = adata if inplace else adata.copy()
    matrix_in = _get_matrix(target, layer)
    is_sparse = sp.issparse(matrix_in)

    _snapshot(target, matrix_in, snapshot_layer)

    if is_sparse:
        transformed = _to_csr_float64(matrix_in)
        transformed.data = np.log1p(transformed.data)
        if base is not None and base != np.e:
            transformed.data /= np.log(base)
    else:
        transformed = np.array(matrix_in, copy=True, dtype=np.float64)
        np.log1p(transformed, out=transformed)
        if base is not None and base != np.e:
            transformed /= np.log(base)

    _set_matrix(target, layer, transformed)
    return None if inplace else target


# ---------------------------------------------------------------------------
# _clr
# ---------------------------------------------------------------------------


def _clr(
    adata: ad.AnnData,
    *,
    layer: str | None = None,
    snapshot_layer: str | None = None,
    base: float | None = None,
    batch_colname: str | None = None,
    inplace: bool = True,
) -> ad.AnnData | None:
    """Centered Log-Ratio (CLR) transformation.

    Per cell: ``CLR(x) = log(x + ε) − mean(log(x + ε))`` where ε = 1e-8.
    When ``batch_colname`` is given, CLR is computed independently per batch
    so that the geometric mean is taken only over features in the same batch.

    Always converts to dense because CLR requires all values (including zeros)
    to compute the geometric mean.
    """
    if base is not None and base <= 0:
        raise ValueError("Logarithm base must be positive.")

    target = adata if inplace else adata.copy()
    matrix_in = _get_matrix(target, layer)

    _snapshot(target, matrix_in, snapshot_layer)

    # CLR requires dense arithmetic (geometric mean needs all values)
    if sp.issparse(matrix_in):
        matrix = np.array(matrix_in.toarray(), copy=True, dtype=np.float64)
    else:
        matrix = np.array(matrix_in, copy=True, dtype=np.float64)

    matrix += 1e-8  # pseudocount to avoid log(0)

    if batch_colname is not None and batch_colname in target.var.columns:
        transformed = np.zeros_like(matrix, dtype=np.float64)
        batch_series = _var_col_to_series(target, batch_colname)
        for batch in batch_series.unique():
            batch_cols = np.where(batch_series == batch)[0]
            if batch_cols.size == 0:
                continue
            sub = matrix[:, batch_cols]
            log_sub = np.log(sub)
            mean_log = np.nanmean(log_sub, axis=1, keepdims=True)
            transformed[:, batch_cols] = log_sub - mean_log
    else:
        log_matrix = np.log(matrix)
        mean_log = np.nanmean(log_matrix, axis=1, keepdims=True)
        transformed = log_matrix - mean_log

    if base is not None and base != np.e:
        transformed /= np.log(base)

    _set_matrix(target, layer, transformed)
    return None if inplace else target
