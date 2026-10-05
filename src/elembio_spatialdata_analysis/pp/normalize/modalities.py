"""Modality-specific normalization pipelines.

Public functions:

* ``normalize_transcripts``  — agg-counts → median-of-ratios → log1p
* ``normalize_proteins``     — agg-counts (optional) → log1p
* ``normalize_imaging``      — **not yet implemented** (stub; ``normalize_cell_paint`` alias)

Each function:
1. Snapshots the raw matrix to ``layers["raw"]`` (configurable).
2. Applies the pipeline steps, each snapshotting its input to the next layer.
3. Writes the final normalized result back to ``adata.X`` (or the specified layer).
"""

from __future__ import annotations

import anndata as ad

from elembio_spatialdata_analysis.pp.normalize._constants import (
    RAW_LAYER,
)
from elembio_spatialdata_analysis.pp.normalize._core import (
    _log1p,
    _normalize_agg_counts,
    _normalize_median_of_ratios,
)


def normalize_transcripts(
    adata: ad.AnnData,
    *,
    batch_colname: str = "batch_name",
    well_colname: str = "Well",
    run_colname: str | None = None,
    layer: str | None = None,
    log_base: float | None = None,
    raw_layer_name: str | None = RAW_LAYER,
    agg_norm_layer_name: str | None = None,
    mor_norm_layer_name: str | None = None,
    target_counts: float | None = None,
    inplace: bool = True,
) -> ad.AnnData | None:
    """Normalize Transcript (or ThreePrimeUntargeted) AnnData.

    Pipeline
    --------
    ``X (raw) → agg-count scaling → median-of-ratios → log1p → X (final)``

    By default only the raw matrix is snapshotted (``layers["raw"]``).
    Intermediate layers are written only when their corresponding
    ``*_layer_name`` argument is set to a non-``None`` string.

    Parameters
    ----------
    adata:
        Input AnnData object.
    batch_colname:
        ``adata.var`` column identifying feature batches.
    well_colname:
        ``adata.obs`` column identifying wells.
    run_colname:
        ``adata.obs`` column identifying the run/experiment.
        ``None`` (default) → auto-detect: tries ``"run_name"`` (written by
        ``Loader``) then ``"run_id"`` (legacy); silently skips run
        disambiguation when neither is present (single-run data).
    layer:
        Layer to read from and write to.  ``None`` → operate on ``adata.X``.
    log_base:
        Logarithm base for the final log1p step.  ``None`` → natural log.
    raw_layer_name:
        Layer name for the raw matrix snapshot taken before any transformation.
        Set to ``None`` to skip (saves one matrix copy).
        Default: ``"raw"``.
    agg_norm_layer_name:
        Layer name for the snapshot taken after agg-count scaling and before
        median-of-ratios.  ``None`` (default) → not written.
        Pass a layer name (conventionally ``"norm_agg_counts"``) to enable.
    mor_norm_layer_name:
        Layer name for the snapshot taken after median-of-ratios and before
        log1p.  ``None`` (default) → not written.
        Pass ``MOR_NORM_LAYER`` (``"norm_med_of_ratios"``) to enable.
    target_counts:
        Fixed normalization target for the agg-count step.  ``None`` →
        data-driven (median total counts per well per batch).
    inplace:
        Modify ``adata`` in place when ``True``; return a copy when ``False``.

    Returns
    -------
    ``None`` when *inplace* is ``True``, otherwise a normalized copy.
    """
    target = adata if inplace else adata.copy()

    _normalize_agg_counts(
        target,
        batch_colname=batch_colname,
        well_colname=well_colname,
        run_colname=run_colname,
        layer=layer,
        snapshot_layer=raw_layer_name,
        target_counts=target_counts,
        inplace=True,
    )

    _normalize_median_of_ratios(
        target,
        batch_colname=batch_colname,
        well_colname=well_colname,
        run_colname=run_colname,
        layer=layer,
        snapshot_layer=agg_norm_layer_name,
        inplace=True,
    )

    _log1p(
        target,
        layer=layer,
        snapshot_layer=mor_norm_layer_name,
        base=log_base,
        inplace=True,
    )

    return None if inplace else target


def normalize_proteins(
    adata: ad.AnnData,
    *,
    batch_colname: str = "batch_name",
    well_colname: str = "Well",
    run_colname: str | None = None,
    layer: str | None = None,
    log_base: float | None = None,
    apply_agg_norm: bool | None = None,
    raw_layer_name: str | None = RAW_LAYER,
    agg_norm_layer_name: str | None = None,
    target_counts: float | None = None,
    inplace: bool = True,
) -> ad.AnnData | None:
    """Normalize Protein AnnData.

    Proteins skip the median-of-ratios step — their variance characteristics
    differ from sequencing data and the assumption that most features are
    unchanged does not hold reliably.

    Pipeline
    --------
    ``X (raw) → agg-count scaling (optional) → log1p → X (final)``

    By default only the raw matrix is snapshotted (``layers["raw"]``).
    The intermediate agg-norm layer is written only when
    ``agg_norm_layer_name`` is set to a non-``None`` string.

    Parameters
    ----------
    adata:
        Input AnnData object.
    batch_colname:
        ``adata.var`` column identifying feature batches.
    well_colname:
        ``adata.obs`` column identifying wells.
    run_colname:
        ``adata.obs`` column identifying the run/experiment.
        ``None`` (default) → auto-detect: tries ``"run_name"`` (written by
        ``Loader``) then ``"run_id"`` (legacy); silently skips run
        disambiguation when neither is present (single-run data).
    layer:
        Layer to read from and write to.  ``None`` → operate on ``adata.X``.
    log_base:
        Logarithm base for log1p.  ``None`` → natural log.
    apply_agg_norm:
        Whether to apply aggregated-count scaling before log1p.
        ``None`` (default) → ``True``.
    raw_layer_name:
        Layer name for the raw matrix snapshot.  ``None`` → not written.
        Default: ``"raw"``.
    agg_norm_layer_name:
        Layer name for the snapshot after agg-count scaling.
        ``None`` (default) → not written.
        Pass a layer name (conventionally ``"norm_agg_counts"``) to enable.
    target_counts:
        Fixed normalization target for the agg-count step.  ``None`` →
        data-driven.  Ignored when ``apply_agg_norm`` is ``False``.
    inplace:
        Modify ``adata`` in place when ``True``; return a copy when ``False``.

    Returns
    -------
    ``None`` when *inplace* is ``True``, otherwise a normalized copy.
    """
    target = adata if inplace else adata.copy()

    if apply_agg_norm is None:
        apply_agg_norm = True

    if apply_agg_norm:
        _normalize_agg_counts(
            target,
            batch_colname=batch_colname,
            well_colname=well_colname,
            run_colname=run_colname,
            layer=layer,
            snapshot_layer=raw_layer_name,
            target_counts=target_counts,
            inplace=True,
        )
        pre_log_snapshot = agg_norm_layer_name
    else:
        from elembio_spatialdata_analysis.pp.normalize._core import (
            _get_matrix,
            _snapshot,
        )

        _snapshot(target, _get_matrix(target, layer), raw_layer_name)
        pre_log_snapshot = None

    _log1p(
        target,
        layer=layer,
        snapshot_layer=pre_log_snapshot,
        base=log_base,
        inplace=True,
    )

    return None if inplace else target


def normalize_imaging(
    adata: ad.AnnData,
    *,
    inplace: bool = True,
) -> ad.AnnData | None:
    """Normalize an imaging (CellPaint / IF Protein) AnnData.

    .. note::
        **Not yet implemented.**  Fluorescence intensities are not counts, so
        the count-based pipelines (aggregated counts, median-of-ratios, log1p)
        do not apply.

    Raises
    ------
    NotImplementedError
        Always; imaging normalization is not implemented yet.
    """
    raise NotImplementedError(
        "Imaging (CellPaint / IF Protein) normalization is not yet implemented: "
        "fluorescence intensities are not counts, so the count-based pipelines don't apply."
    )


# Backward-compat alias (CellPaint -> Imaging rename).
normalize_cell_paint = normalize_imaging
