"""Unified normalization entry point: ``normalize_adata``."""

from __future__ import annotations

import anndata as ad

from elembio_spatialdata_analysis.pp.normalize._constants import (
    RAW_LAYER,
)
from elembio_spatialdata_analysis.pp.normalize.modalities import (
    normalize_imaging,
    normalize_proteins,
    normalize_transcripts,
)
from elembio_spatialdata_analysis.tl.constants import normalize_target_type


def normalize_adata(
    adata: ad.AnnData,
    target_type: str,
    *,
    batch_colname: str = "batch_name",
    well_colname: str = "Well",
    run_colname: str | None = None,
    layer: str | None = None,
    log_base: float | None = None,
    raw_layer_name: str | None = RAW_LAYER,
    agg_norm_layer_name: str | None = None,
    mor_norm_layer_name: str | None = None,
    inplace: bool = True,
    # Protein-specific
    apply_agg_norm: bool | None = None,
    target_counts: float | None = None,
) -> ad.AnnData | None:
    """Normalize an AnnData object based on its modality (target type).

    This is the primary entry point, analogous to ``ebsd.pp.filter_cells``.
    Exposed publicly as ``ebsd.pp.normalize``.
    It dispatches to the appropriate modality-specific pipeline:

    * **Transcript / ThreePrimeUntargeted** →
      :func:`normalize_transcripts` (agg-counts → median-of-ratios → log1p)
    * **Protein** →
      :func:`normalize_proteins` (agg-counts + log1p)
    * **Imaging (CellPaint / IF Protein)** →
      :func:`normalize_imaging` (not yet implemented; raises
      ``NotImplementedError`` with guidance)

    Layer snapshots
    ---------------
    By default only ``adata.layers["raw"]`` is written (the raw matrix
    before any transformation), mirroring Scanpy's ``adata.raw`` convention.
    Intermediate snapshots are opt-in: pass ``agg_norm_layer_name`` and/or
    ``mor_norm_layer_name`` explicitly to enable them.  The final normalized
    result always lands in ``adata.X`` (or the specified ``layer``).

    Parameters
    ----------
    adata:
        AnnData object to normalize.
    target_type:
        Modality identifier.  One of ``"Transcript"``, ``"Protein"``,
        ``"Imaging"`` (or ``"CellPaint"``), ``"ThreePrimeUntargeted"`` (or the
        equivalent lowercase form, e.g. ``"transcript"``, ``"imaging"``).
    batch_colname:
        ``adata.var`` column identifying feature batches.
    well_colname:
        ``adata.obs`` column identifying wells.
    run_colname:
        ``adata.obs`` column identifying the run/experiment.
        ``None`` (default) → auto-detect: tries ``"run_name"`` (written by
        ``Loader.load(..., concatenated=True)``) then ``"run_id"`` (legacy).
        Silently skips run disambiguation when neither is present (single-run
        data).  Pass an explicit column name to override.
    layer:
        Layer to read from and write to.  ``None`` → operate on ``adata.X``.
    log_base:
        Logarithm base for log-transform steps.  ``None`` → natural log.
    raw_layer_name:
        Layer name for the pre-normalization snapshot.
        Default: ``"raw"``.
    agg_norm_layer_name:
        Layer name for the post-agg-norm snapshot.
        ``None`` (default) → not written.
        Pass a layer name (conventionally ``"norm_agg_counts"``) to enable.
    mor_norm_layer_name:
        Layer name for the post-median-of-ratios snapshot (Transcript only).
        ``None`` (default) → not written.
        Pass a layer name (conventionally ``"norm_med_of_ratios"``) to enable.
    inplace:
        Modify ``adata`` in place when ``True``; return a normalized copy
        when ``False``.
    apply_agg_norm:
        Protein-only.  Whether to apply aggregated-count scaling before log1p.
        ``None`` → ``True``.  Ignored for non-protein modalities.
    target_counts:
        Fixed normalization target for the agg-count step.  ``None`` →
        data-driven (median total counts per well per batch).

    Returns
    -------
    ``None`` when *inplace* is ``True``, otherwise a normalized copy.

    Examples
    --------
    >>> ebsd.pp.normalize(rna_adata, "Transcript")
    >>> ebsd.pp.normalize(prot_adata, "Protein")
    >>> norm_copy = ebsd.pp.normalize(rna_adata, "Transcript", inplace=False)
    """
    normalized = normalize_target_type(target_type)

    if normalized in ("transcript", "threeprimeuntargeted"):
        return normalize_transcripts(
            adata,
            batch_colname=batch_colname,
            well_colname=well_colname,
            run_colname=run_colname,
            layer=layer,
            log_base=log_base,
            raw_layer_name=raw_layer_name,
            agg_norm_layer_name=agg_norm_layer_name,
            mor_norm_layer_name=mor_norm_layer_name,
            target_counts=target_counts,
            inplace=inplace,
        )

    if normalized == "protein":
        return normalize_proteins(
            adata,
            batch_colname=batch_colname,
            well_colname=well_colname,
            run_colname=run_colname,
            layer=layer,
            log_base=log_base,
            apply_agg_norm=apply_agg_norm,
            raw_layer_name=raw_layer_name,
            agg_norm_layer_name=agg_norm_layer_name,
            target_counts=target_counts,
            inplace=inplace,
        )

    if normalized == "imaging":
        return normalize_imaging(adata, inplace=inplace)

    raise ValueError(
        f"Unrecognized target_type '{target_type}' "
        f"(normalized: '{normalized}'). "
        "Supported types: Transcript, Protein, Imaging (or CellPaint), ThreePrimeUntargeted."
    )
