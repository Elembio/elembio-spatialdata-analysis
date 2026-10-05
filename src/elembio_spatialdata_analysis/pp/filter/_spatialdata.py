"""SpatialData bookkeeping helpers for the filter module."""

from __future__ import annotations

import anndata as ad

from elembio_spatialdata_analysis.pp.filter._qc_metrics import _obs_col_to_series


def _update_spatialdata_attrs(adata: ad.AnnData) -> None:
    """Update ``spatialdata_attrs`` to reflect regions still present after filtering.

    If ``adata.uns["spatialdata_attrs"]`` exists, updates the ``region`` list to
    only include regions that are still present in ``adata.obs[region_key]``.
    No-op when the key is absent or malformed.
    """
    if "spatialdata_attrs" not in adata.uns:
        return
    attrs = adata.uns["spatialdata_attrs"]
    if not isinstance(attrs, dict):
        return
    region_key = attrs.get("region_key")
    if region_key is None or region_key not in adata.obs.columns:
        return
    attrs["region"] = _obs_col_to_series(adata, region_key).unique().tolist()
