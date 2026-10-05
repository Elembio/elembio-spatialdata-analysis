"""Snapshot layer names for the normalize module."""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Layer name constants — opt-in snapshot names for normalization pipelines.
#
# By default only RAW_LAYER is written (before any transformation).
# Intermediate layers AGG_NORM_LAYER and MOR_NORM_LAYER are NOT written by
# default to avoid storing 2–3 extra full matrix copies in memory.
#
# To enable intermediate snapshots, pass layer names explicitly:
#
#   ebsd.pp.normalize(
#       adata, "Transcript",
#       agg_norm_layer_name="norm_agg_counts",
#       mor_norm_layer_name="norm_med_of_ratios",
#   )
# ---------------------------------------------------------------------------

# Copy of X before any normalization step (written by default).
# Named "raw" to match Scanpy's adata.raw convention.
RAW_LAYER = "raw"
# After aggregated-count (library-size) scaling, before median-of-ratios.
AGG_NORM_LAYER = "norm_agg_counts"
# After median-of-ratios, before log1p (transcript pipelines only).
MOR_NORM_LAYER = "norm_med_of_ratios"
