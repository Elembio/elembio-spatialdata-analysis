"""Imaging (CellPaint) morphology utilities.

``select_morphology_channels`` picks the structural CellPaint channels from
the ``CP*`` imaging batch using ``adata.var`` metadata (panel channels are
often reagent-named, e.g. ``EBRedots_RED``, so names alone aren't reliable).
``get_morphology_feature_groups`` groups features into texture, granularity,
radial, intensity, location, and basic-channel groups, and
``aggregate_morphology_features`` collapses related measurements (e.g.
texture orientations).

Cell geometry (``Area``, ``AreaShape_*``, ...) lives in
``adata.obsm['cell_geometry']``; ``Loader.load_tables(...,
include_cell_geometry=True)`` appends it to ``adata.X`` when needed.
"""

from elembio_spatialdata_analysis.tl.morphology._aggregation import aggregate_morphology_features
from elembio_spatialdata_analysis.tl.morphology._feature_groups import get_morphology_feature_groups
from elembio_spatialdata_analysis.tl.morphology._primary_channels import select_morphology_channels

__all__ = [
    "aggregate_morphology_features",
    "get_morphology_feature_groups",
    "select_morphology_channels",
]
