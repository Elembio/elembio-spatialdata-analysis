"""Well-image plotting.

Public API
----------
- ``list_image_channels(sdata, well)`` — channel names on the well image.
- ``plot_image(sdata, well=...)`` — fluorescence image for one well.
  Omit ``well`` (or pass ``wells=``) for a plate grid of the same images.
  Choose one or more channels, tint each, optionally draw segmentation
  outlines, and optionally fill cells by a table column. Set
  ``show_image=False`` for segmentation without fluorescence.
- ``plot_cell(sdata, cell_id)`` — fixed-size crop around one cell, with
  that cell's stored shape drawn as a vector contour.

Images and labels are well-keyed (``A1-images-cell_paint``, ``A1-labels-Cell``).
Pass an opened ``SpatialData`` — ``Loader`` does not expose rasters.

Well images are large (~15k px at scale 0). ``scale=None`` selects the
coarsest pyramid level. Label pixel values match ``obs["Cell"]``.
"""

from elembio_spatialdata_analysis.vis.spatial._access import (
    label_image,
    list_image_channels,
    resolve_image_key,
    resolve_label_key,
    resolve_shapes_key,
)
from elembio_spatialdata_analysis.vis.spatial._raster import image_composite
from elembio_spatialdata_analysis.vis.spatial.plots import plot_cell, plot_image

__all__ = [
    "image_composite",
    "label_image",
    "list_image_channels",
    "plot_cell",
    "plot_image",
    "resolve_image_key",
    "resolve_label_key",
    "resolve_shapes_key",
]
