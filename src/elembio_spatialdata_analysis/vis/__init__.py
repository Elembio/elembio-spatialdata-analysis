"""Plotting: per-well QC distributions and well images.

QC plots (``plot_qc_per_well``, ``plot_violin_per_well``,
``plot_filtered_histograms_per_well``, ``plot_total_counts_by_batch``,
``plot_imaging_*``) come from ``vis.filter``. Well-image plots
(``list_image_channels``, ``plot_image``, ``plot_cell``) come from
``vis.spatial``. See those subpackages for examples.
"""

from elembio_spatialdata_analysis.vis.filter import (
    plot_cell_paint_intensities,
    plot_cell_paint_morphology,
    plot_filtered_histograms_per_well,
    plot_imaging_intensities,
    plot_imaging_morphology,
    plot_qc_per_well,
    plot_total_counts_by_batch,
    plot_violin_per_well,
)
from elembio_spatialdata_analysis.vis.spatial import list_image_channels, plot_cell, plot_image

__all__ = [
    "list_image_channels",
    "plot_cell",
    "plot_cell_paint_intensities",
    "plot_cell_paint_morphology",
    "plot_filtered_histograms_per_well",
    "plot_image",
    "plot_imaging_intensities",
    "plot_imaging_morphology",
    "plot_qc_per_well",
    "plot_total_counts_by_batch",
    "plot_violin_per_well",
]
