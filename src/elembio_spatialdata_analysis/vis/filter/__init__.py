"""Per-well QC plots: histograms, violins, and pre/post-filter comparisons.

All plot functions are also available on ``elembio_spatialdata_analysis.vis``.
``plot_cell_paint_intensities`` / ``plot_cell_paint_morphology`` are aliases of
the ``plot_imaging_*`` functions. See ``SPECIFICATION.md`` in this folder for
every parameter.

## Examples

```python
import elembio_spatialdata_analysis as ebsd

ebsd.vis.plot_qc_per_well(adata, "total_counts")
ebsd.vis.plot_qc_per_well(adata, "Area")
ebsd.vis.plot_violin_per_well(adata, "pct_counts_mt")
ebsd.vis.plot_filtered_histograms_per_well(adata, flags, metric="total_counts")
ebsd.vis.plot_total_counts_by_batch(adata)
ebsd.vis.plot_imaging_intensities(adata)
ebsd.vis.plot_imaging_morphology(adata, "MeanIntensity")
```

## Prerequisite

Call ``ebsd.pp.compute_qc_metrics(adata, target_type)`` before plotting metrics
from ``obsm['qc_metrics']``.
"""

from __future__ import annotations

from elembio_spatialdata_analysis.vis.filter.plots import (
    plot_cell_paint_intensities,
    plot_cell_paint_morphology,
    plot_filtered_histograms_per_well,
    plot_imaging_intensities,
    plot_imaging_morphology,
    plot_qc_per_well,
    plot_total_counts_by_batch,
    plot_violin_per_well,
)

__all__ = [
    "plot_cell_paint_intensities",
    "plot_cell_paint_morphology",
    "plot_filtered_histograms_per_well",
    "plot_imaging_intensities",
    "plot_imaging_morphology",
    "plot_qc_per_well",
    "plot_total_counts_by_batch",
    "plot_violin_per_well",
]
