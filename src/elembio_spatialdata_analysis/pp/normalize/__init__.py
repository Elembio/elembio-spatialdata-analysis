"""Modality-aware normalization for count tables.

``normalize`` is the main entry point and dispatches by modality. The
per-modality pipelines and the ``log1p`` / ``clr`` primitives are available for
custom workflows. Imaging normalization is not implemented yet. See
``SPECIFICATION.md`` in this folder for the pipelines and their defaults.

## Examples

```python
import elembio_spatialdata_analysis as ebsd

ebsd.pp.normalize(adata, "Transcript")
```

Modality-level:

```python
ebsd.pp.normalize_transcripts(adata)
ebsd.pp.normalize_proteins(adata)
```

Primitives:

```python
ebsd.pp.log1p(adata)
ebsd.pp.clr(adata, batch_colname="batch_name")
```

Keep intermediate layers:

```python
ebsd.pp.normalize(
    adata,
    "Transcript",
    agg_norm_layer_name="norm_agg_counts",
    mor_norm_layer_name="norm_med_of_ratios",
)
```
"""

from __future__ import annotations

from elembio_spatialdata_analysis.pp.normalize._constants import (
    AGG_NORM_LAYER,
    MOR_NORM_LAYER,
    RAW_LAYER,
)
from elembio_spatialdata_analysis.pp.normalize._core import _clr as clr
from elembio_spatialdata_analysis.pp.normalize._core import _log1p as log1p
from elembio_spatialdata_analysis.pp.normalize.modalities import (
    normalize_cell_paint,
    normalize_imaging,
    normalize_proteins,
    normalize_transcripts,
)
from elembio_spatialdata_analysis.pp.normalize.normalize import normalize_adata as normalize

__all__ = [
    "AGG_NORM_LAYER",
    "MOR_NORM_LAYER",
    "RAW_LAYER",
    "clr",
    "log1p",
    "normalize",
    "normalize_cell_paint",
    "normalize_imaging",
    "normalize_proteins",
    "normalize_transcripts",
]
