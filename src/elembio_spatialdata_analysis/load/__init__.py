"""Load SpatialData stores (local, ``.zarr.zip``, or S3) into AnnData and Polars.

``Loader`` is the main entry point and handles one or many runs. ``RunInfo``
holds per-run metadata and point selection, and ``inputs`` has typed parsers
for the pipeline QC files under ``metadata/inputs/``.
"""

from elembio_spatialdata_analysis.load import inputs
from elembio_spatialdata_analysis.load.backends import RunInfo
from elembio_spatialdata_analysis.load.loader import Loader

__all__ = [
    "Loader",
    "RunInfo",
    "inputs",
]
