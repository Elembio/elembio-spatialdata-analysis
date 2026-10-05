"""Analysis tools for Element Biosciences SpatialData objects.

Subpackages: ``load`` (read stores into AnnData / Polars), ``pp`` (QC
filtering and normalization), ``qc`` (run-level QC gates and rollups),
``tl`` (analysis tools), and ``vis`` (plotting).
"""

from importlib.metadata import PackageNotFoundError, metadata, version

try:
    __version__ = version("elembio-spatialdata-analysis")
    __package_name__ = metadata("elembio-spatialdata-analysis")["Name"]
except PackageNotFoundError:
    __version__ = "unknown"
    __package_name__ = "elembio-spatialdata-analysis"

# Subpackages are imported on first access so `import elembio_spatialdata_analysis`
# stays fast and doesn't load the heavy scientific stack up front.
__all__ = ["load", "pp", "qc", "tl", "vis"]


def __getattr__(name):
    """Lazy import of submodules to avoid eager loading of analysis dependencies."""
    import importlib
    import sys

    if name in ("load", "pp", "qc", "tl", "vis"):
        module = importlib.import_module(f"{__name__}.{name}")
        setattr(sys.modules[__name__], name, module)
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
