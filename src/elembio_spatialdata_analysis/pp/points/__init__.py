"""Polony-level (points) preprocessing utilities.

Points / polony parquets (one row per read after demultiplexing) carry
sequencing-level QC fields that aren't aggregated into the per-cell
modality tables.  This module provides helpers to read those fields
robustly across upstream schema versions: ``resolve_quality_column`` finds
the mean-base-quality column (``MeanBaseQ``, or ``MeanBaseQuality`` in older
stores), and ``polony_qc_stats`` summarizes polony quality.
"""

from elembio_spatialdata_analysis.pp.points.qc import (
    POLONY_QUALITY_COLUMNS,
    PolonyQcStats,
    polony_qc_stats,
    resolve_quality_column,
)

__all__ = [
    "POLONY_QUALITY_COLUMNS",
    "PolonyQcStats",
    "polony_qc_stats",
    "resolve_quality_column",
]
