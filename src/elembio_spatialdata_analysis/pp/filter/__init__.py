"""Modality-aware cell and gene QC filtering.

Typical use: ``compute_qc_metrics`` → ``filter_cells`` (or ``filter_modalities``
for several aligned tables) → ``generate_filter_report``. ``filter_genes``
filters genes (mainly for 3′ data), and ``compute_filter_diagnostics`` computes
the flags without removing cells. See ``SPECIFICATION.md`` in this folder for
every filter, flag column, and default threshold.

## Lower-level public helpers

The names below are exported as public API so downstream code never needs
to import the underscore submodules (``_constants`` / ``_flags`` /
``_quantile``) directly.  Renaming or dropping one is a breaking change.

- DEFAULT_BOUNDS_PER_MODALITY, QC_METRICS_OBSM_KEY, FILTER_FLAGS_FULL_UNS_KEY: filter constants / obsm+uns keys.
- AREA_OWNER_MODALITIES: modalities that own the cell-area metric.
- first_failing_reason, flag_data_missing, flag_empty_cells, flag_failed_batches: per-cell filter-flag helpers.
- resolve_threshold_defaults: default quantile thresholds for a modality.
"""

from elembio_spatialdata_analysis.pp.filter._constants import (
    DEFAULT_BOUNDS_PER_MODALITY,
    FILTER_FLAGS_FULL_UNS_KEY,
    QC_METRICS_OBSM_KEY,
)
from elembio_spatialdata_analysis.pp.filter._flags import (
    _flag_data_missing as flag_data_missing,
)
from elembio_spatialdata_analysis.pp.filter._flags import (
    _flag_empty_cells as flag_empty_cells,
)
from elembio_spatialdata_analysis.pp.filter._flags import (
    _flag_failed_batches as flag_failed_batches,
)
from elembio_spatialdata_analysis.pp.filter._flags import (
    first_failing_reason,
)
from elembio_spatialdata_analysis.pp.filter._qc_metrics import compute_qc_metrics
from elembio_spatialdata_analysis.pp.filter._quantile import (
    _resolve_threshold_defaults as resolve_threshold_defaults,
)
from elembio_spatialdata_analysis.pp.filter.filter_cells import (
    FilterDiagnostics,
    ModalityFilterResult,
    compute_filter_diagnostics,
    filter_cells,
    filter_genes,
    get_filter_flags,
)
from elembio_spatialdata_analysis.pp.filter.filter_modalities import (
    _AREA_OWNER_MODALITIES as AREA_OWNER_MODALITIES,
)
from elembio_spatialdata_analysis.pp.filter.filter_modalities import filter_modalities
from elembio_spatialdata_analysis.pp.filter.report import generate_filter_report

__all__ = [
    "AREA_OWNER_MODALITIES",
    "DEFAULT_BOUNDS_PER_MODALITY",
    "FILTER_FLAGS_FULL_UNS_KEY",
    "QC_METRICS_OBSM_KEY",
    "FilterDiagnostics",
    "ModalityFilterResult",
    "compute_filter_diagnostics",
    "compute_qc_metrics",
    "filter_cells",
    "filter_genes",
    "filter_modalities",
    "first_failing_reason",
    "flag_data_missing",
    "flag_empty_cells",
    "flag_failed_batches",
    "generate_filter_report",
    "get_filter_flags",
    "resolve_threshold_defaults",
]
