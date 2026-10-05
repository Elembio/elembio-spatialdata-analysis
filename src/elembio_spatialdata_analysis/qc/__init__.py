"""Run-level quality control: pass / warn / fail gates and summary tables.

``QcGates`` evaluates metrics against thresholds (defaults ship in
``qc/config/qc_gates.yml``), and the ``per_*`` helpers build per-batch-type,
per-well, and per-tile summary tables from a run's metadata.
"""

from elembio_spatialdata_analysis.qc.aggregate import (
    b2f_for_batches,
    join_batch_metadata,
    per_batch_type_summary,
    per_tile_cell_count,
    per_tile_total_counts,
    per_well_overview,
)
from elembio_spatialdata_analysis.qc.gates import GateResult, QcGates, load_default_gates

__all__ = [
    "GateResult",
    "QcGates",
    "b2f_for_batches",
    "join_batch_metadata",
    "load_default_gates",
    "per_batch_type_summary",
    "per_tile_cell_count",
    "per_tile_total_counts",
    "per_well_overview",
]
