"""Typed parsers for ``metadata/inputs/*`` files mirrored into the SpatialData store.

Each parser reads one source pipeline's archive (bases2fastq, cells2stats,
fastq2parquet_3prime) or run-level inputs (RunParameters, Panel, RunManifest)
into a small dataclass plus the raw dict.  IO is delegated to
:meth:`RunInfo._read_json` so parsers work transparently for local
``.zarr`` directories, ``.zarr.zip`` archives, and S3 stores.

Parsers return ``None`` when the source file is missing — callers should
treat absence as "this section does not apply".
"""

from elembio_spatialdata_analysis.load.inputs.bases2fastq import (
    Bases2FastqStats,
    parse_bases2fastq_runstats,
)
from elembio_spatialdata_analysis.load.inputs.cells2stats import (
    Cells2StatsStats,
    parse_average_norm_well_stats,
    parse_cells2stats_runstats,
    parse_cells2stats_versions,
    parse_target_cell_assignment,
    parse_target_counts,
)
from elembio_spatialdata_analysis.load.inputs.fastq2parquet_3p import (
    Fastq2Parquet3pStats,
    parse_fastq2parquet_3p_runstats,
    parse_multiqc_summary,
)
from elembio_spatialdata_analysis.load.inputs.pipeline_provenance import (
    PipelineProvenance,
    parse_pipeline_provenance,
)
from elembio_spatialdata_analysis.load.inputs.run_inputs import (
    PanelInfo,
    RunInputs,
    parse_panel,
    parse_run_inputs,
    parse_run_manifest,
    parse_run_parameters,
)
from elembio_spatialdata_analysis.load.inputs.witness import parse_witness_hashes

__all__ = [
    "Bases2FastqStats",
    "Cells2StatsStats",
    "Fastq2Parquet3pStats",
    "PanelInfo",
    "PipelineProvenance",
    "RunInputs",
    "parse_average_norm_well_stats",
    "parse_bases2fastq_runstats",
    "parse_cells2stats_runstats",
    "parse_cells2stats_versions",
    "parse_fastq2parquet_3p_runstats",
    "parse_multiqc_summary",
    "parse_panel",
    "parse_pipeline_provenance",
    "parse_run_inputs",
    "parse_run_manifest",
    "parse_run_parameters",
    "parse_target_cell_assignment",
    "parse_target_counts",
    "parse_witness_hashes",
]
