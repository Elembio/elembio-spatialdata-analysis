"""Parser for the pipeline-level provenance group, ``metadata/inputs/pipeline/Provenance.json``.

The orchestrating pipeline's provenance group (its own ``metadata/inputs/``
source-mirror, alongside the per-tool ``cells2stats`` / ``bases2fastq`` /
``fastq2parquet_3prime`` RunStats). Carries the orchestrator metadata:
``Workflow`` (version / commit / revision / runName / commandLine), ``Tools``
(container image:tags), and ``Params`` (the resolved workflow parameters).
Present only on runs built with the pipeline-provenance input. Stores built
before the group was renamed carry it as ``metadata/inputs/tetonatlas/``;
both names are read, newest first.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from elembio_spatialdata_analysis.load.backends import RunInfo

_PROVENANCE_PATHS = (
    "metadata/inputs/pipeline/Provenance.json",
    "metadata/inputs/tetonatlas/Provenance.json",  # stores built before the rename
)


@dataclass
class PipelineProvenance:
    """Typed view over the pipeline ``Provenance.json``."""

    raw: dict[str, Any]
    workflow: dict[str, Any] = field(default_factory=dict)
    tools: dict[str, Any] = field(default_factory=dict)
    params: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> PipelineProvenance:
        return cls(
            raw=data,
            workflow=dict(data.get("Workflow") or {}),
            tools=dict(data.get("Tools") or {}),
            params=dict(data.get("Params") or {}),
        )


def parse_pipeline_provenance(run_info: RunInfo) -> PipelineProvenance | None:
    """Read and parse the pipeline ``Provenance.json``; ``None`` when absent.

    Tries the current group name first, then the pre-rename name.
    """
    for path in _PROVENANCE_PATHS:
        try:
            data = run_info._read_json(path)
        except Exception:
            continue
        # Absent and empty both fall through. The builder writes exactly one of these
        # groups per store, so the fallback only ever reaches pre-rename stores.
        if data:
            return PipelineProvenance.from_dict(data)
    return None
