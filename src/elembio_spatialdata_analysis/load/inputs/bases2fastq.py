"""Parser for ``metadata/inputs/bases2fastq/RunStats.json``.

bases2fastq emits a flowcell-level summary (root scalars), per-sequencing-batch
records (``Batches[]``), per-well records nested under each batch
(``Batches[].Wells[]``), and per-base-position telemetry (``Bases[]``) within
each well.  The parser exposes typed access to the run-level scalars and a
flat list of per-batch / per-well rows for downstream rendering.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from elembio_spatialdata_analysis.load.backends import RunInfo

_RUNSTATS_PATH = "metadata/inputs/bases2fastq/RunStats.json"


@dataclass
class Bases2FastqStats:
    """Typed view over ``bases2fastq/RunStats.json``."""

    raw: dict[str, Any]
    run_id: str | None = None
    run_name: str | None = None
    flowcell_id: str | None = None
    analysis_version: str | None = None
    num_polonies: int | None = None
    num_polonies_before_trimming: int | None = None
    total_yield: int | None = None
    percent_q30: float | None = None
    quality_score_mean: float | None = None
    batches: list[dict[str, Any]] = field(default_factory=list)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Bases2FastqStats:
        return cls(
            raw=data,
            run_id=data.get("RunID"),
            run_name=data.get("RunName"),
            flowcell_id=data.get("FlowCellID"),
            analysis_version=data.get("AnalysisVersion"),
            num_polonies=data.get("NumPolonies"),
            num_polonies_before_trimming=data.get("NumPoloniesBeforeTrimming"),
            total_yield=data.get("TotalYield"),
            percent_q30=data.get("PercentQ30"),
            quality_score_mean=data.get("QualityScoreMean"),
            batches=list(data.get("Batches") or []),
        )

    def per_batch_summary(self) -> list[dict[str, Any]]:
        """Return a flat list of one row per sequencing batch.

        ``NumCycles`` is derived from the length of the first well's
        ``Bases[]`` array (one entry per cycle position) so consumers don't
        need to drill into per-base-position telemetry just to learn how
        long a batch was sequenced.
        """
        rows: list[dict[str, Any]] = []
        for batch in self.batches:
            wells = batch.get("Wells") or []
            num_cycles = len(wells[0].get("Bases") or []) if wells else 0
            rows.append(
                {
                    "BatchName": batch.get("BatchName"),
                    "NumPolonies": batch.get("NumPolonies"),
                    "NumPoloniesBeforeTrimming": batch.get("NumPoloniesBeforeTrimming"),
                    "TotalYield": batch.get("TotalYield"),
                    "PercentQ30": batch.get("PercentQ30"),
                    "QualityScoreMean": batch.get("QualityScoreMean"),
                    "NumCycles": num_cycles,
                }
            )
        return rows

    def per_well_summary(self) -> list[dict[str, Any]]:
        """Return a flat list of one row per (batch, well)."""
        rows: list[dict[str, Any]] = []
        for batch in self.batches:
            batch_name = batch.get("BatchName")
            for well in batch.get("Wells") or []:
                rows.append(
                    {
                        "BatchName": batch_name,
                        "WellLocation": well.get("WellLocation"),
                        "WellLabel": well.get("WellLabel"),
                        "NumPolonies": well.get("NumPolonies"),
                        "PercentQ30": well.get("PercentQ30"),
                        "QualityScoreMean": well.get("QualityScoreMean"),
                    }
                )
        return rows

    def per_base_for_batch(self, batch_name: str) -> list[dict[str, Any]] | None:
        """Return per-base-position telemetry for a batch (averaged across its wells).

        Each batch's first well usually has a representative ``Bases[]`` array.
        Returns the raw list of cycle dicts or ``None`` if absent.
        """
        for batch in self.batches:
            if batch.get("BatchName") != batch_name:
                continue
            wells = batch.get("Wells") or []
            if not wells:
                return None
            return list(wells[0].get("Bases") or []) or None
        return None

    def per_well_bases(self, batch_name: str) -> list[dict[str, Any]]:
        """Per-well per-cycle telemetry for a batch.

        Returns one row per well: ``{WellLocation, WellLabel, Bases}`` where
        ``Bases`` is the raw per-cycle list (``[{Position, PercentQ30,
        QualityScoreMean, ...}, ...]``).  Unlike :meth:`per_base_for_batch`
        (which returns only the first well's representative array), this keeps
        every well so a consumer can build a wells x cycles matrix — e.g. the
        per-well %Q30 heatmap.  Empty list when the batch is absent / has no
        per-well ``Bases``.
        """
        for batch in self.batches:
            if batch.get("BatchName") != batch_name:
                continue
            rows: list[dict[str, Any]] = []
            for well in batch.get("Wells") or []:
                bases = well.get("Bases") or []
                if not bases:
                    continue
                rows.append(
                    {
                        "WellLocation": well.get("WellLocation"),
                        "WellLabel": well.get("WellLabel"),
                        "Bases": list(bases),
                    }
                )
            return rows
        return []


def parse_bases2fastq_runstats(run_info: RunInfo) -> Bases2FastqStats | None:
    """Read and parse ``bases2fastq/RunStats.json`` from a SpatialData store."""
    try:
        data = run_info._read_json(_RUNSTATS_PATH)
    except Exception:
        return None
    if not data:
        return None
    return Bases2FastqStats.from_dict(data)
