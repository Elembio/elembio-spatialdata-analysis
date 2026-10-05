"""Parser for ``metadata/inputs/witness_hashes.json``.

The witness file mirrors run-level metadata into the per-pipeline staging
folders and records hash equivalence for each copy, so consumers can flag
any copy mismatch.
"""

from __future__ import annotations

from typing import Any

from elembio_spatialdata_analysis.load.backends import RunInfo

_WITNESS_PATH = "metadata/inputs/witness_hashes.json"


def parse_witness_hashes(run_info: RunInfo) -> dict[str, Any] | None:
    """Read raw ``witness_hashes.json``.

    Each top-level key is a metadata filename (``RunParameters.json``,
    ``RunManifest.json``, ``Panel.json``, optionally ``TargetCellAssignmentManifest.*``)
    mapping to ``{canonical_location, copies}``.  Returns ``None`` when absent.
    """
    try:
        return run_info._read_json(_WITNESS_PATH)
    except Exception:
        return None
