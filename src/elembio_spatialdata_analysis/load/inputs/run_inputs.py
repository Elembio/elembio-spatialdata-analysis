"""Parsers for ``metadata/inputs/run/`` — the canonical run recipe.

Includes ``RunParameters.json`` (instrument / panel / batches / wells / tile
geometry), ``Panel.json`` (panel definition: imaging targets, barcoding
targets, DISS primer tubes, assay types per batch), and ``RunManifest.json``
(operator-supplied per-well metadata: cell type, diameter, primer sequences).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from elembio_spatialdata_analysis.load.backends import RunInfo

_RUN_PARAMS_PATH = "metadata/inputs/run/RunParameters.json"
_RUN_MANIFEST_JSON = "metadata/inputs/run/RunManifest.json"
_PANEL_PATH = "metadata/inputs/run/Panel.json"


@dataclass
class PanelInfo:
    """Typed view over ``Panel.json``."""

    raw: dict[str, Any]
    panel_id: str | None = None
    name: str | None = None
    kit_type: str | None = None
    includes_cell_assignment: bool | None = None
    imaging_targets: list[dict[str, Any]] = field(default_factory=list)
    barcoding_targets: list[dict[str, Any]] = field(default_factory=list)
    diss_primer_tubes: list[dict[str, Any]] = field(default_factory=list)
    imaging_primer_tubes: list[dict[str, Any]] = field(default_factory=list)
    barcoding_primer_tubes: list[dict[str, Any]] = field(default_factory=list)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> PanelInfo:
        return cls(
            raw=data,
            panel_id=data.get("ID"),
            name=data.get("Name"),
            kit_type=data.get("KitType"),
            includes_cell_assignment=data.get("IncludesCellAssignment"),
            imaging_targets=list(data.get("ImagingTargets") or []),
            barcoding_targets=list(data.get("BarcodingTargets") or []),
            diss_primer_tubes=list(data.get("DISSPrimerTubes") or []),
            imaging_primer_tubes=list(data.get("ImagingPrimerTubes") or []),
            barcoding_primer_tubes=list(data.get("BarcodingPrimerTubes") or []),
        )

    def batch_assay_types(self) -> dict[str, str]:
        """Return ``{batch_name: assay_type}`` over all primer-tube definitions.

        Useful for joining sequencing batches (B0X) to their declared assay
        type (Transcript / Protein / SpecializedTargeted / ThreePrimeUntargeted)
        and imaging batches (CP0X) to Imaging.
        """
        out: dict[str, str] = {}
        for src_list, default_type in [
            (self.diss_primer_tubes, None),
            (self.barcoding_primer_tubes, None),
            (self.imaging_primer_tubes, "Imaging"),
        ]:
            for tube in src_list:
                name = tube.get("BatchName") or tube.get("Name")
                assay = tube.get("AssayType") or default_type
                if name and assay:
                    out[str(name)] = str(assay)
        return out

    def control_target_sites(self) -> set[str]:
        """Return the set of OPS target-site names flagged ``IsControl=True``.

        Source of truth for which OPS target site is the QC control: reads
        ``DISSPrimerTubes[].DISSTargets[].IsControl`` from ``Panel.json``.
        The control panel is a known reference (e.g. a multiplexed
        ``QC_Panel``) whose chemistry-level performance gates run health,
        as opposed to the customer's own experimental target sites.

        Returns an empty set when no target site is flagged (legacy panels
        that predate the field, or runs with no control) — callers then
        treat every site as a customer/assay site.
        """
        out: set[str] = set()
        for tube in self.diss_primer_tubes:
            for tgt in tube.get("DISSTargets") or []:
                if tgt.get("IsControl") is True and tgt.get("TargetSiteName"):
                    out.add(str(tgt["TargetSiteName"]))
        return out

    def batch_run_orders(self) -> dict[str, int]:
        """Return ``{batch_name: RunOrder}`` across every primer-tube list.

        ``RunOrder`` is the integer recorded by the instrument recipe
        for the execution sequence on the chip (e.g. CP01=1, B01=3,
        B06=8).  Useful for ordering batch-level tables by the order
        they actually ran on the AVITI24, rather than the file-insertion
        / alphabetical order the build pipeline happens to use.
        """
        out: dict[str, int] = {}
        for src_list in (
            self.imaging_primer_tubes,
            self.diss_primer_tubes,
            self.barcoding_primer_tubes,
        ):
            for tube in src_list:
                name = tube.get("BatchName") or tube.get("Name")
                ro = tube.get("RunOrder")
                if name and ro is not None:
                    try:
                        out[str(name)] = int(ro)
                    except (TypeError, ValueError):
                        continue
        return out


@dataclass
class RunInputs:
    """Aggregate of run-level inputs (RunParameters + Panel + RunManifest)."""

    run_parameters: dict[str, Any] | None
    panel: PanelInfo | None
    run_manifest: dict[str, Any] | None

    @property
    def run_name(self) -> str | None:
        """Run name from ``RunParameters.RunName``, or ``None`` if absent."""
        if self.run_parameters:
            return self.run_parameters.get("RunName")
        return None

    @property
    def run_date(self) -> str | None:
        """Run date string from ``RunParameters.Date``, or ``None`` if absent."""
        if self.run_parameters:
            return self.run_parameters.get("Date")
        return None

    @property
    def instrument_name(self) -> str | None:
        """Instrument name from ``RunParameters.InstrumentName``, or ``None``."""
        if self.run_parameters:
            return self.run_parameters.get("InstrumentName")
        return None

    @property
    def operator_name(self) -> str | None:
        """Operator name from ``RunParameters.OperatorName``, or ``None``."""
        if self.run_parameters:
            return self.run_parameters.get("OperatorName")
        return None

    @property
    def kit_type(self) -> str | None:
        """Kit type from ``RunParameters.KitType``, or ``None``."""
        if self.run_parameters:
            return self.run_parameters.get("KitType")
        return None

    @property
    def panel_name(self) -> str | None:
        """Panel name from ``Panel.json``, falling back to ``RunParameters.PanelName``."""
        if self.panel and self.panel.name:
            return self.panel.name
        if self.run_parameters:
            return self.run_parameters.get("PanelName")
        return None

    @property
    def wells(self) -> list[dict[str, Any]]:
        """Per-well entries from ``RunParameters.Wells`` (empty list if absent)."""
        if self.run_parameters:
            return list(self.run_parameters.get("Wells") or [])
        return []

    @property
    def well_layout(self) -> str | None:
        """Configured flowcell well layout, e.g. ``"OneWell"`` (single mega-well).

        Sourced from ``RunParameters.WellLayout`` — the run-config source of
        truth for how many wells the flowcell physically has, independent of how
        many wells produced data.  Lets per-well grids tell a true single-well
        run apart from a sparsely-populated multi-well plate.
        """
        if self.run_parameters:
            return self.run_parameters.get("WellLayout")
        return None

    @property
    def batch_definitions(self) -> list[dict[str, Any]]:
        """Batch entries from ``RunParameters.Batches`` (empty list if absent)."""
        if self.run_parameters:
            return list(self.run_parameters.get("Batches") or [])
        return []

    @property
    def manifest_wells(self) -> list[dict[str, Any]]:
        """Per-well entries from ``RunManifest.json`` ``Wells`` (empty list if absent)."""
        if self.run_manifest:
            return list(self.run_manifest.get("Wells") or [])
        return []


def parse_run_parameters(run_info: RunInfo) -> dict[str, Any] | None:
    """Read raw ``RunParameters.json`` from ``metadata/inputs/run/``."""
    try:
        return run_info._read_json(_RUN_PARAMS_PATH)
    except Exception:
        return None


def parse_panel(run_info: RunInfo) -> PanelInfo | None:
    """Read ``Panel.json`` and wrap in :class:`PanelInfo`."""
    try:
        data = run_info._read_json(_PANEL_PATH)
    except Exception:
        return None
    if not data:
        return None
    return PanelInfo.from_dict(data)


def parse_run_manifest(run_info: RunInfo) -> dict[str, Any] | None:
    """Read ``RunManifest.json`` (operator-supplied per-well metadata)."""
    try:
        return run_info._read_json(_RUN_MANIFEST_JSON)
    except Exception:
        return None


def parse_run_inputs(run_info: RunInfo) -> RunInputs:
    """Read all three run-level inputs into a single aggregate."""
    return RunInputs(
        run_parameters=parse_run_parameters(run_info),
        panel=parse_panel(run_info),
        run_manifest=parse_run_manifest(run_info),
    )
