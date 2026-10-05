"""Parsers for ``metadata/inputs/cells2stats/`` sidecars.

Includes ``RunStats.json`` (cytometry, demux, target, raw-read summaries),
``AverageNormWellStats.csv`` (per-well CellProfiler-style features), and the
optional OPS-only sidecars ``SpecializedTargetedCounts.json`` (``TargetCounts.json``
before cells2stats 1.4.0-rc8) and
``TargetCellAssignmentManifest.json`` / ``.csv``.
"""

from __future__ import annotations

import io
from dataclasses import dataclass, field
from typing import Any

import pandas as pd

from elembio_spatialdata_analysis.load.backends import RunInfo

_RUNSTATS_PATH = "metadata/inputs/cells2stats/RunStats.json"
_AVG_NORM_PATH = "metadata/inputs/cells2stats/AverageNormWellStats.csv"
# cells2stats 1.4.0-rc8 renamed TargetCounts.json → SpecializedTargetedCounts.json.
# SDOs are durable, so try the rc8+ name first and fall back to the legacy one for
# stores built before the rename.
_TARGET_COUNTS_PATHS: tuple[str, ...] = (
    "metadata/inputs/cells2stats/SpecializedTargetedCounts.json",
    "metadata/inputs/cells2stats/TargetCounts.json",
)
_TARGET_ASSIGN_PATH = "metadata/inputs/cells2stats/TargetCellAssignmentManifest.json"
_VERSIONS_PATH = "metadata/inputs/cells2stats/Versions.json"


# Cells2stats writes ``-999.0`` wherever a percent / ratio metric had a
# zero denominator (e.g. a percentage whose total is 0).
# This is the canonical "no data computed" marker; treat it as NaN here
# so it never flows into aggregations as a real value.
_NAN_SENTINEL = -999.0


def _is_nan_sentinel(v: Any) -> bool:
    """Return True iff *v* is the c2s ``-999.0`` no-data sentinel."""
    try:
        return abs(float(v) - _NAN_SENTINEL) < 1e-6
    except (TypeError, ValueError):
        return False


def _norm(v: Any) -> Any:
    """Sentinel-normalising passthrough: returns ``None`` when *v* is the
    c2s ``-999.0`` no-data sentinel, otherwise *v* unchanged.  Use this
    at the row-emitting boundary so downstream consumers never see the
    sentinel as if it were a real value.
    """
    return None if _is_nan_sentinel(v) else v


# cells2stats ``RunStats.json`` schema version at which the OPS
# ``TargetStats`` layout changed: per-target-site polony stats
# moved out of ``Batches[].TargetSites[].PolonyAssignmentStats`` and into
# ``Batches[].PolonyAssignmentStats.TargetSitePolonyAssignmentStats[]`` (plus
# a new batch/well-level dedup'd aggregate).  ``FileVersion`` < this →
# legacy layout; >= this → new layout.  ``CellAssignmentStats`` is unchanged.
_OPS_TARGETSTATS_V2 = (1, 3, 0)


def _version_tuple(v: Any) -> tuple[int, ...]:
    """Parse a dotted version string into an int tuple (``()`` if unparseable).

    Stops at the first non-numeric component so trailing build suffixes are
    ignored.  An empty/``None`` input yields ``()``, which compares below any
    real version — i.e. an absent ``FileVersion`` falls into the legacy path.
    """
    if not v:
        return ()
    parts: list[int] = []
    for chunk in str(v).split("."):
        lead = ""
        for ch in chunk:
            if not ch.isdigit():
                break
            lead += ch
        if lead == "":
            break
        parts.append(int(lead))
    return tuple(parts)


def _version_at_least(v: Any, minimum: tuple[int, ...]) -> bool:
    """True iff parsed *v* is >= *minimum* (a missing/unparseable *v* is False)."""
    t = _version_tuple(v)
    return bool(t) and t >= minimum


def _percent_perfect_from_mismatch(mismatch: Any) -> float | None:
    """``PercentPerfect = 100 - PercentMismatch`` (the upstream c2s definition).

    Back-fills per-site ``PercentPerfect`` for legacy (< 1.3.0) inputs, which
    carry ``PercentMismatch`` but not ``PercentPerfect``. Returns ``None`` when
    mismatch is missing or the ``-999.0`` no-data sentinel.
    """
    if mismatch is None or _is_nan_sentinel(mismatch):
        return None
    try:
        return 100.0 - float(mismatch)
    except (TypeError, ValueError):
        return None


def _lift_legacy_targetstats_container(container: dict[str, Any]) -> dict[str, Any]:
    """Lift one legacy (< 1.3.0) ``TargetStats`` batch *or* well up to the canonical shape.

    Moves each ``TargetSites[].PolonyAssignmentStats`` into a synthesized
    ``PolonyAssignmentStats.TargetSitePolonyAssignmentStats[]`` (the canonical
    location), leaving ``CellAssignmentStats`` in place, and filling fields the old
    layout never carried (batch/well aggregate ``Percent*`` and per-site
    ``PercentPerfect`` / ``IsControl``) with ``None``.
    """
    site_polony: list[dict[str, Any]] = []
    new_sites: list[dict[str, Any]] = []
    for site in container.get("TargetSites") or []:
        name = site.get("TargetSiteName")
        pol = site.get("PolonyAssignmentStats") or {}
        site_polony.append(
            {
                "TargetSiteName": name,
                "PercentAssigned": pol.get("PercentAssigned"),
                "PercentMismatch": pol.get("PercentMismatch"),
                # Derive PercentPerfect from the upstream definition (100 - PercentMismatch);
                # not emitted before 1.3.0.
                "PercentPerfect": _percent_perfect_from_mismatch(pol.get("PercentMismatch")),
                "IsControl": None,
                "UnassignedSequences": pol.get("UnassignedSequences"),
            }
        )
        new_sites.append(
            {
                "TargetSiteName": name,
                "IsControl": None,
                "CellAssignmentStats": site.get("CellAssignmentStats"),
            }
        )
    lifted = dict(container)
    lifted["PolonyAssignmentStats"] = {
        "PercentAssigned": None,
        "PercentMismatch": None,
        "PercentPerfect": None,
        "TargetSitePolonyAssignmentStats": site_polony,
    }
    lifted["TargetSites"] = new_sites
    return lifted


def _upgrade_targetstats(target_stats: dict[str, Any], file_version: Any) -> dict[str, Any]:
    """Normalise a ``TargetStats`` block UP to the canonical (newest) layout.

    ``FileVersion >= 1.3.0`` is already canonical and passes through unchanged;
    older inputs are lifted (batch + nested wells) so the accessors below only ever
    read one shape. A future schema bump chains a second lift step here.
    """
    if not target_stats or _version_at_least(file_version, _OPS_TARGETSTATS_V2):
        return target_stats
    batches: list[dict[str, Any]] = []
    for batch in target_stats.get("Batches") or []:
        lifted = _lift_legacy_targetstats_container(batch)
        if "Wells" in batch:
            lifted["Wells"] = [_lift_legacy_targetstats_container(w) for w in (batch.get("Wells") or [])]
        batches.append(lifted)
    out = dict(target_stats)
    out["Batches"] = batches
    return out


# cells2stats 1.4.0-rc8 renamed top-level RunStats blocks without bumping
# ``FileVersion`` (still 1.4.0). This loader is rc8+: OPS lives under
# ``SpecializedTargetedStats``; raw reads are synthesized from per-modality
# ``*Stats.ReadStats`` (no legacy ``TargetStats`` / ``RawReadStats`` keys).
_RAW_READ_STATS_PARENTS: tuple[str, ...] = (
    "ThreePrimeUntargetedStats",
    "SpecializedTargetedStats",
    "TargetedTranscriptomeProfilingStats",
)


def _positive_read_count(value: Any) -> float | None:
    """Return a positive read count, or ``None`` for missing / c2s ``-999`` sentinels."""
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if number <= 0:
        return None
    return number


# cells2stats 1.4.0-rc9 renamed per-modality ReadStats counters without bumping
# FileVersion (still 1.4.0). Prefer rc9 names; fall back to rc8 ``PF*`` at parse
# only. Downstream always sees the rc9 keys.
_RAWREAD_KEY_MAP: tuple[tuple[str, str], ...] = (
    ("NumPoloniesBeforeTrimming", "PFReadCount"),
    ("NumPoloniesPerCellBeforeTrimming", "PFReadMeanCountPerCell"),
    ("NumPoloniesPerMM2BeforeTrimming", "PFReadCountPerMM2"),
    ("ExtraCellularRatioBeforeTrimming", "PFReadExtraCellularRatio"),
)


def _rawread_get(data: dict[str, Any], new: str, old: str) -> Any:
    """Read ``new`` (rc9), else ``old`` (rc8 ``PF*``)."""
    value = data.get(new)
    if value is None:
        value = data.get(old)
    return value


def _canonicalize_rawread_batch(batch: dict[str, Any]) -> dict[str, Any]:
    """Copy a ReadStats batch dict onto rc9 field names; drop rc8 ``PF*`` keys."""
    out = dict(batch)
    for new, old in _RAWREAD_KEY_MAP:
        out[new] = _rawread_get(out, new, old)
        out.pop(old, None)
    return out


def _select_target_stats(data: dict[str, Any]) -> dict[str, Any]:
    """Return the OPS stats block under ``SpecializedTargetedStats``.

    Drops a sibling ``ReadStats`` key (rc8 nests per-modality raw-read counters
    next to ``Batches``) so ``target_stats`` stays OPS-shaped for the accessors.
    """
    block = dict(data.get("SpecializedTargetedStats") or {})
    block.pop("ReadStats", None)
    return block


def _synthesize_raw_read_stats(data: dict[str, Any]) -> dict[str, Any]:
    """Merge per-modality ``*Stats.ReadStats`` into one Wells/Batches block.

    rc8 removed run-level ``RawReadStats`` and split it into per-modality
    ``ReadStats``. rc9 renamed the counters from ``PF*`` to
    ``NumPolonies*BeforeTrimming`` / ``ExtraCellularRatioBeforeTrimming``
    without bumping ``FileVersion``. Batches are canonicalized onto the rc9
    names at parse; empty / placeholder modalities (count ``-999``) are skipped.
    Batches are deduped by ``BatchName``. Well and run
    ``NumPoloniesBeforeTrimming`` are sums of kept batch / modality counts
    (not first-parent-wins).
    """
    wells_by_loc: dict[str, dict[str, Any]] = {}
    run_sum = 0.0
    run_any = False
    for parent in _RAW_READ_STATS_PARENTS:
        rs = (data.get(parent) or {}).get("ReadStats") or {}
        run_rc = _positive_read_count(_rawread_get(rs, "NumPoloniesBeforeTrimming", "PFReadCount"))
        if run_rc is not None:
            run_sum += run_rc
            run_any = True
        for well in rs.get("Wells") or []:
            loc = well.get("WellLocation")
            if loc is None:
                continue
            loc_s = str(loc)
            entry = wells_by_loc.setdefault(loc_s, {"WellLocation": loc_s, "batches_by_name": {}})
            for batch in well.get("Batches") or []:
                canon = _canonicalize_rawread_batch(batch)
                if _positive_read_count(canon.get("NumPoloniesBeforeTrimming")) is None:
                    continue
                bname = canon.get("BatchName")
                if bname is None:
                    continue
                entry["batches_by_name"][str(bname)] = canon
    wells_out: list[dict[str, Any]] = []
    for entry in wells_by_loc.values():
        batches = list(entry["batches_by_name"].values())
        if not batches:
            continue
        wells_out.append(
            {
                "WellLocation": entry["WellLocation"],
                "NumPoloniesBeforeTrimming": sum(float(b["NumPoloniesBeforeTrimming"]) for b in batches),
                "Batches": batches,
            }
        )
    if not wells_out:
        return {}
    out: dict[str, Any] = {"Wells": wells_out}
    out["NumPoloniesBeforeTrimming"] = (
        run_sum if run_any else sum(float(w["NumPoloniesBeforeTrimming"]) for w in wells_out)
    )
    return out


def _select_raw_read_stats(data: dict[str, Any]) -> dict[str, Any]:
    """Synthesize raw-read stats from per-modality ``*Stats.ReadStats`` (rc8+)."""
    return _synthesize_raw_read_stats(data)


@dataclass
class Cells2StatsStats:
    """Typed view over ``cells2stats/RunStats.json``."""

    raw: dict[str, Any]
    run_id: str | None = None
    run_name: str | None = None
    flowcell_id: str | None = None
    analysis_version: str | None = None
    file_version: str | None = None
    cyto_stats: dict[str, Any] = field(default_factory=dict)
    demux_stats: dict[str, Any] = field(default_factory=dict)
    target_stats: dict[str, Any] = field(default_factory=dict)
    raw_read_stats: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Cells2StatsStats:
        return cls(
            raw=data,
            run_id=data.get("RunID"),
            run_name=data.get("RunName"),
            flowcell_id=data.get("FlowCellID"),
            analysis_version=data.get("AnalysisVersion"),
            file_version=data.get("FileVersion"),
            cyto_stats=dict(data.get("CytoStats") or {}),
            demux_stats=dict(data.get("DemuxStats") or {}),
            # SpecializedTargetedStats (c2s 1.4.0-rc8+), then normalise site layout
            # UP to the canonical shape.
            target_stats=_upgrade_targetstats(_select_target_stats(data), data.get("FileVersion")),
            raw_read_stats=_select_raw_read_stats(data),
        )

    def _ops_site_polony_entries(self, batch: dict[str, Any]) -> list[tuple[str | None, dict[str, Any]]]:
        """Per-target-site polony dicts for one canonical ``TargetStats`` batch (or well).

        Returns ``(target_site_name, pas)`` pairs from
        ``PolonyAssignmentStats.TargetSitePolonyAssignmentStats[]`` — each ``pas``
        carries the per-site ``PercentAssigned`` / ``PercentMismatch`` /
        ``PercentPerfect`` / ``IsControl`` / ``UnassignedSequences``. Older inputs
        were already lifted into this shape by ``_upgrade_targetstats`` at parse time,
        so no version branching is needed here. (``CellAssignmentStats`` lives on
        ``TargetSites[]`` and is read directly by its own accessors.)
        """
        agg = batch.get("PolonyAssignmentStats") or {}
        return [
            (entry.get("TargetSiteName"), entry)
            for entry in (agg.get("TargetSitePolonyAssignmentStats") or [])
        ]

    def ops_site_polony(self, batch_name: str, target_site_name: str) -> dict[str, Any] | None:
        """Per-site polony dict for one (batch, target-site), or ``None`` if absent.

        Reads the canonical ``PolonyAssignmentStats.TargetSitePolonyAssignmentStats[]``
        location (older inputs were already lifted there by ``_upgrade_targetstats``),
        so callers needn't know which RunStats layout produced the store. Carries
        ``PercentAssigned`` / ``PercentMismatch`` / ``PercentPerfect`` / ``IsControl`` /
        ``UnassignedSequences``.
        """
        for batch in self.target_stats.get("Batches") or []:
            if batch.get("BatchName") != batch_name:
                continue
            for name, pas in self._ops_site_polony_entries(batch):
                if name == target_site_name:
                    return pas
            return None
        return None

    def ops_batch_polony_aggregate(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> list[dict[str, Any]]:
        """Per-batch dedup'd OPS polony aggregate.

        Rows carry ``BatchName``, ``PercentAssigned``, ``PercentMismatch`` from
        ``TargetStats.Batches[].PolonyAssignmentStats`` — the cross-target
        de-duplicated batch-level rate the new (>= 1.3.0) c2s layout introduced.
        Legacy (< 1.3.0) stores have no such block, so the upgrade fills it with
        ``None`` and the values come back as ``None`` here.
        """
        accepted = set(batch_names) if batch_names else None
        rows: list[dict[str, Any]] = []
        for batch in self.target_stats.get("Batches") or []:
            bname = batch.get("BatchName")
            if accepted is not None and bname not in accepted:
                continue
            agg = batch.get("PolonyAssignmentStats") or {}
            if not agg:
                continue
            mismatch = _norm(agg.get("PercentMismatch"))
            perfect = _norm(agg.get("PercentPerfect"))
            if perfect is None:
                perfect = _percent_perfect_from_mismatch(mismatch)
            rows.append(
                {
                    "BatchName": bname,
                    "PercentAssigned": _norm(agg.get("PercentAssigned")),
                    "PercentMismatch": mismatch,
                    "PercentPerfect": perfect,
                }
            )
        return rows

    def ops_polony_assignment_rate_total(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> float | None:
        """Run-level OPS polony assignment rate = total assigned / total polonies.

        The cross-target de-duplicated batch aggregate
        (``PolonyAssignmentStats.PercentAssigned``), averaged across batches when
        there is more than one (rare; proper polony-weighting across batches is a
        follow-up).  This is the correct run headline — unlike
        :meth:`ops_assignment_rate`, an unweighted per-site mean that halves a run
        where one target-site is dead.  Falls back to that per-site mean for legacy
        (< 1.3.0) stores that lack the dedup aggregate.
        """
        agg = self.ops_batch_polony_aggregate(batch_names)
        vals = [r["PercentAssigned"] for r in agg if r.get("PercentAssigned") is not None]
        if vals:
            return sum(vals) / len(vals)
        return self.ops_assignment_rate(batch_names)

    def cyto_per_well(self) -> list[dict[str, Any]]:
        """Per-well cytometry rows (cell count, diameter, batches)."""
        return list(self.cyto_stats.get("Wells") or [])

    def cyto_summary(self) -> dict[str, Any]:
        """Flowcell-level cytometry summary scalars (PercentConfluency, etc.)."""
        return {k: v for k, v in self.cyto_stats.items() if k != "Wells" and not isinstance(v, list)}

    def raw_reads_per_well(self) -> list[dict[str, Any]]:
        """Per-well raw read rows, with nested per-batch ``Batches[]``."""
        return list(self.raw_read_stats.get("Wells") or [])

    def raw_read_aggregates_for_batches(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> dict[str, float] | None:
        """Modality-level pre-trim polony aggregates from synthesized ``Wells[].Batches[]``.

        Sums ``NumPoloniesBeforeTrimming``, ``NumPoloniesPerCellBeforeTrimming``
        and ``NumPoloniesPerMM2BeforeTrimming`` across this modality's batches
        within each well, then averages the per-cell / per-mm² rates across
        wells weighted by the per-well polony count.
        ``ExtraCellularRatioBeforeTrimming`` is aggregated as a count-weighted
        mean across (well, batch).  rc8 ``PF*`` inputs are already lifted onto
        these names by :func:`_synthesize_raw_read_stats`.

        Returns ``None`` when synthesized raw-read stats are empty or no
        matching batch carries the counters.
        """
        wells = self.raw_read_stats.get("Wells") or []
        if not wells:
            return None
        accepted = set(batch_names) if batch_names else None

        total_reads = 0.0
        well_reads_per_cell: list[tuple[float, float]] = []
        well_reads_per_mm2: list[tuple[float, float]] = []
        ec_num = 0.0
        ec_den = 0.0

        for well in wells:
            well_reads = 0.0
            well_per_cell = 0.0
            well_per_mm2 = 0.0
            for batch in well.get("Batches") or []:
                bname = batch.get("BatchName")
                if accepted is not None and bname not in accepted:
                    continue
                rc = batch.get("NumPoloniesBeforeTrimming")
                if rc is None:
                    continue
                rc_f = float(rc)
                total_reads += rc_f
                well_reads += rc_f
                mpc = batch.get("NumPoloniesPerCellBeforeTrimming")
                if mpc is not None:
                    well_per_cell += float(mpc)
                rpm = batch.get("NumPoloniesPerMM2BeforeTrimming")
                if rpm is not None:
                    well_per_mm2 += float(rpm)
                er = batch.get("ExtraCellularRatioBeforeTrimming")
                if er is not None:
                    ec_num += rc_f * float(er)
                    ec_den += rc_f
            if well_reads > 0:
                well_reads_per_cell.append((well_reads, well_per_cell))
                well_reads_per_mm2.append((well_reads, well_per_mm2))

        if total_reads <= 0:
            return None

        def _weighted(pairs: list[tuple[float, float]]) -> float | None:
            num = sum(w * v for w, v in pairs)
            den = sum(w for w, _ in pairs)
            return num / den if den > 0 else None

        return {
            "total_reads": total_reads,
            "reads_per_cell": _weighted(well_reads_per_cell),
            "reads_per_mm2": _weighted(well_reads_per_mm2),
            "extracellular_ratio": ec_num / ec_den if ec_den > 0 else None,
        }

    def demux_per_well(self) -> list[dict[str, Any]]:
        """Per-well demux rows from DemuxStats (may be empty for untargeted runs)."""
        return list(self.demux_stats.get("Wells") or [])

    def target_per_well(self) -> list[dict[str, Any]]:
        """Per-well target rows from TargetStats (may be empty for untargeted runs)."""
        return list(self.target_stats.get("Wells") or [])

    def protein_assignment_rate(self, batch_names: list[str] | tuple[str, ...] | None = None) -> float | None:
        """Polony-weighted Protein assignment rate from ``DemuxStats.Batches``.

        Each Protein batch in ``DemuxStats.Batches`` carries
        ``NumPolonies`` (denominator) and ``PercentAssignedReads`` (% of
        polonies matched to a panel barcode — i.e. survived demux into
        the count table).  This aggregates across the supplied batches
        with polony-count weighting so the run-level value is the
        exact polony-weighted ratio (not a mean-of-percentages).

        Pass ``batch_names`` to restrict to the modality's batches
        when other modalities share ``DemuxStats``.  Pass ``None`` to
        include every batch present.

        Returns ``None`` when DemuxStats is empty or no matching batch
        carries both counters.
        """
        batches = self.demux_stats.get("Batches") or []
        if not batches:
            return None
        accepted = set(batch_names) if batch_names else None
        total_polonies = 0.0
        total_assigned = 0.0
        for batch in batches:
            name = batch.get("BatchName")
            if accepted is not None and name not in accepted:
                continue
            n_polonies = batch.get("NumPolonies")
            pct = batch.get("PercentAssignedReads")
            if n_polonies is None or pct is None or _is_nan_sentinel(pct):
                continue
            total_polonies += float(n_polonies)
            total_assigned += float(n_polonies) * float(pct) / 100.0
        if total_polonies <= 0:
            return None
        return total_assigned / total_polonies * 100.0

    def demux_batch_aggregate(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> list[dict[str, Any]]:
        """Per-batch Protein / Transcript demux aggregate from ``DemuxStats.Batches``.

        Rows carry ``BatchName``, ``PercentAssignedReads``, ``PercentMismatch``
        and ``PercentPerfect`` (``= 100 - PercentMismatch``, the upstream c2s
        definition).  The barcoding (DemuxStats) parallel to
        :meth:`ops_batch_polony_aggregate`, used to surface per-batch assignment
        and perfect-match rates for Protein / Transcript batches.  Returns ``[]``
        when ``DemuxStats`` carries no batches.
        """
        accepted = set(batch_names) if batch_names else None
        rows: list[dict[str, Any]] = []
        for batch in self.demux_stats.get("Batches") or []:
            bname = batch.get("BatchName")
            if accepted is not None and bname not in accepted:
                continue
            mismatch = _norm(batch.get("PercentMismatch"))
            rows.append(
                {
                    "BatchName": bname,
                    "PercentAssignedReads": _norm(batch.get("PercentAssignedReads")),
                    "PercentMismatch": mismatch,
                    "PercentPerfect": _percent_perfect_from_mismatch(mismatch),
                }
            )
        return rows

    def demux_per_well_for_batches(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> list[dict[str, Any]]:
        """Per-(batch, well) rows from ``DemuxStats.Batches[i].Wells``.

        Each row carries ``BatchName``, ``WellLocation``, ``NumPolonies``,
        ``PercentAssignedReads``, ``PercentMismatch``.  Used by the
        Sequencing → Alignment leaves for Protein / Transcript
        Barcoding to render per-well rate bars + per-(well, batch)
        polony-fate stacks.  Returns ``[]`` when ``DemuxStats`` is
        empty or no batch carries Wells.

        Pass ``batch_names`` to restrict to a modality's batches.
        """
        accepted = set(batch_names) if batch_names else None
        rows: list[dict[str, Any]] = []
        for batch in self.demux_stats.get("Batches") or []:
            bname = batch.get("BatchName")
            if accepted is not None and bname not in accepted:
                continue
            for well in batch.get("Wells") or []:
                rows.append(
                    {
                        "BatchName": bname,
                        "WellLocation": well.get("WellLocation"),
                        "NumPolonies": well.get("NumPolonies"),
                        "PercentAssignedReads": _norm(well.get("PercentAssignedReads")),
                        "PercentMismatch": _norm(well.get("PercentMismatch")),
                    }
                )
        return rows

    def demux_per_target_for_batches(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> list[dict[str, Any]]:
        """Per-(batch, target) rows from ``DemuxStats.Batches[i].Targets``.

        Each row carries ``BatchName``, ``TargetName``,
        ``ExpectedSequence``, ``NumPolonies``, ``PercentMismatch``.
        Used by Protein / Transcript Barcoding alignment leaves to
        render a per-target table (which targets dominate the run vs
        which have weak match rates).
        """
        accepted = set(batch_names) if batch_names else None
        rows: list[dict[str, Any]] = []
        for batch in self.demux_stats.get("Batches") or []:
            bname = batch.get("BatchName")
            if accepted is not None and bname not in accepted:
                continue
            for target in batch.get("Targets") or []:
                rows.append(
                    {
                        "BatchName": bname,
                        "TargetName": target.get("TargetName"),
                        "ExpectedSequence": target.get("ExpectedSequence"),
                        "NumPolonies": target.get("NumPolonies"),
                        "PercentMismatch": _norm(target.get("PercentMismatch")),
                    }
                )
        return rows

    def demux_top_unassigned_sequences(
        self,
        batch_names: list[str] | tuple[str, ...] | None = None,
        *,
        n: int = 10,
    ) -> list[dict[str, Any]]:
        """Top-N unassigned-sequence rows from ``DemuxStats.Batches``.

        Pools across matching batches, sums per-sequence counts, and
        returns the ``n`` highest-count entries.  Each row is
        ``{Sequence, Count, BatchNames}`` where ``BatchNames`` is a
        comma-separated list of contributing batches.
        """
        accepted = set(batch_names) if batch_names else None
        pool: dict[str, dict[str, Any]] = {}
        for batch in self.demux_stats.get("Batches") or []:
            bname = batch.get("BatchName")
            if accepted is not None and bname not in accepted:
                continue
            for entry in batch.get("UnassignedSequences") or []:
                seq = entry.get("Sequence")
                count = entry.get("Count") or 0
                if not seq:
                    continue
                existing = pool.setdefault(seq, {"Sequence": seq, "Count": 0, "_batches": set()})
                existing["Count"] += int(count)
                if bname:
                    existing["_batches"].add(str(bname))
        rows = sorted(pool.values(), key=lambda r: r["Count"], reverse=True)[:n]
        for r in rows:
            r["BatchNames"] = ", ".join(sorted(r.pop("_batches")))
        return rows

    def ops_per_site_polony_stats(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> list[dict[str, Any]]:
        """Per-(batch, target-site) rows from ``TargetStats.Batches``.

        Each row carries ``BatchName``, ``TargetSiteName``,
        ``PercentAssigned``, ``PercentMismatch``, ``PercentPerfect``.  Used by
        the OPS alignment leaf for per-site rate bars + the polony-fate stack.

        ``PercentPerfect`` (% of assigned polonies that matched with zero
        mismatches) is read directly when present (cells2stats >= 1.3.0) and
        otherwise derived as ``100 - PercentMismatch`` — the upstream
        definition — so it's available on every RunStats version.
        """
        accepted = set(batch_names) if batch_names else None
        rows: list[dict[str, Any]] = []
        for batch in self.target_stats.get("Batches") or []:
            bname = batch.get("BatchName")
            if accepted is not None and bname not in accepted:
                continue
            for site_name, pas in self._ops_site_polony_entries(batch):
                mismatch = _norm(pas.get("PercentMismatch"))
                perfect = _norm(pas.get("PercentPerfect"))
                if perfect is None:
                    perfect = _percent_perfect_from_mismatch(mismatch)
                rows.append(
                    {
                        "BatchName": bname,
                        "TargetSiteName": site_name,
                        "PercentAssigned": _norm(pas.get("PercentAssigned")),
                        "PercentMismatch": mismatch,
                        "PercentPerfect": perfect,
                    }
                )
        return rows

    def ops_per_site_cell_stats(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> list[dict[str, Any]]:
        """Per-(batch, target-site) rows from
        ``TargetStats.Batches[i].TargetSites[j].CellAssignmentStats``.

        Each row carries the four cell-fate percentages
        (``PercentAssignedPureCells``, ``PercentAssignedMixedCells``,
        ``PercentUnassignedMixedCells``, ``PercentUnassignedLowCountCells``)
        plus identifying fields.  Used by the OPS alignment leaf for
        the per-site cell-assignment-fate stacked bar.
        """
        accepted = set(batch_names) if batch_names else None
        rows: list[dict[str, Any]] = []
        for batch in self.target_stats.get("Batches") or []:
            bname = batch.get("BatchName")
            if accepted is not None and bname not in accepted:
                continue
            for site in batch.get("TargetSites") or []:
                cas = site.get("CellAssignmentStats") or {}
                rows.append(
                    {
                        "BatchName": bname,
                        "TargetSiteName": site.get("TargetSiteName"),
                        "PercentAssignedPureCells": _norm(cas.get("PercentAssignedPureCells")),
                        "PercentAssignedMixedCells": _norm(cas.get("PercentAssignedMixedCells")),
                        "PercentUnassignedMixedCells": _norm(cas.get("PercentUnassignedMixedCells")),
                        "PercentUnassignedLowCountCells": _norm(cas.get("PercentUnassignedLowCountCells")),
                    }
                )
        return rows

    def ops_per_site_cell_assignment(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> list[dict[str, Any]]:
        """Per-(batch, target-site) cell-assignment metrics.

        One row per target site, for **per-target-site** breakdowns instead of
        the unweighted-mean-across-sites the scalar ``ops_*`` aggregators
        return (a mean can hide a single failing site).  Fields per row:

        - ``PercentAssignedCells`` — ``PercentAssignedPureCells +
          PercentAssignedMixedCells`` (the same "assigned" definition
          :meth:`ops_assigned_cells_pct` aggregates).
        - ``MeanUniqueTargetsPerCell``, ``MeanAssignedCountPerCell``,
          ``AssignedCountsPerMM2`` — passthrough from the site's
          ``CellAssignmentStats``.

        Values are sentinel-normalised (``-999.0`` -> ``None``).  Returns an
        empty list when ``TargetStats`` is empty.
        """
        accepted = set(batch_names) if batch_names else None
        rows: list[dict[str, Any]] = []
        for batch in self.target_stats.get("Batches") or []:
            bname = batch.get("BatchName")
            if accepted is not None and bname not in accepted:
                continue
            for site in batch.get("TargetSites") or []:
                cas = site.get("CellAssignmentStats") or {}
                pure = _norm(cas.get("PercentAssignedPureCells"))
                mixed = _norm(cas.get("PercentAssignedMixedCells"))
                pct_assigned = (
                    None if pure is None and mixed is None else float(pure or 0.0) + float(mixed or 0.0)
                )
                rows.append(
                    {
                        "BatchName": bname,
                        "TargetSiteName": site.get("TargetSiteName"),
                        "IsControl": site.get("IsControl"),
                        "PercentAssignedCells": pct_assigned,
                        "MeanUniqueTargetsPerCell": _norm(cas.get("MeanUniqueTargetsPerCell")),
                        "MeanAssignedCountPerCell": _norm(cas.get("MeanAssignedCountPerCell")),
                        "AssignedCountsPerMM2": _norm(cas.get("AssignedCountsPerMM2")),
                    }
                )
        return rows

    def ops_top_unassigned_sequences(
        self,
        batch_names: list[str] | tuple[str, ...] | None = None,
        *,
        n: int = 10,
    ) -> list[dict[str, Any]]:
        """Top-N unassigned-sequence rows from OPS ``PolonyAssignmentStats``.

        Pools across all (batch, target-site) entries, sums per-sequence
        counts, returns the ``n`` highest.  Same shape as
        :meth:`demux_top_unassigned_sequences`.
        """
        accepted = set(batch_names) if batch_names else None
        pool: dict[str, dict[str, Any]] = {}
        for batch in self.target_stats.get("Batches") or []:
            bname = batch.get("BatchName")
            if accepted is not None and bname not in accepted:
                continue
            for _site_name, pas in self._ops_site_polony_entries(batch):
                for entry in pas.get("UnassignedSequences") or []:
                    seq = entry.get("Sequence")
                    count = entry.get("Count") or 0
                    if not seq:
                        continue
                    existing = pool.setdefault(seq, {"Sequence": seq, "Count": 0, "_batches": set()})
                    existing["Count"] += int(count)
                    if bname:
                        existing["_batches"].add(str(bname))
        rows = sorted(pool.values(), key=lambda r: r["Count"], reverse=True)[:n]
        for r in rows:
            r["BatchNames"] = ", ".join(sorted(r.pop("_batches")))
        return rows

    def ops_per_site_unassigned_sequences(
        self,
        batch_names: list[str] | tuple[str, ...] | None = None,
    ) -> list[dict[str, Any]]:
        """Full per-(batch, target-site) unassigned-sequence rows.

        Unlike :meth:`ops_top_unassigned_sequences`, this does NOT pool
        across sites and does NOT top-N truncate.  Every entry that
        upstream ships under
        ``TargetStats.Batches[].TargetSites[].PolonyAssignmentStats.UnassignedSequences[]``
        is returned, tagged with its originating ``BatchName`` and
        ``TargetSiteName``.  Use this for per-site analyses such as
        reverse-complement fractions against each site's own panel
        (the pooled view dilutes the per-site signal — a site with
        0% assigned can be hidden by other sites' contaminants).
        """
        accepted = set(batch_names) if batch_names else None
        rows: list[dict[str, Any]] = []
        for batch in self.target_stats.get("Batches") or []:
            bname = batch.get("BatchName")
            if accepted is not None and bname not in accepted:
                continue
            for sname, pas in self._ops_site_polony_entries(batch):
                for entry in pas.get("UnassignedSequences") or []:
                    seq = entry.get("Sequence")
                    count = entry.get("Count") or 0
                    if not seq:
                        continue
                    rows.append(
                        {
                            "BatchName": bname,
                            "TargetSiteName": sname,
                            "Sequence": str(seq),
                            "Count": int(count),
                        }
                    )
        return rows

    def ops_assigned_cells_pct(self, batch_names: list[str] | tuple[str, ...] | None = None) -> float | None:
        """% of OPS cells confidently assigned to a target — pre-computed
        upstream, no AnnData re-read.

        Numerator: ``PercentAssignedPureCells + PercentAssignedMixedCells``
        from each ``TargetStats.Batches[i].TargetSites[j].CellAssignmentStats``
        block.  ``Pure`` are cells with a clear dominant target call;
        ``Mixed`` are cells assigned despite some target ambiguity.
        Their sum = "cells the upstream pipeline considers assigned".
        ``PercentUnassignedMixedCells + PercentUnassignedLowCountCells``
        is the implicit complement.

        Aggregates across target sites with an unweighted mean (the
        manifest doesn't carry a per-site cell denominator; for the
        common N=1 site case mean = single value exactly).  Pass
        ``batch_names`` to restrict to a subset of batches.

        Returns ``None`` when TargetStats is empty or no site carries
        either ``PercentAssigned*Cells`` field.
        """
        batches = self.target_stats.get("Batches") or []
        if not batches:
            return None
        accepted = set(batch_names) if batch_names else None
        rates: list[float] = []
        for batch in batches:
            name = batch.get("BatchName")
            if accepted is not None and name not in accepted:
                continue
            for site in batch.get("TargetSites") or []:
                cas = site.get("CellAssignmentStats") or {}
                pure = cas.get("PercentAssignedPureCells")
                mixed = cas.get("PercentAssignedMixedCells")
                if _is_nan_sentinel(pure):
                    pure = None
                if _is_nan_sentinel(mixed):
                    mixed = None
                if pure is None and mixed is None:
                    continue
                rates.append(float(pure or 0.0) + float(mixed or 0.0))
        if not rates:
            return None
        return sum(rates) / len(rates)

    def ops_cell_assignment_fates(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> dict[str, float] | None:
        """Per-fate cell breakdown from OPS ``CellAssignmentStats`` blocks.

        Decomposes the headline ``ops_assigned_cells_pct`` value into
        the four mutually-exclusive fates the upstream pipeline ships:

        - ``pure_pct``                — ``PercentAssignedPureCells``,
          cells with a clear single-target call.
        - ``mixed_pct``               — ``PercentAssignedMixedCells``,
          assigned despite secondary target signal (borderline).
        - ``unassigned_mixed_pct``    — ``PercentUnassignedMixedCells``,
          mixed signal but no confident call — doublet candidates.
        - ``unassigned_lowcount_pct`` — ``PercentUnassignedLowCountCells``,
          below detection threshold — undersequenced cells.

        Aggregates across target sites with an unweighted mean (same
        convention as :meth:`ops_assigned_cells_pct`).  Returns ``None``
        when ``TargetStats`` is empty or no qualifying site carries any
        of the four fields.
        """
        batches = self.target_stats.get("Batches") or []
        if not batches:
            return None
        accepted = set(batch_names) if batch_names else None
        keys: tuple[tuple[str, str], ...] = (
            ("pure_pct", "PercentAssignedPureCells"),
            ("mixed_pct", "PercentAssignedMixedCells"),
            ("unassigned_mixed_pct", "PercentUnassignedMixedCells"),
            ("unassigned_lowcount_pct", "PercentUnassignedLowCountCells"),
        )
        accum: dict[str, list[float]] = {out: [] for out, _ in keys}
        for batch in batches:
            name = batch.get("BatchName")
            if accepted is not None and name not in accepted:
                continue
            for site in batch.get("TargetSites") or []:
                cas = site.get("CellAssignmentStats") or {}
                if all(cas.get(src) is None or _is_nan_sentinel(cas.get(src)) for _, src in keys):
                    continue
                for out, src in keys:
                    val = cas.get(src)
                    if val is None or _is_nan_sentinel(val):
                        continue
                    accum[out].append(float(val))
        if not any(accum.values()):
            return None
        return {out: (sum(vals) / len(vals)) if vals else 0.0 for out, vals in accum.items()}

    def ops_mean_unique_targets_per_cell(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> float | None:
        """Mean # distinct OPS targets detected per cell — pre-computed
        upstream as ``CellAssignmentStats.MeanUniqueTargetsPerCell``.

        Together with :meth:`ops_mean_assigned_counts_per_cell`, this
        disambiguates the "mixed-signal" cells in
        :meth:`ops_cell_assignment_fates`: a high counts-per-cell with
        a similarly high targets-per-cell means the dominant call is
        diluted (true doublet/multi-cell suspicion); a low
        counts-per-cell with comparable targets means polonies are
        simply scattered across many targets and no target ever
        accumulates enough signal to clear the confidence threshold.

        Aggregates unweighted-mean across target sites (same convention
        as the other ``ops_*`` aggregators).  Returns ``None`` when
        ``TargetStats`` is empty.
        """
        return self._ops_cas_aggregate("MeanUniqueTargetsPerCell", batch_names)

    def ops_mean_assigned_counts_per_cell(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> float | None:
        """Mean # assigned polonies per cell — pre-computed upstream as
        ``CellAssignmentStats.MeanAssignedCountPerCell``.  See
        :meth:`ops_mean_unique_targets_per_cell` for the diagnostic
        framing the two metrics together provide.
        """
        return self._ops_cas_aggregate("MeanAssignedCountPerCell", batch_names)

    def ops_assigned_counts_per_mm2(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> float | None:
        """Assigned counts per mm² of cell area — pre-computed upstream as
        ``CellAssignmentStats.AssignedCountsPerMM2``.  OPS-side equivalent
        of the modality-AnnData-derived ``Counts / mm²`` that
        :func:`_counts_per_mm2` provides for the other modalities; this
        accessor is the right source for the OPS row because the OPS
        modality table is per-target-site rather than aggregated.
        """
        return self._ops_cas_aggregate("AssignedCountsPerMM2", batch_names)

    def _ops_cas_aggregate(
        self,
        field: str,
        batch_names: list[str] | tuple[str, ...] | None,
    ) -> float | None:
        """Shared unweighted-mean aggregator over
        ``TargetStats.Batches[*].TargetSites[*].CellAssignmentStats[field]``.
        Returns ``None`` when no qualifying site carries the field.
        """
        batches = self.target_stats.get("Batches") or []
        if not batches:
            return None
        accepted = set(batch_names) if batch_names else None
        vals: list[float] = []
        for batch in batches:
            if accepted is not None and batch.get("BatchName") not in accepted:
                continue
            for site in batch.get("TargetSites") or []:
                cas = site.get("CellAssignmentStats") or {}
                v = cas.get(field)
                if v is None or _is_nan_sentinel(v):
                    continue
                try:
                    vals.append(float(v))
                except (TypeError, ValueError):
                    continue
        if not vals:
            return None
        return sum(vals) / len(vals)

    def ops_assignment_rate(self, batch_names: list[str] | tuple[str, ...] | None = None) -> float | None:
        """Mean polony-level OPS assignment rate across target sites.

        OPS records assignment rate in ``TargetStats.Batches[i].
        PolonyAssignmentStats.TargetSitePolonyAssignmentStats[j].PercentAssigned``
        — % of OPS polonies in this target site that mapped to a panel
        target (and thus contributed a count to ``adata.X``).  This
        aggregates across all target sites under matching batches with
        an unweighted mean: the manifest doesn't carry a per-site
        polony denominator, so weighting isn't directly possible — but
        the common case is N=1 site where mean == single value, so the
        approximation is exact for typical runs and reasonable for
        the rare multi-site case.

        Pass ``batch_names`` to restrict to the modality's batches.
        Returns ``None`` when TargetStats is empty or no qualifying
        ``PercentAssigned`` value is present.
        """
        batches = self.target_stats.get("Batches") or []
        if not batches:
            return None
        accepted = set(batch_names) if batch_names else None
        rates: list[float] = []
        for batch in batches:
            name = batch.get("BatchName")
            if accepted is not None and name not in accepted:
                continue
            for _site_name, pas in self._ops_site_polony_entries(batch):
                pct = pas.get("PercentAssigned")
                if pct is None or _is_nan_sentinel(pct):
                    continue
                rates.append(float(pct))
        if not rates:
            return None
        return sum(rates) / len(rates)


def parse_cells2stats_runstats(run_info: RunInfo) -> Cells2StatsStats | None:
    """Read and parse ``cells2stats/RunStats.json``."""
    try:
        data = run_info._read_json(_RUNSTATS_PATH)
    except Exception:
        return None
    if not data:
        return None
    return Cells2StatsStats.from_dict(data)


def parse_average_norm_well_stats(run_info: RunInfo) -> pd.DataFrame | None:
    """Read ``AverageNormWellStats.csv`` as a DataFrame indexed by Well.

    Returns ``None`` when the file is absent or unreadable.  Columns include
    per-channel intensity means (e.g. ``Cell-Membrane.CP01``) plus
    morphology features (``AreaShape_*``, ``Granularity_*``, etc.) with
    ``.std`` siblings for within-well dispersion.
    """
    raw_text = _read_text(run_info, _AVG_NORM_PATH)
    if not raw_text:
        return None
    try:
        df = pd.read_csv(io.StringIO(raw_text))
    except Exception:
        return None
    if df.empty:
        return None
    if "Well" in df.columns:
        df = df.set_index("Well")
    return df


def parse_target_counts(run_info: RunInfo) -> dict[str, Any] | None:
    """Read OPS ``cells2stats/SpecializedTargetedCounts.json``.

    Falls back to the pre-rc8 ``TargetCounts.json`` name. Structure:
    ``Batches[].Wells[].Cells[].Targets[]``. Returns the raw dict or ``None``
    when neither file is present.
    """
    for path in _TARGET_COUNTS_PATHS:
        try:
            data = run_info._read_json(path)
        except Exception:
            continue
        if data:
            return data
    return None


def parse_target_cell_assignment(run_info: RunInfo) -> dict[str, Any] | None:
    """Read OPS ``cells2stats/TargetCellAssignmentManifest.json`` (panel definition)."""
    try:
        data = run_info._read_json(_TARGET_ASSIGN_PATH)
    except Exception:
        return None
    return data or None


def parse_cells2stats_versions(run_info: RunInfo) -> dict[str, Any] | None:
    """Read ``cells2stats/Versions.json`` — upstream tool / file versions.

    Schema (observed across example runs):

    .. code-block:: json

       {
           "FileVersion": "1.0.0",
           "FileVersions": {"AverageNormWellStats.csv": "1.0.1", ...},
           "ProgramVersions": {"CellProfiler": "4.2.6", ...}
       }

    Returns ``None`` when the file is absent so callers can branch
    cleanly.
    """
    try:
        data = run_info._read_json(_VERSIONS_PATH)
    except Exception:
        return None
    return data or None


def panel_target_sequences(manifest: dict[str, Any] | None) -> dict[str, str]:
    """Return ``{ExpectedSequence: TargetName}`` from a panel manifest.

    Pools across every ``TargetSites[].Targets[]`` entry in the
    ``TargetCellAssignmentManifest.json`` payload returned by
    :func:`parse_target_cell_assignment`.  Returns an empty dict when
    the manifest is ``None`` or carries no sequence data — callers
    can branch cleanly on truthiness.

    When two targets in the manifest share a sequence (rare; usually
    a panel-design bug), the later entry wins — duplicates are not
    flagged here.
    """
    if not manifest:
        return {}
    out: dict[str, str] = {}
    for site in manifest.get("TargetSites") or []:
        for target in site.get("Targets") or []:
            seq = target.get("Sequence")
            name = target.get("Name")
            if seq and name:
                out[str(seq).upper()] = str(name)
    return out


def panel_target_sequences_by_site(
    manifest: dict[str, Any] | None,
) -> dict[str, dict[str, str]]:
    """Return ``{site_name: {ExpectedSequence: TargetName}}`` from a panel manifest.

    Per-site variant of :func:`panel_target_sequences` — needed when
    the caller wants to look up reverse-complement matches against a
    specific site's panel, scoped to one site at a time so one site's
    polony matching another site's target isn't misattributed.
    """
    if not manifest:
        return {}
    out: dict[str, dict[str, str]] = {}
    for site in manifest.get("TargetSites") or []:
        sname = site.get("Name") or site.get("TargetSiteName")
        if not sname:
            continue
        site_map: dict[str, str] = {}
        for target in site.get("Targets") or []:
            seq = target.get("Sequence")
            name = target.get("Name")
            if seq and name:
                site_map[str(seq).upper()] = str(name)
        if site_map:
            out[str(sname)] = site_map
    return out


def _read_text(run_info: RunInfo, path: str) -> str | None:
    """Helper: each backend has ``_read_json`` but no generic text helper.

    For CSVs we read via the backend's underlying open paths.  Local + zip
    backends expose ``_root`` / ``_zip_path``; S3 uses fsspec.  We fall back
    to reading the file directly when possible to avoid duplicating logic.
    """
    root = getattr(run_info, "_root", None)
    if root is not None:
        fp = root / path
        if fp.is_file():
            return fp.read_text()
        return None

    zip_path = getattr(run_info, "_zip_path", None)
    if zip_path is not None:
        import zipfile

        zarr_prefix = getattr(run_info, "_zarr_prefix", "")
        full = f"{zarr_prefix}/{path}".lstrip("/")
        try:
            with zipfile.ZipFile(zip_path, "r") as zf:
                if full in zf.namelist():
                    return zf.read(full).decode()
        except Exception:
            return None
        return None

    fs = getattr(run_info, "_fs", None)
    prefix = getattr(run_info, "_prefix", None)
    if fs is not None and prefix:
        key = f"{prefix}/{path}"
        try:
            if fs.exists(key):
                with fs.open(key, "r") as f:
                    return f.read()
        except Exception:
            return None

    return None
