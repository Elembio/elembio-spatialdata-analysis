"""Parsers for ``metadata/inputs/fastq2parquet_3prime/`` sidecars.

Includes ``RunStats.json`` (per-sample alignment / assignment / cell
metrics) and a thin summary view over ``multiqc/multiqc_data.json`` that
exposes the canonical "general stats" rows without expanding individual
modules.
"""

from __future__ import annotations

import contextlib
from dataclasses import dataclass, field
from typing import Any

from elembio_spatialdata_analysis.load.backends import RunInfo

_RUNSTATS_PATH = "metadata/inputs/fastq2parquet_3prime/RunStats.json"
_MULTIQC_PATH = "metadata/inputs/fastq2parquet_3prime/multiqc/multiqc_data.json"

_PROXIMITY_GAP_BIN_ORDER: tuple[str, ...] = (
    "0",
    "1_50",
    "51_100",
    "101_250",
    "251_500",
    "501_1k",
    "1k_5k",
    "5k_plus",
    "no_chrom",
)
_PROXIMITY_GAP_LABELS: dict[str, str] = {
    "0": "0",
    "1_50": "1–50",
    "51_100": "51–100",
    "101_250": "101–250",
    "251_500": "251–500",
    "501_1k": "501–1k",
    "1k_5k": "1k–5k",
    "5k_plus": "5k+",
    "no_chrom": "no chrom",
}

# Newer RunStats renamed the read-filter counters (``*PF`` → ``*FilterPass`` /
# ``*RemovedByFilter``). We read the new key first, falling back to the old ``PF`` key so
# older RunStats.json still parse, and expose the value under the stable old key name.
_PF_KEY_RENAMES: dict[str, str] = {
    "TotalReadsCountablePF": "TotalReadsCountableFilterPass",
    "TotalReadsCountablePFInCell": "TotalReadsCountableFilterPassInCell",
    "TotalReadsLostToPF": "TotalReadsRemovedByFilter",
    "TotalReadsLostToPFInCell": "TotalReadsRemovedByFilterInCell",
    "TotalReadsPF": "TotalReadsFilterPass",
    "TotalReadsPFInCell": "TotalReadsFilterPassInCell",
    "LostToPF_low_quality": "RemovedByFilter_low_quality",
    "LostToPF_high_softclip": "RemovedByFilter_high_softclip",
    "TotalReadsUnassignedUnmappedPF": "TotalReadsUnassignedUnmappedFilterPass",
    "TotalReadsUnassignedNoFeaturePF": "TotalReadsUnassignedNoFeatureFilterPass",
    # Exact no-feature band split: these never had a PF-era name; the entry exists so
    # _rf resolves the FilterPass suffix like every other read-fate counter.
    "TotalReadsUnassignedUnselectedPF": "TotalReadsUnassignedUnselectedFilterPass",
    "TotalReadsUnassignedNoOverlapPF": "TotalReadsUnassignedNoOverlapFilterPass",
    "TotalReadsAssignedOutOfCellPF": "TotalReadsAssignedOutOfCellFilterPass",
    "TotalReadsAssignedUniqueInCellPF": "TotalReadsAssignedUniqueInCellFilterPass",
    "TotalReadsAssignedMultiSimpleInCellPF": "TotalReadsAssignedMultiSimpleInCellFilterPass",
    "TotalReadsAssignedMultiComplexInCellPF": "TotalReadsAssignedMultiComplexInCellFilterPass",
    "TotalReadsMappedPerfectPF": "TotalReadsMappedPerfectFilterPass",
    "TotalReadsMappedSoftclippedPF": "TotalReadsMappedSoftclippedFilterPass",
    "TotalReadsMappedWithEditsPF": "TotalReadsMappedWithEditsFilterPass",
    "MeanReadBaseQualityPF": "MeanReadBaseQualityFilterPass",
    "MeanReadBaseQualityPFfail": "MeanReadBaseQualityFilterFail",
    "MeanReadErrorRatePF": "MeanReadErrorRateFilterPass",
    "MeanReadErrorRatePFfail": "MeanReadErrorRateFilterFail",
}


def _rf(section: dict[str, Any], stable_key: str) -> Any:
    """Read a read-filter counter by its new RunStats key, falling back to the old PF key
    for older RunStats. ``stable_key`` is the old name, which the value is exposed under."""
    new_key = _PF_KEY_RENAMES.get(stable_key, stable_key)
    val = section.get(new_key)
    return val if val is not None else section.get(stable_key)


# ── Feature-tier vocabulary crosswalk ────────────────────────────────────────────────
# Newer RunStats (exon-union + 3'UTR-priority model) renamed the 3' feature tiers. Tier
# tokens are enumerated out of the RunStats keys — ``CountedReads_tier_<tier>``,
# ``AssignedReads_tier_<tier>_<mtype>``, ``SoftclipPct_tier_<tier>``,
# ``ErrorRateMean_tier_<tier>``, ``CountableReads_tier_<tier>``,
# ``UnselectedTierPct_<tier>`` — plus the nested ``TierQuality`` dict. Those keys are
# parsed by PREFIX-STRIP (``key[len(prefix):]``) + known-mtype SUFFIX-match, never a naive
# ``split("_")``, so the multi-underscore current names (``three_prime_utr``,
# ``three_prime_utr_extension``, ``chrM_heavy``, ``chrM_light``) round-trip intact.
#
# This crosswalk normalizes old tier tokens to the new vocabulary so older RunStats.json
# expose the current tier names and cross-run aggregation merges aliases correctly. It mirrors ``_PF_KEY_RENAMES``: a
# read-time fallback, never a hard error. Several old tiers MERGE into one new tier (their
# per-sample counts are summed); the retired ``synthetic_3utr`` tier folds into the
# ``uncategorized`` catch-all so per-tier counts still sum to the counted total.
_TIER_RENAMES: dict[str, str] = {
    "canonical_3utr": "three_prime_utr",
    "alt_isoform_3utr": "three_prime_utr",
    "cds_edge_3utr": "three_prime_utr_extension",
    "polya_extension": "three_prime_utr_extension",
    "exon_body": "exon",
    "chrM_umbrella_heavy": "chrM_heavy",
    "chrM_umbrella_light": "chrM_light",
    "synthetic_3utr": "uncategorized",  # retired tier → catch-all (preserves totals)
    # NOTE: intronic_sense (an older name for `intron`) is intentionally NOT
    # crosswalked — no current run emits it, and the breakdown tests use it as a
    # retained example tier. Left as passthrough.
}


def canonical_tier(tier: str) -> str:
    """Map an old feature-tier token to its current exon-union name.

    Unknown / already-current tokens pass through unchanged (fallback, never raises) —
    so new-vocabulary RunStats are untouched and any not-yet-crosswalked token still
    parses. See :data:`_TIER_RENAMES` for the mapping."""
    return _TIER_RENAMES.get(tier, tier)


@dataclass
class Fastq2Parquet3pStats:
    """Typed view over ``fastq2parquet_3prime/RunStats.json``."""

    raw: dict[str, Any]
    run_id: str | None = None
    run_name: str | None = None
    flowcell_id: str | None = None
    analysis_version: str | None = None
    run_summary: dict[str, Any] = field(default_factory=dict)
    samples: list[dict[str, Any]] = field(default_factory=list)
    # Resolved workflow parameters from RunStats.json ``Params`` (read-filter
    # K / Q / soft-clip, count policy, features, biotypes, slop); empty on
    # older runs that don't record them.
    params: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> Fastq2Parquet3pStats:
        return cls(
            raw=data,
            run_id=data.get("RunID"),
            run_name=data.get("RunName"),
            flowcell_id=data.get("FlowCellID"),
            analysis_version=data.get("AnalysisVersion"),
            run_summary=dict(data.get("RunSummary") or {}),
            samples=list(data.get("Samples") or []),
            params=dict(data.get("Params") or {}),
        )

    def per_sample_summary(self) -> list[dict[str, Any]]:
        """Flat list of one row per (well, batch) sample with key metrics."""
        rows: list[dict[str, Any]] = []
        for sample in self.samples:
            alignment = sample.get("Alignment") or {}
            on_target = sample.get("OnTarget") or {}
            cell_metrics = sample.get("CellMetrics") or {}
            in_cell = cell_metrics.get("InCell") or {}
            rows.append(
                {
                    "SampleID": sample.get("SampleID"),
                    "WellLocation": sample.get("WellLocation"),
                    "WellLabel": sample.get("WellLabel"),
                    "BatchName": sample.get("BatchName"),
                    "AssayType": sample.get("AssayType"),
                    "TotalReads": (sample.get("Reads") or {}).get("TotalReads"),
                    "AlignmentRate": alignment.get("AlignmentRate"),
                    "OnTargetRate": on_target.get("OnTargetRate"),
                    "GenesDetected": (sample.get("Genes") or {}).get("GenesDetected"),
                    "CellsDetected": cell_metrics.get("CellsDetected"),
                    "InCellReadCount": in_cell.get("ReadCount"),
                    "MeanReadCountPerCell": in_cell.get("MeanReadCountPerCell"),
                    "ReadCountPerMM2": in_cell.get("ReadCountPerMM2"),
                }
            )
        return rows

    def per_sample_alignment(self) -> list[dict[str, Any]]:
        """Flat list of one row per sample with read-fate + alignment fields.

        Returns the per-(well, batch) alignment and read-fate numbers.
        ``RunSummary``-derived rates aren't covered here (they're scalars on
        ``run_summary``). ``per_sample_summary`` is the shallow overview;
        this is the deeper accessor.
        """
        rows: list[dict[str, Any]] = []
        for sample in self.samples:
            alignment = sample.get("Alignment") or {}
            on_target = sample.get("OnTarget") or {}
            mapping = sample.get("Mapping") or {}
            assignment = sample.get("Assignment") or {}
            filtered = sample.get("FilteredPolonies") or {}
            strand = sample.get("StrandSpecificity") or {}
            # ReadQC section holds the read-quality counts + the PF-sliced per-read
            # means (TotalReadsMapped* / MeanRead* live there, NOT in Assignment).
            readqc = sample.get("ReadQC") or {}
            rows.append(
                {
                    "SampleID": sample.get("SampleID"),
                    "WellLocation": sample.get("WellLocation"),
                    "WellLabel": sample.get("WellLabel"),
                    "BatchName": sample.get("BatchName"),
                    "AssayType": sample.get("AssayType"),
                    # Rates (Alignment / OnTarget)
                    "PctReadsMapped": alignment.get("PctReadsMapped"),
                    "UniquelyMappedReadsPct": alignment.get("UniquelyMappedReadsPct"),
                    "MismatchRatePerBasePct": alignment.get("MismatchRatePerBasePct"),
                    # On-target is SAF-derived and MAPPED-denominated in newer RunStats:
                    # the fraction of mapped reads that resolved a winner in the
                    # SELECTED collapse. Older RunStats used Picard CollectHsMetrics
                    # base-level coverage on a GTF 3'UTR BED — same key, different
                    # meaning, so a mixed-version cross-run comparison is not valid.
                    # Older RunStats call PctOffTargetMapped ``PctOffBait``; both are read.
                    "PctOnTarget": on_target.get("PctOnTarget"),
                    "PctOffTargetMapped": on_target.get("PctOffTargetMapped", on_target.get("PctOffBait")),
                    # Read-fate counters.  ``TotalFilteredPolonies`` is the
                    # canonical filtered-OUT count (anti-join of instrument
                    # rawreads against pipeline ``read_id``s).  Some older runs
                    # report ``TotalFilteredPolonies == TotalInstrumentReads``
                    # (i.e. 100%) due to an upstream bug.
                    # ``TotalDistinctPipelineReads`` lets consumers sanity-check
                    # the upstream uniqueness contract (one row per ``read_id``
                    # in the pipeline parquet).
                    "TotalInstrumentReads": filtered.get("TotalInstrumentReads"),
                    "TotalPipelineReads": filtered.get("TotalPipelineReads"),
                    "TotalDistinctPipelineReads": filtered.get("TotalDistinctPipelineReads"),
                    "TotalFilteredPolonies": filtered.get("TotalFilteredPolonies"),
                    "TotalReadsUnmapped": mapping.get("TotalReadsUnmapped"),
                    "TotalReadsAssigned": assignment.get("TotalReadsAssigned"),
                    "TotalReadsAssignedInCell": assignment.get("TotalReadsAssignedInCell"),
                    "TotalReadsAssignedOutOfCell": assignment.get("TotalReadsAssignedOutOfCell"),
                    # Assigned sub-categories — pipeline default counts
                    # ``Unique + MultiSimple`` into the gene matrix and drops
                    # ``MultiComplex`` (cross-gene multi-mappers).  The exact
                    # ``IsCounted`` formula is configurable per workflow run.
                    "TotalReadsAssignedUniqueInCell": assignment.get("TotalReadsAssignedUniqueInCell"),
                    "TotalReadsAssignedMultiSimpleInCell": assignment.get(
                        "TotalReadsAssignedMultiSimpleInCell"
                    ),
                    "TotalReadsAssignedMultiComplexInCell": assignment.get(
                        "TotalReadsAssignedMultiComplexInCell"
                    ),
                    "TotalReadsIsCounted": assignment.get("TotalReadsIsCounted"),
                    "TotalReadsIsCountedInCell": assignment.get("TotalReadsIsCountedInCell"),
                    "TotalReadsUnassignedNoFeature": assignment.get("TotalReadsUnassignedNoFeature"),
                    "TotalReadsUnassignedNoFeatureInCell": assignment.get(
                        "TotalReadsUnassignedNoFeatureInCell"
                    ),
                    "TotalReadsUnassignedUnmapped": assignment.get("TotalReadsUnassignedUnmapped"),
                    # EXACT split of the no-feature band (newer RunStats, annotate-then-
                    # select). Unselected = hit a real gene outside the run's selection;
                    # NoOverlap = no sense-strand feature anywhere. Counts, not percentages,
                    # and they partition TotalReadsUnassignedNoFeature exactly.
                    #
                    # NoOverlap is NOT "intergenic": the counting pass is sense-stranded, so
                    # it still holds reads antisense to a gene. It replaces the older
                    # UnassignedIntergenicPct, which was measured on both strands on a
                    # 10% downsample and therefore is NOT comparable to it directly — the
                    # sampled unselected share must be multiplied by its own sense fraction
                    # first.
                    #
                    # Absent on older SDOs; the older percentages are read alongside so
                    # those runs still parse.
                    "TotalReadsUnassignedUnselected": assignment.get("TotalReadsUnassignedUnselected"),
                    "TotalReadsUnassignedUnselectedInCell": assignment.get(
                        "TotalReadsUnassignedUnselectedInCell"
                    ),
                    "TotalReadsUnassignedNoOverlap": assignment.get("TotalReadsUnassignedNoOverlap"),
                    "TotalReadsUnassignedNoOverlapInCell": assignment.get(
                        "TotalReadsUnassignedNoOverlapInCell"
                    ),
                    "UnassignedIntergenicPct": assignment.get("UnassignedIntergenicPct"),
                    "UnassignedUnselectedPct": assignment.get("UnassignedUnselectedPct"),
                    "ProximityGap_pct_le_100": (sample.get("Reads") or {}).get("ProximityGap_pct_le_100"),
                    "ProximityGap_pct_le_250": (sample.get("Reads") or {}).get("ProximityGap_pct_le_250"),
                    # Read-filter dimension. IsCountable = engine eligibility (pre-filter,
                    # == IsCounted when the filter did not run); CountableFilterPass =
                    # IsCountable AND read_filter_pass; RemovedByFilter = IsCountable AND NOT
                    # read_filter_pass; FilterPass = the filter-passing total. Present only
                    # when the read filter ran. Read via _rf (new key, old-key fallback);
                    # exposed under the stable internal names below.
                    "TotalReadsIsCountable": assignment.get("TotalReadsIsCountable"),
                    "TotalReadsIsCountableInCell": assignment.get("TotalReadsIsCountableInCell"),
                    "TotalReadsCountablePF": _rf(assignment, "TotalReadsCountablePF"),
                    "TotalReadsCountablePFInCell": _rf(assignment, "TotalReadsCountablePFInCell"),
                    "TotalReadsLostToPF": _rf(assignment, "TotalReadsLostToPF"),
                    "TotalReadsLostToPFInCell": _rf(assignment, "TotalReadsLostToPFInCell"),
                    "TotalReadsPF": _rf(assignment, "TotalReadsPF"),
                    "TotalReadsPFInCell": _rf(assignment, "TotalReadsPFInCell"),
                    # Removal reason partition (sums to TotalReadsLostToPF/RemovedByFilter)
                    "LostToPF_low_quality": _rf(assignment, "LostToPF_low_quality"),
                    "LostToPF_high_softclip": _rf(assignment, "LostToPF_high_softclip"),
                    # Filter-pass slices of each read-fate + quality bucket. {fate}PF = the
                    # filter-passing subset of {fate}; filter-fail is {fate} - {fate}PF.
                    # Present only when the read filter ran.
                    "TotalReadsUnassignedUnmappedPF": _rf(assignment, "TotalReadsUnassignedUnmappedPF"),
                    "TotalReadsUnassignedNoFeaturePF": _rf(assignment, "TotalReadsUnassignedNoFeaturePF"),
                    "TotalReadsUnassignedUnselectedPF": _rf(assignment, "TotalReadsUnassignedUnselectedPF"),
                    "TotalReadsUnassignedNoOverlapPF": _rf(assignment, "TotalReadsUnassignedNoOverlapPF"),
                    "TotalReadsAssignedOutOfCellPF": _rf(assignment, "TotalReadsAssignedOutOfCellPF"),
                    "TotalReadsAssignedUniqueInCellPF": _rf(assignment, "TotalReadsAssignedUniqueInCellPF"),
                    "TotalReadsAssignedMultiSimpleInCellPF": _rf(
                        assignment, "TotalReadsAssignedMultiSimpleInCellPF"
                    ),
                    "TotalReadsAssignedMultiComplexInCellPF": _rf(
                        assignment, "TotalReadsAssignedMultiComplexInCellPF"
                    ),
                    # Quality counts + per-read means live in the ReadQC section.
                    "TotalReadsMappedPerfectPF": _rf(readqc, "TotalReadsMappedPerfectPF"),
                    "TotalReadsMappedSoftclippedPF": _rf(readqc, "TotalReadsMappedSoftclippedPF"),
                    "TotalReadsMappedWithEditsPF": _rf(readqc, "TotalReadsMappedWithEditsPF"),
                    # Filter-sliced per-read means (ALL / filter-pass / filter-fail)
                    "MeanReadBaseQuality": readqc.get("MeanReadBaseQuality"),
                    "MeanReadBaseQualityPF": _rf(readqc, "MeanReadBaseQualityPF"),
                    "MeanReadBaseQualityPFfail": _rf(readqc, "MeanReadBaseQualityPFfail"),
                    "MeanReadErrorRate": readqc.get("MeanReadErrorRate"),
                    "MeanReadErrorRatePF": _rf(readqc, "MeanReadErrorRatePF"),
                    "MeanReadErrorRatePFfail": _rf(readqc, "MeanReadErrorRatePFfail"),
                    # Strand specificity (Picard CollectRnaSeqMetrics) — sense-vs-antisense
                    # to the gene model; PctCorrectStrandReads is the library strand-spec.
                    "PctCorrectStrandReads": strand.get("PctCorrectStrandReads"),
                    "CorrectStrandReads": strand.get("CorrectStrandReads"),
                    "IncorrectStrandReads": strand.get("IncorrectStrandReads"),
                    # Not-counted taxonomy (NotCounted_<reason>) — nested for the breakdown table.
                    "not_counted_breakdown": {
                        k: v for k, v in assignment.items() if k.startswith("NotCounted_")
                    },
                }
            )
        return rows

    def assignment_rate(self, batch_names: list[str] | tuple[str, ...] | None = None) -> float | None:
        """Run-level 3' Untargeted assignment rate.

        ``Σ TotalReadsIsCountedInCell / Σ TotalInstrumentReads × 100``
        across samples — the fraction of raw sequencer polonies that
        survive trim → align → assign → counted-class → in-cell and
        end up as a count in ``adata.X``.

        Pass ``batch_names`` to restrict to a subset (e.g. only the 3p
        batches when the f2p stats happen to carry other batch types).
        Returns ``None`` when no rows qualify or the denominator is
        zero.

        Mirrors the signature of
        :meth:`Cells2StatsStats.protein_assignment_rate` and
        :meth:`Cells2StatsStats.ops_assignment_rate` so callers can
        request "assignment rate for this modality" uniformly.
        """
        counted, instrument = self.in_cell_and_instrument(batch_names)
        if counted is None or not instrument:
            return None
        return counted / instrument * 100.0

    def in_cell_and_instrument(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> tuple[float | None, float | None]:
        """``(Σ TotalReadsIsCountedInCell, Σ TotalInstrumentReads)`` across samples
        (optionally restricted to *batch_names*); ``(None, None)`` when no sample
        carries both counters.

        The parquet ``IsCountedInCell`` numerator here is the read-level
        **cross-check** — the table (``adata.X``) is the truth for in-cell yield.
        """
        accepted = set(batch_names) if batch_names else None
        instrument = 0.0
        counted = 0.0
        any_counter = False
        for sample in self.samples:
            batch = sample.get("BatchName")
            if accepted is not None and batch not in accepted:
                continue
            filtered = sample.get("FilteredPolonies") or {}
            assigned = sample.get("Assignment") or {}
            try:
                inst = filtered.get("TotalInstrumentReads")
                cnt = assigned.get("TotalReadsIsCountedInCell")
            except AttributeError:
                continue
            if inst is None or cnt is None:
                continue
            try:
                instrument += float(inst)
                counted += float(cnt)
                any_counter = True
            except (TypeError, ValueError):
                continue
        if not any_counter:
            return None, None
        return counted, instrument

    def read_quality_summary(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> dict[str, int] | None:
        """Run-level always-on read-quality distribution counts, summed across samples.

        Sums the per-sample ``ReadQC`` ``ReadQuality*`` counters:
        ``ReadQualityNReads`` / ``ReadQualityNMappedReads`` (denominators) +
        ``ReadQualityBaseQMin{K}Below{Q}`` / ``ReadQualitySoftclipOver{T}`` (threshold
        counts). Raw counts — divide by the denominators to estimate what a read filter
        at any (K, Q) / softclip_max would remove. Returns ``None`` when no sample carries
        them (older RunStats).
        """
        accepted = set(batch_names) if batch_names else None
        acc: dict[str, int] = {}
        found = False
        for sample in self.samples:
            if accepted is not None and sample.get("BatchName") not in accepted:
                continue
            readqc = sample.get("ReadQC") or {}
            for key, value in readqc.items():
                if not key.startswith("ReadQuality"):
                    continue
                try:
                    acc[key] = acc.get(key, 0) + int(value)
                except (TypeError, ValueError):
                    continue
                found = True
        return acc if found else None

    def countable_reads_breakdown(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> dict[str, list[dict[str, Any]]] | None:
        """Run-level countable-read breakdown by gene biotype and feature tier.

        Reads the per-sample ``Reads`` section ``CountableReads_<bin>`` /
        ``CountableReads_<bin>_PF`` counters (biotype × tier cross-tab) and
        sums them across samples (optionally restricted to *batch_names*).
        ``CountableReads_tier_<bin>`` keys are the feature-tier dimension; the
        remaining ``CountableReads_<bin>`` keys are gene biotypes. Each bin is
        split ALL (``IsCountable``) vs PF-pass (``IsCountable AND pf_pass``);
        PF-fail is derived (``all - pf_pass``). These are **read-level** counts
        — distinct from the UMI-level biotype composition of the count matrix.

        Returns ``{"biotype": [...], "tier": [...]}`` where each list holds one
        dict per bin ``{"bin", "all", "pf_pass", "pf_fail", "pct_kept"}`` sorted
        by ``all`` descending, or ``None`` when no sample carries the counters
        (older or read-filter-off SDOs).
        """
        accepted = set(batch_names) if batch_names else None
        # dimension -> bin -> {"all": int, "pf_pass": int}
        acc: dict[str, dict[str, dict[str, int]]] = {"biotype": {}, "tier": {}}
        found_any = False
        for sample in self.samples:
            if accepted is not None and sample.get("BatchName") not in accepted:
                continue
            reads = sample.get("Reads") or {}
            for key, value in reads.items():
                if not key.startswith("CountableReads_"):
                    continue
                try:
                    count = int(value)
                except (TypeError, ValueError):
                    continue
                found_any = True
                rest = key[len("CountableReads_") :]
                # Newer RunStats renamed the pass-slice suffix _PF -> _FilterPass
                # (accept both so older and newer RunStats parse).
                if rest.endswith("_FilterPass"):
                    is_pf, rest = True, rest[: -len("_FilterPass")]
                elif rest.endswith("_PF"):
                    is_pf, rest = True, rest[: -len("_PF")]
                else:
                    is_pf = False
                if rest.startswith("tier_"):
                    # Normalize legacy tier tokens to the current vocabulary so old +
                    # new RunStats aggregate under one bin (aliases merge / sum).
                    dim, bin_name = "tier", canonical_tier(rest[len("tier_") :])
                else:
                    dim, bin_name = "biotype", rest
                slot = acc[dim].setdefault(bin_name, {"all": 0, "pf_pass": 0})
                slot["pf_pass" if is_pf else "all"] += count
        if not found_any:
            return None

        def _rows(bins: dict[str, dict[str, int]]) -> list[dict[str, Any]]:
            out: list[dict[str, Any]] = []
            for bin_name, sums in bins.items():
                total, pf = sums["all"], sums["pf_pass"]
                out.append(
                    {
                        "bin": bin_name,
                        "all": total,
                        "pf_pass": pf,
                        "pf_fail": total - pf,
                        "pct_kept": (pf / total * 100.0) if total else None,
                    }
                )
            out.sort(key=lambda r: r["all"], reverse=True)
            return out

        return {"biotype": _rows(acc["biotype"]), "tier": _rows(acc["tier"])}

    def tier_mapping_type_breakdown(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> list[dict[str, Any]] | None:
        """Per-tier mapping-type breakdown of all assigned reads.

        Reads ``AssignedReads_tier_{tier}_{mtype}`` keys (unique / multi_simple /
        multi_complex) from the ``Reads`` section and sums across samples.
        Covers ALL assigned reads regardless of IsCounted — multi_complex reads are
        never counted by policy, so they appear here but not in the count matrix.
        Returns a list of dicts ``{"tier", "unique", "multi_simple", "multi_complex",
        "total"}`` sorted by total descending, or ``None`` when the keys are absent
        (older RunStats).
        """
        accepted = set(batch_names) if batch_names else None
        _MTYPES = ("unique", "multi_simple", "multi_complex")
        acc: dict[str, dict[str, int]] = {}
        found_any = False
        for sample in self.samples:
            if accepted is not None and sample.get("BatchName") not in accepted:
                continue
            reads = sample.get("Reads") or {}
            for key, value in reads.items():
                if not key.startswith("AssignedReads_tier_"):
                    continue
                rest = key[len("AssignedReads_tier_") :]
                # expect "{tier}_{mtype}" — split on the last underscore-delimited mtype suffix
                matched = next((m for m in _MTYPES if rest.endswith(f"_{m}")), None)
                if matched is None:
                    continue
                tier = canonical_tier(rest[: -(len(matched) + 1)])
                try:
                    count = int(value)
                except (TypeError, ValueError):
                    continue
                found_any = True
                acc.setdefault(tier, {"unique": 0, "multi_simple": 0, "multi_complex": 0})[matched] += count
        if not found_any:
            return None
        out = []
        for tier, counts in acc.items():
            total = counts["unique"] + counts["multi_simple"] + counts["multi_complex"]
            out.append({"tier": tier, **counts, "total": total})
        out.sort(key=lambda r: r["total"], reverse=True)
        return out

    def tier_quality_breakdown(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> list[dict[str, Any]] | None:
        """Per-tier mean softclip % and error rate for counted reads.

        Reads ``SoftclipPct_tier_{tier}`` and ``ErrorRateMean_tier_{tier}`` from the
        ``ReadQC`` section (weighted by ``CountedReads_tier_{tier}`` from ``Reads`` so the
        run-level mean reflects actual read-count proportions). Returns a list of dicts
        ``{"tier", "softclip_pct", "error_rate"}`` sorted by ``softclip_pct`` descending, or
        ``None`` when the keys are absent (older RunStats).
        """
        accepted = set(batch_names) if batch_names else None
        # Weighted accumulation: {tier: [sum_sc_weighted, sum_er_weighted, total_reads]}
        acc: dict[str, list[float]] = {}
        found_any = False
        for sample in self.samples:
            if accepted is not None and sample.get("BatchName") not in accepted:
                continue
            readqc = sample.get("ReadQC") or {}
            reads = sample.get("Reads") or {}
            for key, sc_val in readqc.items():
                if not key.startswith("SoftclipPct_tier_"):
                    continue
                # Correlate ErrorRateMean_/CountedReads_ by the RAW token (the sibling
                # keys carry the same legacy token), then bin under the canonical name.
                tier_raw = key[len("SoftclipPct_tier_") :]
                er_val = readqc.get(f"ErrorRateMean_tier_{tier_raw}")
                n_reads = reads.get(f"CountedReads_tier_{tier_raw}", 0)
                try:
                    sc = float(sc_val)
                    er = float(er_val) if er_val is not None else 0.0
                    n = float(n_reads)
                except (TypeError, ValueError):
                    continue
                found_any = True
                tier = canonical_tier(tier_raw)
                slot = acc.setdefault(tier, [0.0, 0.0, 0.0])
                slot[0] += sc * n
                slot[1] += er * n
                slot[2] += n
        if not found_any:
            return None
        out = []
        for tier, (sc_w, er_w, n) in acc.items():
            out.append(
                {
                    "tier": tier,
                    "softclip_pct": sc_w / n if n else 0.0,
                    "error_rate": er_w / n if n else 0.0,
                    "n_reads": n,
                }
            )
        out.sort(key=lambda r: r["softclip_pct"], reverse=True)
        return out

    def uncounted_quality_breakdown(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> list[dict[str, Any]] | None:
        """Softclip % and read count for non-counted mapped reads by ``why_not_counted`` reason.

        Reads ``SoftclipPct_uncounted_{reason}`` and ``ReadCount_uncounted_{reason}`` from
        the ``Assignment`` section. Weighted-averages softclip across
        samples. Returns a list of dicts ``{"reason", "softclip_pct", "n_reads"}`` or
        ``None`` when absent (older RunStats).
        """
        accepted = set(batch_names) if batch_names else None
        acc: dict[str, list[float]] = {}
        found_any = False
        for sample in self.samples:
            if accepted is not None and sample.get("BatchName") not in accepted:
                continue
            asg = sample.get("Assignment") or {}
            for key, sc_val in asg.items():
                if not key.startswith("SoftclipPct_uncounted_"):
                    continue
                reason = key[len("SoftclipPct_uncounted_") :]
                n_val = asg.get(f"ReadCount_uncounted_{reason}", 0)
                try:
                    sc = float(sc_val)
                    n = float(n_val)
                except (TypeError, ValueError):
                    continue
                found_any = True
                slot = acc.setdefault(reason, [0.0, 0.0])
                slot[0] += sc * n
                slot[1] += n
        if not found_any:
            return None
        out = []
        for reason, (sc_w, n) in acc.items():
            out.append(
                {
                    "reason": reason,
                    "softclip_pct": sc_w / n if n else 0.0,
                    "n_reads": n,
                }
            )
        out.sort(key=lambda r: r["n_reads"], reverse=True)
        return out

    def tier_quality_table(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> list[dict[str, Any]] | None:
        """Per-tier × mapping-type quality diagnostics.

        Reads the ``TierQuality`` top-level dict from each sample
        (``{tier: {mtype: {softclip_pct, error_rate_mean, mean_read_length,
        mean_aln_length, n}}}``). Aggregates across samples using read-count
        weighting (``n`` stored per entry). Returns a list of dicts
        ``{"tier", "mtype", "softclip_pct", "error_rate_mean",
        "mean_read_length", "mean_aln_length", "n"}`` sorted by tier then
        mtype, or ``None`` when ``TierQuality`` is absent (older runs).
        """
        accepted = set(batch_names) if batch_names else None
        # acc: {tier: {mtype: {metric: weighted_sum, "n": total_n}}}
        acc: dict[str, dict[str, dict[str, float]]] = {}
        found_any = False
        for sample in self.samples:
            if accepted is not None and sample.get("BatchName") not in accepted:
                continue
            tq = sample.get("TierQuality")
            if not isinstance(tq, dict):
                continue
            for tier_raw, mtypes in tq.items():
                if not isinstance(mtypes, dict):
                    continue
                tier = canonical_tier(tier_raw)
                for mtype, entry in mtypes.items():
                    if not isinstance(entry, dict):
                        continue
                    try:
                        n = float(entry["n"])
                    except (KeyError, TypeError, ValueError):
                        continue
                    found_any = True
                    slot = acc.setdefault(tier, {}).setdefault(mtype, {"n": 0.0})
                    slot["n"] += n
                    for k, v in entry.items():
                        if k == "n":
                            continue
                        with contextlib.suppress(TypeError, ValueError):
                            slot[k] = slot.get(k, 0.0) + float(v) * n
        if not found_any:
            return None
        out = []
        for tier, mtypes in acc.items():
            for mtype, slot in mtypes.items():
                total_n = slot["n"]
                row: dict[str, Any] = {"tier": tier, "mtype": mtype, "n": total_n}
                for k, v in slot.items():
                    if k == "n":
                        continue
                    row[k] = v / total_n if total_n else 0.0
                out.append(row)
        out.sort(key=lambda r: (r["tier"], r["mtype"]))
        return out

    def uncounted_quality_table(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> list[dict[str, Any]] | None:
        """Per-reason quality diagnostics for non-counted reads.

        Reads the ``UncountedQuality`` top-level dict from each sample
        (``{reason: {softclip_pct, error_rate_mean, mean_read_length,
        mean_aln_length, n}}``). Aggregates across samples using read-count
        weighting. Returns a list of dicts ``{"reason", "softclip_pct",
        "error_rate_mean", "mean_read_length", "mean_aln_length", "n"}``
        sorted by ``n`` descending, or ``None`` when ``UncountedQuality`` is
        absent (older runs).
        """
        accepted = set(batch_names) if batch_names else None
        acc: dict[str, dict[str, float]] = {}
        found_any = False
        for sample in self.samples:
            if accepted is not None and sample.get("BatchName") not in accepted:
                continue
            uq = sample.get("UncountedQuality")
            if not isinstance(uq, dict):
                continue
            for reason, entry in uq.items():
                if not isinstance(entry, dict):
                    continue
                try:
                    n = float(entry["n"])
                except (KeyError, TypeError, ValueError):
                    continue
                found_any = True
                slot = acc.setdefault(reason, {"n": 0.0})
                slot["n"] += n
                for k, v in entry.items():
                    if k == "n":
                        continue
                    with contextlib.suppress(TypeError, ValueError):
                        slot[k] = slot.get(k, 0.0) + float(v) * n
        if not found_any:
            return None
        out = []
        for reason, slot in acc.items():
            total_n = slot["n"]
            row: dict[str, Any] = {"reason": reason, "n": total_n}
            for k, v in slot.items():
                if k == "n":
                    continue
                row[k] = v / total_n if total_n else 0.0
            out.append(row)
        out.sort(key=lambda r: r["n"], reverse=True)
        return out

    def read_accounting_totals(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> dict[str, float] | None:
        """Run-level read-fate totals (for a read-accounting waterfall).

        Sums the read-fate RunStats keys across samples. Returns a flat dict of read-fate counters, or ``None`` when the base counters
        are absent (extremely old runs with no Assignment section).
        """
        accepted = set(batch_names) if batch_names else None
        _KEYS = [
            # Reads section
            ("Reads", "TotalReads"),
            ("Reads", "TotalReadsInCell"),
            ("Reads", "TotalReadsOutOfCell"),
            # Assignment section
            ("Assignment", "TotalReadsIsCountedInCell"),
            ("Assignment", "TotalReadsIsCounted"),
            ("Assignment", "TotalReadsAssignedOutOfCell"),
            ("Assignment", "TotalReadsAssignedMultiComplexInCell"),
            ("Assignment", "TotalReadsUnassignedNoFeatureInCell"),
            ("Assignment", "TotalReadsUnassignedUnmappedInCell"),
            ("Assignment", "TotalReadsUnassignedNoFeatureOutOfCell"),
            ("Assignment", "TotalReadsUnassignedUnmappedOutOfCell"),
            ("Assignment", "NotCounted_multimap_complex"),
            ("Assignment", "NotCounted_unmapped"),
            ("Assignment", "NotCounted_unassigned_no_feature"),
            # Newer PolonyCategory split; absent on older runs, where
            # NotCounted_unassigned_no_feature still carries both populations.
            ("Assignment", "NotCounted_unassigned_unselected"),
        ]
        acc: dict[str, float] = {}
        found_any = False
        for sample in self.samples:
            if accepted is not None and sample.get("BatchName") not in accepted:
                continue
            for section, key in _KEYS:
                val = (sample.get(section) or {}).get(key)
                if val is None:
                    continue
                try:
                    acc[key] = acc.get(key, 0.0) + float(val)
                    found_any = True
                except (TypeError, ValueError):
                    pass
        return acc if found_any else None

    def read_fate_hierarchy(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> dict[str, Any] | None:
        """Hierarchical read-fate breakdown.

        Location (In-cell / Out-of-cell) → Assigned / Unassigned → mapping type
        (unique / multi_simple / multi_complex) or unassigned reason
        (no-feature / unmapped).  Every node sums to its parent and
        In-cell + Out-of-cell = TotalReads, all from RunStats counters.  Returns ``{"total", "counted_in_cell",
        "counted_total", "rows"}`` where each row is ``{"label", "depth", "reads",
        "parent_reads"}`` in display order, or ``None`` when the base location
        counters are absent.
        """
        accepted = set(batch_names) if batch_names else None
        acc: dict[str, float] = {}
        for sample in self.samples:
            if accepted is not None and sample.get("BatchName") not in accepted:
                continue
            for section in ("Reads", "Assignment"):
                for k, v in (sample.get(section) or {}).items():
                    if isinstance(v, (int, float)):
                        acc[k] = acc.get(k, 0.0) + float(v)

        total = acc.get("TotalReads")
        if not total or acc.get("TotalReadsInCell") is None:
            return None

        def g(key: str) -> float:
            return acc.get(key, 0.0)

        rows: list[dict[str, Any]] = []
        for loc, loc_label in (("InCell", "In-cell"), ("OutOfCell", "Out of cell")):
            loc_reads = g(f"TotalReads{loc}")
            rows.append({"label": loc_label, "depth": 0, "reads": loc_reads, "parent_reads": total})
            assigned = g(f"TotalReadsAssigned{loc}")
            rows.append({"label": "Assigned", "depth": 1, "reads": assigned, "parent_reads": loc_reads})
            for mt, mt_label in (
                ("Unique", "unique"),
                ("MultiSimple", "multi-simple"),
                ("MultiComplex", "multi-complex"),
            ):
                rows.append(
                    {
                        "label": mt_label,
                        "depth": 2,
                        "reads": g(f"TotalReadsAssigned{mt}{loc}"),
                        "parent_reads": assigned,
                    }
                )
            unassigned = g(f"TotalReadsUnassigned{loc}")
            rows.append({"label": "Unassigned", "depth": 1, "reads": unassigned, "parent_reads": loc_reads})
            no_feature = g(f"TotalReadsUnassignedNoFeature{loc}")
            rows.append(
                {
                    "label": "no-feature",
                    "depth": 2,
                    "reads": no_feature,
                    "parent_reads": unassigned,
                }
            )
            # Exact split of the no-feature band (newer RunStats). Emitted only when
            # present, so older runs keep the two-row shape rather than showing zeros.
            unselected = g(f"TotalReadsUnassignedUnselected{loc}")
            no_overlap = g(f"TotalReadsUnassignedNoOverlap{loc}")
            if unselected or no_overlap:
                rows.append(
                    {
                        "label": "unselected gene",
                        "depth": 3,
                        "reads": unselected,
                        "parent_reads": no_feature,
                    }
                )
                rows.append(
                    {
                        "label": "no overlap",
                        "depth": 3,
                        "reads": no_overlap,
                        "parent_reads": no_feature,
                    }
                )
            rows.append(
                {
                    "label": "unmapped",
                    "depth": 2,
                    "reads": g(f"TotalReadsUnassignedUnmapped{loc}"),
                    "parent_reads": unassigned,
                }
            )
        return {
            "total": total,
            "counted_in_cell": g("TotalReadsIsCountedInCell"),
            "counted_total": g("TotalReadsIsCounted"),
            "rows": rows,
        }

    def unselected_reads_breakdown(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> dict[str, list[dict[str, Any]]] | None:
        """Composition of the *unselected* mapped reads — what we are NOT counting.

        Prefers the EXACT per-sample read counts in newer RunStats
        (``UnselectedReads_tier_<tier>`` / ``UnselectedReads_biotype_<bin>``, summed straight
        across samples). Falls back to the older sampled percentages
        (``UnselectedTierPct_`` / ``UnselectedBiotypePct_``), aggregated band-weighted by each
        sample's estimated unselected-read count so wells contribute in proportion to size.

        Returns ``{"tier": [...], "biotype": [...], "exact": bool}`` where each list is one dict
        per bin ``{"bin", "reads", "pct"}`` sorted by ``reads`` descending, or ``None`` when
        neither source is present.

        **No strand dimension.** ``UnselectedStrandPct_*`` came from a both-strand
        full-annotation pass that newer runs don't perform, so it has no exact equivalent.
        """
        accepted = set(batch_names) if batch_names else None
        exact: dict[str, dict[str, float]] = {"tier": {}, "biotype": {}}
        sampled: dict[str, dict[str, float]] = {"tier": {}, "biotype": {}}
        _EXACT_PREFIXES = (
            ("UnselectedReads_tier_", "tier"),
            ("UnselectedReads_biotype_", "biotype"),
        )
        _SAMPLED_PREFIXES = (
            ("UnselectedTierPct_", "tier"),
            ("UnselectedBiotypePct_", "biotype"),
        )
        for sample in self.samples:
            if accepted is not None and sample.get("BatchName") not in accepted:
                continue
            reads = sample.get("Reads") or {}
            got_exact = False
            for key, value in reads.items():
                for prefix, dim in _EXACT_PREFIXES:
                    if not key.startswith(prefix):
                        continue
                    bin_name = key[len(prefix) :]
                    if dim == "tier":
                        bin_name = canonical_tier(bin_name)
                    try:
                        n = float(value)
                    except (TypeError, ValueError):
                        break
                    got_exact = True
                    exact[dim][bin_name] = exact[dim].get(bin_name, 0.0) + n
                    break
            if got_exact:
                continue
            # Sampled fallback: percentages of THIS sample's unselected reads, so they
            # need band weighting before they can be summed across wells.
            asg = sample.get("Assignment") or {}
            try:
                unsel_reads = (
                    float(asg.get("TotalReadsUnassignedNoFeature"))
                    * float(asg.get("UnassignedUnselectedPct"))
                    / 100.0
                )
            except (TypeError, ValueError):
                continue
            if unsel_reads <= 0:
                continue
            for key, value in reads.items():
                for prefix, dim in _SAMPLED_PREFIXES:
                    if not key.startswith(prefix):
                        continue
                    bin_name = key[len(prefix) :]
                    if dim == "tier":
                        bin_name = canonical_tier(bin_name)
                    try:
                        pct = float(value)
                    except (TypeError, ValueError):
                        break
                    sampled[dim][bin_name] = sampled[dim].get(bin_name, 0.0) + unsel_reads * pct / 100.0
                    break

        # Exact wins outright — never blended, same reasoning as unselected_split.
        found_exact = any(exact[d] for d in exact)
        acc = exact if found_exact else sampled
        if not any(acc[d] for d in acc):
            return None

        def _rows(bins: dict[str, float]) -> list[dict[str, Any]]:
            total = sum(bins.values()) or 1.0
            out = [{"bin": b, "reads": v, "pct": v / total * 100.0} for b, v in bins.items()]
            out.sort(key=lambda r: r["reads"], reverse=True)
            return out

        return {"tier": _rows(acc["tier"]), "biotype": _rows(acc["biotype"]), "exact": found_exact}

    def proximity_gap_breakdown(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> dict[str, Any] | None:
        """Run-level histogram of distance from unassigned mapped reads to selected SAF.

        Sums ``ProximityGap_bin_*`` / ``ProximityGap_{Unselected,NoOverlap}_bin_*``
        across samples. Returns ``None`` when no sample carries the keys.
        Does not change assignment counts.
        """
        accepted = set(batch_names) if batch_names else None
        total = dict.fromkeys(_PROXIMITY_GAP_BIN_ORDER, 0.0)
        unselected = dict.fromkeys(_PROXIMITY_GAP_BIN_ORDER, 0.0)
        no_overlap = dict.fromkeys(_PROXIMITY_GAP_BIN_ORDER, 0.0)
        n_samples = 0
        pct_le_100_w = pct_le_250_w = weight = 0.0
        for sample in self.samples:
            if accepted is not None and sample.get("BatchName") not in accepted:
                continue
            reads = sample.get("Reads") or {}
            n = reads.get("ProximityGap_n")
            if n is None:
                continue
            try:
                n_f = float(n)
            except (TypeError, ValueError):
                continue
            if n_f <= 0:
                continue
            got = False
            for key in _PROXIMITY_GAP_BIN_ORDER:
                raw = reads.get(f"ProximityGap_bin_{key}")
                if raw is None:
                    continue
                try:
                    total[key] += float(raw)
                    got = True
                except (TypeError, ValueError):
                    continue
                try:
                    unselected[key] += float(reads.get(f"ProximityGap_Unselected_bin_{key}") or 0)
                    no_overlap[key] += float(reads.get(f"ProximityGap_NoOverlap_bin_{key}") or 0)
                except (TypeError, ValueError):
                    pass
            if not got:
                continue
            n_samples += 1
            weight += n_f
            try:
                pct_le_100_w += float(reads.get("ProximityGap_pct_le_100") or 0) * n_f
                pct_le_250_w += float(reads.get("ProximityGap_pct_le_250") or 0) * n_f
            except (TypeError, ValueError):
                pass

        n_total = sum(total.values())
        n_children = sum(unselected.values()) + sum(no_overlap.values())
        if n_samples == 0 or n_total <= 0 or n_children <= 0:
            return None

        def _rows(bins: dict[str, float]) -> list[dict[str, Any]]:
            return [
                {
                    "bin": key,
                    "label": _PROXIMITY_GAP_LABELS[key],
                    "reads": bins[key],
                    "pct": bins[key] / n_total * 100.0,
                }
                for key in _PROXIMITY_GAP_BIN_ORDER
            ]

        return {
            "n": n_total,
            "n_samples": n_samples,
            "pct_le_100": (pct_le_100_w / weight) if weight else 0.0,
            "pct_le_250": (pct_le_250_w / weight) if weight else 0.0,
            "bins": _rows(total),
            "unselected": _rows(unselected),
            "no_overlap": _rows(no_overlap),
        }

    def unselected_split(
        self, batch_names: list[str] | tuple[str, ...] | None = None
    ) -> dict[str, float | None] | None:
        """Run-level split of the no-feature band into unselected vs no-overlap.

        Prefers the EXACT per-sample counts in newer RunStats
        (``TotalReadsUnassignedUnselected`` / ``TotalReadsUnassignedNoOverlap``), which
        partition ``TotalReadsUnassignedNoFeature`` with no sampling. Falls back to the
        older sampled percentages (``UnassignedIntergenicPct``), aggregated
        band-weighted so the run-level % reflects the true per-sample band sizes.

        Returns ``{"unselected_pct", "no_overlap_pct", "unselected_reads",
        "no_overlap_reads", "no_feature_reads", "exact", "downsample_frac"}``, or ``None``
        when neither source is present. ``exact`` says which source was used;
        ``downsample_frac`` is ``None`` on the exact path.

        The two paths are NOT interchangeable at the same numeric value. The sampled
        diagnostic measured overlap at ``-s 0`` (both strands) while counting is ``-s 2``
        (sense only), so its unselected share runs ~25 points high — reconciling only
        after multiplying by its own reported sense fraction. Never mix them in one series.
        """
        accepted = set(batch_names) if batch_names else None
        # Separate accumulators per source: a mixed-version run must not divide an
        # exact numerator by a band that also includes sampled-only samples.
        band_exact = band_sampled = 0.0
        unselected_exact = no_overlap_exact = 0.0
        unselected_est = 0.0
        downsample_frac: float | None = None
        for sample in self.samples:
            if accepted is not None and sample.get("BatchName") not in accepted:
                continue
            asg = sample.get("Assignment") or {}
            try:
                band = float(asg.get("TotalReadsUnassignedNoFeature"))
            except (TypeError, ValueError):
                continue
            unsel, no_ovl = (
                asg.get("TotalReadsUnassignedUnselected"),
                asg.get("TotalReadsUnassignedNoOverlap"),
            )
            if unsel is not None and no_ovl is not None:
                try:
                    unselected_exact += float(unsel)
                    no_overlap_exact += float(no_ovl)
                except (TypeError, ValueError):
                    continue
                band_exact += band
                continue
            try:
                upct = float(asg.get("UnassignedUnselectedPct"))
            except (TypeError, ValueError):
                continue
            band_sampled += band
            unselected_est += band * upct / 100.0
            if downsample_frac is None:
                raw = (sample.get("Reads") or {}).get("UnselectedDownsampleFrac")
                downsample_frac = float(raw) if raw is not None else None

        # Never blend the two sources — their unselected shares differ by ~25 points
        # for strandedness reasons, so a blended run-level number would be meaningless.
        # Exact wins outright, and the sampled-only samples are then excluded entirely.
        found_exact = band_exact > 0
        if found_exact:
            band_total = band_exact
            unselected_reads, no_overlap_reads = unselected_exact, no_overlap_exact
        elif band_sampled > 0:
            band_total = band_sampled
            unselected_reads = unselected_est
            no_overlap_reads = band_sampled - unselected_est
        else:
            return None
        return {
            "unselected_pct": unselected_reads / band_total * 100.0,
            "no_overlap_pct": no_overlap_reads / band_total * 100.0,
            "unselected_reads": unselected_reads,
            "no_overlap_reads": no_overlap_reads,
            "no_feature_reads": band_total,
            "exact": found_exact,
            "downsample_frac": None if found_exact else downsample_frac,
        }


def parse_fastq2parquet_3p_runstats(run_info: RunInfo) -> Fastq2Parquet3pStats | None:
    """Read and parse ``fastq2parquet_3prime/RunStats.json``."""
    try:
        data = run_info._read_json(_RUNSTATS_PATH)
    except Exception:
        return None
    if not data:
        return None
    return Fastq2Parquet3pStats.from_dict(data)


def parse_multiqc_summary(run_info: RunInfo) -> dict[str, Any] | None:
    """Extract a thin summary from ``multiqc_data.json``.

    Returns a dict with:
    - ``general_stats``: rows from ``report_general_stats_data`` (per-sample)
    - ``modules_present``: list of module short names from ``report_data_sources``
    - ``multiqc_version``: from ``config_short_version``

    Full plot blobs and saved raw data are intentionally not surfaced — the
    multiqc HTML is also archived alongside the JSON for deeper inspection.
    """
    try:
        data = run_info._read_json(_MULTIQC_PATH)
    except Exception:
        return None
    if not data:
        return None

    general_stats: list[dict[str, Any]] = []
    raw_general = data.get("report_general_stats_data")
    if isinstance(raw_general, list):
        for block in raw_general:
            if isinstance(block, dict):
                for sample, metrics in block.items():
                    if isinstance(metrics, dict):
                        row = {"sample": sample}
                        row.update(metrics)
                        general_stats.append(row)

    modules_present = sorted((data.get("report_data_sources") or {}).keys())

    return {
        "general_stats": general_stats,
        "modules_present": modules_present,
        "multiqc_version": data.get("config_short_version"),
    }
