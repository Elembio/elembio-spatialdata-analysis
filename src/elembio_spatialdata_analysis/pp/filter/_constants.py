"""Constants and default thresholds for the filter module."""

from __future__ import annotations

from typing import Any

# Sentinel for "caller did not supply this argument" (distinct from None which means "disabled").
_UNSET = object()

QC_METRICS_OBSM_KEY = "qc_metrics"
FILTER_FLAGS_OBSM_KEY = "filter_flags"
FILTER_PARAMS_UNS_KEY = "filter_params"
# Full pre-subset flags stored in uns so generate_filter_report works after inplace filtering.
FILTER_FLAGS_FULL_UNS_KEY = "filter_flags_full"

# obsm keys in the AnnData tables
OBSM_TOTAL_COUNTS = "total_counts"
OBSM_UNASSIGNED_COUNTS = "unassigned_counts"
OBSM_CONTROL_VARIABLES = "control_variables"
OBSM_CELL_GEOMETRY = "cell_geometry"
OBSM_DATA_MISSING = "data_missing"
OBSM_SPATIAL_ARTIFACTS = "spatial_artifacts"

# Modality rollup columns in obsm['spatial_artifacts'].
# Per-unit columns (imaging-Actin, ThreePrimeUntargeted-B01, OPS-B03) exist
# but filter_cells gates only on these rollups. Cross-modality OR is "any".
SPATIAL_ARTIFACT_ROLLUP: dict[str, str] = {
    "imaging": "any-imaging",
    "threeprimeuntargeted": "any-ThreePrimeUntargeted",
    "ops": "any-OPS",
}

# var column used to identify GAPDH positive controls
VAR_CONTROL_TYPE_COL = "control_type"
VAR_GAPDH_CONTROL_TYPE = "Positive Control 1"

# Default per-modality filter bounds.
# None means the bound is disabled by default.  Every modality dict
# declares every key explicitly so the global fallback in
# _FILTER_PARAM_DEFAULTS never silently substitutes a value the modality
# wouldn't expect (e.g. NSB filtering on 3' Untargeted, which has no NSB
# control probes).  max_quant_gapdh is None (opt-in only) because GAPDH
# counts live in .X and require explicit opt-in to avoid accidental
# filtering on positive controls.
DEFAULT_BOUNDS_PER_MODALITY: dict[str, dict[str, Any]] = {
    "transcript": {
        "min_quant_count": 0.025,
        "max_quant_count": 0.99,
        "min_quant_area": 0.01,
        "max_quant_area": 0.99,
        "max_quant_nsb": 0.99,
        "max_quant_gapdh": None,
        # Targeted panel, no MT genes by design.
        "max_pct_mt": None,
        # Disabled by default: unassigned counts are a well-level background metric,
        # not per-cell quality. Low assignment rate cells are already caught by
        # min_quant_count. Pass min_assignment_rate=0.5 explicitly to opt in.
        "min_assignment_rate": None,
        "min_cells_per_gene": None,
        "min_counts_per_gene": None,
        # Opt-in: drop cells in obsm['spatial_artifacts']['any-*'] regions.
        "filter_spatial_artifacts": False,
    },
    "threeprimeuntargeted": {
        "min_quant_count": 0.025,
        "max_quant_count": 0.99,
        "min_quant_area": 0.01,
        "max_quant_area": 0.99,
        # PolyA capture, no probe-based NSB or GAPDH controls.
        "max_quant_nsb": None,
        "max_quant_gapdh": None,
        "max_pct_mt": 0.10,
        "min_assignment_rate": None,
        # Whole-transcriptome capture: drop genes seen in fewer than 10 cells
        # (in at least one well) as noise. Other modalities are targeted panels
        # and disable gene filtering by default.
        "min_cells_per_gene": 10,
        "min_counts_per_gene": None,
        "filter_spatial_artifacts": False,
    },
    "protein": {
        "min_quant_count": 0.025,
        "max_quant_count": 0.99,
        "min_quant_area": 0.01,
        "max_quant_area": 0.99,
        # NSB filtering is Transcript-only; the Protein path does not apply it.
        "max_quant_nsb": None,
        "max_quant_gapdh": None,
        "max_pct_mt": None,
        "min_assignment_rate": None,
        "min_cells_per_gene": None,
        "min_counts_per_gene": None,
        "filter_spatial_artifacts": False,
    },
    # OPS is a classification readout (one-hot cell-line identity per
    # cell), not an expression measurement.  Per-cell total OPS counts
    # are not a quality metric: a low-count cell isn't "noisy", it's a
    # cell where the assignment chemistry under-detected, which is
    # already captured by the per-cell OPS Confidence metric.
    # Count-quantile filtering would drop cells already visible in the
    # confidence distribution and isn't aligned with OPS's QC model.
    # Area filtering stays on because cell segmentation outliers are
    # genuine segmentation artifacts shared across modalities.
    "ops": {
        "min_quant_count": None,
        "max_quant_count": None,
        "min_quant_area": 0.01,
        "max_quant_area": 0.99,
        # Pooled-screen barcode panel, no NSB / GAPDH / MT in scope.
        "max_quant_nsb": None,
        "max_quant_gapdh": None,
        "max_pct_mt": None,
        "min_assignment_rate": None,
        "min_cells_per_gene": None,
        "min_counts_per_gene": None,
        "filter_spatial_artifacts": False,
    },
    "imaging": {
        "min_quant_count": 0.025,
        "max_quant_count": 0.99,
        "min_quant_area": 0.01,
        "max_quant_area": 0.99,
        # Imaging modality, no sequencing-based control probes.
        "max_quant_nsb": None,
        "max_quant_gapdh": None,
        "max_pct_mt": None,
        "min_assignment_rate": None,
        "min_cells_per_gene": None,
        "min_counts_per_gene": None,
        "filter_spatial_artifacts": False,
    },
}

# Global fallback defaults when a modality does not specify a key.
_FILTER_PARAM_DEFAULTS: dict[str, Any] = {
    "min_quant_count": 0.025,
    "max_quant_count": None,
    "min_quant_area": 0.01,
    "max_quant_area": 0.99,
    "max_quant_nsb": 0.99,
    "max_quant_gapdh": None,
    "max_pct_mt": None,
    "min_assignment_rate": None,
    "min_cells_per_gene": None,
    "min_counts_per_gene": None,
    "filter_spatial_artifacts": False,
}
