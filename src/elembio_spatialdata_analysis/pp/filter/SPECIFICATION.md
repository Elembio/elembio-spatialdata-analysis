# pp/filter — Filtering Module Specification

## Overview

The `pp/filter` package provides cell-level and gene-level quality-control filtering for
AnnData objects extracted from a SpatialData store. It is modality-aware and supports
Transcript, Protein, Imaging (CellPaint + IF Protein), ThreePrimeUntargeted, and OPS assays.

---

## Package Structure

```
pp/filter/
├── __init__.py            # Public API re-exports
├── _constants.py          # Sentinel, obsm keys, per-modality defaults
├── _quantile.py           # QuantileSpec type, quantile resolution helpers, _resolve_threshold_defaults
├── _spatialdata.py        # _update_spatialdata_attrs
├── _qc_metrics.py         # Per-metric computation (_compute_*, _add_obsm_qc_matrix)
├── _flags.py              # Per-filter flag functions (_flag_*)
├── filter_cells.py        # Public API: filter_cells, filter_genes, compute_filter_diagnostics, get_filter_flags
├── filter_modalities.py   # Public API: filter_modalities
├── report.py              # Public API: generate_filter_report
└── SPECIFICATION.md       # This file
```

---

## Data Contracts

All AnnData objects are read from `sdata.tables["global-table-<Modality>"]`.

### Key obsm / layers / var fields

| Field | Location | Description |
|---|---|---|
| Assigned counts per batch | `adata.obsm['total_counts']` | DataFrame, columns = batch names |
| Unassigned counts per batch | `adata.obsm['unassigned_counts']` | DataFrame, columns = batch names |
| NSB + control probes | `adata.obsm['control_variables']` | DataFrame, columns = probe names |
| GAPDH positive controls | `adata.X` rows where `adata.var['control_type'] == 'Positive Control 1'` | Sparse/dense matrix |
| Nuclear counts | `adata.layers['nuclear']` | Same shape as `adata.X` |
| Cell geometry (Area, etc.) | `adata.obsm['cell_geometry']` | DataFrame; fallback to `adata.obs` |
| Data missing flags | `adata.obsm['data_missing']` | Boolean DataFrame |
| Spatial artifacts | `adata.obsm['spatial_artifacts']` | Boolean DataFrame written by the pipeline. Shared across every cell-aligned table. Columns: per-unit (`imaging-{target}`, `ThreePrimeUntargeted-{batch}`, `OPS-{batch}`), modality rollups (`any-imaging`, `any-ThreePrimeUntargeted`, `any-OPS`), and `any` (OR of rollups). Absent on stores built without spatial-artifact detection and on older stores. When detection was disabled for a run, the full column set is present but all `False`. |
| Artifact params | `adata.uns['spatial_artifact_params']` | Build metadata for the artifact frame (`version`, `spatial_method`, `detection_enabled`, `detection_disabled_reason`, `rollup_columns`, per-modality params). |
| QC metrics (computed) | `adata.obsm['qc_metrics']` | Written by `_add_obsm_qc_matrix` |
| Filter flags (computed) | `adata.obsm['filter_flags']` | Written by `filter_cells` / `compute_filter_diagnostics` (subset to kept cells) |
| Filter flags (full) | `adata.uns['filter_flags_full']` | Pre-subset flag table written by `filter_cells` (primary source for `generate_filter_report`) |
| Filter params (computed) | `adata.uns['filter_params']` | Resolved threshold dict written by `filter_cells` / `compute_filter_diagnostics` |

### adata.obs columns

| Column | Description |
|---|---|
| `Cell` | Cell identifier |
| `Well` | Well label (e.g. `"A1"`) |
| `batch_name` | Feature batch identifier |
| `X_tile`, `Y_tile` | Cell centroid, tile-local pixels (well-canvas coordinates are in `obsm['spatial']`) |

### adata.var columns

| Column | Description |
|---|---|
| `target_type` | Feature type (`"Transcript"`, `"Protein"`, etc.) |
| `control_type` | Control classification (`"Positive Control 1"` = GAPDH, `"Negative Control"` = NSB) |
| `batch_name` | Batch the feature belongs to |

---

## Public API

### `compute_qc_metrics`

```python
compute_qc_metrics(
    adata: AnnData,
    target_type: str,
    *,
    batch_colname: str = "batch_name",
) -> None
```

**Purpose:** Compute QC metrics and store them in ``adata.obsm['qc_metrics']``. Always inplace, no return value — mirrors ``sc.pp.calculate_qc_metrics`` style.

This is the **explicit first step** before plotting or filtering:

```python
ebsd.pp.compute_qc_metrics(rna_adata, "Transcript")  # writes obsm['qc_metrics']
ebsd.vis.plot_total_counts_per_well(rna_adata)  # pre-filter plots work
ebsd.vis.plot_nsb_counts_per_well(rna_adata)

ebsd.pp.filter_cells(rna_adata, "Transcript")  # skips recomputing qc_metrics
ebsd.vis.plot_total_counts_per_well(rna_adata)  # post-filter plots work
ebsd.pp.generate_filter_report(rna_adata)
```

``filter_cells`` will skip recomputing if ``obsm['qc_metrics']`` is already present. No threshold parameters — this function only computes, never filters.

**Metrics written per modality** (all to ``adata.obsm['qc_metrics']``):

| Metric | Modalities | Source |
|---|---|---|
| `total_counts` | All | `obsm['total_counts']` sum, fallback `X` row sum. Imaging: row sum of the primary morphology channels (see Imaging below) |
| `n_genes_by_counts` | Transcript, ThreePrimeUntargeted | `adata.X > 0` count per cell |
| `total_assignment_rate` | Transcript, Protein, ThreePrimeUntargeted | `obsm['total_counts']` / (`obsm['total_counts']` + `obsm['unassigned_counts']`) |
| `total_nsb_counts` | Transcript | `obsm['control_variables']` NSB columns sum |
| `total_gapdh_counts` | Transcript | `adata.X` for `var['control_type'] == 'Positive Control 1'` |
| `nuclear_frac` | Transcript, Protein, ThreePrimeUntargeted | `layers['nuclear']` / `X` row sum (written only when `layers['nuclear']` is present) |
| `pct_counts_mt` | ThreePrimeUntargeted | scanpy `calculate_qc_metrics` (`var["mt"]` = `MT-*` genes) or pre-computed in `obs` |
| `pct_counts_chrm` | ThreePrimeUntargeted | scanpy `calculate_qc_metrics` (`var["chrm"]` = `MT-*` ∪ matrix gene `chrM_umbrella`) or pre-computed in `obs` |
| `pct_counts_ribo`, `pct_counts_hb` | ThreePrimeUntargeted | scanpy `calculate_qc_metrics` (`var["ribo"]` = `RPS*`/`RPL*`; `var["hb"]` = `HB*` excluding `HBP*`) or pre-computed in `obs` |
| `total_counts_mt`, `total_counts_chrm`, `total_counts_ribo`, `total_counts_hb` | ThreePrimeUntargeted | Same scanpy call; absolute counts for each gene set |
| `pct_counts_in_top_10_genes`, `pct_counts_in_top_20_genes` | ThreePrimeUntargeted | Taken from `obs` when pre-computed (scanpy is called with `percent_top=None`) |

---

### `filter_cells`

```python
filter_cells(
    adata: AnnData,
    target_type: str | None = None,
    *,
    batch_colname: str = "batch_name",
    well_colname: str = "Well",
    area_colname: str = "Area",
    min_quant_count: QuantileSpec = _UNSET,
    max_quant_count: QuantileSpec = _UNSET,
    min_quant_area: QuantileSpec = _UNSET,
    max_quant_area: QuantileSpec = _UNSET,
    max_quant_nsb: QuantileSpec = _UNSET,
    max_quant_gapdh: QuantileSpec = _UNSET,
    max_pct_mt: float | None = _UNSET,
    min_assignment_rate: float | None = _UNSET,
    filter_spatial_artifacts: bool = _UNSET,
    inplace: bool = True,
    diagnostics: bool = False,
) -> AnnData | tuple[AnnData, FilterDiagnostics]
```

**Purpose:** Filter cells in an AnnData object based on QC metrics. Gene filtering is handled separately via `filter_genes`.

**Parameters:**
- `target_type` — Required. Modality string. Accepted values:
  - `"Transcript"` or `"transcript"`
  - `"Protein"` or `"protein"`
  - `"Imaging"` or `"imaging"` (also accepts `"CellPaint"` / `"cell_paint"` for backward compat)
  - `"ThreePrimeUntargeted"` or `"threeprimeuntargeted"`
  - `"OPS"` or `"ops"` (also accepts `"SpecializedTargeted"`)
- `batch_colname` — Column in `adata.obs` identifying feature batches. Default `"batch_name"`.
- `well_colname` — Column in `adata.obs` identifying wells. Default `"Well"`.
- `area_colname` — Column name for cell area. Default `"Area"`. Resolved from `obsm['cell_geometry']` first, then `obs`.
- `min_quant_count` — Lower quantile threshold for per-cell counts per batch per well. `None` disables. Default: modality-specific (0.025 for Transcript/Protein/ThreePrimeUntargeted; `None` for OPS).
- `max_quant_count` — Upper quantile threshold for per-cell counts. `None` or `>= 1.0` disables. Default: modality-specific (0.99 for Transcript/Protein/ThreePrimeUntargeted, i.e. the top 1% of cells by count are dropped; `None` for OPS).
- `min_quant_area` — Lower quantile threshold for cell area per well. Default: 0.01.
- `max_quant_area` — Upper quantile threshold for cell area per well. Default: 0.99.
- `max_quant_nsb` — Upper quantile threshold for NSB control counts. Reads from `obsm['control_variables']`. Default: 0.99 (Transcript only).
- `max_quant_gapdh` — Upper quantile threshold for GAPDH positive-control counts. Reads from `adata.X` via `var['control_type']`. **Default: `None` (opt-in only)**. Set to e.g. `0.99` to enable.
- `max_pct_mt` — Maximum mitochondrial gene percentage (0–1). Reads from `obsm['qc_metrics']['pct_counts_mt']`. Default: `0.10` for ThreePrimeUntargeted; `None` (disabled) for every other modality.
- `min_assignment_rate` — Minimum fraction of assigned counts per cell. Reads from `obsm['total_counts']` + `obsm['unassigned_counts']`. **Default: `None` (disabled for all modalities).** `obsm['unassigned_counts']` is a well-level background count — transcripts not assigned to any cell — not a per-cell quality metric. Low-count cells are already removed by `min_quant_count`; enabling this filter additionally penalises cells in high-background wells regardless of their own count quality. Pass `min_assignment_rate=0.5` explicitly to opt in.
- `filter_spatial_artifacts` — If `True`, drop cells flagged in the **modality rollup** column of `obsm['spatial_artifacts']`. Default: `False` (flag is written as all-True; `generate_filter_report` hides the row). Column used: `any-imaging` (Imaging), `any-ThreePrimeUntargeted` (3′), `any-OPS` (OPS). Transcript / Protein have no detector — opting in warns and is a no-op. Missing slot or missing rollup column (store built without spatial-artifact detection, or an older store) warns and is a no-op. Per-unit columns (`imaging-Actin`, `ThreePrimeUntargeted-B01`, …) are left for inspection; `filter_cells` does not gate on them. `filter_modalities` forwards this kwarg per table; when every present artifact-bearing modality opts in, the intersection equals `~any`.
- `inplace` — If `True`, subset `adata` in place (default). Pass `False` to return a filtered copy without modifying the original.
- `diagnostics` — If `True`, return `(adata, FilterDiagnostics)`. Default: `False`.

**Returns:** `AnnData` or `(AnnData, FilterDiagnostics)`.

---

### `filter_genes`

```python
filter_genes(
    adata: AnnData,
    target_type: str | None = None,
    *,
    well_colname: str = "Well",
    min_cells_per_gene: int | None = _UNSET,
    min_counts_per_gene: float | None = _UNSET,
    inplace: bool = True,
) -> AnnData
```

**Purpose:** Filter genes from an AnnData object based on per-well expression thresholds. Typically called after `filter_cells`. Primarily intended for `ThreePrimeUntargeted` data but works for any modality.

**Parameters:**
- `target_type` — Optional. When provided, default thresholds are looked up from `DEFAULT_BOUNDS_PER_MODALITY`. For `ThreePrimeUntargeted` the defaults are non-`None`; for all other modalities they default to `None` (no filtering unless explicitly passed).
- `well_colname` — Column in `adata.obs` identifying wells. Default `"Well"`.
- `min_cells_per_gene` — Minimum cells expressing each gene within at least one well. Default: modality-specific (`None` for most modalities; `10` for ThreePrimeUntargeted).
- `min_counts_per_gene` — Minimum total counts per gene within at least one well. Default: modality-specific.
- `inplace` — If `True`, subset `adata` in place (default). Pass `False` to return a filtered copy.

**Returns:** `AnnData` with low-quality genes removed.

**QuantileSpec type:** `float | Mapping[str, float | Mapping[str, float]] | None`

Override schema for per-well / per-batch quantiles:
```python
# Global scalar
min_quant_count = 0.025

# Per-batch override (applies across all wells)
min_quant_count = {"B01": 0.01, "__default__": 0.025}

# Per-well override
min_quant_count = {"A1": {"B01": 0.01, "__default__": 0.02}, "__default__": 0.025}

# Disable for a specific batch
min_quant_count = {"B01": None, "__default__": 0.025}
```

---

### `compute_filter_diagnostics`

```python
compute_filter_diagnostics(
    adata: AnnData,
    target_type: str | None = None,
    *,
    ...same threshold params as filter_cells...,
    min_cells_per_gene: int | None = None,
    min_counts_per_gene: float | None = None,
    inplace: bool = True,
    return_diagnostics: bool = False,
) -> FilterDiagnostics | None
```

**Purpose:** Compute QC metrics and filter flags without subsetting `adata`. Always writes `adata.obsm['filter_flags']` and `adata.uns['filter_params']`, and writes `adata.obsm['qc_metrics']` if it is not already present. `inplace` does not change this; it is only recorded in the thresholds. `min_cells_per_gene` / `min_counts_per_gene` are recorded in the thresholds but do not produce flags.

**Returns:** `FilterDiagnostics` if `return_diagnostics=True`, else `None`.

**Namespace:** `compute_filter_diagnostics` and `get_filter_flags` are exported from `ebsd.pp.filter` only, not from `ebsd.pp`.

---

### `get_filter_flags`

```python
get_filter_flags(
    adata: AnnData,
    *,
    target_type: str,
    batch_colname: str = "batch_name",
    well_colname: str = "Well",
    area_colname: str = "Area",
    min_total_count: float = 1.0,
    min_quant_count: QuantileSpec = 0.025,
    max_quant_count: QuantileSpec = None,
    min_quant_area: QuantileSpec = 0.01,
    max_quant_area: QuantileSpec = 0.99,
    max_quant_nsb: QuantileSpec = 0.99,
    max_quant_gapdh: QuantileSpec = None,
    max_pct_mt: float | None = None,
    min_assignment_rate: float | None = None,
    filter_spatial_artifacts: bool = False,
) -> pd.DataFrame
```

**Purpose:** Low-level dispatcher that returns a boolean DataFrame of filter flags. Each column is a `passes_*` flag; `passes_all_filters` is the conjunction.

---

### `filter_modalities`

```python
filter_modalities(
    adata_by_modality: Mapping[str, AnnData],
    *,
    modality_order: Sequence[str] | None = None,
    modality_kwargs: Mapping[str, Mapping[str, Any]] | None = None,
    return_pass_summary: bool = True,
    disable_area_for_non_owners: bool = True,
    **common_kwargs: Any,
) -> tuple[dict[str, ModalityFilterResult], pd.DataFrame | None]
```

**Purpose:** Filter multiple modalities and intersect passing cells across all modalities.

> **Note on modality identifiers**: CellProfiler morphology features live in the `global-table-imaging` table (or legacy `global-table-CellPaint`) and are filtered through the imaging path. Use `"Imaging"` (or `"CellPaint"` for backward compat) as the `target_type`.

**Area ownership in multimodal runs:**
In a multimodal experiment all modalities share the same cell segmentation. Cell area
quality is therefore enforced exactly once — by the Transcript (or ThreePrimeUntargeted)
modality. When `disable_area_for_non_owners=True` (default), `filter_modalities`
automatically sets `min_quant_area=None` and `max_quant_area=None` for every
modality outside `AREA_OWNER_MODALITIES` (i.e. Protein, Imaging, and OPS) before
forwarding to `filter_cells`. Individual overrides via `modality_kwargs` still take
precedence.

Set `disable_area_for_non_owners=False` to apply the same area defaults to every modality
(useful if you are filtering only Protein or Imaging in isolation via `filter_modalities`
without a Transcript table in the same call).

**Returns:**
- `results_by_modality` — dict mapping modality label → `ModalityFilterResult` (filtered AnnData + diagnostics).
- `pass_summary` — DataFrame with columns `passes_<modality>`, `n_modalities_passed`, `passes_all_modalities`, `passes_any_modalities`.

---

### `generate_filter_report`

```python
generate_filter_report(
    filter_results: AnnData | FilterDiagnostics | pd.DataFrame,
    *,
    filter_order: list[str] | None = None,
    include_total: bool = True,
    exclude_filters: Sequence[str] | None = None,
) -> pd.DataFrame
```

**Purpose:** Summarize filtering outcomes with unique vs shared failure counts per filter.

**Accepts:**
- `AnnData` — reads `uns['filter_flags_full']` (written by `filter_cells` before inplace subsetting) with fallback to `obsm['filter_flags']`, plus `uns['filter_params']`. This is the primary usage pattern. Disabled gates (threshold `None`, or `filter_spatial_artifacts=False`) are omitted from the table.
- `FilterDiagnostics` — uses `flags` and `thresholds` fields directly (still supported).
- `pd.DataFrame` — raw boolean flags only (no threshold info in the report).

**Typical usage:**
```python
ebsd.pp.filter_cells(rna_adata, "Transcript")  # writes to obsm + uns in place
ebsd.pp.generate_filter_report(rna_adata)  # reads from obsm + uns
```

**Returns:** DataFrame with columns: `filter`, `cells_passed`, `cells_failed`, `cells_failed_unique`, `cells_failed_overlap`, `pct_failed_total`, `pct_failed_unique`, `pct_failed_overlap`, `pct_passed_total`, `thresholds`. With `include_total=True` (default), a final `TOTAL` row summarizes failures across all filters.

---

### Lower-level helpers

Also exported from `ebsd.pp.filter` as stable public API:

| Name | Purpose |
|---|---|
| `first_failing_reason(flags, *, order=FIRST_FAILING_ORDER)` | Count removed cells by their first failing `passes_*` flag (each cell counted once) |
| `flag_data_missing(adata)` | Pass flag: no indicator set in `obsm['data_missing']` (all-True when absent) |
| `flag_empty_cells(adata, min_total_counts=1.0, feature_names=None)` | Pass flag: row sum of `X` (optionally restricted to `feature_names`) ≥ `min_total_counts` |
| `flag_failed_batches(adata, batch_colname="batch_name", feature_names=None)` | Pass flag: no NaN values in `X` (failed batches) |
| `resolve_threshold_defaults(normalized_target_type, overrides)` | Merge `DEFAULT_BOUNDS_PER_MODALITY` for a modality with explicit overrides |
| `DEFAULT_BOUNDS_PER_MODALITY` | Per-modality default thresholds (see below) |
| `AREA_OWNER_MODALITIES` | Modalities that enforce area QC in `filter_modalities`: `{"transcript", "threeprimeuntargeted"}` |
| `QC_METRICS_OBSM_KEY` | `"qc_metrics"` |
| `FILTER_FLAGS_FULL_UNS_KEY` | `"filter_flags_full"` |

---

## Dataclasses

### `FilterDiagnostics`

```python
@dataclass(slots=True)
class FilterDiagnostics:
    flags: pd.DataFrame  # Boolean per-cell filter flags (obs × filters)
    thresholds: dict[str, Any]  # Resolved threshold values used
    metadata: dict[str, Any]  # Optional extra metadata
```

### `ModalityFilterResult`

```python
@dataclass(slots=True)
class ModalityFilterResult:
    adata: ad.AnnData  # Filtered AnnData (cells passing all modalities)
    diagnostics: FilterDiagnostics
```

---

## Filter Criteria by Modality

> **`passes_counts`** is the shared count-depth quantile gate. Transcript and
> Protein compute it **per batch** (`obsm['total_counts']`); ThreePrimeUntargeted
> computes it from the **summed `adata.X`, per well** (its batches are summed into
> one expression matrix, so per-batch filtering would wrongly drop cells that
> legitimately lack counts in some batches). Same column, method chosen per
> modality.

### Transcript

| Filter | Flag column | Data source | Default threshold |
|---|---|---|---|
| NaN in X | `passes_nan_filter` | `adata.X` | Always applied |
| Empty cells | `passes_empty_filter` | `adata.X` row sum | `min_total_count=1.0` |
| Total counts (per batch) | `passes_counts` | `obsm['total_counts']` | `min_quant_count=0.025, max_quant_count=0.99` |
| Assignment rate | `passes_assignment_rate` | `obsm['total_counts']` + `obsm['unassigned_counts']` | `min_assignment_rate=None` (disabled; opt-in) |
| Cell area | `passes_area` | `obsm['cell_geometry']['Area']` | `min=0.01, max=0.99` |
| NSB outlier | `passes_nsb_outlier` | `obsm['control_variables']` (NSB* cols) | `max_quant_nsb=0.99` |
| GAPDH outlier | `passes_gapdh_outlier` | `adata.X` (var `control_type=='Positive Control 1'`) | `max_quant_gapdh=None` (opt-in) |
| Spatial artifacts | `passes_spatial_artifacts` | `obsm['spatial_artifacts']` | `filter_spatial_artifacts=False` (opt-in; no Transcript rollup — no-op if enabled) |

### Protein

| Filter | Flag column | Data source | Standalone default | Multimodal default |
|---|---|---|---|---|
| NaN in X | `passes_nan_filter` | `adata.X` | Always applied | Always applied |
| Empty cells | `passes_empty_filter` | `adata.X` row sum | `min_total_count=1.0` | `min_total_count=1.0` |
| Total counts (per batch) | `passes_counts` | `obsm['total_counts']` | `min_quant_count=0.025, max_quant_count=0.99` | same |
| Assignment rate | `passes_assignment_rate` | `obsm['total_counts']` + `obsm['unassigned_counts']` | `None` (opt-in) | `None` (opt-in) |
| Cell area | `passes_area` | `obsm['cell_geometry']['Area']` | `min=0.01, max=0.99` | `None` (disabled — Transcript owns area QC) |
| Spatial artifacts | `passes_spatial_artifacts` | `obsm['spatial_artifacts']` | `False` (opt-in; no Protein rollup — no-op if enabled) | same |

### Imaging (CellPaint + IF Protein)

| Filter | Flag column | Data source | Standalone default | Multimodal default |
|---|---|---|---|---|
| NaN in primary morphology | `passes_nan_filter` | `adata.X` (primary features only) | Always applied | Always applied |
| Empty cells (primary morphology) | `passes_empty_filter` | `adata.X` row sum (primary features only) | `min_total_count=1.0` | `min_total_count=1.0` |
| Primary feature intensities | `passes_imaging_intensity` | `adata.X` (primary features) | `min=0.025, max=0.99` | `min=0.025, max=0.99` |
| Cell area | `passes_area` | `obsm['cell_geometry']['Area']` | `min=0.01, max=0.99` | `None` (disabled — Transcript owns area QC) |
| Spatial artifacts | `passes_spatial_artifacts` | `obsm['spatial_artifacts']['any-imaging']` | `filter_spatial_artifacts=False` (opt-in) | same |

Primary features for the **NaN filter, empty-cell filter, and intensity gate** are the **structural morphology channels**, selected by `tl.morphology.select_morphology_channels` (batch-anchored):

1. **CellPaint morphology batch(es) only** — `batch_name` starts with `"CP"`. This excludes the protein / antibody (`P*`) batch, whose markers and IgG controls must never gate cells (real IF panels label every imaging target `target_type == "CellPaint"`, so the batch, not `target_type`, is the discriminator).
2. **Whole-cell base intensities** — `feature_class == "base"` (current `global-table-imaging` format) or `target_type == "CellPaint"` (legacy `global-table-CellPaint` format), with `suffix == ""` and `is_background == False`.
3. **Fail-safe no-op** — if no CellPaint (`CP*`) batch is present, the NaN filter, empty-cell filter, and intensity gate are skipped and all cells pass (a warning is emitted). None of these gates fall back to the protein batch.

The selected channels are logged (info level) each run. IF antibody markers (`target_type == "ProteinIF"`, in the `P*` batch), CellProfiler descriptors (`feature_class == "cellprofiler"`), compartment variants (`_Nuclear` / `_Cytoplasmic`), and background channels are **never** included in NaN, empty, or intensity gating; they remain in the table for downstream analysis. A NaN confined to an IF / background / CellProfiler / compartment column does **not** fail `passes_nan_filter` or `passes_empty_filter` (the empty gate must not use a full-matrix `.sum()`, which would turn those NaNs into a failed row). `compute_qc_metrics` for Imaging likewise writes `obsm['qc_metrics']['total_counts']` from the same morphology channel sum. The gates only answer "is this a well-segmented, well-stained cell?"

Use `"Imaging"` (or `"CellPaint"` for backward compat) for all imaging-based assays.

### ThreePrimeUntargeted

| Filter | Flag column | Data source | Default threshold |
|---|---|---|---|
| NaN in X | `passes_nan_filter` | `adata.X` | Always applied |
| Empty cells | `passes_empty_filter` | `adata.X` row sum | `min_total_count=1.0` |
| Total counts (summed, per well) | `passes_counts` | `adata.X` row sum | `min_quant_count=0.025, max_quant_count=0.99` |
| Assignment rate | `passes_assignment_rate` | `obsm['total_counts']` + `obsm['unassigned_counts']` | `min_assignment_rate=None` (no default; disabled, opt-in) |
| Cell area | `passes_area` | `obsm['cell_geometry']['Area']` | `min=0.01, max=0.99` |
| MT outlier | `passes_mt_outlier` | `obsm['qc_metrics']['pct_counts_mt']` | `max_pct_mt=0.10` |
| Spatial artifacts | `passes_spatial_artifacts` | `obsm['spatial_artifacts']['any-ThreePrimeUntargeted']` | `filter_spatial_artifacts=False` (opt-in) |

Note: NSB and GAPDH filters are **not applied** for ThreePrimeUntargeted.

### OPS

OPS tables share Protein's count/area shape but have no `unassigned_counts`. Count-quantile gates are **off** by default (identity chemistry, not expression depth). Confidence is not this gate.

| Filter | Flag column | Data source | Default threshold |
|---|---|---|---|
| NaN in X | `passes_nan_filter` | `adata.X` | Always applied |
| Empty cells | `passes_empty_filter` | `adata.X` row sum | `min_total_count=1.0` |
| Total counts | `passes_counts` | `obsm['total_counts']` | `min_quant_count=None, max_quant_count=None` (off) |
| Assignment rate | `passes_assignment_rate` | — | Always `True` (the OPS path forces `min_assignment_rate=None`) |
| Cell area | `passes_area` | `obsm['cell_geometry']['Area']` | `min=0.01, max=0.99` |
| Spatial artifacts | `passes_spatial_artifacts` | `obsm['spatial_artifacts']['any-OPS']` | `filter_spatial_artifacts=False` (opt-in) |

---

## Default Thresholds (`DEFAULT_BOUNDS_PER_MODALITY`)

Defined in `_constants.py` and applied when a threshold argument is not passed. Values below are the effective defaults per modality (`None` = disabled):

| Threshold | Transcript | Protein | Imaging | ThreePrimeUntargeted | OPS |
|---|---|---|---|---|---|
| `min_quant_count` | 0.025 | 0.025 | 0.025 | 0.025 | `None` |
| `max_quant_count` | 0.99 | 0.99 | 0.99 | 0.99 | `None` |
| `min_quant_area` | 0.01 | 0.01 | 0.01 | 0.01 | 0.01 |
| `max_quant_area` | 0.99 | 0.99 | 0.99 | 0.99 | 0.99 |
| `max_quant_nsb` | 0.99 | — | `None` | `None` | `None` |
| `max_quant_gapdh` | `None` | `None` | `None` | `None` | `None` |
| `max_pct_mt` | `None` | `None` | `None` | 0.10 | `None` |
| `min_assignment_rate` | `None` | `None` | `None` | `None` | `None` |
| `filter_spatial_artifacts` | `False` | `False` | `False` | `False` | `False` |
| `min_cells_per_gene` (`filter_genes`) | `None` | `None` | `None` | 10 | `None` |
| `min_counts_per_gene` (`filter_genes`) | `None` | `None` | `None` | `None` | `None` |

For Imaging, `min_quant_count` / `max_quant_count` drive the primary-feature intensity gate (`passes_imaging_intensity`). `—` = not applied for that modality. In `filter_modalities`, area bounds are additionally disabled for non-owner modalities (see above).
