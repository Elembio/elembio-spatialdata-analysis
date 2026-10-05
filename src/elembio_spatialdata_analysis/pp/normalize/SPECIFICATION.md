# pp/normalize — Normalization Module Specification

## Overview

The `pp/normalize` package provides modality-aware normalization for AnnData objects
extracted from a SpatialData store.  It mirrors the structure and conventions of
`pp/filter/`:

- A single entry point (`ebsd.pp.normalize`) dispatches to modality-specific pipelines.
- The raw matrix is snapshotted to `layers["raw"]` by default; intermediate snapshots
  are opt-in.
- Defaults are defined in the function signatures.

---

## Package Structure

```
pp/normalize/
├── __init__.py       # Public re-exports
├── _constants.py     # Layer name constants
├── _resolve.py       # _resolve_well_groups helper
├── _core.py          # Private primitives: _normalize_agg_counts,
│                     #   _normalize_median_of_ratios, _log1p, _clr
├── modalities.py     # normalize_transcripts, normalize_proteins,
│                     #   normalize_imaging (stub; normalize_cell_paint alias)
├── normalize.py      # normalize_adata, exported as ebsd.pp.normalize
└── SPECIFICATION.md  # This file
```

---

## Public API

### Primary entry point

```python
ebsd.pp.normalize(
    adata,
    target_type,   # "Transcript", "Protein", "Imaging" (or "CellPaint"), "ThreePrimeUntargeted"
                   # (or lowercase equivalents: "transcript", "protein", "imaging", "threeprimeuntargeted")
    *,
    batch_colname="batch_name",   # adata.var column identifying feature batches
    well_colname="Well",          # adata.obs column identifying wells
    run_colname=None,             # adata.obs run column; None → auto-detect (see below)
    layer=None,                   # read from / write to this layer; None → adata.X
    log_base=None,                # None → natural log
    raw_layer_name="raw",         # pre-normalization snapshot; None → not written
    agg_norm_layer_name=None,     # post-agg-norm snapshot; None → not written
    mor_norm_layer_name=None,     # post-median-of-ratios snapshot (Transcript only); None → not written
    inplace=True,
    apply_agg_norm=None,          # Protein only; None → True
    target_counts=None,           # fixed agg-norm target; None → data-driven
)
```

Analogous to `ebsd.pp.filter_cells`.  Most users should call only this function.
Returns `None` when `inplace=True`, otherwise a normalized copy.

Any other modality (for example `"OPS"`) raises `ValueError`.

**Run disambiguation.** Wells are grouped per run so that identically named wells
from different runs aren't pooled.  With `run_colname=None` the column is
auto-detected: `"run_name"` (written by the `Loader` when concatenating runs), then
`"run_id"`.  If neither exists (single-run data), wells are grouped by label alone.

### Modality-level functions

```python
ebsd.pp.normalize_transcripts(adata, *, ...)
ebsd.pp.normalize_proteins(adata, *, ...)
ebsd.pp.normalize_imaging(adata, *, inplace=True)   # raises NotImplementedError (normalize_cell_paint alias)
```

These take the same keyword arguments as `ebsd.pp.normalize` for their modality
(`normalize_transcripts` has no `apply_agg_norm`; `normalize_proteins` has no
`mor_norm_layer_name`).  Use them directly when you want to skip the dispatch layer.

### Primitive transforms (power-user API)

```python
ebsd.pp.log1p(adata, *, layer=None, snapshot_layer=None, base=None, inplace=True)
ebsd.pp.clr(adata, *, layer=None, snapshot_layer=None, base=None, batch_colname=None, inplace=True)
```

Analogous to `scanpy.pp.log1p`.  These operate directly on any layer and do
not enforce any modality contract.  `snapshot_layer` copies the input matrix to
that layer before transforming.  `clr` computes the geometric mean per batch when
`batch_colname` is given, otherwise across all features.

---

## Layer-Tracking Convention

Only the raw snapshot is written by default.  Intermediate snapshots are opt-in:
pass a layer name to `agg_norm_layer_name` / `mor_norm_layer_name`.  Pass
`raw_layer_name=None` to skip even the raw snapshot.  The default names are
shown below (the constants in `_constants.py` hold the same strings):

```python
ebsd.pp.normalize(adata, "Transcript",
                  agg_norm_layer_name="norm_agg_counts",
                  mor_norm_layer_name="norm_med_of_ratios")
```

| Constant | Default name | Written when |
|---|---|---|
| `RAW_LAYER` | `"raw"` | By default, before any transformation (`raw_layer_name`) |
| `AGG_NORM_LAYER` | `"norm_agg_counts"` | Opt-in via `agg_norm_layer_name`: after agg-count scaling, before the next step |
| `MOR_NORM_LAYER` | `"norm_med_of_ratios"` | Opt-in via `mor_norm_layer_name`: after median-of-ratios, before log1p (Transcript only) |

The final normalized result **always lands in `adata.X`** (or the `layer` kwarg
if provided).  With the default `raw_layer_name`, raw counts are recoverable via
`adata.layers["raw"]`.

---

## Pipelines by Modality

### Transcript / ThreePrimeUntargeted

Both modalities share the same pipeline and are routed through `normalize_transcripts`.
Their biological differences are handled correctly because the pipeline is batch-aware:

| Property | Targeted Transcript | 3' Untargeted |
|---|---|---|
| Gene panel | Fixed, curated (500–2000 genes) | Whole transcriptome (~20,000 genes) |
| Batch structure | Multiple reagent batches per cell (B01, B02 …) | Typically **one batch** per run |
| Sparsity | Moderate | Very high (90–95%+ zeros) |
| MT genes | Not present | Present (MT% is a real QC metric) |
| NSB / GAPDH probes | Present | Not present |

For untargeted data, raw sequencing batches are **summed before analysis** so the
AnnData already holds a single unified count matrix with no `batch_name` column in
`adata.var`.  Both `_normalize_agg_counts` and `_normalize_median_of_ratios` detect
the missing column and treat all genes as one batch, making agg-norm equivalent to
`scanpy.pp.normalize_total`.  The MOR step then computes one global well-level
scaling factor — consistent with DESeq2-style practice for scRNA-seq.

```
X (raw)
 │
 ├─► layers["raw"]                      ← snapshot (default)
 ▼
Aggregated-count scaling (per batch, per well)
 │
 ├─► layers["norm_agg_counts"]          ← snapshot (opt-in)
 ▼
Median-of-ratios correction (per batch, per well)
 │
 ├─► layers["norm_med_of_ratios"]       ← snapshot (opt-in)
 ▼
log1p
 │
 ▼
X (final)
```

**Aggregated-count scaling** divides each cell's counts for a given batch by
the cell's total counts for that batch, then multiplies by a target value.
The target is either a fixed number (`target_counts=10_000`) or, by default,
data-driven: for each batch, the median across wells of each well's median
total counts.  This normalizes for variation in sequencing depth / capture
efficiency within a batch.

**Median-of-ratios** computes a single scaling factor per well (per batch)
by comparing each well's mean expression profile to a pseudo-reference built
from the geometric mean across wells.  This corrects for systematic
well-to-well differences in signal intensity that remain after depth
normalization.

**log1p** stabilizes variance and compresses dynamic range.

### Protein

```
X (raw) → layers["raw"] (default)
 → Aggregated-count scaling → layers["norm_agg_counts"] (opt-in snapshot)
 → log1p → X (final)
```

The median-of-ratios step is **skipped** for proteins because the assumption
that most features are unchanged across wells does not hold reliably for
antibody-derived data.

Use `apply_agg_norm=False` to skip aggregated-count scaling before log1p (raw
snapshot only).  For CLR on protein data, call `ebsd.pp.clr` as a primitive
after preparing counts as needed.

### Imaging (CellPaint + IF Protein)

**Not yet implemented.**  `normalize_imaging` raises `NotImplementedError`.

Imaging tables hold **fluorescence intensities**, not counts, so the count-based pipelines
(aggregated counts → median-of-ratios → log1p) do not apply. The table mixes several feature
types, distinguishable by `adata.var`:

| feature | `var` markers | Examples |
|---|---|---|
| CellPaint structural stains | `target_type == "CellPaint"`, `CP*` batch | `Mito_RED`, `WGA_GRN`, `EBRedots_RED` |
| IF Protein markers | `target_type == "ProteinIF"`, `P*` batch | `CD8A_Rb_RED`, `phosT180-p38_Rb_GRN` |
| Background channels | `is_background == True` (`Bkg_*`), for subtraction | `Bkg_RED_CP01` |
| CellProfiler features | `feature_class == "cellprofiler"` | `Granularity_3_Nucleus`, `Texture_Entropy_..._Actin` |

An intensity-appropriate transform (and possibly background subtraction and
per-well correction) is needed; it is not implemented yet.
