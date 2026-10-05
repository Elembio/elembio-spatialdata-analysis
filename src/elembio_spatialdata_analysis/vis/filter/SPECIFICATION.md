# `vis.filter` — Filtering QC Visualization

## Overview

The `vis.filter` subpackage provides all plotting functions for per-well QC
diagnostics and filter assessment.  It mirrors the structure of `pp.filter`:
private engine modules contain no public API; all user-facing functions live in
`plots.py` and are re-exported at both `ebsd.vis.filter.*` and `ebsd.vis.*`.

## Package structure

```
vis/filter/
├── __init__.py       Public re-exports (6 functions + 2 aliases)
├── _engine.py        Histogram engine, layout helpers, metric registry
├── _violin.py        Violin rendering engine (shared by batch + CellPaint plots)
├── plots.py          All public plot functions: plot_qc_per_well, plot_violin_per_well,
│                     plot_filtered_histograms_per_well, plot_total_counts_by_batch,
│                     plot_imaging_intensities, plot_imaging_morphology
└── SPECIFICATION.md  This file
```

## Prerequisite

Most functions require QC metrics to be pre-computed:

```python
ebsd.pp.compute_qc_metrics(adata, "Transcript")  # or "Protein", "Imaging", etc.
```

This stores per-cell metrics in `adata.obsm['qc_metrics']`.  Attempting to plot
a metric that is not there raises a clear `KeyError` with a helpful message.

## Public API

### `plot_qc_per_well(adata, metric, ...)`

Single entry point for all per-well QC histogram plots.

```python
ebsd.vis.plot_qc_per_well(adata, "total_counts")
ebsd.vis.plot_qc_per_well(adata, "Area")  # from obsm['cell_geometry']
ebsd.vis.plot_qc_per_well(adata, "nuclear_frac")  # any qc_metrics column
ebsd.vis.plot_qc_per_well(adata, "my_custom_col")  # falls back to obs lookup
```

**`metric` resolution order:**

1. `"total_counts"` — row sums of `layer` when `layer` is given; otherwise
   `adata.obsm['qc_metrics']['total_counts']` if present (avoids touching a
   lazy `X`), else row sums of `adata.X`.
2. Column in `adata.obsm['qc_metrics']` (requires `compute_qc_metrics` first).
3. Column in `adata.obs`.
4. Column in `adata.obsm['cell_geometry']` (e.g. `"Area"`).
5. `KeyError` with a clear message.

**Built-in recognised metric names** (auto x-label + display options):

| `metric` | X-axis label | clip [0,1] |
|---|---|---|
| `"total_counts"` | Total counts per cell | no |
| `"total_assignment_rate"` | Assignment rate | yes |
| `"total_nsb_counts"` | NSB counts per cell | no |
| `"total_gapdh_counts"` | GAPDH counts per cell | no |
| `"nuclear_frac"` | Nuclear assignment rate | yes |
| `"n_genes_by_counts"` | Number of genes detected per cell | no |
| `"pct_counts_mt"` | Mitochondrial gene percentage | yes |
| `"Area"` | Cell area (px²) | no |

Any other string falls through the resolution order and uses the metric name
itself as the x-axis label.

**Key parameters:**

| Parameter | Default | Description |
|---|---|---|
| `metric` | `"total_counts"` | What to plot (see table above) |
| `layer` | `None` | Expression layer to sum for `metric="total_counts"` |
| `well_colname` | `"Well"` | `adata.obs` column for well labels |
| `wells` | `None` | Explicit well order |
| `orient_flowcell` | `None` | Flowcell layout (`"12_well"` / `"6_well"`); auto-detected from the well set when possible |
| `bins` | `40` | Histogram bins; `"fd"` for Freedman-Diaconis |
| `cols` | `4` | Max subplot columns |
| `logx` | `False` | Log-scale x axis |
| `show_threshold` | `False` | Draw a fixed threshold line |
| `threshold` | `None` | Threshold value |
| `show_percentiles` | `False` | Draw lower/upper percentile lines |
| `lower_pct` / `upper_pct` | `0.025` / `0.975` | Percentile bounds |
| `show_percentile` / `max_pct` | `False` / `0.99` | Draw a single upper-tail percentile line (distinct from `show_percentiles`) |

---

### `plot_violin_per_well(adata, metric, ...)`

One violin per well for any QC metric.  `metric` uses the same resolution
order as `plot_qc_per_well` (`obsm['qc_metrics']` → `obs` → `obsm['cell_geometry']`).

```python
ebsd.vis.plot_violin_per_well(adata, "pct_counts_mt")
ebsd.vis.plot_violin_per_well(adata, "Area", scale=pixel_size_um**2, ylabel="Cell area (µm²)")
```

**Key parameters:**

| Parameter | Default | Description |
|---|---|---|
| `metric` | *(required)* | What to plot (same resolution as `plot_qc_per_well`) |
| `well_colname` | `"Well"` | `adata.obs` column for well labels |
| `wells` | `None` | Explicit well order |
| `scale` | `1.0` | Multiply values before plotting (e.g. px² → µm²); set `ylabel` to match |
| `ylabel` / `title` | `None` | Axis label / title overrides |
| `figsize` | `(5.5, 3.4)` | Figure size in inches |
| `color` / `alpha` | `"#4C78A8"` / `0.55` | Violin fill |
| `show_median` | `True` | Draw each well's median |
| `hline` / `hline_label` / `hline_color` | `None` / `None` / `"#e53e3e"` | Optional absolute reference line |
| `dpi` | `150` | Figure resolution |

---

### `plot_filtered_histograms_per_well(adata, flags, ...)`

Overlays raw (all cells) and filtered (passing cells) distributions per well.

```python
diag = ebsd.pp.filter.compute_filter_diagnostics(adata, "Transcript", return_diagnostics=True)
flags = diag.flags
ebsd.vis.plot_filtered_histograms_per_well(adata, flags, metric="total_counts")
ebsd.vis.plot_filtered_histograms_per_well(adata, flags, metric="Area")
```

**Key parameters:**

| Parameter | Default | Description |
|---|---|---|
| `metric` | `"total_counts"` | What to plot (same resolution as `plot_qc_per_well`) |
| `obs_col` | `None` | Explicit `adata.obs` column when `metric` isn't found in `obsm['qc_metrics']` |
| `layer` | `None` | Expression layer for `metric="total_counts"` |
| `flag_col` | `"passes_all_filters"` | Boolean column in `flags` |
| `well_col` | `"Well"` | `adata.obs` column for well labels |
| `wells` | `None` | Explicit well order |
| `bins` / `cols` | `40` / `6` | Histogram bins / max subplot columns |
| `figsize` | `(18.0, 6.0)` | Figure size in inches |
| `colors` | `None` | `{"original": ..., "filtered": ...}` colour overrides |
| `show_raw` | `True` | Show raw (unfiltered) distribution as background |

---

### `plot_total_counts_by_batch(adata, ...)`

Violin plots of per-cell total counts, one violin per batch per well.  Reads
pre-aggregated counts from `adata.obsm['total_counts']` when that key exists
(avoids full X access); a batch with no column there is skipped.  Only when
`obsm['total_counts']` is absent entirely does it sum `adata.X` per batch.

```python
ebsd.vis.plot_total_counts_by_batch(adata)
ebsd.vis.plot_total_counts_by_batch(adata, show_percentiles=True, lower_pct=0.025, upper_pct=0.975)
```

**Key parameters:**

| Parameter | Default | Description |
|---|---|---|
| `batch_colname` | `"batch_name"` | `adata.var` column for batch labels |
| `well_colname` | `"Well"` | `adata.obs` column for well labels |
| `show_percentiles` | `False` | Draw per-violin threshold lines |
| `lower_pct` / `upper_pct` | `0.025` / `0.975` | `QuantileSpec` (scalar or per-well/batch dict) |
| `ylim_percentile` | `99.5` | Upper y-axis percentile cap (`None` = no cap) |
| `logy` | `False` | Log-scale y axis |
| `orient_flowcell` | `None` | Flowcell layout key (e.g. `"12_well"`) |
| `max_violin_points` | `2000` | Downsample for KDE performance |

---

### `plot_imaging_intensities(adata, ...)`

Violin plots of the structural imaging (CellPaint) channel intensities per well.
Channels are selected by `tl.morphology.select_morphology_channels` (the CellPaint
`CP*` batch), so the protein/antibody panel, CellProfiler features, compartment
variants, and backgrounds are excluded automatically. Pass `channels=[...]` to
plot an explicit set.

Alias: `plot_cell_paint_intensities`.

```python
ebsd.vis.plot_imaging_intensities(adata)
```

**Key parameters:**

| Parameter | Default | Description |
|---|---|---|
| `channels` | `None` | Explicit `var_names` to plot; overrides batch-anchored selection |
| `target_type_value` | `None` | Deprecated: ignored, and passing a value warns |
| `show_percentiles` | `False` | Draw per-violin threshold lines |

---

### `plot_imaging_morphology(adata, metric, ...)`

Violin plots of a single imaging (CellPaint) morphology metric, one violin per
channel per well.  Feature-to-channel mapping uses the last `_`-delimited token of
each matching `var_name` as the channel label.

Alias: `plot_cell_paint_morphology`.

```python
ebsd.vis.plot_imaging_morphology(adata, "MeanIntensity")
ebsd.vis.plot_imaging_morphology(adata, "Granularity_3", channels=["Actin", "Nucleus"])
```

**Key parameters:**

| Parameter | Default | Description |
|---|---|---|
| `metric` | *(required)* | Case-insensitive substring match against `var_names` |
| `channels` | `None` | Restrict to specific channel suffixes |
| `max_violin_points` | `2000` | Downsample for KDE performance |

---

## `QuantileSpec` format

`lower_pct` and `upper_pct` in the violin functions accept a `QuantileSpec`:

```python
# Scalar — same threshold for all wells and batches
lower_pct = 0.025

# Per-well dict — different thresholds per well
lower_pct = {"A1": 0.01, "B1": 0.025, "__default__": 0.025}

# Per-well-per-batch nested dict (plot_total_counts_by_batch)
lower_pct = {
    "A1": {"total_counts_B02": 0.01, "__default__": 0.025},
    "__default__": 0.025,
}
```

Inner keys are `total_counts_{batch}` for `plot_total_counts_by_batch` and the
channel `var_names` for `plot_imaging_intensities`.

## Private module responsibilities

### `_engine.py`

- `FLOWCELL_ORIENTATIONS` — registry of named flowcell layouts
- `_resolve_well_layout` / `_prepare_well_grid` — well ordering and grid sizing
- `_plot_per_well_histogram` — core histogram rendering engine
- `_resolve_qc_values` — metric resolution (obsm → obs → cell_geometry)
- `_QC_METRIC_META` — display metadata registry; extend here to add new metrics

**To add a new built-in metric:**

```python
# in _engine.py
_QC_METRIC_META["my_new_metric"] = ("My Metric Label", False, False)  # (label, clip_0_1, default_logx)
```

### `_violin.py`

- `_render_violin_grid` — shared violin rendering for batch and CellPaint plots
- `_normalize_quantile_value` / `_normalize_quantile_spec` — `QuantileSpec` normalisation
