# SpatialData Object Overview

---

## What is a SpatialData object?

A SpatialData object (`.zarr` store) is the primary output of the Element Biosciences
spatial profiling pipeline. It is a multi-modal container that holds all assay data,
segmentation masks, fluorescence images, and run metadata for one experiment in a single
file. It is built on the [SpatialData](https://spatialdata.scverse.org/) framework and
can be opened with:

```python
import spatialdata as sd
sdata = sd.read_zarr("/path/to/run.zarr")
```

---

## Top-Level Structure

A SpatialData object has five element namespaces. Every element name follows the pattern
`{well_id}-{element_type}-{descriptor}`.

| Namespace | Key pattern | Description |
|-----------|-------------|-------------|
| `sdata.tables` | `global-table-{Modality}` | Cell × feature count matrices (AnnData). One per assay modality plus a metadata table. |
| `sdata.points` | `{well_id}-points-{batch_name}` | Per-molecule decoded reads (Dask DataFrame). One per well × batch combination. |
| `sdata.images` | `{well_id}-images-cell_paint` | Multi-resolution fluorescence images (DataTree, `c×y×x`). One per well. |
| `sdata.labels` | `{well_id}-labels-{Cell\|Nuclear}` | Integer segmentation masks (DataTree, `y×x`). Cell and Nuclear masks per well. |
| `sdata.shapes` | `{well_id}-shapes-{Cell\|Nuclear}` | Cell boundary polygons (GeoDataFrame). Cell and Nuclear outlines per well. |

### Example run: imaging + OPS

```
SpatialData object
├── Images     ({well}-images-cell_paint)
├── Labels     ({well}-labels-Cell, {well}-labels-Nuclear)
├── Points     ({well}-points-{batch})
├── Shapes     ({well}-shapes-Cell, {well}-shapes-Nuclear)
└── Tables     global-table-imaging
               global-table-B01-Cellecta_OPS_A549
               global-table-B01-QC_Panel
               global-table-metadata
```

### Example run: legacy barcoding (12 wells)

```
SpatialData object
├── Images     12  ({well}-images-cell_paint)
├── Labels     24  ({well}-labels-Cell, {well}-labels-Nuclear)
├── Points     84  ({well}-points-{B01..B07})
├── Shapes     24  ({well}-shapes-Cell, {well}-shapes-Nuclear)
└── Tables      4  (global-table-Transcript, -Protein, -CellPaint, -metadata)
```

Legacy barcoding SDOs store imaging under `global-table-CellPaint`; the analysis
library loads it as the `imaging` modality transparently.

---

## Coordinate Systems

Every spatial element is registered in one or more coordinate systems:

| System | Elements | Description |
|--------|----------|-------------|
| `global` | All elements | Shared experiment-wide pixel space. All wells are stitched into a single coordinate frame. |
| `well_{id}` | All elements for that well | Per-well local pixel space (e.g. `well_A1`, `well_B2`). |

Coordinates in `global` space use pixel units. Physical units (µm) are stored as
separate columns: `Xum`, `Yum` in points, and tile-local `Xum_tile`, `Yum_tile` in
cell geometry. Per-cell centroids on the well canvas are in each table's
`obsm["spatial"]` (see the AnnData overview, *Cell coordinates*).

---

## Tables

Tables are AnnData objects stored at `sdata.tables[name]`. Every cell appears in every
modality table: all tables share the same `n_obs` and the same cell index
(`cell_id_global`).

Which tables are present depends on how the run was configured. A run may carry any
combination of the assay tables below, alongside a metadata table:

| Table name | Modality (target type) | Description |
|------------|------------------------|-------------|
| `global-table-ThreePrimeUntargeted` | `ThreePrimeUntargeted` (3′ transcriptome) | Cell × gene whole-transcriptome count matrix |
| `global-table-{batch}-{panel}` | `OPS` / `SpecializedTargeted` | Cell × target count matrix, one table per specialized panel (e.g. `global-table-B01-Cellecta_OPS_A549`) |
| `global-table-imaging` | `imaging` (CellPaint morphology + IF Protein) | Cell × intensity + morphology feature matrix |
| `global-table-metadata` | (none) | Run-level metadata only (no cell data); all content is in `.uns` |

The 3′ transcriptome and OPS / specialized panels are dissociation-sequencing (`DISSBatch`)
assays; each OPS panel gets its own table, and all OPS panels resolve to the `OPS` target
type. The imaging table carries CellPaint morphology and, where the panel includes them,
IF protein markers.

### Legacy barcoding tables

Older barcoding runs additionally emit `Transcript` and `Protein` tables. These modalities
are deprecated and current runs do not produce them, but the analysis library loads them
transparently for backward compatibility with older SDOs.

| Table name | Modality (target type) | Description |
|------------|------------------------|-------------|
| `global-table-Transcript` | `Transcript` (legacy) | Cell × gene count matrix |
| `global-table-Protein` | `Protein` (legacy) | Cell × antibody count matrix |

See [anndata_object_overview.md](anndata_object_overview.md) for the full field
reference for each table.

---

## Points

Points are per-molecule decoded reads stored as lazy Dask DataFrames.

**Key pattern**: `{well_id}-points-{batch_name}` (e.g. `A1-points-B02`)

**Columns** (17 total):

| Column | dtype | Description |
|--------|-------|-------------|
| `x` | int | Pixel x-coordinate in well space |
| `y` | int | Pixel y-coordinate in well space |
| `z` | int | Z-plane index (0 for 2D runs) |
| `Xum` | float32 | Physical x-coordinate in micrometres |
| `Yum` | float32 | Physical y-coordinate in micrometres |
| `Zum` | float32 | Physical z-coordinate in micrometres |
| `target_index` | int | Index into the batch's target list |
| `target_name` | str | Human-readable target name (e.g. `"AKT1"`) |
| `cell` | int | Local cell label within the tile (0 = unassigned) |
| `label_Cell` | int | Cell segmentation label within the tile |
| `label_Nuclear` | int | Nuclear segmentation label within the tile |
| `label_Cell_global` | int | Global cell ID across the whole experiment (links to `cell_id_global` in tables) |
| `is_nuclear` | bool | Whether the molecule was decoded inside the nuclear mask |
| `tile_id` | str | Tile identifier (e.g. `"L2R01C01S1"`) |
| `well_id` | str | Well identifier (e.g. `"A1"`) |
| `x_tile` | int | Pixel x-coordinate within the tile |
| `y_tile` | int | Pixel y-coordinate within the tile |

**Loading**:
```python
# Lazy (default) — recommended for large datasets
pts = sdata.points["A1-points-B02"]          # Dask DataFrame
df  = pts.compute()                           # materialise to pandas

# Via the loader API
from elembio_spatialdata_analysis.load import load_points
loaded = load_points(sdata, target_type="Transcript")
```

---

## Images

Images are multi-resolution fluorescence stacks stored as `xarray.DataTree` objects.

**Key pattern**: `{well_id}-images-cell_paint`

**Structure**: 6 pyramid levels (`scale0`–`scale5`), each with dims `(c, y, x)`.

| Property | Value |
|----------|-------|
| Channels (`c`) | `Cell-Membrane`, `Nucleus`, `Actin`, `Cell-Membrane_Nuclear`, `Nucleus_Nuclear`, `Mitochondria` |
| Full resolution (scale0) | ~14 688 × 16 488 pixels per well |
| dtype | `uint16` |
| Chunking | `(1, 918, 916)` per channel slice |

**Loading**:
```python
img = sdata.images["A1-images-cell_paint"]
full_res = img["scale0"].ds["image"]   # xarray DataArray (c, y, x)
```

---

## Labels

Labels are integer segmentation masks stored as `xarray.DataTree` objects.

**Key pattern**: `{well_id}-labels-{Cell|Nuclear}`

| Property | Value |
|----------|-------|
| Dims | `(y, x)` |
| dtype | integer |
| Pyramid levels | 6 (`scale0`–`scale5`) |
| Values | 0 = background; positive integers = cell/nucleus IDs (match `cell_id_local` in tables) |

---

## Shapes

Shapes are cell boundary polygons stored as GeoDataFrames.

**Key pattern**: `{well_id}-shapes-{Cell|Nuclear}`

| Column | Description |
|--------|-------------|
| `geometry` | Shapely `Polygon` or `MultiPolygon` in pixel coordinates (global space) |

The GeoDataFrame index corresponds to the local cell label within the well (matches
`cell_id_local` in tables).

---

## Metadata Table (`global-table-metadata`)

This table has no cell rows (`shape = 0 × 0`). All content lives in `.uns`. It is the
canonical source of run-level metadata.

### `.uns` keys

| Key | Type | Description |
|-----|------|-------------|
| `batches` | DataFrame | One row per batch. Columns: `batch_name`, `batch_class`, `batch_type`, `min_cycles`, `run_order`, `num_targets`, `pmg_mask`, `barcode_mask`, `fastq_mask`, `default_mismatch`, `num_diss_targets` |
| `wells` | DataFrame | One row per well. Index = `well_name`. Columns: `well_index`, `well_label`, `well_layout`, `cell_type`, `cell_diameter`, `target_site_names` |
| `tiles` | DataFrame | One row per tile. Index = `tile_name`. Columns: `tile_index`, `x_mm`, `y_mm`, `well_name`, `well_label` |
| `targets_barcoding` | DataFrame | **Legacy barcoding runs only.** All barcoding targets (Transcript + Protein). Columns: `target_name`, `target_id`, `target_type`, `target_sequence`, `control_type`, `batch_name`, `probe_concentration` |
| `diss_targets_metadata` | DataFrame | **Current (3′ / OPS) runs.** Dissociation-sequencing target metadata. Columns: `TargetSiteName`, `TargetMask`, `CellAssignment`, `RTOligos`, `SequencingPrimer`, `batch_name` |
| `targets_imaging` | DataFrame | All imaging targets (CellPaint + IF Protein). Columns: `target_name`, `target_id`, `target_type`, `target_sequence`, `control_type`, `batch_name`, `probe_concentration` |
| `build_provenance` | dict | Build command, pipeline version, git hash, timestamp, and build arguments |
| `missing_batches` | list[str] | Present only when a count batch declared in the panel had no points data at all, so no table exists for it. Lists those batch names. |
| `runstats_run_metadata` | DataFrame | Run IDs, analysis version, flow cell ID, run name |
| `runstats_rawread_summary` | DataFrame | Total raw read count |
| `runstats_rawread_wells` | DataFrame | Raw read counts per well |
| `runstats_demux_summary` | DataFrame | Overall demux stats: polonies, % assigned, % mismatch |
| `runstats_demux_batches` | DataFrame | Demux stats per batch |
| `runstats_demux_wells` | DataFrame | Demux stats per well × batch |
| `runstats_demux_targets` | DataFrame | Demux stats per target: polony count, % mismatch |
| `runstats_demux_unassigned_sequences` | DataFrame | Top unassigned sequences per batch |
| `runstats_cyto_summary` | DataFrame | Cytometry summary: confluency, nucleated %, cell count, counts/mm² |
| `runstats_cyto_wells` | DataFrame | Cytometry stats per well: confluency, cell count, mean counts/cell, median diameter |
| `runstats_cyto_batches` | DataFrame | Cytometry stats per well × batch: mean counts/cell, extracellular ratio, % cellular assigned |
| `runstats_cyto_control_targets` | DataFrame | Control target stats per well × batch: mean counts/cell, counts/mm² |

**Access example**:
```python
meta = sdata.tables["global-table-metadata"]
wells_df   = meta.uns["wells"]
batches_df = meta.uns["batches"]
targets_df = meta.uns["targets_barcoding"]
```

---

## Naming Conventions

| Pattern | Example | Meaning |
|---------|---------|---------|
| `{well_id}` | `A1`, `B2`, `F1` | Well position on the plate |
| `{batch_name}` | `B01`, `B02`, `CP01` | Sequencing batch identifier |
| `global-table-{Modality}` | `global-table-Transcript` | Experiment-wide cell × feature table |
| `{well}-points-{batch}` | `A1-points-B02` | Per-molecule reads for one well + batch |
| `{well}-images-cell_paint` | `A1-images-cell_paint` | Fluorescence image for one well |
| `{well}-labels-Cell` | `A1-labels-Cell` | Cell segmentation mask for one well |
| `{well}-shapes-Cell` | `A1-shapes-Cell` | Cell boundary polygons for one well |

---

## Key Linkages Between Elements

```
sdata.tables["global-table-Transcript"].obs["cell_id_global"]
    ↕  links to
sdata.points["A1-points-B02"]["label_Cell_global"]

sdata.tables["global-table-Transcript"].obs["cell_id_local"]
    ↕  links to
sdata.labels["A1-labels-Cell"]  (pixel value)
sdata.shapes["A1-shapes-Cell"]  (GeoDataFrame index)

sdata.tables["global-table-Transcript"].obs["region_key"]
    ↕  identifies which label layer the cell belongs to
    e.g. "A1-labels-Cell"

sdata.tables["global-table-metadata"].uns["wells"]
sdata.tables["global-table-metadata"].uns["batches"]
sdata.tables["global-table-metadata"].uns["targets_barcoding"]
    ↕  describe the experimental design for all tables
```

---

## Loading via the API

`Loader` is the entry point; it accepts one store or many (a path, a list of paths, or an
in-memory `SpatialData`).

```python
from elembio_spatialdata_analysis.load import Loader

loader = Loader("/path/to/run.zarr")           # or a list of stores

# Discovery
print(loader.available_tables)         # table suffixes, e.g. {'imaging', 'B01-QC_Panel'}
print(loader.available_target_types)   # e.g. {'imaging', 'OPS'}
                                       #   or {'imaging', 'ThreePrimeUntargeted'}
                                       # legacy SDOs: {'Transcript', 'Protein', 'imaging'}

# Load one table by suffix (lazy by default; pass lazy=False to materialise).
# Returns a single AnnData for one run, or concatenated across runs by default.
adata = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)

# imaging <-> CellPaint is aliased: this finds global-table-imaging on current
# stores and global-table-CellPaint on legacy stores.
imaging = loader.load_tables(table="imaging", lazy=False)

# Per-run instead of concatenated -> dict keyed by run
per_run = loader.load_tables(table="imaging", concatenated=False, lazy=False)

# Load points, filtered by target type / batch / well / tile
points = loader.load_points(target_type="ThreePrimeUntargeted")

# Write per-run tables to disk (h5ad) with provenance
paths = loader.export_tables("/out/dir", fmt="h5ad")
```

OPS/specialized panels are one table per batch × panel; discover their suffixes with
`run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))` on a `RunInfo`
(see the load [SPECIFICATION.md](../src/elembio_spatialdata_analysis/load/SPECIFICATION.md)).
