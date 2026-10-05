# AnnData Object Overview

The tables in a SpatialData object (`sdata.tables`) are standard
[AnnData](https://anndata.readthedocs.io/) objects. Each table is a **cell × feature**
matrix where every row is one segmented cell and every column is one measured feature
(gene, protein, or morphology measurement).

---

## Modalities

Which tables a run contains depends on how it was configured. Current runs draw from three
assays: the 3′ transcriptome (`ThreePrimeUntargeted`), OPS / specialized panels
(`SpecializedTargeted`), and imaging (CellPaint morphology, optionally with IF protein).
The 3′ and OPS tables are dissociation-sequencing (`DISSBatch`) count matrices. Every run
also has a `metadata` table.

`Transcript` and `Protein` barcoding are deprecated modalities. They share the cell × gene
count-matrix layout of the 3′ transcriptome table and remain supported for reading older
SDOs, but current runs do not emit them.

## Shared Structure (All Modality Tables)

All modality tables share the same `.obs` columns and the same `.uns` keys for
SpatialData linkage. Differences are in `.var`, `.obsm`, `.varm`, and `.layers`.

### `.obs` — Cell Annotations (shared across all tables)

Every row is one cell. The same cells appear in every table in a run.


| Column           | dtype    | Description                                                                                                    |
| ---------------- | -------- | -------------------------------------------------------------------------------------------------------------- |
| `Cell`           | uint32   | Local cell label within the tile (matches pixel value in `{well}-labels-Cell`)                                 |
| `Well`           | str      | Well identifier (e.g. `"A1"`, `"B2"`)                                                                          |
| `WellLabel`      | str      | Human-readable well label / condition name (e.g. `"Hela"`)                                                     |
| `tile_id`        | str      | Tile identifier where this cell was segmented (e.g. `"L2R01C01S1"`)                                            |
| `X_tile`         | float32  | Cell centroid x-coordinate in pixels, **tile-local** (relative to the cell's tile)                             |
| `Y_tile`         | float32  | Cell centroid y-coordinate in pixels, **tile-local** (relative to the cell's tile)                             |
| `cell_id_local`  | int32    | Cell ID unique within its well (matches GeoDataFrame index in `{well}-shapes-Cell`)                            |
| `cell_id_global` | int32    | Cell ID unique across the entire experiment (links to `label_Cell_global` in points)                           |
| `region_key`     | category | Name of the label layer this cell belongs to (e.g. `"A1-labels-Cell"`) — used by SpatialData for spatial joins |
| `well_id`        | str      | Well identifier (duplicate of `Well`; kept for SpatialData compatibility)                                      |

#### Cell coordinates

Each cell's centroid is stored in three places, in different frames:

| Location | Columns | Frame |
|---|---|---|
| `obs` | `X_tile`, `Y_tile` | Pixels, relative to the cell's tile |
| `obsm["cell_geometry"]` | `Xum_tile`, `Yum_tile` | Micrometres, relative to the cell's tile |
| `obsm["spatial"]` | `X`, `Y` | Pixels on the stitched **well canvas**, the same frame as the well's images, labels, shapes, and points. In tile-mode stores these equal `X_tile` / `Y_tile`. |

Use `obsm["spatial"]` for spatial plots and neighborhood analysis (it is what
Squidpy and other AnnData spatial tools read). It is a DataFrame, not an array;
call `.to_numpy()` if a tool needs an `ndarray`. Cells whose tile could not be
placed on the well canvas have `NaN` coordinates.

Older stores use `obs["X"]` / `obs["Y"]` and `cell_geometry["Xum"]` / `["Yum"]`
for the same tile-local values, and have no `obsm["spatial"]`.


### `.uns` — Unstructured Metadata (shared across all tables)


| Key                 | Type | Description                                                                                                                                                                                               |
| ------------------- | ---- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `spatialdata_attrs` | dict | SpatialData linkage: `region_key` (column name), `region` (list of label layer names), `instance_key` (`"cell_id_global"`)                                                                                |
| `metadata_attrs`    | dict | Modality metadata: `batch_name` (array of batch names in this table), `batch_class` (`"DISSBatch"`, `"ImagingBatch"`, or legacy `"BarcodingBatch"`), `target_type` (e.g. `"ThreePrimeUntargeted"`), `merge_strategy` (`"expand_features"`) |
| `spatial_artifact_params` | dict | Present when `obsm["spatial_artifacts"]` was written: `version`, `spatial_method`, `detection_enabled`, `detection_disabled_reason`, `rollup_columns` (`imaging` → `any-imaging`, `ThreePrimeUntargeted` → `any-ThreePrimeUntargeted`, `OPS` → `any-OPS`, `all` → `any`), plus per-modality circle/flag params. |

### Spatial artifacts

`obsm["spatial_artifacts"]` is a boolean DataFrame (one row per cell) written by
the pipeline onto **every** cell-aligned table. `True` means the cell
sits in a detected artifact region. Column names:

| Pattern | Meaning | Example |
|---|---|---|
| `imaging-{target}` | Imaging target unit | `imaging-Actin` |
| `ThreePrimeUntargeted-{batch}` | 3′ batch unit | `ThreePrimeUntargeted-B01` |
| `OPS-{batch}` | OPS batch unit | `OPS-B03` |
| `any-imaging` / `any-ThreePrimeUntargeted` / `any-OPS` | OR of that modality's units | |
| `any` | OR of the modality rollups | |

Absent on stores built without spatial-artifact detection and on older stores.
When detection was disabled for a run, the full column set is present but every
cell is `False` (see `uns["spatial_artifact_params"]["detection_enabled"]`).

`ebsd.pp.filter_cells` does not apply this gate by default. Pass
`filter_spatial_artifacts=True` to drop cells in the **modality rollup**
(`any-imaging` / `any-ThreePrimeUntargeted` / `any-OPS`).
`filter_modalities(..., filter_spatial_artifacts=True)` applies each table's
own rollup, then intersects; that intersection equals `~any` when every
present artifact-bearing modality is in the call.

---

## `global-table-ThreePrimeUntargeted` (3′ transcriptome)

**Shape**: `n_obs × n_vars` = cells × genes (whole transcriptome).

The 3′ transcriptome is a dissociation-sequencing (`DISSBatch`) assay. `.X` holds
whole-cell transcript counts per gene, with gene-level annotations in `.var`, matching the
cell × gene count layout of the legacy
[`global-table-Transcript`](#global-table-transcript-legacy-barcoding) table. Its per-cell
`.obsm` follows the DISS convention shared with the OPS panels below, and DISS target
design is recorded in `metadata.uns["diss_targets_metadata"]`.

---

## OPS and specialized panels

OPS / `SpecializedTargeted` panels are targeted dissociation-sequencing (`DISSBatch`)
count matrices. Each panel is a separate table named `global-table-{batch}-{panel}` (e.g.
`global-table-B01-A549`, `global-table-B01-QC_Panel`), and all normalize to
the `OPS` target type. Each shares the common `.obs` and `.uns` linkage structure above.

### `.var`: OPS Target Annotations

| Column             | dtype | Description                                             |
| ------------------ | ----- | ------------------------------------------------------- |
| `target_id`        | int   | Numeric target ID from the panel design                 |
| `target_gene`      | str   | Gene the target maps to                                 |
| `target_sequence`  | str   | Barcode sequence used to identify this target           |
| `control_type`     | str   | Empty for real targets; populated for controls          |
| `batch_name`       | str   | Specialized batch this target belongs to (e.g. `"B01"`) |
| `target_site_name` | str   | Target site / probe-site identifier                     |

### `.obsm`

| Key                 | Description                                                                     |
| ------------------- | ------------------------------------------------------------------------------- |
| `total_counts`      | Total assigned counts per cell per batch                                        |
| `target_assignment` | Per-cell target-assignment detail                                               |
| `data_missing`      | Boolean flag per batch: whether data is missing for this cell in that batch     |
| `cell_geometry`     | Cell morphology measurements (same 23 fields as every table)                    |
| `spatial`           | Cell centroids (`X`, `Y`) in well-canvas pixels (see [Cell coordinates](#cell-coordinates)) |
| `spatial_artifacts` | Boolean artifact membership (same frame on every cell-aligned table; see below) |

OPS tables have no `nuclear` layer. DISS target design lives in
`metadata.uns["diss_targets_metadata"]` (see the
[SpatialData overview](spatialdata_object_overview.md#metadata-table-global-table-metadata)).

---

## `global-table-imaging`

**Shape**: `n_obs × n_vars` = cells × imaging features (e.g. 153 770 × 558).
Legacy stores may name this `global-table-CellPaint`; both load transparently.

### `.X` — Main Count Matrix


| Property | Value                                                         |
| -------- | ------------------------------------------------------------- |
| dtype    | float32                                                       |
| Contents | Per-channel intensities (CellPaint stains + IF Protein markers) and CellProfiler-derived morphology features |
| Layers   | None (no `nuclear` layer)                                     |


### `.var` — Feature Annotations


| Column                | dtype   | Description                             |
| --------------------- | ------- | --------------------------------------- |
| `target_id`           | float64 | Numeric channel ID                      |
| `target_type`         | str     | Per feature: `"CellPaint"` (stains) or `"ProteinIF"` (IF markers); NaN for background / CellProfiler |
| `feature_class`       | str     | `"base"` (raw intensity) or `"cellprofiler"` (derived metric) |
| `suffix`              | str     | `""` (whole-cell) or `_Nuclear` / `_Cytoplasmic` / threshold variant |
| `is_background`       | bool    | `True` for per-color background channels (`Bkg_*`) |
| `batch_name`          | str     | Imaging batch: `CP*` = CellPaint stains, `P*` = IF Protein |
| `target_sequence`     | str     | Empty (not applicable for imaging)      |
| `control_type`        | str     | Empty                                   |
| `probe_concentration` | float64 | Stain concentration                     |


`**.var_names`**: feature names combining channel and measurement type
(e.g. `"Cell-Membrane"`, `"Nucleus"`, `"Actin"`, `"Cell-Membrane_Nuclear"`)

### `.obsm`


| Key             | Shape       | Description                                                |
| --------------- | ----------- | ---------------------------------------------------------- |
| `cell_geometry` | (n_obs, 23) | Cell morphology measurements. Same fields as other tables. |
| `spatial` | (n_obs, 2) | Cell centroids (`X`, `Y`) in well-canvas pixels (see [Cell coordinates](#cell-coordinates)). |
| `spatial_artifacts` | (n_obs, run-specific) | Boolean artifact membership. Same columns on every table in the run. |


### `.varm`

None.

---

## `global-table-Transcript` (legacy barcoding)

> **Legacy / deprecated.** Emitted only by older barcoding runs. The
> current 3′ transcriptome (`ThreePrimeUntargeted`) table follows the same cell × gene
> count-matrix layout; see [`global-table-ThreePrimeUntargeted`](#global-table-threeprimeuntargeted-3-transcriptome).

**Shape**: `n_obs × n_vars` = cells × genes (e.g. 153 770 × 344)

### `.X` — Main Count Matrix


| Property | Value                                                                |
| -------- | -------------------------------------------------------------------- |
| dtype    | float64                                                              |
| Contents | Whole-cell transcript counts per gene (cytoplasm + nucleus combined) |
| Sparsity | Sparse (CSR)                                                         |


### `.var` — Gene Annotations


| Column                | dtype   | Description                                                             |
| --------------------- | ------- | ----------------------------------------------------------------------- |
| `target_id`           | int64   | Numeric target ID from the panel design                                 |
| `target_type`         | str     | Always `"Transcript"` for this table                                    |
| `target_sequence`     | str     | Barcode sequence used to identify this target                           |
| `control_type`        | str     | Empty string for real genes; e.g. `"Negative Control 1"` for NSB probes |
| `batch_name`          | str     | Sequencing batch this target belongs to (e.g. `"B02"`, `"B03"`)         |
| `probe_concentration` | float64 | Probe concentration in nM                                               |


`**.var_names`** (index): gene/target names (e.g. `"AIP"`, `"ATF2"`, `"AURKB"`)

### `.layers`


| Layer     | dtype   | Description                                                                                                                            |
| --------- | ------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| `nuclear` | float64 | Nuclear-only transcript counts (reads decoded inside the nuclear mask). Complement of cytoplasmic counts: `cytoplasmic = X - nuclear`. |


### `.obsm` — Per-Cell DataFrames


| Key                 | Shape                  | Description                                                                                                                                                                  |
| ------------------- | ---------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `total_counts`      | (n_obs, 2 × n_batches) | Total assigned counts per cell per batch. Columns: `{batch_name}` (whole-cell) and `{batch_name}_nuclear` (nuclear only). e.g. `B02`, `B02_nuclear`, `B03`, `B03_nuclear`, … |
| `unassigned_counts` | (n_obs, 2 × n_batches) | Counts that decoded to a valid barcode but were not assigned to any target. Columns: `{batch_name}_Unassigned`, `{batch_name}_Unassigned_nuclear`.                           |
| `control_variables` | (n_obs, 4 × n_batches) | Per-cell counts for NSB (non-specific binding) negative control probes. Columns: `{batch_name}_NSB_1` … `{batch_name}_NSB_4`. Used as a proxy for background signal.         |
| `data_missing`      | (n_obs, n_batches)     | Boolean flag per batch indicating whether data is missing for this cell in that batch (e.g. cell falls outside the imaged area for that batch).                              |
| `cell_geometry`     | (n_obs, 23)            | Cell morphology measurements. See [Cell Geometry Fields](#cell-geometry-fields) below.                                                                                       |
| `spatial`           | (n_obs, 2)             | Cell centroids (`X`, `Y`) in well-canvas pixels. See [Cell coordinates](#cell-coordinates).                                                                                  |
| `spatial_artifacts` | (n_obs, run-specific)  | Boolean artifact membership. Same frame on every cell-aligned table; see [Spatial artifacts](#spatial-artifacts).                                                            |


### `.varm` — Per-Gene DataFrames


| Key                    | Shape             | Description                                                                                                                                                 |
| ---------------------- | ----------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `extracellular_counts` | (n_vars, n_wells) | Total counts for each gene that decoded outside any cell boundary, per well. Columns are well IDs (e.g. `"A1"`, `"B2"`). Useful for estimating ambient RNA. |


### `.uns` (Transcript-specific)


| Key                          | Type | Description                                                                                                                                        |
| ---------------------------- | ---- | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| `control_variables_metadata` | dict | Metadata for the NSB control probes: arrays of `target_id`, `target_name`, `target_sequence`, `control_type`, `batch_name`, `probe_concentration`. |


---

## `global-table-Protein` (legacy barcoding)

> **Legacy / deprecated.** Emitted only by older barcoding runs. IF protein
> markers instead live in the `imaging` table (`target_type == "ProteinIF"`).

**Shape**: `n_obs × n_vars` = cells × antibody targets (e.g. 153 770 × 50)

### `.X` — Main Count Matrix


| Property | Value                                                                |
| -------- | -------------------------------------------------------------------- |
| dtype    | float64                                                              |
| Contents | Whole-cell antibody counts per target (cytoplasm + nucleus combined) |


### `.var` — Target Annotations


| Column                | dtype   | Description                                    |
| --------------------- | ------- | ---------------------------------------------- |
| `target_id`           | int64   | Numeric target ID                              |
| `target_type`         | str     | Always `"Protein"`                             |
| `target_sequence`     | str     | 9-mer barcode sequence                         |
| `control_type`        | str     | Empty for real targets; populated for controls |
| `batch_name`          | str     | Batch this target belongs to (e.g. `"B01"`)    |
| `probe_concentration` | float64 | Antibody probe concentration in nM             |


`**.var_names`**: antibody/protein names (e.g. `"AKT1"`, `"phosS473-AKT1"`, `"CREB"`)

### `.layers`


| Layer     | dtype   | Description                  |
| --------- | ------- | ---------------------------- |
| `nuclear` | float64 | Nuclear-only antibody counts |


### `.obsm`


| Key                 | Shape       | Description                                                                                                    |
| ------------------- | ----------- | -------------------------------------------------------------------------------------------------------------- |
| `total_counts`      | (n_obs, 2)  | Columns: `{batch_name}` (whole-cell total), `{batch_name}_nuclear` (nuclear total). e.g. `B01`, `B01_nuclear`. |
| `unassigned_counts` | (n_obs, 2)  | Columns: `Unassigned` (whole-cell), `Unassigned_nuclear`.                                                      |
| `data_missing`      | (n_obs, 1)  | Boolean flag: whether data is missing for this cell in the protein batch.                                      |
| `cell_geometry`     | (n_obs, 23) | Cell morphology measurements. Same fields as Transcript table.                                                 |
| `spatial`           | (n_obs, 2)  | Cell centroids (`X`, `Y`) in well-canvas pixels. See [Cell coordinates](#cell-coordinates).                    |
| `spatial_artifacts` | (n_obs, *)  | Boolean artifact membership. Same frame as other tables.                                                       |


### `.varm`


| Key                    | Shape             | Description                                        |
| ---------------------- | ----------------- | -------------------------------------------------- |
| `extracellular_counts` | (n_vars, n_wells) | Extracellular antibody counts per target per well. |


---

## Cell Geometry Fields

`obsm["cell_geometry"]` is a DataFrame present in **every modality table** with
identical values. It contains CellProfiler-derived shape measurements for each cell.


| Column                         | dtype   | Description                                                      |
| ------------------------------ | ------- | ---------------------------------------------------------------- |
| `Area`                         | uint32  | Cell area in pixels                                              |
| `NuclearArea`                  | uint32  | Nuclear area in pixels (0 if no nucleus detected)                |
| `Xum_tile`                     | float32 | Cell centroid x-coordinate in micrometres, tile-local            |
| `Yum_tile`                     | float32 | Cell centroid y-coordinate in micrometres, tile-local            |
| `AreaUm`                       | float32 | Cell area in µm²                                                 |
| `NuclearAreaUm`                | float32 | Nuclear area in µm²                                              |
| `AreaShape_Compactness`        | float32 | Ratio of perimeter² to area (higher = less compact)              |
| `AreaShape_ConvexArea`         | float32 | Area of the convex hull in pixels                                |
| `AreaShape_Eccentricity`       | float32 | Eccentricity of the best-fit ellipse (0 = circle, 1 = line)      |
| `AreaShape_EquivalentDiameter` | float32 | Diameter of a circle with the same area                          |
| `AreaShape_EulerNumber`        | float32 | Topological measure (1 = no holes)                               |
| `AreaShape_Extent`             | float32 | Ratio of cell area to bounding box area                          |
| `AreaShape_FormFactor`         | float32 | 4π × area / perimeter² (1 = perfect circle)                      |
| `AreaShape_MajorAxisLength`    | float32 | Length of the major axis of the best-fit ellipse (pixels)        |
| `AreaShape_MaxFeretDiameter`   | float32 | Maximum caliper diameter                                         |
| `AreaShape_MaximumRadius`      | float32 | Maximum distance from centroid to boundary                       |
| `AreaShape_MeanRadius`         | float32 | Mean distance from centroid to boundary                          |
| `AreaShape_MedianRadius`       | float32 | Median distance from centroid to boundary                        |
| `AreaShape_MinFeretDiameter`   | float32 | Minimum caliper diameter                                         |
| `AreaShape_MinorAxisLength`    | float32 | Length of the minor axis of the best-fit ellipse (pixels)        |
| `AreaShape_Orientation`        | float32 | Angle of the major axis relative to x-axis (degrees, −90 to +90) |
| `AreaShape_Perimeter`          | float32 | Perimeter of the cell boundary in pixels                         |
| `AreaShape_Solidity`           | float32 | Ratio of cell area to convex hull area (1 = fully convex)        |


**Access**:

```python
adata = sdata.tables["global-table-Transcript"]
geom  = adata.obsm["cell_geometry"]   # pandas DataFrame, shape (n_obs, 23)
```

To treat geometry as expression-like variables (same columns appended to `X` /
`var`), load with the flag:

```python
from elembio_spatialdata_analysis.load import Loader

loader = Loader("/path/to/run.zarr")
adata = loader.load_tables(table="Transcript", include_cell_geometry=True, lazy=False)
# Geometry columns (Area, AreaShape_*, …) appear in adata.var alongside genes
```

Geometry is **not** merged into `.obs` by the loader. For a flat `obs` column,
copy or join from `obsm["cell_geometry"]` (e.g. `adata.obs["Area"] = adata.obsm["cell_geometry"]["Area"]`).

---

## `global-table-metadata`

**Shape**: `0 × 0` — no cell rows, no feature columns. All content is in `.uns`.

See the [SpatialData Object Overview](spatialdata_object_overview.md#metadata-table-global-table-metadata)
for the full `.uns` key reference.

**Quick access**:

```python
meta = sdata.tables["global-table-metadata"]

wells    = meta.uns["wells"]               # DataFrame: well layout, cell type, etc.
batches  = meta.uns["batches"]             # DataFrame: batch class, type, cycles, etc.
targets  = meta.uns["targets_imaging"]     # DataFrame: imaging targets (current + legacy runs)
# DISS runs: meta.uns["diss_targets_metadata"];  legacy barcoding: meta.uns["targets_barcoding"]
tiles    = meta.uns["tiles"]               # DataFrame: tile positions in mm
```

---

## Common Access Patterns

```python
import spatialdata as sd

sdata = sd.read_zarr("/path/to/run.zarr")
adata = sdata.tables["global-table-Transcript"]

# Cell × gene matrix
X = adata.X                          # sparse float64

# Nuclear counts only
nuclear = adata.layers["nuclear"]    # sparse float64

# Cytoplasmic counts
cytoplasmic = adata.X - adata.layers["nuclear"]

# Cell metadata
obs = adata.obs                      # pandas DataFrame
wells_in_data = obs["Well"].unique() # e.g. ['A1', 'A2', ...]

# Gene metadata
var = adata.var                      # pandas DataFrame
batches = var["batch_name"].unique() # e.g. ['B02', 'B03', ...]

# Total counts per cell per batch
total = adata.obsm["total_counts"]   # DataFrame, cols: B02, B02_nuclear, ...

# NSB control counts (background proxy)
nsb = adata.obsm["control_variables"]  # DataFrame, cols: B02_NSB_1, ...

# Extracellular (ambient) counts per gene per well
ambient = adata.varm["extracellular_counts"]  # DataFrame, cols: A1, A2, ...

# Cell geometry
geom = adata.obsm["cell_geometry"]   # DataFrame, 23 morphology columns
```

