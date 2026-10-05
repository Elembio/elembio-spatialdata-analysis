# `load/` Module Specification

This document specifies the architecture, data model, and behavior of the
`load/` module -- the unified loading API for SpatialData zarr stores.

---

## 1. Module Structure

```
load/
+-- __init__.py        # Public exports: Loader, RunInfo, inputs
+-- loader.py          # Loader class, multi-store orchestration, table/points loading, export
+-- backends.py        # RunInfo ABC, concrete backends, factory function
+-- merge.py           # AnnData and Polars concatenation utilities
+-- inputs/            # Typed parsers for metadata/inputs/* sidecars
|   +-- __init__.py
|   +-- bases2fastq.py            # Bases2FastqStats from RunStats.json
|   +-- cells2stats.py            # Cells2StatsStats + AverageNormWellStats CSV +
|   |                             # OPS-only SpecializedTargetedCounts (pre-rc8:
|   |                             # TargetCounts) / TargetCellAssignmentManifest
|   +-- fastq2parquet_3p.py       # Fastq2Parquet3pStats + thin multiqc summary
|   +-- pipeline_provenance.py    # PipelineProvenance from tetonatlas/Provenance.json
|   +-- run_inputs.py             # RunParameters / Panel / RunManifest aggregate
|   +-- witness.py                # witness_hashes.json
+-- SPECIFICATION.md   # This file
```

### Public API

| Symbol                          | Type     | Purpose                                                     |
|---------------------------------|----------|-------------------------------------------------------------|
| `Loader`                        | class    | Multi-run data loader (main entrypoint).                    |
| `RunInfo`                       | ABC      | Per-run metadata, point selection, and IO.                  |
| `load.inputs.parse_*`           | functions| Read parsed metadata sidecars from a `RunInfo` backend.     |
| `load.inputs.{Bases2FastqStats, Cells2StatsStats, Fastq2Parquet3pStats, PipelineProvenance, RunInputs, PanelInfo}` | dataclasses | Typed views with raw-dict fallbacks. |

The `inputs/` subpackage delegates IO to `RunInfo._read_json` so it works transparently for local directories, `.zarr.zip` archives, and S3 stores.  Parsers return `None` when the source file is absent — callers should treat absence as "this section does not apply".

All other symbols are private (`_`-prefixed).

---

## 2. Supported Store Inputs

The `Loader` accepts any of the following as a store input (`StoreInput`):

| Input type                        | Backend selected    | Example                                           |
|-----------------------------------|---------------------|---------------------------------------------------|
| Local `.zarr` directory           | `_LocalDirBackend`  | `"/data/RUN-0001.zarr"`                           |
| Local `.zarr.zip` file            | `_LocalZipBackend`  | `"/data/RUN-0001.zarr.zip"`                       |
| S3 `.zarr` directory              | `_S3DirBackend`     | `"s3://bucket/spatialdata_runs/RUN-0001.zarr"`    |
| S3 `.zarr.zip` file               | `_S3ZipBackend`     | `"s3://bucket/spatialdata_runs/RUN-0001.zarr.zip"` |
| In-memory `SpatialData` object    | `_SDataBackend`     | `sd.read_zarr("RUN-0001.zarr")`                   |

Backend selection is determined by `make_backend()` based on the input type
(SpatialData object), string prefix (`s3://`), and suffix (`.zarr.zip`).

### Storage Options

For S3 backends, optional `storage_options` (a `dict`) are forwarded to
`s3fs.S3FileSystem` and `fsspec`. Common keys:

| Key       | Purpose                          | Example    |
|-----------|----------------------------------|------------|
| `profile` | AWS CLI named profile            | `"my-profile"` |
| `region`  | AWS region override              | `"us-west-2"` |

For Polars parquet reads, s3fs-style keys are translated automatically
(e.g. `profile` -> `aws_profile`) by `_S3DirBackend._polars_storage_options()`.

Usage:
```python
loader = Loader("s3://bucket/run.zarr", storage_options={"profile": "my-profile"})
```

Storage options are ignored for local and in-memory backends.

---

## 3. `RunInfo` ABC

All backends inherit from the `RunInfo` abstract base class, which is also
the public-facing handle users interact with via `loader._runs`. It provides:

1. **Abstract IO methods**: `_list_tables`, `_read_table`, `_has_table`,
   `_list_points`, `_read_points`, `_read_table_metadata_attrs`,
   `_read_json`, `_read_sdata_attrs`, `_is_local` (property), `_source`
   (property). `_is_zip` (property) defaults to `False`; zip backends
   override it.
2. **Metadata caching**: reads `global-table-metadata.uns` once on first
   access (lazy property `_metadata`). Emits a `UserWarning` if the
   metadata table is missing. All derived lookups are also lazily
   computed and cached.
3. **Run-name derivation**: `_derive_run_name()` is called during `__init__`.
   The base implementation strips `.zarr`/`.zarr.zip` suffixes from
   `_source`. `_LocalDirBackend` and `_SDataBackend` override it to first
   try `RunName` from `metadata/inputs/run/RunParameters.json`.
4. **Source derivation**: each backend implements `_source` as a property
   derived from its own state (path, URL, or object id).
5. **Point selection**: `select_points()` and lookup helpers.
6. **Summary**: `summary()` returns a human-readable multi-line string
   with run name, source, tables, batches, wells, tiles, points,
   and build version.

```python
class RunInfo(abc.ABC):
    def __init__(self) -> None: ...

    # Public attributes and properties
    run_name: str                           # derived by _derive_run_name()
    run_date: str | None                    # lazy; sdata attrs run_metadata.run_parameters.Date
    _source: str                            # abstract property, per-backend

    batches: pd.DataFrame                   # lazy, cached
    wells: pd.DataFrame                     # lazy, cached
    tiles: pd.DataFrame                     # lazy, cached
    targets_barcoding: pd.DataFrame
    targets_imaging: pd.DataFrame
    diss_targets_metadata: pd.DataFrame
    build_provenance: dict
    table_names: list[str]                  # delegates to _list_tables()
    point_names: list[str]                  # delegates to _list_points()
    batch_names: list[str]                  # derived from batches
    batch_types: dict[str, str]
    well_names: list[str]
    available_tables: set[str]
    available_target_types: set[str]

    # Summary / reachability
    def summary(self) -> str: ...
    def source_exists(self) -> bool: ...

    # Point selection
    def select_points(...) -> list[str]: ...
    def batches_for_target_type(target_type) -> list[str]: ...
    def table_suffixes_for_target_type(target_type) -> list[str]: ...  # OPS panel discovery
    def points_for_target_type(target_type) -> list[str]: ...
    def points_for_batch(batch_name) -> list[str]: ...
    def points_for_well(well_id) -> list[str]: ...
    def points_for_tile(tile_name) -> list[str]: ...

    # Abstract IO methods every backend implements
    def _list_tables(self) -> list[str]: ...
    def _read_table(self, name: str, *, lazy: bool = True) -> ad.AnnData: ...
    def _has_table(self, name: str) -> bool: ...
    def _list_points(self) -> list[str]: ...
    def _read_points(self, name: str, *, lazy: bool = True) -> pl.LazyFrame: ...
    def _read_table_metadata_attrs(self, name: str) -> dict: ...
    def _read_json(self, path: str) -> dict | None: ...
    def _read_sdata_attrs(self) -> dict: ...
    @property
    def _is_local(self) -> bool: ...
    @property
    def _source(self) -> str: ...
```

### Init Order

Concrete backends set their own attributes **before** calling
`super().__init__()`, so that `_source`, `_read_json`, and other abstract
methods are available during base-class initialization (e.g. for
`_derive_run_name()`).

### Backend-Specific Behavior

| Backend           | Tables read via         | Points read via              | lazy support          |
|-------------------|------------------------|------------------------------|-----------------------|
| `_LocalDirBackend`| `read_lazy(path)` when `lazy=True`, else `ad.read_zarr(path)` | `pl.scan_parquet(path)` | Full |
| `_LocalZipBackend`| `read_lazy` (lazy) or `ad.read_zarr` (keeps `ZipStore` open) | Extract element to the persistent cache, then `pl.scan_parquet` | Full (points lazy after extract) |
| `_S3DirBackend`   | `read_lazy(mapper)` when `lazy=True`, else `ad.read_zarr(mapper)` | `pl.scan_parquet(s3://...)` | Full |
| `_S3ZipBackend`   | Per-table extract to persistent ETag cache; metadata via remote zarr | Per-element extract to the same ETag cache, then `pl.scan_parquet` | Lazy after extract |
| `_SDataBackend`   | `sdata.tables[name]`   | `sdata.points[name]` -> PyArrow -> Polars | No (in-memory) |

With `lazy=True`, tables are opened with `anndata.experimental.read_lazy`
and `obs`/`var` are materialised to pandas; if `read_lazy` is unavailable or
fails, the backend falls back to `ad.read_zarr`.

### Partitioned Parquet Handling

Points parquet may be a single file or a partitioned directory containing
`part.N.parquet` files. The backends handle both layouts:

- **Local dir**: returns a `**/*.parquet` glob for directories, plain path otherwise.
- **Local zip / S3 zip**: extract the single file or every `*.parquet` part
  (preserving the layout) into the cache, then scan the extracted files.
- **S3 dir**: checks `isdir()` and uses `**/*.parquet` glob if true.

### `_LocalZipBackend` Points Cache

Points elements are extracted once into
`$ELEMBIO_SPATIALDATA_CACHE/local_zip/<sha1(path)[:12]>_<mtime>_<size>/<name>/`
(a `.{name}.extracted` sentinel marks completion). The key includes the
zip's mtime and size, so replacing the zip in place triggers re-extraction.

### `_SDataBackend` Points Conversion

In-memory `SpatialData` points (dask/pandas DataFrames) are converted to
Polars via PyArrow (`pa.Table.from_pandas`) to preserve nullable integer
types that would otherwise be lost in a direct pandas -> Polars conversion.

### `_LocalZipBackend` Table Caching

`_list_tables()` caches its result after the first call to avoid
repeatedly opening the zip archive for table discovery.

### S3 Zip Selective Extraction

`_S3ZipBackend` reads the zip central directory remotely via `fsspec`'s
`ZipFileSystem` (one range request to the end of the file). Table
discovery, root attrs (`zarr.json` / `.zattrs`), and per-table
`uns/metadata_attrs` are read without extracting full table bodies (only
the small ``metadata_attrs`` subtree is cached when needed).
When a table body is requested, only that table's zarr subtree is
extracted sequentially into a persistent on-disk cache keyed by S3 URL
hash and object ETag (default root: `~/.cache/elembio-spatialdata`,
override with `ELEMBIO_SPATIALDATA_CACHE`). Points elements are extracted
the same way, one element at a time, under `points/<name>/` in the same
cache, and then scanned lazily.

Per-element extraction (tables and points) uses a lock directory
(`.{name}.lock`) and a sentinel file (`.{name}.extracted`). The lock holder
heartbeats `owner.json` every `ELEMBIO_S3ZIP_LOCK_HEARTBEAT_SEC` (default
30s) during extraction; stale locks are reclaimed when `owner.json` has not
been touched for `ELEMBIO_S3ZIP_LOCK_STALE_SEC` (default 120s). A process
waiting on another's lock raises `TimeoutError` after
`ELEMBIO_S3ZIP_LOCK_MAX_WAIT_SEC` (default 1800s). Missing
`metadata_attrs` subtrees are cached via the sentinel so later reads do
not re-hit S3. After a successful extract, older ETag sibling
directories for the same URL are pruned only when they have no active
`.lock` directories and are older than `ELEMBIO_S3ZIP_PRUNE_MIN_AGE_SEC`
(default 3600s).

The in-memory `_S3ZipBackend` caches (`_extract_cache`,
`_metadata_attrs_cache`, `_points_extract_cache`) map
`(url, etag, name)` to local paths. An `atexit` handler clears these dicts;
persistent cache files remain on disk for reuse across sessions.

### Factory Function (`make_backend`)

`make_backend(store_input, *, storage_options=None)` selects and instantiates the
appropriate backend based on the input type. Each backend derives its own
`run_name` and `_source` during construction -- no external ID or source
arguments are needed.

### Multi-Store Orchestration (`_open_stores` in `loader.py`)

`_open_stores(inputs, storage_options)` in `loader.py` normalizes inputs
into a deduplicated `dict[str, RunInfo]`. Duplicate detection is based on
string-normalized paths. Collision on `run_name` is resolved by appending
`__N` suffixes (e.g. `"RUN-0001"`, `"RUN-0001__2"`).

---

## 4. Data Model

### 4.1 Zarr Store Layout

A SpatialData zarr store has this structure:

```
run.zarr/
+-- metadata/
|   +-- inputs/
|       +-- run/                     # Source-mirror folder for run-dir files
|       |   +-- RunParameters.json   # Contains RunName
|       |   +-- Panel.json
|       |   +-- RunManifest.json     # (only the .json form is read)
|       +-- cells2stats/             # Source-mirror folder for c2s files (incl. Versions.json)
|       +-- bases2fastq/             # Source-mirror folder for b2f files
|       +-- fastq2parquet_3prime/    # Source-mirror folder for f2p_3p files (incl. multiqc/)
|       +-- pipeline/Provenance.json    # Pipeline-level provenance (optional; older stores: tetonatlas/)
|       +-- witness_hashes.json      # Copy-equivalence hashes (optional)
+-- tables/
|   +-- global-table-metadata/       # Metadata-only AnnData (no X)
|   +-- global-table-imaging/        # AnnData: cells x imaging features (CellPaint + IF Protein)
|   +-- global-table-ThreePrimeUntargeted/  # AnnData: cells x genes (3′ transcriptome)
|   +-- global-table-{batch}-{site}/ # AnnData per OPS batch × target site
|   +-- global-table-Transcript/     # (legacy) cells x transcript targets
|   +-- global-table-Protein/        # (legacy) cells x protein targets
+-- points/
    +-- {prefix}-points-{batch}/     # One directory per (well-or-tile, batch)
    |   +-- points.parquet           # Single file or partitioned directory
    +-- ...
```

### 4.2 Table Naming Convention

Tables follow the pattern `global-table-{name}`:

| Table name                  | Target type            | Batch class      |
|-----------------------------|------------------------|------------------|
| `global-table-imaging`      | `imaging`              | `ImagingBatch` / `PreAmpImagingBatch` |
| `global-table-ThreePrimeUntargeted` | `ThreePrimeUntargeted` | `DISSBatch` |
| `global-table-{batch}-{panel}` | `OPS` / `SpecializedTargeted` | `DISSBatch` |
| `global-table-Transcript`   | `Transcript` *(legacy)*  | `BarcodingBatch` |
| `global-table-Protein`      | `Protein` *(legacy)*     | `BarcodingBatch` |
| `global-table-metadata`     | *(special)*            | N/A              |

> **OPS has no aggregated table.** OPS / specialized panels are stored one
> table per batch × panel (e.g. `global-table-B01-A549`), not a
> single `global-table-OPS`. Discover the suffixes with
> `RunInfo.table_suffixes_for_target_type(("OPS", "SpecializedTargeted"))`
> (§5.2), then pass each to `Loader.load_tables(table=...)`.

> **Modalities:** Current runs contain `imaging`, `ThreePrimeUntargeted`
> (3′ transcriptome), and `OPS` / `SpecializedTargeted` tables. `Transcript`
> and `Protein` barcoding tables are deprecated but still readable. Older
> SDOs may name the imaging table `global-table-CellPaint`; `"CellPaint"`
> normalizes to `"imaging"` throughout.

In addition, the store may contain per-batch tables
(`global-table-{batch_name}`), OPS per-batch-site tables
(`global-table-{batch}-{site}`), and `_per_batch` layer tables
(e.g. `global-table-Transcript_per_batch`). All non-metadata tables are
included in `RunInfo.available_tables`.

### 4.3 Table Structure (AnnData)

Each table is an AnnData object:

**`obs`** (one row per cell):

| Column          | Type   | Description                          |
|-----------------|--------|--------------------------------------|
| `Cell`          | int    | Cell identifier                      |
| `Well`          | str    | Well location (e.g. `"A1"`)          |
| `WellLabel`     | str    | Well label from sample sheet         |
| `tile_id`       | str    | Tile identifier                      |
| `X_tile`, `Y_tile` | float | Cell centroid, tile-local pixels (older stores: `X`, `Y`) |
| `cell_id_local` | int    | Local cell ID within tile            |
| `cell_id_global`| int    | Global cell ID (tile_index << 16 | local_id) |
| `region_key`    | str    | SpatialData region reference         |
| `well_id`       | str    | Well identifier                      |

**`var`** (one row per feature/target):

| Column               | Type   | Description                    |
|----------------------|--------|--------------------------------|
| `target_id`          | int    | Numeric target identifier      |
| `target_type`        | str    | `"Transcript"`, `"Protein"`, etc. |
| `target_sequence`    | str    | Barcode sequence               |
| `control_type`       | str    | Control classification or `""` |
| `batch_name`         | str    | Which batch this target belongs to |
| `probe_concentration`| float  | Probe concentration (nM)       |

**`X`**: count matrix (cells x features), `float32`.

**`layers`**: `{"nuclear": ...}` for barcoding tables (nuclear-only counts).

**`obsm`** (optional per-cell matrices):

| Key                  | Type      | Description                                  |
|----------------------|-----------|----------------------------------------------|
| `total_counts`       | DataFrame | Per-batch total counts per cell               |
| `unassigned_counts`  | DataFrame | Per-batch unassigned counts per cell          |
| `cell_geometry`      | DataFrame | Cell morphology (Area, NuclearArea, Xum_tile, Yum_tile, ...) |
| `control_variables`  | DataFrame | Control variable values per cell              |
| `data_missing`       | DataFrame | Missing data flags per cell                   |
| `target_assignment`  | DataFrame | OPS target assignment (OPS tables only)       |
| `spatial`            | DataFrame | Cell centroids (`X`, `Y`, float32) in well-canvas pixels; tile-local in tile mode; `NaN` for cells whose tile has no canvas placement |
| `spatial_artifacts`  | DataFrame | Boolean per-cell artifact membership (see filter SPECIFICATION) |

**`varm`** (optional per-variable matrices):

| Key                   | Type      | Description                               |
|-----------------------|-----------|-------------------------------------------|
| `extracellular_counts`| DataFrame | Per-well extracellular counts (columns are well IDs) |

**`uns`**:

| Key                | Type | Description                                    |
|--------------------|------|------------------------------------------------|
| `metadata_attrs`   | dict | Batch/target metadata (see Section 4.5)        |
| `control_variables_metadata` | dict | Control variable definitions (if present) |
| `spatial_artifact_params` | dict | Settings and column map for `obsm["spatial_artifacts"]` (`version`, `spatial_method`, `detection_enabled`, `detection_disabled_reason`, `rollup_columns`, per-modality params). Absent when artifacts were not written. |

Note: `spatialdata_attrs` is stripped from `.uns` at load time since it is
a SpatialData framework internal not relevant to analysis.

### 4.4 Points Element Naming

Point element names follow the pattern:

```
{prefix}-points-{batch_name}
```

Where `{prefix}` depends on the output mode:

| Mode      | Prefix value    | Example element name         |
|-----------|-----------------|------------------------------|
| Well mode | Well ID         | `A1-points-B01`              |
| Tile mode | Tile name       | `L1R01C01S1-points-B01`      |

The mode is auto-detected by `RunInfo._point_mode` (a lazy-cached property)
by comparing the set of parsed prefixes against well names and tile names
from the metadata. If prefixes are a subset of tile names (and not of well
names), tile mode is used.

Parsing is done by `split_points_name(name)` which returns
`(prefix, batch_name)` or `(None, None)` if the name doesn't match.

### 4.5 Table `metadata_attrs`

Each table carries a `uns["metadata_attrs"]` dict:

```python
{
    "target_type": "Protein",  # str
    "batch_class": "BarcodingBatch",  # str
    "batch_name": np.array(["B01"], dtype=StringDType()),  # ndarray of str
    "merge_strategy": "expand_features",  # str
}
```

This is the **canonical source** for the batch -> target_type mapping. The
loader reads `metadata_attrs` from each table via lightweight zarr reads
(`_read_table_metadata_attrs`) on first access -- avoiding full X matrix
loads -- to build and cache `RunInfo._batch_to_target_types`.

### 4.6 Points Parquet Schema

| Column              | Polars dtype | Description                              |
|---------------------|-------------|------------------------------------------|
| `x`                 | `Int64`     | X pixel coordinate (global)              |
| `y`                 | `Int64`     | Y pixel coordinate (global)              |
| `z`                 | `Int64`     | Z coordinate                             |
| `Xum`               | `Float32`   | X coordinate in micrometers              |
| `Yum`               | `Float32`   | Y coordinate in micrometers              |
| `Zum`               | `Float32`   | Z coordinate in micrometers              |
| `target_index`      | `UInt16`    | Numeric target identifier                |
| `cell`              | varies      | Local cell assignment (e.g. `Int32` in OPS points, zero-padded `String` in 3′ Untargeted points) |
| `is_nuclear`        | `Boolean`   | Whether the point is in the nucleus      |
| `target_name`       | `String`    | Target gene/protein name                 |
| `tile_id`           | `String`    | Tile identifier                          |
| `well_id`           | `String`    | Well identifier                          |
| `label_{Mask}`      | `Int64`     | Cell label per segmentation mask (local, within tile). Column name varies by mask (e.g. `label_Cell`, `label_Nuclear`). |
| `label_{Mask}_global` | `Int32`   | Global cell label per mask (`tile_index << 16 \| local_id`). Column name varies (e.g. `label_Cell_global`). |
| `x_tile`            | `Int64`     | X pixel coordinate (within tile)         |
| `y_tile`            | `Int64`     | Y pixel coordinate (within tile)         |

**Note on `label_` columns**: The `label_Cell` / `label_Cell_global` names shown
above are for the default cell segmentation mask. Different segmentation masks
produce different column names (e.g. `label_Nuclei`, `label_Nuclei_global`).
These columns are not validated.

**Required columns** validated when `Loader(..., validate_points=True)`:

- **Always required** (`_ALWAYS_REQUIRED_POINTS_COLUMNS`): `well_id`,
  `tile_id`, `x`, `y`. These are the structural columns needed to identify
  a point within a run.
- **Conditionally required**: `target_name` may be absent for ThreePrime/OPS
  batches where the upstream alignment pipeline was not run.

Validation samples one point element per batch (not all elements) and emits
separate warnings for missing structural columns vs missing `target_name`.

### 4.7 Metadata Table (`global-table-metadata`)

The `global-table-metadata` table is a metadata-only AnnData whose `uns` dict
contains run-level DataFrames:

| `uns` key              | Type      | Description                                         |
|------------------------|-----------|-----------------------------------------------------|
| `batches`              | DataFrame | One row per batch: `batch_name`, `batch_class`, `batch_type`, `min_cycles`, `run_order`, `num_targets`, `pmg_mask`, `barcode_mask`, `fastq_mask`, `default_mismatch`, `num_diss_targets` |
| `wells`                | DataFrame | One row per well (indexed by `well_name`): `well_index`, `well_label`, `well_layout`, `cell_type`, `cell_diameter`, `target_site_names` |
| `tiles`                | DataFrame | One row per tile (indexed by `tile_name`): `tile_index`, `x_mm`, `y_mm`, `well_name`, `well_label` |
| `targets_barcoding`    | DataFrame | Per-target metadata for barcoding batches: `target_name`, `target_id`, `target_type`, `target_sequence`, `control_type`, `batch_name`, `probe_concentration` |
| `targets_imaging`      | DataFrame | Per-target metadata for imaging batches (same columns as above) |
| `diss_targets_metadata`| DataFrame | Per-target metadata for DISS/OPS batches (variable columns from Panel.json DISSTargets schema) |
| `build_provenance`     | dict      | Build pipeline version and parameters                |
| `missing_batches`      | list[str] | Panel-declared count batches with no points data (no table built); present only when non-empty |
| `c2s_runstats_*`       | DataFrame | cells2stats RunStats: `run_metadata` (with folded `Demux_*`/`Cyto_*`/`RawRead_*`/`Versions_*` scalars), `cyto_wells`, `cyto_batches`, `rawread_wells`, `rawread_batches`, `demux_*`, `targetstats_*` |
| `b2f_runstats_*`       | DataFrame | bases2fastq RunStats: `run_metadata`, `batches`, `wells` |
| `f2p_3p_runstats_*`    | DataFrame | fastq2parquet_3prime RunStats: `run_metadata` (with folded `RunSummary_*` scalars), `samples`, `sample_cell_metrics`, `sample_cell_metric_totals` |

A missing `global-table-metadata` table emits a `UserWarning` on first
metadata access, and all metadata properties return empty defaults.

---

## 5. `RunInfo` -- Properties and Selection

`RunInfo` instances are created by `make_backend()` via `_open_stores()`
and stored as `dict[str, RunInfo]` in `Loader._runs`. All metadata
properties are lazily computed and cached on first access.

### 5.1 Properties

| Property                | Return type       | Description                                  |
|-------------------------|-------------------|----------------------------------------------|
| `run_name`              | `str`             | Human-readable run identifier                |
| `run_date`              | `str \| None`     | Instrument run start date (ISO string from sdata attrs); lazy, cached |
| `_source`               | `str`             | Abstract property: path, URL, or object id   |
| `table_names`           | `list[str]`       | All table names in this store                |
| `point_names`           | `list[str]`       | All point element names in this store        |
| `batches`               | `DataFrame`       | One row per batch (from metadata)            |
| `wells`                 | `DataFrame`       | One row per well (from metadata)             |
| `tiles`                 | `DataFrame`       | One row per tile (from metadata)             |
| `targets_barcoding`     | `DataFrame`       | Per-target metadata for barcoding batches    |
| `targets_imaging`       | `DataFrame`       | Per-target metadata for imaging batches      |
| `diss_targets_metadata` | `DataFrame`       | Per-target metadata for DISS/OPS batches     |
| `build_provenance`      | `dict`            | Build pipeline version and parameters        |
| `batch_names`           | `list[str]`       | Sorted batch names from `batches` DataFrame  |
| `batch_types`           | `dict[str, str]`  | `{batch_name: batch_type}` from `batches`    |
| `well_names`            | `list[str]`       | Sorted well IDs from `wells` index           |
| `available_tables`      | `set[str]`        | Table suffixes from `global-table-*` names (excluding metadata) |
| `available_target_types`| `set[str]`        | Target types derived from table metadata     |

### 5.2 Methods

| Method                                  | Returns                | Description |
|-----------------------------------------|------------------------|-------------|
| `summary()`                             | `str`                  | Multi-line human-readable summary (run name, source, tables, batches with types, wells, tiles, points with mode, build version) |
| `source_exists()`                       | `bool`                 | Cheap reachability check of the store path/URL (not a validity check); always `True` for in-memory backends |
| `batches_for_target_type(target_type)`  | `list[str]`            | Batch names whose targets include the given type |
| `table_suffixes_for_target_type(target_type)` | `list[str]`      | Table suffixes whose `uns['metadata_attrs']['target_type']` matches. Accepts a string or a tuple of aliases (e.g. `("OPS", "SpecializedTargeted")`). The correct way to discover OPS/specialized panel tables (returns per-site suffixes like `["B01-A549"]`) |
| `points_for_target_type(target_type)`   | `list[str]`            | Point elements for batches of the given target type |
| `points_for_batch(batch_name)`          | `list[str]`            | Point elements for a specific batch |
| `points_for_well(well_id)`              | `list[str]`            | Point elements for a well (resolves through tiles in tile mode) |
| `points_for_tile(tile_name)`            | `list[str]`            | Point elements for a tile (direct lookup in tile mode; `[]` in well mode) |
| `select_points(**filters)`              | `list[str]`            | Composable intersection of filters (see below) |

### 5.3 `select_points()` -- Composable Filters

```python
run.select_points(
    target_type: str | list[str] | None = None,
    batch: str | list[str] | None = None,
    well: str | list[str] | None = None,
    tile: str | list[str] | None = None,
    point_name: str | list[str] | None = None,
)
```

- All filters compose as **intersection** (AND logic).
- `point_name` acts as an **explicit override** -- when provided, other filters
  are ignored.
- With no filters, returns all point names in the run.
- `well` works in both modes: in tile mode it resolves through `_tile_to_well`.
- `tile` works in tile mode only; in well mode it emits a `UserWarning` and
  the tile filter is ignored (the other filters still apply).
- Unknown `target_type`, `batch`, `well`, `tile`, or `point_name` values emit
  a `UserWarning` listing the available values.

### 5.4 `available_tables`

Returns all `global-table-{name}` suffixes except the metadata table.
These are the valid values for the `table` parameter in `Loader.load_tables()`.

### 5.5 `available_target_types`

Returns the set of `target_type` values found in table metadata
(`uns['metadata_attrs']['target_type']`). Derived from
`_batch_to_target_types`. These are the valid values for the
`target_type` parameter in `Loader.load_points()`.

---

## 6. `Loader` -- Multi-Run Entrypoint

### 6.1 Construction

```python
Loader(
    stores: str | Path | SpatialData | list[...],
    storage_options: dict | None = None,
    *,
    validate_points: bool = False,
)
```

1. Normalizes inputs, rejects duplicates, and opens each store via
   `make_backend()` (`_open_stores()` in `loader.py`).
2. Each backend is a `RunInfo` instance with lazily cached metadata.
3. If `validate_points` is `True`, validates points schema
   (`_validate_points_schema()`): samples one point element per batch per run
   and checks for required columns. Defaults to `False` because this is slow
   for stores with many batches (e.g. 3′ untargeted runs) or remote/zip-backed
   stores. Set `validate_points=True` to surface schema issues at construction.

### 6.2 Discovery

| Property / Method       | Return type  | Description                              |
|-------------------------|-------------|------------------------------------------|
| `available_tables`      | `set[str]`  | Union of table suffixes across all runs  |
| `available_target_types`| `set[str]`  | Union of target types across all runs    |
| `summary()`             | `str`       | Multi-line summary of all runs           |
| `export_tables(out_dir, ...)` | `list[Path]` | Write tables to local files (see Section 6.5) |
| `_runs`                 | `dict[str, RunInfo]` | Keyed by run name (with `__N` dedup) |

### 6.3 `load_tables()`

```python
loader.load_tables(
    *,
    table: str,                              # Required, table name suffix
    concatenated: bool = True,
    lazy: bool = True,
    include_target_assignment: bool = False,
    include_cell_geometry: bool = False,
)
```

**Behavior:**

1. Resolves `table` to `global-table-{table}`.
2. Reads the table from each run via `RunInfo._read_table()`.
3. Strips `spatialdata_attrs` from `.uns` (framework internal, not relevant
   to analysis).
4. Applies enrichments (if requested):
   - `include_cell_geometry` (default `False`): when `False`, leaves
     `obsm["cell_geometry"]` unchanged from the store; when `True`, merges
     those columns into the main `X` matrix (appended as new variables in
     `.var`).
   - `include_target_assignment`: merges `obsm["target_assignment"]` into
     `obs`. Checks for column name conflicts with existing `obs` columns
     and drops duplicate columns before concatenation. Specific to OPS tables.
5. If only one run produced a table, returns that `AnnData` directly
   (regardless of `concatenated`).
6. If multiple runs produced tables and `concatenated=True`, concatenates
   them via `_concat_tables_rowwise` (see Section 7.1).

**Return types:**

| Runs loaded | `concatenated` | Return type              |
|-------------|----------------|--------------------------|
| 1           | *(ignored)*    | `AnnData`                |
| >1          | `True` (default) | `AnnData`              |
| >1          | `False`        | `dict[str, AnnData]`     |

### 6.4 `load_points()`

```python
loader.load_points(
    *,
    target_type: str | list[str] | None = None,
    batch: str | list[str] | None = None,
    well: str | list[str] | None = None,
    tile: str | list[str] | None = None,
    point_name: str | list[str] | None = None,
    concatenated: bool = True,
    lazy: bool = True,
)
```

**Behavior:**

1. For each run, calls `RunInfo.select_points()` with the given filters.
2. Reads each selected point element via `RunInfo._read_points()`.
3. Adds `batch_name` and `run_name` columns to each LazyFrame.
4. If `concatenated=True`, concatenates all frames with a `source_key`
   column (`{run_key}/{point_name}`).

**Return types:**

| `concatenated` | Return type                          |
|-----------------|--------------------------------------|
| `True` (default) | `LazyFrame`                        |
| `False`         | `dict[str, dict[str, LazyFrame]]`    |

### 6.5 `export_tables()`

```python
loader.export_tables(
    out_dir: str | Path,
    *,
    tables: Sequence[str] | None = None,   # default: all available_tables
    fmt: str = "h5ad",                      # or "zarr"
    overwrite: bool = False,
    include_metadata: bool = True,          # also export global-table-metadata
) -> list[Path]
```

Writes one file per `(run, table)` named `{run}__{table}.{ext}`, reusing
`load_tables()` (one table materialised at a time). Source run, store, and
build version are stamped into `uns['export_provenance']`. Spatial
coordinates are kept; images are not. Existing files are skipped with a
warning unless `overwrite=True`. Useful for pulling just the tables from an
S3 store without copying images and points.

---

## 7. Concatenation (`merge.py`)

### 7.1 Row-wise Table Concatenation (`_concat_tables_rowwise`)

Concatenates AnnData objects across runs (axis=0). The merge strategy is
hardcoded (not user-configurable) with semantics chosen for correctness:

**`.var_names`** (target/feature axis):
- Outer-joined (union of all targets). Cells missing a target get `NaN` in `X`.

**`.var`** (per-target metadata columns):
- Kept only when `.var_names` are identical across all runs.
- Cleared when var_names differ (values like `target_id`, `batch_name` no
  longer match across runs). Emits a `UserWarning` listing the differences.

**`.varm`** (per-variable matrices, e.g. `extracellular_counts`):
- When var_names are identical: DataFrames are column-stacked with
  `{run_name}_{column}` prefixes (e.g. `RUN-0001_A1`, `RUN-0002_A1`).
- When var_names differ: dropped entirely.
- Handles both eager (`pd.DataFrame`) and lazy (`Dataset2D`) varm values
  by materializing lazy values via `.to_memory()`.

**`.obs_names`**: made unique across runs with `ad.concat(..., index_unique="-")`,
i.e. each cell name gets a `-{source_key}` suffix.

**`.obsm`** (per-cell matrices):
- Outer-merged across runs. Missing columns are filled with `NaN`.

**`.uns`** (unstructured metadata):
- `spatialdata_attrs` is always stripped (framework internal).
- Remaining entries are kept only if identical across all runs (recursive
  deep comparison handles nested dicts, numpy arrays, and DataFrames).
- Differing entries are dropped silently.

**`.obs`** provenance columns added:

| Column         | Description                                    |
|----------------|------------------------------------------------|
| `source_key`   | Run key from the `Loader` (e.g. `"RUN-0001"`) |
| `run_name`     | `source_key` with `__N` dedup suffix stripped  |
| `source_store` | Original path/URI of the store                 |

**StringDType handling**: numpy `StringDType` columns in `.obs` and `.var`
are downcast to `object` before concatenation to avoid `DTypePromotionError`
when the same column has mixed dtypes across runs.

### 7.2 Points Frame Concatenation (`_concat_frames`)

Concatenates Polars LazyFrames with `how="diagonal_relaxed"` (tolerates
schema differences). Adds a `source_key` column for provenance.

---

## 8. Error Handling and Warnings

| Condition                                    | Behavior                        |
|----------------------------------------------|---------------------------------|
| Missing `global-table-metadata`              | `UserWarning` on first metadata access (metadata will be empty) |
| Missing table for requested name             | `UserWarning` with available tables listed   |
| Missing structural points columns            | `UserWarning` at construction (only with `validate_points=True`) |
| Missing `target_name` in points              | `UserWarning` (only with `validate_points=True`; batch may not be loadable by type) |
| No matching points for a filter              | `UserWarning` (empty result)    |
| Unknown `target_type` / `batch` / `well` / `tile` / `point_name` in `select_points()` | `UserWarning` (unknown values skipped) |
| `tile=` filter on a well-mode run            | `UserWarning`; tile filter ignored |
| S3 zip extraction lock not released in time  | `TimeoutError` after `ELEMBIO_S3ZIP_LOCK_MAX_WAIT_SEC` |
| Column conflicts in `_enrich_target_assignment` | Drops duplicates silently     |
| `.var_names` mismatch during concatenation   | `UserWarning`; `.var` and `.varm` cleared |
| `.obsm` schema mismatch during concatenation | `UserWarning`                   |
| No tables loaded                             | `ValueError`                    |
| No points loaded for concatenation           | `ValueError`                    |
| Duplicate store inputs                       | `ValueError` at construction    |

---

## 9. Usage Examples

### Single run, local

```python
from elembio_spatialdata_analysis.load import Loader

loader = Loader("/data/RUN-0001.zarr")

# Summary
print(loader.summary())

# Inspect run metadata
run = next(iter(loader._runs.values()))
run.well_names  # ['A1', 'A2', 'B1', ...]
run.batch_names  # ['B01', 'B02', ..., 'CP01', 'CP02']
run.available_tables  # {'imaging', 'ThreePrimeUntargeted', 'B01-A549'}
run.available_target_types  # {'imaging', 'ThreePrimeUntargeted', 'OPS'}
run.build_provenance  # {'pipeline_version': '1.0.0', ...}

# Load the 3′ transcriptome table (concatenated by default)
adata = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)

# Load as dict (one AnnData per run)
tables = loader.load_tables(table="ThreePrimeUntargeted", concatenated=False, lazy=False)
adata = tables["RUN-0001"]

# OPS: one table per batch × target site -- discover the suffixes first
for suffix in run.table_suffixes_for_target_type(("OPS", "SpecializedTargeted")):
    ops = loader.load_tables(table=suffix, include_target_assignment=True, lazy=False)

# Load points for batch B01 in well A1 (concatenated LazyFrame by default)
points = loader.load_points(batch="B01", well="A1")
```

### Single run, S3

```python
loader = Loader(
    "s3://your-bucket/spatialdata_runs/RUN-0001.zarr",
    storage_options={"profile": "my-profile"},
)
```

### Multi-run concatenation

```python
loader = Loader(
    [
        "/data/RUN-0001.zarr",
        "/data/RUN-0002.zarr",
        "/data/RUN-0003.zarr",
    ]
)

# Tables concatenated across runs by default
adata = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
# adata.obs["run_name"] identifies the source run
# adata.varm["extracellular_counts"] (if present) columns prefixed with run name

# Points concatenated across runs by default
lf = loader.load_points(batch="B01", well="A1")
# lf has "source_key" and "run_name" columns
```

### Mixed S3 + local

```python
loader = Loader(
    [
        "s3://your-bucket/spatialdata_runs/RUN-0001.zarr",
        "/data/RUN-0002.zarr",
    ],
    storage_options={"profile": "my-profile"},
)
```

### Tile-mode run

```python
loader = Loader("/data/tile-mode-run.zarr")
run = next(iter(loader._runs.values()))

run._point_mode  # "tile"
run.points_for_tile("L1R01C01S1")  # ['L1R01C01S1-points-B01', ...]
run.points_for_well("A1")  # resolves through tiles metadata

# Filter by tile
points = loader.load_points(tile="L1R01C01S1", batch="B01")

# Filter by well (resolves to tiles automatically)
points = loader.load_points(well="A1", batch="B01")
```

### Imaging with cell geometry

```python
adata = loader.load_tables(
    table="imaging",  # also finds a legacy global-table-CellPaint
    include_cell_geometry=True,
    lazy=False,
)
# adata.X now includes Area, NuclearArea, Xum, Yum, AreaUm, ... as extra variables
# adata.var lists both the original variables and the geometry columns
```

### Cross-run concatenation with different feature sets

```python
loader = Loader(["/data/RUN-0001.zarr", "/data/RUN-0004.zarr"])

adata = loader.load_tables(table="ThreePrimeUntargeted", lazy=False)
# adata.var_names: union of all features across runs
# adata.var: empty (columns cleared -- metadata no longer matches)
# adata.varm: empty (dropped -- var_names differ)
# adata.X: NaN for targets missing in a given run
# Warning emitted listing the var_names differences
```
