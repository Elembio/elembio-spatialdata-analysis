"""High-level, multi-run loading API for SpatialData stores.

This module defines :class:`Loader`, the main user-facing entrypoint for:

- discovery over run-level metadata (batches, wells, tiles, targets), and
- loading either tables (:class:`anndata.AnnData`) or points (:class:`polars.LazyFrame`)
  from one or many stores with optional concatenation.

Each store is opened at construction; run metadata is read on first access
and cached in :class:`RunInfo`.
"""

from __future__ import annotations

import warnings
from collections import Counter
from collections.abc import Sequence
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl
import spatialdata as sd

from elembio_spatialdata_analysis.load.backends import (
    RunInfo,
    StoreInput,
    make_backend,
    split_points_name,
)
from elembio_spatialdata_analysis.load.merge import (
    _concat_frames,
    _concat_tables_rowwise,
)

# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------

_LoadedTables = dict[str, ad.AnnData]
"""``{source_key: adata}``."""

_LoadedPoints = dict[str, dict[str, pl.LazyFrame]]
"""``{run_key: {points_name: LazyFrame}}``."""

_LoadedTablesResult = _LoadedTables | ad.AnnData
_LoadedPointsResult = _LoadedPoints | pl.LazyFrame

# Structural columns that must always be present in points.
# ``target_name`` may be absent for ThreePrime/OPS batches when the upstream
# alignment pipeline was not run (spec §4.5, §10).
_ALWAYS_REQUIRED_POINTS_COLUMNS = frozenset({"well_id", "tile_id", "x", "y"})


# ---------------------------------------------------------------------------
# Table loading
# ---------------------------------------------------------------------------

# Imaging is stored as ``global-table-imaging`` on current SDOs and
# ``global-table-CellPaint`` on legacy stores.  Either name resolves to
# whichever key a given store actually contains, so old and new stores load
# through the same ``load_tables(table=...)`` call.
_IMAGING_TABLE_ALIASES = ("imaging", "CellPaint")


def _resolve_table_key(table: str, table_names: Sequence[str]) -> str | None:
    """Resolve *table* to an existing ``global-table-{...}`` key in one store.

    Returns the matching key, or ``None`` if the store has no such table. The
    exact name is preferred; for the imaging pair (``"imaging"`` /
    ``"CellPaint"``) the alias is tried when the requested name is absent, so a
    request for ``"imaging"`` still finds a legacy ``global-table-CellPaint``
    (and vice versa).
    """
    requested = f"global-table-{table}"
    if requested in table_names:
        return requested
    if table in _IMAGING_TABLE_ALIASES:
        for alias in _IMAGING_TABLE_ALIASES:
            key = f"global-table-{alias}"
            if key in table_names:
                return key
    return None


def _load_tables_core(
    keyed_stores: list[tuple[str, RunInfo]],
    *,
    table: str,
    lazy: bool,
    include_target_assignment: bool,
    include_cell_geometry: bool,
) -> tuple[_LoadedTables, dict[str, str]]:
    """Load a single table from each keyed store.

    Parameters
    ----------
    keyed_stores
        Ordered ``(source_key, run_info)`` pairs to load from.
    table
        Table to load (e.g. ``"Transcript"``).  Resolved to
        ``global-table-{table}`` within each store.
    lazy
        Forwarded to :meth:`RunInfo._read_table`.
    include_target_assignment
        If ``True``, merge ``obsm['target_assignment']`` into ``obs``.
    include_cell_geometry
        If ``False`` (default), ``obsm['cell_geometry']`` is left as stored in
        the table (no merge). If ``True``, merge those columns into the main
        ``X`` matrix (appended as new variables in ``.var``). Emits a
        ``UserWarning`` if ``True`` for non-imaging tables (the merge still
        runs).

    Returns
    -------
    tuple[_LoadedTables, dict[str, str]]
        ``(loaded, source_store_lookup)`` where:

        - ``loaded`` is ``{source_key: AnnData}``
        - ``source_store_lookup`` maps ``source_key -> store source URI/path``
    """
    if include_cell_geometry and table not in _IMAGING_TABLE_ALIASES:
        warnings.warn(
            "`include_cell_geometry=True` is intended for imaging tables.",
            UserWarning,
            stacklevel=2,
        )

    loaded: _LoadedTables = {}
    source_store_lookup: dict[str, str] = {}

    for source_key, run_info in keyed_stores:
        source_store_lookup[source_key] = run_info._source
        table_key = _resolve_table_key(table, run_info.table_names)
        if table_key is None:
            warnings.warn(
                f"Table {table!r} not found in store {source_key!r}; skipping. "
                f"Available tables: {sorted(run_info.available_tables)}",
                UserWarning,
                stacklevel=2,
            )
            continue

        adata = run_info._read_table(table_key, lazy=lazy)
        needs_copy = include_target_assignment or include_cell_geometry
        if needs_copy:
            adata = adata.copy()
        if include_target_assignment:
            _enrich_target_assignment(adata)
        if include_cell_geometry:
            adata = _enrich_cell_geometry(adata)
        adata.uns.pop("spatialdata_attrs", None)
        loaded[source_key] = adata

    return loaded, source_store_lookup


# ---------------------------------------------------------------------------
# Enrichment helpers
# ---------------------------------------------------------------------------


def _enrich_target_assignment(adata: ad.AnnData) -> ad.AnnData:
    """Merge ``.obsm['target_assignment']`` into ``.obs`` (in-place).

    Parameters
    ----------
    adata
        Table to enrich.

    Returns
    -------
    ad.AnnData
        The same object passed in, optionally modified in-place.

    Notes
    -----
    This enrichment is specific to OPS perturbation tables. It is a no-op if
    the key is missing or if ``obsm['target_assignment']`` is not a DataFrame.
    """
    obs = adata.obs
    if not isinstance(obs, pd.DataFrame):
        return adata
    if "target_assignment" in adata.obsm:
        ta = adata.obsm["target_assignment"]
        if isinstance(ta, pd.DataFrame):
            # Drop any columns already present in obs to avoid duplicates.
            overlap = set(ta.columns) & set(obs.columns)
            if overlap:
                ta = ta.drop(columns=list(overlap))
            if not ta.empty:
                adata.obs = pd.concat([obs, ta], axis="columns")
    return adata


def _enrich_cell_geometry(adata: ad.AnnData) -> ad.AnnData:
    """Merge ``.obsm['cell_geometry']`` columns into the main ``X`` matrix (in-place).

    The geometry columns (e.g. ``Area``, ``NuclearArea``, ``Xum_tile``, ``Yum_tile``,
    ``AreaUm``) are appended as new variables to ``X`` and registered in
    ``.var``.

    Parameters
    ----------
    adata
        Table to enrich.

    Returns
    -------
    ad.AnnData
        The same object passed in, optionally modified in-place.

    Notes
    -----
    No-op if ``cell_geometry`` is missing from ``obsm`` or not a DataFrame.
    """
    if "cell_geometry" not in adata.obsm:
        return adata
    cg = adata.obsm["cell_geometry"]
    if not isinstance(cg, pd.DataFrame) or cg.empty:
        return adata

    cg_array = cg.to_numpy(dtype=np.float32)

    # Build a geometry-only AnnData with matching var columns so that
    # ad.concat preserves the original var metadata instead of dropping it.
    geom_var = pd.DataFrame(index=cg.columns)
    for col in adata.var.columns:
        dtype = adata.var[col].dtype
        if pd.api.types.is_bool_dtype(dtype):
            geom_var[col] = False
        else:
            geom_var[col] = adata.var[col].dtype.type() if pd.api.types.is_numeric_dtype(dtype) else ""
    if "feature_class" in geom_var.columns:
        geom_var["feature_class"] = "geometry"

    geom_adata = ad.AnnData(
        X=cg_array,
        obs=adata.obs[[]],
        var=geom_var,
    )
    merged = ad.concat([adata, geom_adata], axis=1, merge="first")
    merged.obs = adata.obs
    merged.uns = adata.uns
    merged.obsm = adata.obsm
    return merged


# ---------------------------------------------------------------------------
# Multi-store orchestration
# ---------------------------------------------------------------------------


def _open_stores(
    inputs: StoreInput | Sequence[StoreInput],
    storage_options: dict[str, object] | None = None,
) -> dict[str, RunInfo]:
    """Open one or many store inputs, returning ``{key: RunInfo}``.

    Normalizes the input to a list, rejects duplicates, selects the
    appropriate backend for each input, derives a human-readable run ID,
    and assigns unique keys (appending ``__N`` on collision).
    """
    # Normalize to list.
    if isinstance(inputs, (str, Path, sd.SpatialData)):
        input_list: list[StoreInput] = [inputs]
    else:
        input_list = list(inputs)
    if not input_list:
        raise ValueError("No stores provided.")

    # Reject duplicate inputs.
    seen_inputs: dict[str, int] = {}
    for idx, inp in enumerate(input_list):
        dup_key = f"sdata@{id(inp)}" if isinstance(inp, sd.SpatialData) else str(inp).rstrip("/")
        if dup_key in seen_inputs:
            raise ValueError(
                f"Duplicate store at positions {seen_inputs[dup_key]} and {idx}: {inp!r}. "
                "Concatenating the same data twice is not allowed."
            )
        seen_inputs[dup_key] = idx

    # Open each store and assign unique keys.
    seen_names: Counter[str] = Counter()
    result: dict[str, RunInfo] = {}
    for inp in input_list:
        store = make_backend(inp, storage_options=storage_options)
        seen_names[store.run_name] += 1
        key = (
            store.run_name
            if seen_names[store.run_name] == 1
            else f"{store.run_name}__{seen_names[store.run_name]}"
        )
        result[key] = store
    return result


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------


class Loader:
    """Multi-run data loader with metadata discovery and explicit read APIs.

    Opens each store on construction and exposes it as a :class:`RunInfo`,
    which reads ``global-table-metadata`` on first access and caches it.
    Optionally validates the points schema up front (``validate_points=True``).

    Parameters
    ----------
    stores
        One or many store inputs (local path, S3 URI, ``.zarr.zip``, or
        ``SpatialData`` object).  Duplicates are rejected.
    storage_options
        Extra keyword arguments forwarded to ``s3fs.S3FileSystem`` and
        ``fsspec`` for S3 backends (e.g. ``{"profile": "my-profile"}``).
        Ignored for local backends.
    validate_points
        If ``True``, sample one points element per batch per run at construction
        and check required columns. Defaults to ``False`` because this can be
        slow for stores with many batches (e.g. 3′ runs) or remote / zip-backed
        paths. Set ``True`` to surface schema issues eagerly; otherwise invalid
        points will only surface when calling :meth:`load_points`.
    """

    def __init__(
        self,
        stores: StoreInput | Sequence[StoreInput],
        storage_options: dict[str, object] | None = None,
        *,
        validate_points: bool = False,
    ) -> None:
        self._runs = _open_stores(stores, storage_options=storage_options)
        if validate_points:
            self._validate_points_schema()

    def __repr__(self) -> str:
        runs = list(self._runs.keys())
        tables = sorted(self.available_tables)
        return f"Loader(runs={runs}, tables={tables})"

    def summary(self) -> str:
        """Return a multi-line human-readable summary of all runs."""
        blocks: list[str] = []
        for key, run_info in self._runs.items():
            header = f"--- {key} ---"
            body = "\n".join(f"  {line}" for line in run_info.summary().splitlines())
            blocks.append(f"{header}\n{body}")
        return "\n\n".join(blocks)

    # -------------------------------------------------------------------
    # Discovery
    # -------------------------------------------------------------------

    @property
    def available_tables(self) -> set[str]:
        """Table names loadable via :meth:`load_tables` (union across all runs)."""
        tables: set[str] = set()
        for run in self._runs.values():
            tables |= run.available_tables
        return tables

    @property
    def available_target_types(self) -> set[str]:
        """Target types from table metadata (union across all runs).

        Useful for discovering valid values for the ``target_type`` parameter
        in :meth:`load_points`.
        """
        types: set[str] = set()
        for run in self._runs.values():
            types |= run.available_target_types
        return types

    # -------------------------------------------------------------------
    # Validation
    # -------------------------------------------------------------------

    def _validate_points_schema(self) -> None:
        """Check points elements per run have the required columns.

        Samples one element per batch to avoid checking every element.
        Structural columns (``x``, ``y``, ``tile_id``, ``well_id``) are
        always required.  ``target_name`` may be absent for ThreePrime/OPS
        batches where the upstream alignment pipeline was not run.
        """
        for run_key, run_info in self._runs.items():
            # Sample one element per batch.
            checked: set[str | None] = set()
            for pname in run_info.point_names:
                _, batch = split_points_name(pname)
                if batch in checked:
                    continue
                checked.add(batch)
                lf = run_info._read_points(pname)
                columns = set(lf.collect_schema().names())
                missing_structural = _ALWAYS_REQUIRED_POINTS_COLUMNS - columns
                if missing_structural:
                    warnings.warn(
                        f"Run {run_key!r}: points element {pname!r} is missing "
                        f"structural columns {sorted(missing_structural)}. "
                        "Points loading will fail.",
                        UserWarning,
                        stacklevel=2,
                    )
                elif "target_name" not in columns:
                    warnings.warn(
                        f"Run {run_key!r}: points element {pname!r} (batch "
                        f"{batch!r}) is missing 'target_name'. This batch's "
                        "points may not be loadable by target type.",
                        UserWarning,
                        stacklevel=2,
                    )

    # -------------------------------------------------------------------
    # Load tables
    # -------------------------------------------------------------------

    def load_tables(
        self,
        *,
        table: str,
        concatenated: bool = True,
        lazy: bool = True,
        include_target_assignment: bool = False,
        include_cell_geometry: bool = False,
    ) -> _LoadedTablesResult:
        """Load a table from one or many runs.

        Parameters
        ----------
        table
            Table to load (e.g. ``"Transcript"``).  Resolved to
            ``global-table-{table}`` within each store.  Imaging is aliased:
            ``"imaging"`` and ``"CellPaint"`` each resolve to whichever key the
            store actually contains, so legacy (``global-table-CellPaint``) and
            current (``global-table-imaging``) stores load through one call.
        concatenated
            If ``True``, concatenate loaded tables across runs.
        lazy
            If ``True``, use lazy backend reads where available.
        include_target_assignment
            Merge ``obsm['target_assignment']`` into ``obs`` (OPS tables).
        include_cell_geometry
            If ``False`` (default), keep ``obsm['cell_geometry']`` as read from
            the store. If ``True``, merge those columns into ``X`` (new
            ``.var`` entries). Emits a ``UserWarning`` when ``True`` for
            non-imaging tables.

        Returns
        -------
        _LoadedTablesResult
            - ``AnnData`` when a single run is loaded (regardless of
              *concatenated*) or when ``concatenated=True``,
            - ``dict[str, AnnData]`` (keyed by run) when multiple runs are
              loaded and ``concatenated=False``.
        """
        loaded, source_store_lookup = _load_tables_core(
            list(self._runs.items()),
            table=table,
            lazy=lazy,
            include_target_assignment=include_target_assignment,
            include_cell_geometry=include_cell_geometry,
        )
        if not loaded:
            raise ValueError(
                f"No tables were loaded for table={table!r}. "
                f"Available tables: {sorted(self.available_tables)}"
            )

        if len(loaded) == 1:
            return next(iter(loaded.values()))

        if not concatenated:
            return loaded

        return _concat_tables_rowwise(loaded, source_stores=source_store_lookup)

    # -------------------------------------------------------------------
    # Export tables
    # -------------------------------------------------------------------

    def export_tables(
        self,
        out_dir: str | Path,
        *,
        tables: Sequence[str] | None = None,
        fmt: str = "h5ad",
        overwrite: bool = False,
        include_metadata: bool = True,
    ) -> list[Path]:
        """Export each table to a local file, one per ``(run, table)``.

        The AnnData tables are a tiny fraction of an SDO (images / points dominate),
        so this pulls just the tables locally, e.g. from S3, without copying the whole
        store. Reuses :meth:`load_tables`, so ``obs`` / ``var`` come back as pandas and
        the canonical enrichment applies. Tables are materialised one at a time.

        Files are named ``{run}__{table}.{ext}``; the source run / store / build version
        are stamped into ``uns['export_provenance']``. Spatial coordinates are preserved,
        imagery is not, so table-level analysis works locally but image overlays need the
        full store.

        Parameters
        ----------
        out_dir
            Directory to write into (created if absent).
        tables
            Table names to export; default is all of :attr:`available_tables`.
        fmt
            ``"h5ad"`` (default) or ``"zarr"`` (streamable for very large tables).
        overwrite
            Overwrite existing files; otherwise an existing file is skipped (with a warning).
        include_metadata
            Also export the small ``global-table-metadata`` table (runstats, provenance,
            targets). Default ``True``.

        Returns
        -------
        list[pathlib.Path]
            The files written (or skipped).
        """
        if fmt not in ("h5ad", "zarr"):
            raise ValueError(f"fmt must be 'h5ad' or 'zarr', got {fmt!r}")
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        ext = "h5ad" if fmt == "h5ad" else "zarr"

        names = sorted(self.available_tables) if tables is None else list(tables)

        written: list[Path] = []
        for name in names:
            result = self.load_tables(table=name, concatenated=False, lazy=False)
            if isinstance(result, dict):
                per_run = result
            else:
                run_key = next(
                    (k for k, ri in self._runs.items() if name in ri.available_tables),
                    next(iter(self._runs)),
                )
                per_run = {run_key: result}
            for run_key, adata in per_run.items():
                written.append(self._export_write(out, run_key, name, adata, ext, overwrite))

        if include_metadata:
            from elembio_spatialdata_analysis.load.backends import _METADATA_TABLE

            for run_key, run_info in self._runs.items():
                try:
                    meta = run_info._read_table(_METADATA_TABLE, lazy=False)
                except Exception:
                    continue
                written.append(self._export_write(out, run_key, "metadata", meta, ext, overwrite))

        return written

    def _export_write(
        self, out: Path, run_key: str, name: str, adata: ad.AnnData, ext: str, overwrite: bool
    ) -> Path:
        """Stamp provenance and write one table; skip (warn) if it exists and not overwriting."""
        run_info = self._runs[run_key]
        try:
            build_version = (run_info.build_provenance or {}).get("pipeline_version")
        except Exception:
            build_version = None
        adata.uns["export_provenance"] = {
            "run": run_key,
            "source": str(run_info._source),
            "table": name,
            "build_pipeline_version": build_version,
        }
        path = out / f"{run_key}__{name}.{ext}"
        if path.exists() and not overwrite:
            warnings.warn(
                f"{path.name} already exists; skipping (pass overwrite=True to replace).",
                stacklevel=2,
            )
            return path
        if ext == "h5ad":
            adata.write_h5ad(path)
        else:
            adata.write_zarr(path)
        return path

    # -------------------------------------------------------------------
    # Load points
    # -------------------------------------------------------------------

    def load_points(
        self,
        *,
        target_type: str | Sequence[str] | None = None,
        batch: str | Sequence[str] | None = None,
        well: str | Sequence[str] | None = None,
        tile: str | Sequence[str] | None = None,
        point_name: str | Sequence[str] | None = None,
        concatenated: bool = True,
        lazy: bool = True,
    ) -> _LoadedPointsResult:
        """Load points from one or many runs with composable filters.

        Parameters
        ----------
        target_type
            Select points whose batch maps to this target type.
        batch
            Select points from specific batch(es).
        well
            Select points from specific well(s).  In tile mode, resolves
            wells to their constituent tiles via the tiles metadata.
        tile
            Select points from specific tile(s).  Meaningful in tile mode;
            in well mode returns no matches.
        point_name
            Explicit element name(s); overrides other filters.
        concatenated
            If ``True``, concatenate all loaded frames into one.
        lazy
            If ``True``, use lazy backend reads where available.

        Returns
        -------
        _LoadedPointsResult
            - ``dict[str, dict[str, LazyFrame]]`` when ``concatenated=False``,
            - ``LazyFrame`` when ``concatenated=True``.
        """
        loaded: _LoadedPoints = {}

        for run_key, run_info in self._runs.items():
            selected = run_info.select_points(
                target_type=target_type,
                batch=batch,
                well=well,
                tile=tile,
                point_name=point_name,
            )
            if not selected:
                hints: list[str] = []
                if target_type is not None:
                    hints.append(
                        f"target_type={target_type!r} (available: {sorted(run_info.available_target_types)})"
                    )
                if batch is not None:
                    hints.append(f"batch={batch!r} (available: {run_info.batch_names})")
                if well is not None:
                    hints.append(f"well={well!r} (available: {run_info.well_names})")
                if tile is not None:
                    hints.append(f"tile={tile!r}")
                if point_name is not None:
                    hints.append(f"point_name={point_name!r}")
                filter_detail = "; ".join(hints) if hints else "no filters"
                warnings.warn(
                    f"No matching points in run {run_key!r} for {filter_detail}.",
                    UserWarning,
                    stacklevel=2,
                )
                loaded[run_key] = {}
                continue

            run_points: dict[str, pl.LazyFrame] = {}
            for pname in selected:
                lf = run_info._read_points(pname, lazy=lazy)
                _, batch_name = split_points_name(pname)
                lf = lf.with_columns(
                    pl.lit(batch_name or "unknown").alias("batch_name"),
                    pl.lit(run_info.run_name).alias("run_name"),
                )
                run_points[pname] = lf
            loaded[run_key] = run_points

        if not concatenated:
            return loaded

        all_frames: dict[str, pl.LazyFrame] = {}
        for run_key, points_by_name in loaded.items():
            for pname, lf in points_by_name.items():
                all_frames[f"{run_key}/{pname}"] = lf

        if not all_frames:
            raise ValueError("No points were loaded for concatenation.")

        return _concat_frames(all_frames)
