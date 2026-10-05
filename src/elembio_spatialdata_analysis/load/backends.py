"""RunInfo ABC, concrete storage backends, and multi-store factory.

:class:`RunInfo` is the public-facing handle for a single run.  It defines
the abstract IO interface that concrete backends implement, and provides
lazily cached metadata, point selection, and discovery properties on top.

Concrete backends:
- :class:`_LocalDirBackend` — ``.zarr`` directory on local filesystem
- :class:`_LocalZipBackend` — ``.zarr.zip`` file on local filesystem
- :class:`_S3DirBackend` — ``.zarr`` directory on S3
- :class:`_S3ZipBackend` — ``.zarr.zip`` file on S3 (lazy per-table cache)
- :class:`_SDataBackend` — in-memory ``SpatialData`` object
"""

from __future__ import annotations

import abc
import atexit
import hashlib
import json
import logging
import os
import shutil
import time
import warnings
import zipfile
from collections.abc import Sequence
from functools import cached_property
from pathlib import Path
from typing import Any, ClassVar, Literal

import anndata as ad
import fsspec
import numpy as np
import pandas as pd
import polars as pl
import pyarrow as pa
import s3fs
import spatialdata as sd
import zarr
from zarr.storage import ZipStore

_logger = logging.getLogger(__name__)

try:
    from anndata.experimental import read_lazy as _read_lazy
except Exception:  # pragma: no cover
    _read_lazy: Any = None


def _materialize_obs_var(adata: ad.AnnData) -> ad.AnnData:
    """Ensure ``obs`` / ``var`` are pandas DataFrames, keeping ``X`` lazy.

    ``anndata.experimental.read_lazy`` returns ``obs`` / ``var`` as xarray
    ``Dataset2D`` objects (StringDType-backed). Every downstream pandas consumer
    (``pp.filter``, ``vis``, ``tl.morphology``) breaks on the first
    ``.astype(str)`` / ``.str`` / ``.fillna`` op against those. ``obs`` / ``var``
    are small, so materialise them here while leaving ``X`` and layers lazy.
    """
    if not isinstance(adata.var, pd.DataFrame):
        adata.var = pd.DataFrame(
            {c: np.asarray(adata.var[c]) for c in list(adata.var.columns)},
            index=pd.Index(np.asarray(adata.var_names)),
        )
    if not isinstance(adata.obs, pd.DataFrame):
        adata.obs = pd.DataFrame(
            {c: np.asarray(adata.obs[c]) for c in list(adata.obs.columns)},
            index=pd.Index(np.asarray(adata.obs_names)),
        )
    return adata


StoreInput = sd.SpatialData | str | Path
"""A single data source accepted by the public API."""

_StorageOptions = dict[str, Any]
"""Kwargs forwarded to ``s3fs.S3FileSystem`` / ``fsspec`` for S3 access."""

_METADATA_TABLE = "global-table-metadata"
_RUN_PARAMS_PATH = "metadata/inputs/run/RunParameters.json"
_PERSISTENT_CACHE_ROOT = Path(
    os.environ.get(
        "ELEMBIO_SPATIALDATA_CACHE",
        str(Path.home() / ".cache" / "elembio-spatialdata"),
    )
)
_S3ZIP_LOCK_POLL_SEC = 0.1


def _parse_float_env(name: str, default: float) -> float:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return float(raw)
    except ValueError:
        return default


_S3ZIP_LOCK_STALE_SEC = _parse_float_env("ELEMBIO_S3ZIP_LOCK_STALE_SEC", 120.0)
_S3ZIP_LOCK_MAX_WAIT_SEC = _parse_float_env("ELEMBIO_S3ZIP_LOCK_MAX_WAIT_SEC", 1800.0)
_S3ZIP_LOCK_HEARTBEAT_SEC = _parse_float_env("ELEMBIO_S3ZIP_LOCK_HEARTBEAT_SEC", 30.0)
_S3ZIP_PRUNE_MIN_AGE_SEC = _parse_float_env("ELEMBIO_S3ZIP_PRUNE_MIN_AGE_SEC", 3600.0)


def _read_zarr_metadata_attrs_group(ma_group: zarr.Group) -> dict:
    """Read scalar values from an open ``metadata_attrs`` zarr group."""
    result: dict = {}
    for key in list(ma_group):  # ty: ignore[invalid-argument-type]
        val = ma_group[key][()]  # ty: ignore[invalid-argument-type,not-subscriptable]
        if hasattr(val, "item") and getattr(val, "ndim", None) == 0:
            val = val.item()  # ty: ignore[call-non-callable]
        result[key] = val
    return result


def _read_zarr_metadata_attrs(root: zarr.Group) -> dict:
    """Read ``uns['metadata_attrs']`` scalar values from an open zarr group."""
    if "uns" not in root or "metadata_attrs" not in root["uns"]:
        return {}
    ma_group = root["uns"]["metadata_attrs"]  # ty: ignore[invalid-argument-type]
    return _read_zarr_metadata_attrs_group(ma_group)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def split_points_name(points_name: str) -> tuple[str | None, str | None]:
    """Split ``'{prefix}-points-{batch_name}'`` into ``(prefix, batch_name)``.

    Returns ``(None, None)`` when the name does not match the pattern.
    """
    marker = "-points-"
    if marker not in points_name:
        return None, None
    prefix, batch = points_name.split(marker, 1)
    return prefix, batch


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------


def make_backend(
    store_input: StoreInput,
    *,
    storage_options: _StorageOptions | None = None,
) -> RunInfo:
    """Create the appropriate backend for a single store input."""
    if isinstance(store_input, sd.SpatialData):
        return _SDataBackend(store_input)

    path_str = str(store_input)
    is_s3 = path_str.startswith("s3://")
    is_zip = path_str.endswith(".zarr.zip")

    if is_s3 and is_zip:
        return _S3ZipBackend(path_str, storage_options=storage_options)
    if is_s3:
        return _S3DirBackend(path_str, storage_options=storage_options)
    if is_zip:
        return _LocalZipBackend(Path(path_str))
    return _LocalDirBackend(Path(path_str))


# ---------------------------------------------------------------------------
# ABC
# ---------------------------------------------------------------------------


class RunInfo(abc.ABC):
    """Per-run metadata, point selection, and IO — the public-facing handle.

    Each concrete backend (local dir, zip, S3, in-memory) inherits from this
    ABC and implements the abstract IO methods.  All derived metadata
    (DataFrames, lookup dicts, point mode) is lazily computed and cached.

    Users obtain ``RunInfo`` instances from :class:`Loader` via
    ``loader._runs`` and interact with the discovery properties and
    ``select_points()`` method.
    """

    def __init__(self) -> None:
        self.run_name: str = self._derive_run_name()
        # Lazy caches (populated on first access).
        self._metadata_cache: dict | None = None
        self._metadata_frames_cache: dict[str, pd.DataFrame] | None = None
        self._points_by_prefix_cache: dict[str, list[str]] | None = None
        self._points_by_batch_cache: dict[str, list[str]] | None = None
        self._point_mode_cache: Literal["well", "tile"] | None = None
        self._tile_to_well_cache: dict[str, str] | None = None
        self._batch_to_target_types_cache: dict[str, set[str]] | None = None

    def __repr__(self) -> str:
        return f"{type(self).__name__}(run_name={self.run_name!r}, source={self._source!r})"

    def summary(self) -> str:
        """Return a multi-line human-readable summary of this run."""
        lines = [
            f"run_name:       {self.run_name}",
            f"source:         {self._source}",
            f"tables:         {sorted(self.available_tables)}",
        ]

        bt = self.batch_types
        if bt:
            batch_strs = [f"{name} ({btype})" for name, btype in sorted(bt.items())]
            lines.append(f"batches ({len(bt)}):   {', '.join(batch_strs)}")
        else:
            lines.append("batches:        (none)")

        wells = self.well_names
        if wells:
            preview = ", ".join(wells[:8])
            suffix = f", ... ({len(wells)} total)" if len(wells) > 8 else ""
            lines.append(f"wells ({len(wells)}):     {preview}{suffix}")
        else:
            lines.append("wells:          (none)")

        tiles = self.tiles
        lines.append(f"tiles:          {len(tiles)}")

        points = self.point_names
        lines.append(f"points:         {len(points)} elements ({self._point_mode} mode)")

        prov = self.build_provenance
        if prov:
            version = prov.get("pipeline_version", "unknown")
            lines.append(f"build:          {version}")

        return "\n".join(lines)

    # ------------------------------------------------------------------
    # Source and run-name derivation
    # ------------------------------------------------------------------

    @property
    @abc.abstractmethod
    def _source(self) -> str:
        """Human-readable source label (path, URL, or object id)."""
        ...

    def _derive_run_name(self) -> str:
        """Derive a human-readable run identifier from :attr:`_source`.

        Strips ``.zarr`` / ``.zarr.zip`` suffixes and path components.
        Subclasses that can cheaply read ``RunParameters.json`` override
        this to try ``RunName`` first.
        """
        source = self._source.rstrip("/")
        name = source.rsplit("/", 1)[-1]
        if name.endswith(".zarr.zip"):
            return name[: -len(".zarr.zip")]
        if name.endswith(".zarr"):
            return name[: -len(".zarr")]
        return name

    def source_exists(self) -> bool:
        """Cheap reachability check for the underlying store.

        Returns ``True`` when :attr:`_source` resolves to something readable
        (an S3 object/prefix or a local path that exists).  It does **not**
        assert the store holds valid tables — only that the source is there.

        Callers use this to fail fast on a missing/unreadable source (a
        typo'd path, wrong bucket or project folder) instead of silently
        producing empty results.  A present-but-empty store still returns
        ``True`` and is handled downstream by the metadata/table warnings.
        In-memory backends have no path and override this to ``True``.
        """
        try:
            fs, path = fsspec.core.url_to_fs(
                self._source, **dict(getattr(self, "_storage_options", {}) or {})
            )
            return bool(fs.exists(path))
        except Exception:
            return False

    @cached_property
    def run_date(self) -> str | None:
        """Instrument run start date, lazily read from sdata attrs on first access.

        Kept lazy so that constructing a remote backend does not trigger
        network I/O. The underlying read is implemented by
        :meth:`_derive_run_date`.
        """
        return self._derive_run_date()

    def _derive_run_date(self) -> str | None:
        """Read instrument run start date from ``SpatialData.attrs['run_metadata']``.

        Returns the raw ISO-8601 string from
        ``run_metadata['run_parameters']['Date']``, or ``None`` when the
        attribute is missing or unreadable. Every backend implements
        ``_read_sdata_attrs``, so the base implementation works for all.
        """
        try:
            attrs = self._read_sdata_attrs()
            rm_blob = attrs.get("run_metadata")
            if not rm_blob:
                return None
            rm = json.loads(rm_blob) if isinstance(rm_blob, str) else rm_blob
            rp = (rm.get("run_parameters") if isinstance(rm, dict) else None) or {}
            date_val = rp.get("Date")
            return str(date_val) if date_val else None
        except Exception:
            return None

    # ------------------------------------------------------------------
    # Raw metadata (lazy, cached)
    # ------------------------------------------------------------------

    @property
    def _metadata(self) -> dict:
        """``global-table-metadata.uns`` dict, read once and cached."""
        if self._metadata_cache is None:
            self._metadata_cache = self._read_metadata()
        return self._metadata_cache

    def _read_metadata(self) -> dict:
        if not self._has_table(_METADATA_TABLE):
            warnings.warn(
                f"Store {self._source!r} has no '{_METADATA_TABLE}' table. "
                "Batch, well, tile, and target metadata will be empty.",
                UserWarning,
                stacklevel=2,
            )
            return {}
        try:
            adata = self._read_table(_METADATA_TABLE, lazy=False)
            return dict(adata.uns)
        except Exception:
            _logger.warning("Failed to read metadata from %s", self._source, exc_info=True)
            return {}

    _METADATA_FRAME_KEYS = (
        "batches",
        "wells",
        "tiles",
        "targets_barcoding",
        "targets_imaging",
        "diss_targets_metadata",
    )

    # ------------------------------------------------------------------
    # Metadata frames (lazy, cached)
    # ------------------------------------------------------------------

    def _metadata_frames(self) -> dict[str, pd.DataFrame]:
        """Extract DataFrames from ``global-table-metadata.uns``."""
        if self._metadata_frames_cache is None:
            meta = self._metadata
            self._metadata_frames_cache = {
                name: value.copy() if isinstance((value := meta.get(name)), pd.DataFrame) else pd.DataFrame()
                for name in self._METADATA_FRAME_KEYS
            }
        return self._metadata_frames_cache

    @property
    def batches(self) -> pd.DataFrame:
        """One row per batch from ``global-table-metadata``."""
        return self._metadata_frames()["batches"]

    @property
    def wells(self) -> pd.DataFrame:
        """One row per well from ``global-table-metadata``."""
        return self._metadata_frames()["wells"]

    @property
    def tiles(self) -> pd.DataFrame:
        """One row per tile from ``global-table-metadata``."""
        return self._metadata_frames()["tiles"]

    @property
    def targets_barcoding(self) -> pd.DataFrame:
        """Per-target metadata for barcoding batches."""
        return self._metadata_frames()["targets_barcoding"]

    @property
    def targets_imaging(self) -> pd.DataFrame:
        """Per-target metadata for imaging batches."""
        return self._metadata_frames()["targets_imaging"]

    @property
    def diss_targets_metadata(self) -> pd.DataFrame:
        """Per-target metadata for DISS/OPS batches."""
        return self._metadata_frames()["diss_targets_metadata"]

    @property
    def build_provenance(self) -> dict:
        """Build pipeline version and parameters."""
        return dict(self._metadata.get("build_provenance", {}))

    # ------------------------------------------------------------------
    # Derived lookups (lazy, cached)
    # ------------------------------------------------------------------

    def _ensure_point_lookups(self) -> None:
        """Parse point names and cache prefix/batch dicts (once)."""
        if self._points_by_prefix_cache is not None:
            return
        by_prefix: dict[str, list[str]] = {}
        by_batch: dict[str, list[str]] = {}
        for pname in self._list_points():
            prefix, batch_name = split_points_name(pname)
            if prefix is not None:
                by_prefix.setdefault(prefix, []).append(pname)
            if batch_name is not None:
                by_batch.setdefault(batch_name, []).append(pname)
        self._points_by_prefix_cache = by_prefix
        self._points_by_batch_cache = by_batch

    @property
    def _points_by_prefix(self) -> dict[str, list[str]]:
        self._ensure_point_lookups()
        assert self._points_by_prefix_cache is not None
        return self._points_by_prefix_cache

    @property
    def _points_by_batch(self) -> dict[str, list[str]]:
        self._ensure_point_lookups()
        assert self._points_by_batch_cache is not None
        return self._points_by_batch_cache

    @property
    def _point_mode(self) -> Literal["well", "tile"]:
        """Auto-detected: ``"well"`` or ``"tile"``."""
        if self._point_mode_cache is None:
            prefixes = set(self._points_by_prefix.keys())
            wells_df = self.wells
            tiles_df = self.tiles
            well_names = set(wells_df.index.tolist()) if not wells_df.empty else set()
            tile_names = set(tiles_df.index.tolist()) if not tiles_df.empty else set()
            if prefixes and prefixes <= tile_names and not prefixes <= well_names:
                self._point_mode_cache = "tile"
            else:
                self._point_mode_cache = "well"
        return self._point_mode_cache

    @property
    def _tile_to_well(self) -> dict[str, str]:
        if self._tile_to_well_cache is None:
            tiles_df = self.tiles
            if tiles_df.empty or "well_name" not in tiles_df.columns:
                self._tile_to_well_cache = {}
            else:
                self._tile_to_well_cache = {
                    str(tile_name): str(row["well_name"]) for tile_name, row in tiles_df.iterrows()
                }
        return self._tile_to_well_cache

    @property
    def _batch_to_target_types(self) -> dict[str, set[str]]:
        if self._batch_to_target_types_cache is None:
            result: dict[str, set[str]] = {}
            for tname in self._list_tables():
                if tname == _METADATA_TABLE:
                    continue
                try:
                    attrs = self._read_table_metadata_attrs(tname)
                except Exception:
                    continue
                target_type = attrs.get("target_type")
                batch_names_val = attrs.get("batch_name")
                if target_type is None or batch_names_val is None:
                    continue
                if isinstance(batch_names_val, str):
                    batch_names_val = [batch_names_val]
                for bn in batch_names_val:
                    result.setdefault(str(bn), set()).add(str(target_type))
            self._batch_to_target_types_cache = result
        return self._batch_to_target_types_cache

    # ------------------------------------------------------------------
    # Public properties
    # ------------------------------------------------------------------

    @property
    def table_names(self) -> list[str]:
        """All table names in this store."""
        return self._list_tables()

    @property
    def point_names(self) -> list[str]:
        """All point element names in this store."""
        return self._list_points()

    @property
    def batch_names(self) -> list[str]:
        """Sorted batch names from :attr:`batches`."""
        if self.batches.empty:
            return []
        return sorted(self.batches["batch_name"].tolist())

    @property
    def batch_types(self) -> dict[str, str]:
        """``{batch_name: batch_type}`` from :attr:`batches`."""
        if self.batches.empty:
            return {}
        return dict(zip(self.batches["batch_name"], self.batches["batch_type"], strict=False))

    @property
    def well_names(self) -> list[str]:
        """Sorted well identifiers from :attr:`wells`."""
        if self.wells.empty:
            return []
        return sorted(self.wells.index.tolist())

    @property
    def available_tables(self) -> set[str]:
        """Table suffixes present in this store.

        Extracts suffixes from ``global-table-{name}`` entries,
        filtering out only the metadata table.
        """
        prefix = "global-table-"
        return {
            t.removeprefix(prefix)
            for t in self._list_tables()
            if t.startswith(prefix) and t != "global-table-metadata"
        }

    @property
    def available_target_types(self) -> set[str]:
        """Target types derived from table metadata.

        Unlike :attr:`available_tables` (which returns all table suffixes),
        this returns only values recorded as ``target_type`` in each table's
        ``uns['metadata_attrs']``.  Useful for filtering points by target type.
        """
        types: set[str] = set()
        for batch_types in self._batch_to_target_types.values():
            types |= batch_types
        return types

    # ------------------------------------------------------------------
    # Point selection
    # ------------------------------------------------------------------

    def batches_for_target_type(self, target_type: str) -> list[str]:
        """Batch names whose targets include *target_type*."""
        return sorted(b for b, types in self._batch_to_target_types.items() if target_type in types)

    def table_suffixes_for_target_type(self, target_type: str | tuple[str, ...] | list[str]) -> list[str]:
        """Table suffixes whose ``uns['metadata_attrs']['target_type']`` matches.

        Used by modalities like OPS that emit one ``global-table-{batch}-{sample}``
        table per target-site instead of a single aggregated
        ``global-table-{TargetType}``.  Returns the suffixes (e.g.
        ``"B01-A549"``) in sorted order — ``[]`` when no matching
        tables exist.

        ``target_type`` accepts a single string or an iterable of
        aliases.  OPS in particular shows up in metadata as either
        ``"OPS"`` or ``"SpecializedTargeted"`` depending on the
        pipeline version — passing ``("OPS", "SpecializedTargeted")``
        catches both.

        Discovery reads ``uns['metadata_attrs']`` from every
        ``global-table-*`` entry (skipping the metadata table itself),
        so the result is independent of name shape — works equally well
        for per-site (``B01-A549``) and aggregated (``Transcript``)
        layouts when the metadata is set.
        """
        accepted = {target_type} if isinstance(target_type, str) else set(target_type)
        suffixes: list[str] = []
        for tname in self._list_tables():
            if tname == _METADATA_TABLE or not tname.startswith("global-table-"):
                continue
            try:
                attrs = self._read_table_metadata_attrs(tname)
            except Exception:
                continue
            if attrs.get("target_type") in accepted:
                suffixes.append(tname.removeprefix("global-table-"))
        return sorted(suffixes)

    def points_for_target_type(self, target_type: str) -> list[str]:
        """Point element names whose batch maps to *target_type*."""
        result: list[str] = []
        for b in self.batches_for_target_type(target_type):
            result.extend(self._points_by_batch.get(b, []))
        return sorted(result)

    def points_for_batch(self, batch_name: str) -> list[str]:
        """Point element names for a given batch."""
        return list(self._points_by_batch.get(batch_name, []))

    def points_for_well(self, well_id: str) -> list[str]:
        """Point element names for a given well.

        In well mode the lookup is direct.  In tile mode the tiles metadata
        is used to resolve which tiles belong to the well.
        """
        if self._point_mode == "well":
            return list(self._points_by_prefix.get(well_id, []))
        tile_names = [t for t, w in self._tile_to_well.items() if w == well_id]
        result: list[str] = []
        for t in tile_names:
            result.extend(self._points_by_prefix.get(t, []))
        return sorted(result)

    def points_for_tile(self, tile_name: str) -> list[str]:
        """Point element names for a given tile.

        In tile mode the lookup is direct.  In well mode this returns ``[]``.
        """
        if self._point_mode == "tile":
            return list(self._points_by_prefix.get(tile_name, []))
        return []

    def select_points(
        self,
        *,
        target_type: str | Sequence[str] | None = None,
        batch: str | Sequence[str] | None = None,
        well: str | Sequence[str] | None = None,
        tile: str | Sequence[str] | None = None,
        point_name: str | Sequence[str] | None = None,
    ) -> list[str]:
        """Select point element names matching the given filters.

        Filters compose as intersection.  When *point_name* is given it acts
        as an explicit override and the other filters are ignored.  With no
        filters at all, returns every point name in the run.
        """
        if point_name is not None:
            names = [point_name] if isinstance(point_name, str) else list(point_name)
            available = set(self.point_names)
            unknown = set(names) - available
            if unknown:
                warnings.warn(
                    f"Point(s) {sorted(unknown)!r} not found in run {self.run_name!r}; skipping.",
                    UserWarning,
                    stacklevel=2,
                )
            return sorted(n for n in names if n in available)

        candidates = set(self.point_names)

        if target_type is not None:
            types = {target_type} if isinstance(target_type, str) else set(target_type)
            unknown_types = types - self.available_target_types
            if unknown_types:
                warnings.warn(
                    f"Run {self.run_name!r}: target type(s) {sorted(unknown_types)!r} "
                    f"not recognised. Available: {sorted(self.available_target_types)}",
                    UserWarning,
                    stacklevel=2,
                )
            matching_batches = {b for b, tt in self._batch_to_target_types.items() if tt & types}
            type_points: set[str] = set()
            for b in matching_batches:
                type_points.update(self._points_by_batch.get(b, []))
            candidates &= type_points

        if batch is not None:
            batches = [batch] if isinstance(batch, str) else list(batch)
            known_batches = set(self.batch_names)
            unknown_batches = set(batches) - known_batches
            if unknown_batches:
                warnings.warn(
                    f"Run {self.run_name!r}: batch(es) {sorted(unknown_batches)!r} "
                    f"not recognised. Available: {sorted(known_batches)}",
                    UserWarning,
                    stacklevel=2,
                )
            batch_points: set[str] = set()
            for b in batches:
                batch_points.update(self._points_by_batch.get(b, []))
            candidates &= batch_points

        if well is not None:
            wells = [well] if isinstance(well, str) else list(well)
            known_wells = set(self.well_names)
            unknown_wells = set(wells) - known_wells
            if unknown_wells:
                warnings.warn(
                    f"Run {self.run_name!r}: well(s) {sorted(unknown_wells)!r} "
                    f"not recognised. Available: {sorted(known_wells)}",
                    UserWarning,
                    stacklevel=2,
                )
            well_points: set[str] = set()
            for w in wells:
                well_points.update(self.points_for_well(w))
            candidates &= well_points

        if tile is not None:
            if self._point_mode == "well":
                warnings.warn(
                    f"Run {self.run_name!r} uses well mode; tile filter has no effect. Use well= instead.",
                    UserWarning,
                    stacklevel=2,
                )
            else:
                tiles = [tile] if isinstance(tile, str) else list(tile)
                known_tiles = set(self._points_by_prefix.keys())
                unknown_tiles = set(tiles) - known_tiles
                if unknown_tiles:
                    warnings.warn(
                        f"Run {self.run_name!r}: tile(s) {sorted(unknown_tiles)!r} "
                        f"not recognised. Available: {sorted(known_tiles)}",
                        UserWarning,
                        stacklevel=2,
                    )
                tile_points: set[str] = set()
                for t in tiles:
                    tile_points.update(self.points_for_tile(t))
                candidates &= tile_points

        return sorted(candidates)

    # ------------------------------------------------------------------
    # Abstract IO methods
    # ------------------------------------------------------------------

    @abc.abstractmethod
    def _list_tables(self) -> list[str]: ...

    @abc.abstractmethod
    def _read_table(self, name: str, *, lazy: bool = True) -> ad.AnnData: ...

    @abc.abstractmethod
    def _has_table(self, name: str) -> bool: ...

    @abc.abstractmethod
    def _list_points(self) -> list[str]: ...

    @abc.abstractmethod
    def _read_points(self, name: str, *, lazy: bool = True) -> pl.LazyFrame: ...

    @abc.abstractmethod
    def _read_table_metadata_attrs(self, name: str) -> dict:
        """Read only ``uns['metadata_attrs']`` from a table without loading X."""
        ...

    @abc.abstractmethod
    def _read_json(self, path: str) -> dict | None:
        """Read a JSON file from the store.  *path* is relative to the zarr root."""
        ...

    @abc.abstractmethod
    def _read_sdata_attrs(self) -> dict:
        """Read the SpatialData zarr-root ``attrs`` as a plain dict.

        Returns ``{}`` when the attrs are missing or unreadable.
        """
        ...

    @property
    @abc.abstractmethod
    def _is_local(self) -> bool: ...

    @property
    def _is_zip(self) -> bool:
        """True when the backing store is a ``.zarr.zip`` archive.

        Defaults to ``False``; zip subclasses override to ``True``.
        Callers can use this to skip per-polony aggregations that
        require seekable random access on the points parquets (zip
        archives only expose streaming reads, so the points layer would
        have to extract to disk and read back, which is the source of
        the multi-well multi-batch OOM behavior).
        """
        return False


# ---------------------------------------------------------------------------
# LocalDirBackend
# ---------------------------------------------------------------------------


class _LocalDirBackend(RunInfo):
    """Backend for a ``.zarr`` directory on the local filesystem."""

    def __init__(self, root: Path) -> None:
        self._root = root
        super().__init__()

    @property
    def _source(self) -> str:
        return str(self._root)

    def _derive_run_name(self) -> str:
        try:
            data = self._read_json(_RUN_PARAMS_PATH)
            if data and "RunName" in data:
                return data["RunName"]
        except Exception:
            pass
        return super()._derive_run_name()

    def _list_tables(self) -> list[str]:
        tables_dir = self._root / "tables"
        if not tables_dir.is_dir():
            return []
        return sorted(p.name for p in tables_dir.iterdir() if p.is_dir())

    def _read_table(self, name: str, *, lazy: bool = True) -> ad.AnnData:
        table_path = str(self._root / "tables" / name)
        if lazy and _read_lazy is not None:
            try:
                return _materialize_obs_var(_read_lazy(table_path))
            except Exception:
                pass
        return ad.read_zarr(table_path)

    def _has_table(self, name: str) -> bool:
        return (self._root / "tables" / name).is_dir()

    def _read_table_metadata_attrs(self, name: str) -> dict:
        table_path = self._root / "tables" / name
        try:
            root = zarr.open_group(table_path, mode="r")
            if "uns" not in root or "metadata_attrs" not in root["uns"]:
                return {}
            ma_group = root["uns"]["metadata_attrs"]  # ty: ignore[invalid-argument-type]
            result: dict = {}
            for key in list(ma_group):  # ty: ignore[invalid-argument-type]
                val = ma_group[key][()]  # ty: ignore[invalid-argument-type,not-subscriptable]
                if hasattr(val, "item") and val.ndim == 0:
                    val = val.item()  # ty: ignore[call-non-callable]
                result[key] = val
            return result
        except Exception:
            return {}

    def _list_points(self) -> list[str]:
        points_dir = self._root / "points"
        if not points_dir.is_dir():
            return []
        return sorted(p.name for p in points_dir.iterdir() if p.is_dir() and (p / "points.parquet").exists())

    def _read_points(self, name: str, *, lazy: bool = True) -> pl.LazyFrame:
        parquet_path = self._root / "points" / name / "points.parquet"
        source = str(parquet_path / "**/*.parquet") if parquet_path.is_dir() else str(parquet_path)
        if lazy:
            return pl.scan_parquet(source)
        return pl.read_parquet(source).lazy()

    def _read_json(self, path: str) -> dict | None:
        fp = self._root / path
        if not fp.is_file():
            return None
        return json.loads(fp.read_text())

    def _read_sdata_attrs(self) -> dict:
        try:
            return dict(zarr.open_group(self._root, mode="r").attrs)
        except Exception:
            return {}

    @property
    def _is_local(self) -> bool:
        return True


# ---------------------------------------------------------------------------
# LocalZipBackend
# ---------------------------------------------------------------------------


class _LocalZipBackend(RunInfo):
    """Backend for a ``.zarr.zip`` file on the local filesystem."""

    _points_extract_cache: ClassVar[dict[tuple[str, str, str], Path]] = {}

    def __init__(self, zip_path: Path) -> None:
        self._zip_path = zip_path
        self._zarr_prefix = self._detect_zarr_prefix()
        self._tables_cache: list[str] | None = None
        super().__init__()

    @property
    def _source(self) -> str:
        return str(self._zip_path)

    def _detect_zarr_prefix(self) -> str:
        """Return the required top-level ``*.zarr`` prefix within the zip."""
        with zipfile.ZipFile(self._zip_path, "r") as zf:
            entries = [name for name in zf.namelist() if name and not name.endswith("/")]

        root_candidates = {name.split("/", 1)[0] for name in entries if "/" in name}
        zarr_roots = sorted(root for root in root_candidates if root.endswith(".zarr"))
        if len(zarr_roots) != 1:
            raise ValueError(
                f"Expected exactly one top-level '*.zarr' directory in {self._zip_path!s}, "
                f"found {zarr_roots or 'none'}."
            )

        prefix = zarr_roots[0]
        expected_prefix = f"{prefix}/"
        unexpected = [name for name in entries if not name.startswith(expected_prefix)]
        if unexpected:
            raise ValueError(
                f"Expected all archive entries to be under {expected_prefix!r}, "
                f"found entries like {unexpected[:3]}."
            )
        return prefix

    def _list_tables(self) -> list[str]:
        if self._tables_cache is not None:
            return list(self._tables_cache)
        zs = ZipStore(self._zip_path, mode="r")
        try:
            root = zarr.open_group(zs, mode="r", path=self._zarr_prefix)
            if "tables" not in root:
                self._tables_cache = []
                return []
            self._tables_cache = sorted(root["tables"].group_keys())
            return list(self._tables_cache)
        finally:
            zs.close()

    def _read_table(self, name: str, *, lazy: bool = True) -> ad.AnnData:
        zs = ZipStore(self._zip_path, mode="r")
        root = zarr.open_group(zs, mode="r", path=self._zarr_prefix)
        table_group = root["tables"][name]
        if lazy and _read_lazy is not None:
            try:
                # The lazy AnnData reads X from the ZipStore on demand, so the
                # store must stay open; it is kept alive by references from the
                # returned object (do not close it here). obs/var are
                # materialized to pandas so the result matches the other
                # backends' lazy contract.
                return _materialize_obs_var(_read_lazy(table_group))
            except Exception:
                pass
        try:
            return ad.read_zarr(table_group)
        finally:
            zs.close()

    def _has_table(self, name: str) -> bool:
        return name in self._list_tables()

    def _read_table_metadata_attrs(self, name: str) -> dict:
        zs = ZipStore(self._zip_path, mode="r")
        try:
            root = zarr.open_group(zs, mode="r", path=self._zarr_prefix)
            table_group = root["tables"][name]  # ty: ignore[invalid-argument-type]
            if "uns" not in table_group or "metadata_attrs" not in table_group["uns"]:  # ty: ignore[unsupported-operator,invalid-argument-type,not-subscriptable]
                return {}
            ma_group = table_group["uns"]["metadata_attrs"]  # ty: ignore[invalid-argument-type,not-subscriptable]
            result: dict = {}
            for key in list(ma_group):  # ty: ignore[invalid-argument-type]
                val = ma_group[key][()]  # ty: ignore[invalid-argument-type,not-subscriptable]
                if hasattr(val, "item") and val.ndim == 0:
                    val = val.item()  # ty: ignore[call-non-callable]
                result[key] = val
            return result
        except Exception:
            return {}
        finally:
            zs.close()

    def _list_points(self) -> list[str]:
        names: set[str] = set()
        points_prefix = f"{self._zarr_prefix}/points/"
        with zipfile.ZipFile(self._zip_path, "r") as zf:
            for entry in zf.namelist():
                if not entry.startswith(points_prefix):
                    continue
                relative = entry.removeprefix(points_prefix)
                parts = relative.split("/")
                if len(parts) >= 3 and parts[1] == "points.parquet" and parts[-1].endswith(".parquet"):
                    names.add(parts[0])
                    continue
                if len(parts) == 2 and parts[1] == "points.parquet":
                    names.add(parts[0])
        return sorted(names)

    def _points_cache_discriminator(self) -> str:
        """Stable cache key piece for this zip: ``sha1(path)`` + mtime + size.

        Re-extracts when the file is replaced in place; the path hash
        keeps zips with the same basename in different dirs separate.
        """
        st = self._zip_path.stat()
        path_hash = hashlib.sha1(str(self._zip_path).encode()).hexdigest()[:12]
        return f"{path_hash}_{int(st.st_mtime)}_{st.st_size}"

    def _persistent_points_root(self) -> Path:
        """Persistent cache directory for points extracted from this zip."""
        _PERSISTENT_CACHE_ROOT.mkdir(parents=True, exist_ok=True)
        return _PERSISTENT_CACHE_ROOT / "local_zip" / self._points_cache_discriminator()

    def _ensure_local_points(self, name: str) -> Path:
        """Extract a points element from the local zip into the persistent cache.

        Returns a directory containing either ``points.parquet`` (single
        file) or ``points.parquet/<part>.parquet`` (partitioned), mirroring
        the layout inside the zip so ``pl.scan_parquet`` can pushdown.
        """
        cache_key = (str(self._zip_path), self._points_cache_discriminator(), name)
        cached = self._points_extract_cache.get(cache_key)
        if cached is not None and cached.exists():
            return cached

        local_dir = self._persistent_points_root() / name
        sentinel = local_dir.parent / f".{name}.extracted"
        if sentinel.is_file():
            self._points_extract_cache[cache_key] = local_dir
            return local_dir

        point_root = f"{self._zarr_prefix}/points/{name}/points.parquet"
        point_prefix = f"{point_root}/"
        local_dir.mkdir(parents=True, exist_ok=True)
        total_size = 0
        with zipfile.ZipFile(self._zip_path, "r") as zf:
            all_names = zf.namelist()
            if point_root in all_names:
                dest = local_dir / "points.parquet"
                with zf.open(point_root, "r") as src, dest.open("wb") as out:
                    shutil.copyfileobj(src, out, length=1 << 22)
                total_size = dest.stat().st_size
            else:
                part_files = sorted(
                    n for n in all_names if n.startswith(point_prefix) and n.endswith(".parquet")
                )
                if not part_files:
                    raise KeyError(f"There is no item named '{point_root}' in the archive")
                for entry in part_files:
                    rel = entry.removeprefix(f"{self._zarr_prefix}/points/{name}/")
                    dest = local_dir / rel
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    with zf.open(entry, "r") as src, dest.open("wb") as out:
                        shutil.copyfileobj(src, out, length=1 << 22)
                    total_size += dest.stat().st_size
        sentinel.touch()
        _logger.info(
            "Extracted points %r (%.1f MB) from %s to %s",
            name,
            total_size / 1e6,
            self._zip_path,
            local_dir,
        )
        self._points_extract_cache[cache_key] = local_dir
        return local_dir

    def _read_points(self, name: str, *, lazy: bool = True) -> pl.LazyFrame:
        local_dir = self._ensure_local_points(name)
        parquet_root = local_dir / "points.parquet"
        if parquet_root.is_file():
            sources: Path | list[Path] = parquet_root
        elif parquet_root.is_dir():
            files = sorted(parquet_root.rglob("*.parquet"))
            if not files:
                raise KeyError(f"No parquet files found in {parquet_root}")
            sources = files
        else:
            raise KeyError(f"Points element {name!r} not found at {parquet_root}")
        if lazy:
            return pl.scan_parquet(sources)
        return pl.read_parquet(sources).lazy()

    def _read_json(self, path: str) -> dict | None:
        normalized = path.strip("/")
        archive_path = f"{self._zarr_prefix}/{normalized}"
        with zipfile.ZipFile(self._zip_path, "r") as zf:
            if archive_path not in {name.strip("/") for name in zf.namelist()}:
                return None
            return json.loads(zf.read(archive_path))

    def _read_sdata_attrs(self) -> dict:
        zs = ZipStore(self._zip_path, mode="r")
        try:
            root = zarr.open_group(zs, mode="r", path=self._zarr_prefix)
            return dict(root.attrs)
        except Exception:
            return {}
        finally:
            zs.close()

    @property
    def _is_local(self) -> bool:
        return True

    @property
    def _is_zip(self) -> bool:
        return True


# ---------------------------------------------------------------------------
# S3DirBackend
# ---------------------------------------------------------------------------


class _S3DirBackend(RunInfo):
    """Backend for a ``.zarr`` directory on S3."""

    def __init__(
        self,
        s3_url: str,
        *,
        storage_options: _StorageOptions | None = None,
    ) -> None:
        self._url = s3_url.rstrip("/")
        self._storage_options = storage_options or {}
        self._fs = s3fs.S3FileSystem(**self._storage_options)
        self._prefix = self._url.removeprefix("s3://")
        super().__init__()

    @property
    def _source(self) -> str:
        return self._url

    def _list_tables(self) -> list[str]:
        tables_prefix = f"{self._prefix}/tables"
        if not self._fs.isdir(tables_prefix):
            return []
        entries = self._fs.ls(tables_prefix, detail=False)
        return sorted(e.rsplit("/", 1)[-1] for e in entries if self._fs.isdir(e))

    def _read_table(self, name: str, *, lazy: bool = True) -> ad.AnnData:
        table_url = f"{self._url}/tables/{name}"
        mapper = fsspec.get_mapper(table_url, **self._storage_options)
        if lazy and _read_lazy is not None:
            try:
                return _materialize_obs_var(_read_lazy(mapper))
            except Exception:
                pass
        return ad.read_zarr(mapper)

    def _has_table(self, name: str) -> bool:
        return self._fs.isdir(f"{self._prefix}/tables/{name}")

    def _read_table_metadata_attrs(self, name: str) -> dict:
        table_url = f"{self._url}/tables/{name}"
        try:
            mapper = fsspec.get_mapper(table_url, **self._storage_options)
            root = zarr.open_group(mapper, mode="r")
            if "uns" not in root or "metadata_attrs" not in root["uns"]:
                return {}
            ma_group = root["uns"]["metadata_attrs"]  # ty: ignore[invalid-argument-type]
            result: dict = {}
            for key in list(ma_group):  # ty: ignore[invalid-argument-type]
                val = ma_group[key][()]  # ty: ignore[invalid-argument-type,not-subscriptable]
                if hasattr(val, "item") and val.ndim == 0:
                    val = val.item()  # ty: ignore[call-non-callable]
                result[key] = val
            return result
        except Exception:
            return {}

    def _list_points(self) -> list[str]:
        points_prefix = f"{self._prefix}/points"
        if not self._fs.isdir(points_prefix):
            return []
        entries = self._fs.ls(points_prefix, detail=False)
        names: list[str] = []
        for e in entries:
            if self._fs.isdir(e) and self._fs.exists(f"{e}/points.parquet"):
                names.append(e.rsplit("/", 1)[-1])
        return sorted(names)

    @staticmethod
    def _polars_storage_options(opts: _StorageOptions) -> dict[str, Any]:
        """Translate s3fs-style storage options to Polars object-store options."""
        mapping = {"profile": "aws_profile", "region": "aws_region"}
        return {mapping.get(k, k): v for k, v in opts.items()}

    def _read_points(self, name: str, *, lazy: bool = True) -> pl.LazyFrame:
        parquet_url = f"{self._url}/points/{name}/points.parquet"
        parquet_key = f"{self._prefix}/points/{name}/points.parquet"
        source = f"{parquet_url}/**/*.parquet" if self._fs.isdir(parquet_key) else parquet_url
        pl_storage = (
            self._polars_storage_options(self._storage_options) if self._storage_options else {"anon": False}
        )
        if lazy:
            return pl.scan_parquet(source, storage_options=pl_storage)
        return pl.read_parquet(source, storage_options=pl_storage).lazy()

    def _read_json(self, path: str) -> dict | None:
        key = f"{self._prefix}/{path}"
        if not self._fs.exists(key):
            return None
        with self._fs.open(key, "r") as f:
            return json.load(f)

    def _read_sdata_attrs(self) -> dict:
        try:
            mapper = fsspec.get_mapper(self._url, **self._storage_options)
            return dict(zarr.open_group(mapper, mode="r").attrs)
        except Exception:
            return {}

    @property
    def _is_local(self) -> bool:
        return False


# ---------------------------------------------------------------------------
# S3ZipBackend  (lazy per-table extraction into persistent ETag cache)
# ---------------------------------------------------------------------------


class _S3ZipBackend(RunInfo):
    """Backend for a ``.zarr.zip`` file on S3.

    Discovery reads only the zip central directory. Tables are extracted
    lazily, one table at a time, into an ETag-keyed on-disk cache so repeat
    sessions can reuse unchanged S3 objects.
    """

    _extract_cache: ClassVar[dict[tuple[str, str, str], Path]] = {}
    _metadata_attrs_cache: ClassVar[dict[tuple[str, str, str], Path]] = {}
    _points_extract_cache: ClassVar[dict[tuple[str, str, str], Path]] = {}

    @classmethod
    def _cleanup_extract_cache(cls) -> None:
        """Clear the in-memory cache; persistent table extracts stay on disk."""
        cls._extract_cache.clear()
        cls._metadata_attrs_cache.clear()
        cls._points_extract_cache.clear()

    def __init__(
        self,
        s3_url: str,
        *,
        storage_options: _StorageOptions | None = None,
    ) -> None:
        self._s3_url = s3_url
        self._storage_options = storage_options or {}
        self._remote_fs: fsspec.AbstractFileSystem | None = None
        self._zarr_prefix: str | None = None
        self._etag_cached: str | None = None
        self._tables_cache: list[str] | None = None
        super().__init__()

    @property
    def _source(self) -> str:
        return self._s3_url

    def _remote(self) -> fsspec.AbstractFileSystem:
        """Lazily open the remote zip via fsspec (reads central directory only)."""
        if self._remote_fs is None:
            s3_opts = dict(self._storage_options)
            # ZipFileSystem uses target_options for remote ``fo``, not ``s3=`` (ignored).
            self._remote_fs = fsspec.filesystem(
                "zip",
                fo=self._s3_url,
                target_protocol="s3",
                target_options=s3_opts,
            )
            entries = self._remote_fs.ls("/")
            zarr_roots = [(e["name"] if isinstance(e, dict) else e).rstrip("/") for e in entries]
            zarr_roots = [r for r in zarr_roots if r.endswith(".zarr")]
            if len(zarr_roots) != 1:
                raise ValueError(
                    f"Expected exactly one top-level '*.zarr' directory in {self._s3_url!s}, "
                    f"found {zarr_roots or 'none'}."
                )
            self._zarr_prefix = zarr_roots[0]
        return self._remote_fs

    def _prefix(self) -> str:
        """Zarr prefix inside the zip archive."""
        self._remote()
        assert self._zarr_prefix is not None
        return self._zarr_prefix

    def _etag(self) -> str:
        """Stable cache discriminator for the current S3 object."""
        if self._etag_cached is None:
            s3 = s3fs.S3FileSystem(**self._storage_options)
            info = s3.info(self._s3_url)
            etag = (info.get("ETag") or info.get("etag") or "").strip('"')
            self._etag_cached = etag or hashlib.sha1(self._s3_url.encode()).hexdigest()
        return self._etag_cached

    def _ensure_cache_root(self) -> Path:
        """Create and validate the persistent cache root directory."""
        root = _PERSISTENT_CACHE_ROOT
        try:
            root.mkdir(parents=True, exist_ok=True)
            probe = root / ".write_probe"
            probe.write_text("")
            probe.unlink(missing_ok=True)
        except OSError as exc:
            raise OSError(f"SpatialData cache root is not writable: {root}") from exc
        return root

    def _persistent_root(self) -> Path:
        self._ensure_cache_root()
        url_hash = hashlib.sha1(self._s3_url.encode()).hexdigest()[:12]
        return _PERSISTENT_CACHE_ROOT / url_hash / self._etag()

    def _url_cache_parent(self) -> Path:
        self._ensure_cache_root()
        url_hash = hashlib.sha1(self._s3_url.encode()).hexdigest()[:12]
        return _PERSISTENT_CACHE_ROOT / url_hash

    def _prune_stale_etag_caches(self) -> None:
        """Remove old sibling ETag cache directories for the same S3 URL."""
        parent = self._url_cache_parent()
        current = self._etag()
        if not parent.is_dir():
            return
        now = time.time()
        for child in parent.iterdir():
            if not child.is_dir() or child.name == current:
                continue
            if self._etag_cache_has_active_locks(child):
                continue
            if (now - child.stat().st_mtime) < _S3ZIP_PRUNE_MIN_AGE_SEC:
                continue
            shutil.rmtree(child, ignore_errors=True)

    @staticmethod
    def _etag_cache_has_active_locks(etag_dir: Path) -> bool:
        """True when another process may be extracting into this ETag directory."""
        if not etag_dir.is_dir():
            return False
        return any(path.is_dir() and path.name.endswith(".lock") for path in etag_dir.rglob("*"))

    @staticmethod
    def _touch_lock_heartbeat(lock_dir: Path) -> None:
        """Refresh lock freshness while a long extraction is in progress."""
        owner_file = lock_dir / "owner.json"
        if owner_file.is_file():
            owner_file.touch()

    @staticmethod
    def _lock_is_stale(lock_dir: Path) -> bool:
        if not lock_dir.is_dir():
            return False
        owner_file = lock_dir / "owner.json"
        ref = owner_file if owner_file.is_file() else lock_dir
        return (time.time() - ref.stat().st_mtime) > _S3ZIP_LOCK_STALE_SEC

    def _maybe_heartbeat_lock(self, lock_dir: Path, last_heartbeat: float) -> float:
        """Touch the lock if the heartbeat interval has elapsed; return updated timestamp."""
        now = time.monotonic()
        if now - last_heartbeat >= _S3ZIP_LOCK_HEARTBEAT_SEC:
            self._touch_lock_heartbeat(lock_dir)
            return now
        return last_heartbeat

    def _acquire_table_lock(self, lock_dir: Path, sentinel: Path) -> None:
        """Acquire a per-table extraction lock with stale-lock recovery."""
        deadline = time.monotonic() + _S3ZIP_LOCK_MAX_WAIT_SEC
        while True:
            if sentinel.is_file():
                return
            try:
                lock_dir.mkdir(parents=True)
                owner = {
                    "pid": os.getpid(),
                    "hostname": os.uname().nodename,
                    "started": time.time(),
                }
                (lock_dir / "owner.json").write_text(json.dumps(owner))
                return
            except FileExistsError:
                if sentinel.is_file():
                    return
                if self._lock_is_stale(lock_dir):
                    shutil.rmtree(lock_dir, ignore_errors=True)
                    continue
                if time.monotonic() >= deadline:
                    raise TimeoutError(
                        f"Timed out waiting for extraction lock {lock_dir!s} "
                        f"after {_S3ZIP_LOCK_MAX_WAIT_SEC:.0f}s"
                    ) from None
                time.sleep(_S3ZIP_LOCK_POLL_SEC)

    def _ensure_local_table(self, name: str) -> Path:
        """Extract one table from the remote zip into the persistent cache."""
        etag = self._etag()
        cache_key = (self._s3_url, etag, name)
        if cache_key in self._extract_cache:
            return self._extract_cache[cache_key]

        prefix = self._prefix()
        local_dir = self._persistent_root() / prefix / "tables" / name
        sentinel = local_dir.parent / f".{name}.extracted"
        if sentinel.is_file():
            self._extract_cache[cache_key] = local_dir
            return local_dir

        lock_dir = local_dir.parent / f".{name}.lock"
        self._acquire_table_lock(lock_dir, sentinel)
        if sentinel.is_file():
            self._extract_cache[cache_key] = local_dir
            return local_dir

        try:
            if sentinel.is_file():
                self._extract_cache[cache_key] = local_dir
                return local_dir

            fs = self._remote()
            remote_root = f"{prefix}/tables/{name}"
            files = fs.find(remote_root)
            if not files:
                raise KeyError(f"No files found under {remote_root}")

            local_dir.mkdir(parents=True, exist_ok=True)

            total_size = 0
            last_heartbeat = time.monotonic()
            for f_path in files:
                last_heartbeat = self._maybe_heartbeat_lock(lock_dir, last_heartbeat)
                rel = f_path.removeprefix(f"{remote_root}/")
                dest = local_dir / rel
                dest.parent.mkdir(parents=True, exist_ok=True)
                with fs.open(f_path, "rb") as src:
                    data = src.read()
                dest.write_bytes(data)
                total_size += len(data)

            sentinel.touch()
            self._prune_stale_etag_caches()
            _logger.info(
                "Extracted table %r (%d files, %.1f MB) from %s to %s",
                name,
                len(files),
                total_size / 1e6,
                self._s3_url,
                local_dir,
            )
        finally:
            shutil.rmtree(lock_dir, ignore_errors=True)

        self._extract_cache[cache_key] = local_dir
        return local_dir

    def _ensure_local_metadata_attrs(self, name: str) -> Path:
        """Extract only ``uns/metadata_attrs`` for a table (small remote read)."""
        etag = self._etag()
        cache_key = (self._s3_url, etag, name)
        if cache_key in self._metadata_attrs_cache:
            return self._metadata_attrs_cache[cache_key]

        prefix = self._prefix()
        local_dir = self._persistent_root() / prefix / "tables" / name / "uns" / "metadata_attrs"
        sentinel = local_dir.parent.parent / f".{name}.metadata_attrs.extracted"
        if sentinel.is_file():
            self._metadata_attrs_cache[cache_key] = local_dir
            return local_dir

        lock_dir = local_dir.parent.parent / f".{name}.metadata_attrs.lock"
        self._acquire_table_lock(lock_dir, sentinel)
        if sentinel.is_file():
            self._metadata_attrs_cache[cache_key] = local_dir
            return local_dir

        try:
            if sentinel.is_file():
                self._metadata_attrs_cache[cache_key] = local_dir
                return local_dir

            fs = self._remote()
            remote_root = f"{prefix}/tables/{name}/uns/metadata_attrs"
            try:
                files = fs.find(remote_root)
            except FileNotFoundError:
                files = []

            if not files:
                sentinel.parent.mkdir(parents=True, exist_ok=True)
                sentinel.touch()
                self._metadata_attrs_cache[cache_key] = local_dir
                return local_dir

            local_dir.mkdir(parents=True, exist_ok=True)
            last_heartbeat = time.monotonic()
            for f_path in files:
                last_heartbeat = self._maybe_heartbeat_lock(lock_dir, last_heartbeat)
                rel = f_path.removeprefix(f"{remote_root}/")
                if not rel:
                    continue
                dest = local_dir / rel
                dest.parent.mkdir(parents=True, exist_ok=True)
                with fs.open(f_path, "rb") as src:
                    dest.write_bytes(src.read())

            sentinel.touch()
        finally:
            shutil.rmtree(lock_dir, ignore_errors=True)

        self._metadata_attrs_cache[cache_key] = local_dir
        return local_dir

    def _list_tables(self) -> list[str]:
        if self._tables_cache is not None:
            return list(self._tables_cache)
        fs = self._remote()
        prefix = self._prefix()
        try:
            entries = fs.ls(f"{prefix}/tables", detail=True)
        except FileNotFoundError:
            self._tables_cache = []
            return []
        names: list[str] = []
        for entry in entries:
            name = entry["name"] if isinstance(entry, dict) else entry
            etype = entry.get("type") if isinstance(entry, dict) else None
            if etype in ("directory", None):
                names.append(name.rstrip("/").rsplit("/", 1)[-1])
        self._tables_cache = sorted(names)
        return list(self._tables_cache)

    def _read_table(self, name: str, *, lazy: bool = True) -> ad.AnnData:
        table_path = str(self._ensure_local_table(name))
        if lazy and _read_lazy is not None:
            try:
                return _materialize_obs_var(_read_lazy(table_path))
            except Exception:
                pass
        return ad.read_zarr(table_path)

    def _has_table(self, name: str) -> bool:
        return name in self._list_tables()

    def _read_table_metadata_attrs(self, name: str) -> dict:
        """Read ``uns/metadata_attrs`` without extracting the full table body."""
        try:
            meta_dir = self._ensure_local_metadata_attrs(name)
            if not meta_dir.is_dir():
                return {}
            ma_group = zarr.open_group(meta_dir, mode="r")
            return _read_zarr_metadata_attrs_group(ma_group)
        except Exception:
            return {}

    def _list_points(self) -> list[str]:
        fs = self._remote()
        prefix = self._prefix()
        points_path = f"{prefix}/points"
        try:
            entries = fs.ls(points_path)
        except FileNotFoundError:
            return []
        names: list[str] = []
        for entry in entries:
            name = entry["name"] if isinstance(entry, dict) else entry
            info = fs.info(name)
            if info.get("type") == "directory":
                point_name = name.rstrip("/").rsplit("/", 1)[-1]
                names.append(point_name)
        return sorted(names)

    def _ensure_local_points(self, name: str) -> Path:
        """Extract a points element from the remote zip into the persistent cache.

        Mirrors :meth:`_ensure_local_table` but lands the parquet bytes
        under ``points/<name>/`` so ``pl.scan_parquet`` can pushdown
        column/row-group reads instead of materialising the whole file.
        """
        etag = self._etag()
        cache_key = (self._s3_url, etag, name)
        if cache_key in self._points_extract_cache:
            return self._points_extract_cache[cache_key]

        prefix = self._prefix()
        local_dir = self._persistent_root() / prefix / "points" / name
        sentinel = local_dir.parent / f".{name}.extracted"
        if sentinel.is_file():
            self._points_extract_cache[cache_key] = local_dir
            return local_dir

        lock_dir = local_dir.parent / f".{name}.lock"
        self._acquire_table_lock(lock_dir, sentinel)
        if sentinel.is_file():
            self._points_extract_cache[cache_key] = local_dir
            return local_dir

        try:
            if sentinel.is_file():
                self._points_extract_cache[cache_key] = local_dir
                return local_dir

            fs = self._remote()
            point_root = f"{prefix}/points/{name}/points.parquet"
            try:
                info = fs.info(point_root)
            except (FileNotFoundError, KeyError):
                info = {}

            local_dir.mkdir(parents=True, exist_ok=True)
            total_size = 0
            last_heartbeat = time.monotonic()

            if info.get("type") == "file":
                dest = local_dir / "points.parquet"
                with fs.open(point_root, "rb") as src, dest.open("wb") as out:
                    shutil.copyfileobj(src, out, length=1 << 22)
                total_size = dest.stat().st_size
            else:
                part_files = sorted(p for p in fs.find(point_root) if p.endswith(".parquet"))
                if not part_files:
                    raise KeyError(f"No parquet files found under {point_root}")
                for f_path in part_files:
                    last_heartbeat = self._maybe_heartbeat_lock(lock_dir, last_heartbeat)
                    rel = f_path.removeprefix(f"{prefix}/points/{name}/")
                    dest = local_dir / rel
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    with fs.open(f_path, "rb") as src, dest.open("wb") as out:
                        shutil.copyfileobj(src, out, length=1 << 22)
                    total_size += dest.stat().st_size

            sentinel.touch()
            _logger.info(
                "Extracted points %r (%.1f MB) from %s to %s",
                name,
                total_size / 1e6,
                self._s3_url,
                local_dir,
            )
        finally:
            shutil.rmtree(lock_dir, ignore_errors=True)

        self._points_extract_cache[cache_key] = local_dir
        return local_dir

    def _read_points(self, name: str, *, lazy: bool = True) -> pl.LazyFrame:
        local_dir = self._ensure_local_points(name)
        parquet_root = local_dir / "points.parquet"
        if parquet_root.is_file():
            sources: Path | list[Path] = parquet_root
        elif parquet_root.is_dir():
            files = sorted(parquet_root.rglob("*.parquet"))
            if not files:
                raise KeyError(f"No parquet files found in {parquet_root}")
            sources = files
        else:
            raise KeyError(f"Points element {name!r} not found at {parquet_root}")
        if lazy:
            return pl.scan_parquet(sources)
        return pl.read_parquet(sources).lazy()

    def _read_json(self, path: str) -> dict | None:
        fs = self._remote()
        prefix = self._prefix()
        full_path = f"{prefix}/{path.strip('/')}"
        try:
            with fs.open(full_path, "rb") as f:
                return json.loads(f.read())
        except (FileNotFoundError, KeyError):
            return None

    def _read_sdata_attrs(self) -> dict:
        for attrs_path in ("zarr.json", ".zattrs"):
            data = self._read_json(attrs_path)
            if not data:
                continue
            attrs = data.get("attributes") if attrs_path == "zarr.json" else data
            if isinstance(attrs, dict):
                return dict(attrs)
        return {}

    @property
    def _is_local(self) -> bool:
        return False

    @property
    def _is_zip(self) -> bool:
        return True


atexit.register(_S3ZipBackend._cleanup_extract_cache)


# ---------------------------------------------------------------------------
# SDataBackend  (wraps an in-memory SpatialData object)
# ---------------------------------------------------------------------------


class _SDataBackend(RunInfo):
    """Backend wrapping an in-memory ``SpatialData`` object."""

    def __init__(self, sdata: sd.SpatialData) -> None:
        self._sdata = sdata
        super().__init__()

    @property
    def _source(self) -> str:
        sdata_path = getattr(self._sdata, "path", None)
        return str(sdata_path) if sdata_path else f"sdata@{id(self._sdata)}"

    def source_exists(self) -> bool:
        # In-memory object is always reachable; ``_source`` may be a synthetic id.
        return True

    def _derive_run_name(self) -> str:
        try:
            data = self._read_json(_RUN_PARAMS_PATH)
            if data and "RunName" in data:
                return data["RunName"]
        except Exception:
            pass
        return super()._derive_run_name() if not self._source.startswith("sdata@") else "sdata"

    def _list_tables(self) -> list[str]:
        return sorted(self._sdata.tables.keys())

    def _read_table(self, name: str, *, lazy: bool = True) -> ad.AnnData:
        # In-memory tables are already loaded, so `lazy` is a no-op here.
        # Materialize obs/var for parity with the other backends (a no-op unless
        # the in-memory SpatialData happens to hold lazy Dataset2D frames).
        return _materialize_obs_var(self._sdata.tables[name])

    def _has_table(self, name: str) -> bool:
        return name in self._sdata.tables

    def _read_table_metadata_attrs(self, name: str) -> dict:
        if name not in self._sdata.tables:
            return {}
        adata = self._sdata.tables[name]
        return dict(adata.uns.get("metadata_attrs", {}))

    def _list_points(self) -> list[str]:
        return sorted(self._sdata.points.keys())

    def _read_points(self, name: str, *, lazy: bool = True) -> pl.LazyFrame:
        points_table = self._sdata.points[name]
        points_pd = points_table.compute() if hasattr(points_table, "compute") else points_table
        if not isinstance(points_pd, pd.DataFrame):
            points_pd = pd.DataFrame(points_pd)
        arrow_table = pa.Table.from_pandas(points_pd, preserve_index=False)
        return pl.from_arrow(arrow_table).lazy()

    def _read_json(self, path: str) -> dict | None:
        sdata_path = getattr(self._sdata, "path", None)
        if sdata_path is None:
            return None
        fp = Path(str(sdata_path)) / path
        if not fp.is_file():
            return None
        return json.loads(fp.read_text())

    def _read_sdata_attrs(self) -> dict:
        try:
            return dict(self._sdata.attrs or {})
        except Exception:
            return {}

    @property
    def _is_local(self) -> bool:
        return True
