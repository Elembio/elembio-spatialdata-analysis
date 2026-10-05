"""Read well-keyed images and label masks from a SpatialData store."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd
from numpy.typing import NDArray


def _store_keys(store) -> list[str]:
    if store is None:
        return []
    keys = getattr(store, "keys", None)
    if callable(keys):
        return [str(k) for k in keys()]
    if isinstance(store, Mapping):
        return [str(k) for k in store]
    return []


def list_image_wells(sdata) -> list[str]:
    """Well ids that have a ``{well}-images-cell_paint`` element."""
    names = _store_keys(getattr(sdata, "images", None))
    suffix = "-images-cell_paint"
    wells = [n[: -len(suffix)] for n in names if n.endswith(suffix)]
    return wells


def has_image(sdata, well: str) -> bool:
    key = f"{well}-images-cell_paint"
    return key in _store_keys(getattr(sdata, "images", None))


def resolve_image_key(sdata, well: str) -> str:
    """Return the image key for ``well`` (``"{well}-images-cell_paint"``); ``KeyError`` if absent."""
    key = f"{well}-images-cell_paint"
    names = _store_keys(getattr(sdata, "images", None))
    if key in names:
        return key
    available = sorted(n for n in names if "-images-" in n)
    raise KeyError(f"No image for well {well!r}. Available: {available}")


def resolve_label_key(sdata, well: str, compartment: str = "Cell") -> str:
    """Return the label key for ``well`` (``"{well}-labels-{compartment}"``); ``KeyError`` if absent."""
    key = f"{well}-labels-{compartment}"
    names = _store_keys(getattr(sdata, "labels", None))
    if key in names:
        return key
    available = sorted(n for n in names if n.endswith(f"-labels-{compartment}"))
    raise KeyError(f"No {compartment} labels for well {well!r}. Available: {available}")


def resolve_shapes_key(sdata, well: str, compartment: str = "Cell") -> str:
    """Return the shapes key for ``well`` (``"{well}-shapes-{compartment}"``); ``KeyError`` if absent."""
    key = f"{well}-shapes-{compartment}"
    names = _store_keys(getattr(sdata, "shapes", None))
    if key in names:
        return key
    available = sorted(n for n in names if n.endswith(f"-shapes-{compartment}"))
    raise KeyError(f"No {compartment} shapes for well {well!r}. Available: {available}")


def find_cell_well(sdata, cell_id: int) -> str:
    """Well of ``cell_id``, read from a table on the store.

    Uses any ``sdata.tables`` entry whose ``obs`` carries both ``Cell`` and
    ``Well``. Cell ids are unique across the store, so the first hit wins.
    """
    tables = getattr(sdata, "tables", None)
    keys = _store_keys(tables)
    searched: list[str] = []
    for name in keys:
        obs = getattr(tables[name], "obs", None)
        if obs is None or not {"Cell", "Well"} <= set(obs.columns):
            continue
        searched.append(name)
        ids = pd.to_numeric(obs["Cell"], errors="coerce")
        hit = obs.loc[ids == int(cell_id), "Well"]
        if len(hit) == 0:
            continue
        wells = {str(w) for w in hit}
        if len(wells) > 1:
            raise ValueError(f"Cell {cell_id} maps to several wells in {name!r}: {sorted(wells)}.")
        return str(hit.iloc[0])
    if not searched:
        raise ValueError(
            f"No table on this store has obs['Cell'] and obs['Well'], so well "
            f"lookup for cell {cell_id} is not possible. Pass well= instead."
        )
    raise KeyError(f"Cell {cell_id} not found in {searched}.")


def cell_shape(sdata, well: str, cell_id: int, *, compartment: str = "Cell"):
    """Return the stored polygon for one cell and the shapes element name."""
    key = resolve_shapes_key(sdata, well, compartment=compartment)
    gdf = sdata.shapes[key]
    ids = pd.to_numeric(pd.Series(np.asarray(gdf.index)), errors="coerce")
    match = np.flatnonzero((ids == int(cell_id)).to_numpy())
    if match.size == 0:
        raise KeyError(f"Cell {cell_id} not in {key!r}.")
    if match.size > 1:
        raise ValueError(f"Cell {cell_id} appears {match.size} times in {key!r}.")
    return gdf.geometry.iloc[int(match[0])], key


def image_level_info(sdata, well: str, scale: int | str | None) -> tuple[str, tuple[int, int], float, float]:
    """``(key, (h, w), y_ratio, x_ratio)`` for one pyramid level of a well image.

    The ratios convert full-resolution intrinsic coordinates — the space the
    stored shapes live in — to this level. Metadata only; no pixels are read.
    """
    key = resolve_image_key(sdata, well)
    element = sdata.images[key]
    level = _select_scale(element, scale)
    full = _select_scale(element, 0)
    h, w = array_hw(level)
    h0, w0 = array_hw(full)
    return key, (h, w), h0 / h, w0 / w


def full_to_level(value, ratio: float):
    """Map full-resolution pixel coordinates onto a coarser level's grid."""
    return (np.asarray(value, dtype=np.float64) + 0.5) / float(ratio) - 0.5


def array_hw(array) -> tuple[int, int]:
    """Spatial ``(height, width)`` of a raster element."""
    shape = tuple(int(s) for s in array.shape)
    return shape[-2], shape[-1]


def window_array(array, y0: int, y1: int, x0: int, x1: int):
    """Slice ``(y, x)``, staying lazy when the backend supports ``isel``."""
    isel = getattr(array, "isel", None)
    if callable(isel):
        try:
            return isel(y=slice(y0, y1), x=slice(x0, x1))
        except (TypeError, ValueError, KeyError):
            pass
    data = _to_numpy(array)
    if data.ndim == 2:
        return data[y0:y1, x0:x1]
    return data[..., y0:y1, x0:x1]


def _to_numpy(obj: Any) -> NDArray:
    if hasattr(obj, "compute") and not isinstance(obj, np.ndarray):
        obj = obj.compute()
    if hasattr(obj, "to_numpy") and not isinstance(obj, np.ndarray):
        return np.asarray(obj.to_numpy())
    if hasattr(obj, "values") and not isinstance(obj, np.ndarray | pd.Series | pd.DataFrame):
        return np.asarray(obj.values)
    return np.asarray(obj)


def _unwrap_image_node(node: Any) -> Any:
    ds = getattr(node, "ds", None)
    if ds is not None and hasattr(ds, "__contains__") and "image" in ds:
        return ds["image"]
    image = getattr(node, "image", None)
    if image is not None and image is not node:
        return image
    return node


def _scale_keys(element: Any) -> list[str] | None:
    keys_fn = getattr(element, "keys", None)
    if not callable(keys_fn):
        return None
    try:
        raw = [str(k) for k in keys_fn()]
    except Exception:
        return None
    scale_like = [k for k in raw if k.isdigit() or k.startswith("scale")]
    if not scale_like:
        return None

    def _sort(k: str) -> int:
        digits = "".join(c for c in k if c.isdigit())
        return int(digits) if digits else 0

    return sorted(scale_like, key=_sort)


def _select_scale(element: Any, scale: int | str | None) -> Any:
    keys = _scale_keys(element)
    if keys is None:
        return _unwrap_image_node(element)
    if scale is None:
        key = keys[-1]
    elif isinstance(scale, int):
        key = keys[scale] if 0 <= scale < len(keys) else keys[-1]
    else:
        key = str(scale)
        if key not in keys:
            raise KeyError(f"Scale {scale!r} not in {keys}.")
    return _unwrap_image_node(element[key])


def _channel_names(array) -> list[str]:
    coords = getattr(array, "coords", None)
    if coords is not None and "c" in coords:
        raw = coords["c"]
        values = getattr(raw, "values", raw)
        return [str(v) for v in np.asarray(values)]
    return []


def _extract_channel(array, name: str | None, available: Sequence[str]) -> NDArray[np.float32]:
    data = _to_numpy(array)
    if data.ndim == 2:
        return data.astype(np.float32, copy=False)
    if name is None:
        return np.zeros(data.shape[-2:], dtype=np.float32)
    if available and name in available:
        if hasattr(array, "sel"):
            try:
                return _to_numpy(array.sel(c=name)).astype(np.float32, copy=False).squeeze()
            except Exception:
                pass
        return data[list(available).index(name)].astype(np.float32, copy=False)
    try:
        return data[int(name)].astype(np.float32, copy=False)
    except (TypeError, ValueError, IndexError):
        return np.zeros(data.shape[-2:], dtype=np.float32)


def _percentile_norm(plane: NDArray, lo: float, hi: float) -> NDArray[np.float32]:
    a = np.asarray(plane, dtype=np.float32)
    lo_v, hi_v = np.percentile(a, [lo, hi])
    if hi_v <= lo_v:
        hi_v = lo_v + 1e-6
    return np.clip((a - lo_v) / (hi_v - lo_v), 0.0, 1.0)


def _name_list(value: str | Sequence[str], *, what: str) -> list[str]:
    if isinstance(value, str):
        return [value]
    names = [str(v) for v in value]
    if not names:
        raise ValueError(f"{what} is empty.")
    return names


def list_image_channels(sdata, well: str, *, scale: int | str | None = None) -> list[str]:
    """Return channel names on ``{well}-images-cell_paint``.

    Reads coordinates only — it does not load pixels. Names are the same at
    every pyramid level. Includes morphology stains and IF targets present
    on that image; it does not invent channels.
    """
    key = resolve_image_key(sdata, well)
    array = _select_scale(sdata.images[key], scale)
    return _channel_names(array)


def label_image(
    sdata,
    well: str,
    *,
    compartment: str = "Cell",
    scale: int | str | None = None,
) -> tuple[NDArray[np.int64], str]:
    """Return an ``(H, W)`` integer label mask and the label element name."""
    key = resolve_label_key(sdata, well, compartment=compartment)
    array = _select_scale(sdata.labels[key], scale)
    lab = _to_numpy(array).squeeze().astype(np.int64, copy=False)
    if lab.ndim != 2:
        raise ValueError(f"Label array for {key!r} has shape {lab.shape}; expected (y, x).")
    return lab, key


def _match_hw(rgb: NDArray, lab: NDArray) -> tuple[NDArray, NDArray]:
    h = min(rgb.shape[0], lab.shape[0])
    w = min(rgb.shape[1], lab.shape[1])
    return rgb[:h, :w], lab[:h, :w]
