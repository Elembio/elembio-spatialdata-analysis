"""Channel selection and tinted RGB composites."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
from matplotlib.colors import to_rgb
from numpy.typing import NDArray

from elembio_spatialdata_analysis.vis.spatial._access import (
    _channel_names,
    _extract_channel,
    _name_list,
    _percentile_norm,
    _select_scale,
    list_image_channels,
    resolve_image_key,
    window_array,
)

_DEFAULT_RGB = ("Nucleus", "Actin", "Cell-Membrane")
_DEFAULT_TINTS = {"Nucleus": "red", "Actin": "green", "Cell-Membrane": "blue"}
_FALLBACK_TINTS = ("red", "green", "blue", "magenta", "cyan", "yellow", "orange", "white")


def _as_rgb(color) -> NDArray[np.float32]:
    return np.asarray(to_rgb(color), dtype=np.float32)


def _resolve_channel_names(sdata, well: str, channels: str | Sequence[str] | None) -> list[str]:
    available = list_image_channels(sdata, well)
    if channels is None:
        names = [c for c in _DEFAULT_RGB if not available or c in available]
        if not names:
            names = list(available[:3])
        if not names:
            raise KeyError(f"No channels on well {well!r}.")
        return names
    names = _name_list(channels, what="channels")
    missing = [n for n in names if available and n not in available]
    if missing:
        raise KeyError(f"Channel(s) {missing} not in {available}.")
    return names


def _resolve_tints(
    names: Sequence[str],
    channel_colors: str | Sequence[str] | Mapping[str, str] | None,
) -> list[NDArray[np.float32]]:
    if channel_colors is None:
        return [
            _as_rgb(_DEFAULT_TINTS.get(name, _FALLBACK_TINTS[i % len(_FALLBACK_TINTS)]))
            for i, name in enumerate(names)
        ]
    if isinstance(channel_colors, Mapping):
        missing = [n for n in names if n not in channel_colors]
        if missing:
            raise KeyError(f"channel_colors missing {missing}.")
        return [_as_rgb(channel_colors[n]) for n in names]
    if isinstance(channel_colors, str):
        if len(names) != 1:
            raise ValueError("A single channel color requires one channel.")
        return [_as_rgb(channel_colors)]
    colors = list(channel_colors)
    if len(colors) != len(names):
        raise ValueError(f"channel_colors has {len(colors)} colors for {len(names)} channels.")
    return [_as_rgb(c) for c in colors]


def _composite(planes: Sequence[NDArray], tints: Sequence[NDArray]) -> NDArray[np.float32]:
    rgb = np.zeros((*planes[0].shape, 3), dtype=np.float32)
    for plane, tint in zip(planes, tints, strict=True):
        rgb += np.asarray(plane, dtype=np.float32)[..., None] * tint
    return np.clip(rgb, 0.0, 1.0)


def image_planes(
    sdata,
    well: str,
    channels: Sequence[str],
    *,
    scale: int | str | None = None,
    pct: tuple[float, float] = (2.0, 99.8),
) -> tuple[list[NDArray[np.float32]], str]:
    """Return percentile-normalized grayscale planes and the image element name."""
    key = resolve_image_key(sdata, well)
    array = _select_scale(sdata.images[key], scale)
    available = _channel_names(array)
    missing = [name for name in channels if available and name not in available]
    if missing:
        raise KeyError(f"Channel(s) {missing} not in {available}.")
    planes = [_percentile_norm(_extract_channel(array, name, available), *pct) for name in channels]
    return planes, key


def image_planes_window(
    sdata,
    well: str,
    channels: Sequence[str],
    *,
    window: tuple[int, int, int, int],
    scale: int | str | None = None,
    pct: tuple[float, float] = (2.0, 99.8),
) -> tuple[list[NDArray[np.float32]], str]:
    """Like ``image_planes``, but reads only ``(y0, y1, x0, x1)`` of the level."""
    key = resolve_image_key(sdata, well)
    array = _select_scale(sdata.images[key], scale)
    available = _channel_names(array)
    missing = [name for name in channels if available and name not in available]
    if missing:
        raise KeyError(f"Channel(s) {missing} not in {available}.")
    y0, y1, x0, x1 = window
    cropped = window_array(array, y0, y1, x0, x1)
    planes = [_percentile_norm(_extract_channel(cropped, name, available), *pct) for name in channels]
    return planes, key


def image_composite(
    sdata,
    well: str,
    *,
    channels: str | Sequence[str] | None = None,
    channel_colors: str | Sequence[str] | Mapping[str, str] | None = None,
    scale: int | str | None = None,
    pct: tuple[float, float] = (2.0, 99.8),
) -> tuple[NDArray[np.float32], str, list[str], list[NDArray[np.float32]]]:
    """Return a tinted RGB composite, the image key, channel names, and tints."""
    names = _resolve_channel_names(sdata, well, channels)
    planes, key = image_planes(sdata, well, names, scale=scale, pct=pct)
    tints = _resolve_tints(names, channel_colors)
    return _composite(planes, tints), key, names, tints
