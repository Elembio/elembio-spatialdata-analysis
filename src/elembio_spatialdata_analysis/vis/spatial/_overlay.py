"""Segmentation strokes and table-value cell fills."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from matplotlib import colormaps as mpl_colormaps
from matplotlib.axes import Axes
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.patches import Patch, PathPatch
from matplotlib.path import Path
from numpy.typing import NDArray
from pandas.api.types import is_bool_dtype, is_numeric_dtype, is_object_dtype

from elembio_spatialdata_analysis.pp.filter._constants import QC_METRICS_OBSM_KEY


def _boundary_segments(lab: NDArray) -> NDArray[np.float64]:
    """Line segments along label boundaries, in imshow (x=col, y=row) coords."""
    lab = np.asarray(lab)
    parts: list[NDArray[np.float64]] = []
    dv = lab[1:, :] != lab[:-1, :]
    rows, cols = np.nonzero(dv)
    if rows.size:
        y = rows.astype(np.float64) + 0.5
        x0 = cols.astype(np.float64) - 0.5
        x1 = cols.astype(np.float64) + 0.5
        parts.append(np.stack([np.column_stack([x0, y]), np.column_stack([x1, y])], axis=1))
    dh = lab[:, 1:] != lab[:, :-1]
    rows, cols = np.nonzero(dh)
    if rows.size:
        x = cols.astype(np.float64) + 0.5
        y0 = rows.astype(np.float64) - 0.5
        y1 = rows.astype(np.float64) + 0.5
        parts.append(np.stack([np.column_stack([x, y0]), np.column_stack([x, y1])], axis=1))
    if not parts:
        return np.zeros((0, 2, 2), dtype=np.float64)
    return np.concatenate(parts, axis=0)


def _draw_outlines(ax: Axes, lab: NDArray, *, width: float = 1.0) -> None:
    """Stroke mask boundaries. ``width`` is a matplotlib line weight, not a pixel fill."""
    w = float(width)
    if w <= 0:
        return
    segs = _boundary_segments(lab)
    if segs.size == 0:
        return
    ax.add_collection(LineCollection(segs, colors="white", linewidths=w, antialiased=True))


def _chaikin(ring: NDArray[np.float64], iterations: int = 2) -> NDArray[np.float64]:
    """Corner-cut a closed ring, turning right angles into rounded corners."""
    pts = ring[:-1] if np.allclose(ring[0], ring[-1]) else ring
    if len(pts) < 3:
        return ring
    for _ in range(iterations):
        nxt = np.roll(pts, -1, axis=0)
        cut = np.empty((len(pts) * 2, 2), dtype=np.float64)
        cut[0::2] = 0.75 * pts + 0.25 * nxt
        cut[1::2] = 0.25 * pts + 0.75 * nxt
        pts = cut
    return np.vstack([pts, pts[:1]])


def _shape_rings(geom, *, smooth: float = 0.0) -> list[NDArray[np.float64]]:
    """Exterior and interior rings of a Polygon or MultiPolygon, as xy arrays.

    Stored shapes trace the segmentation mask pixel by pixel. ``smooth`` is a
    simplification tolerance in the geometry's own units; simplified rings are
    then corner-cut so the outline reads as a cell boundary.
    """
    if smooth > 0:
        simplified = geom.simplify(float(smooth))
        if not simplified.is_empty and simplified.geom_type in {"Polygon", "MultiPolygon"}:
            geom = simplified
    parts = list(getattr(geom, "geoms", [])) or [geom]
    rings: list[NDArray[np.float64]] = []
    for part in parts:
        exterior = getattr(part, "exterior", None)
        if exterior is None:
            continue
        for ring in (exterior, *part.interiors):
            xy = np.asarray(ring.coords, dtype=np.float64)[:, :2]
            if xy.shape[0] < 3:
                continue
            rings.append(_chaikin(xy) if smooth > 0 else xy)
    if not rings:
        raise ValueError(f"Geometry {geom.geom_type!r} has no polygon boundary to draw.")
    return rings


def _draw_shape_outline(
    ax: Axes,
    geom,
    *,
    to_display,
    color: str = "white",
    width: float = 1.5,
    smooth: float = 0.0,
) -> None:
    """Stroke one stored shape as a vector contour. ``to_display`` maps xy."""
    if float(width) <= 0:
        return
    verts: list[NDArray[np.float64]] = []
    codes: list[int] = []
    for ring in _shape_rings(geom, smooth=smooth):
        mapped = to_display(ring)
        verts.append(mapped)
        codes.extend([Path.MOVETO, *[Path.LINETO] * (len(mapped) - 1)])
    path = Path(np.concatenate(verts, axis=0), codes)
    ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, linewidth=float(width), antialiased=True))


def _color_series(adata, color: str) -> pd.Series:
    if color in adata.obs.columns:
        series = adata.obs[color]
    elif color in adata.var_names:
        idx = list(adata.var_names).index(color)
        x = adata.X[:, idx]
        values = np.asarray(x.toarray()).ravel() if hasattr(x, "toarray") else np.asarray(x).ravel()
        series = pd.Series(values, index=adata.obs_names, name=color)
    elif QC_METRICS_OBSM_KEY in adata.obsm:
        qc = adata.obsm[QC_METRICS_OBSM_KEY]
        qc_df = qc if isinstance(qc, pd.DataFrame) else pd.DataFrame(qc, index=adata.obs_names)
        if color not in qc_df.columns:
            raise KeyError(f"{color!r} not found in obs, var_names, or obsm['{QC_METRICS_OBSM_KEY}'].")
        series = qc_df[color]
    else:
        raise KeyError(f"{color!r} not found in obs columns or var_names.")
    if not series.index.equals(adata.obs_names):
        series = series.reindex(adata.obs_names)
    return series


def _is_categorical(series: pd.Series) -> bool:
    return bool(
        isinstance(series.dtype, pd.CategoricalDtype) or is_bool_dtype(series) or is_object_dtype(series)
    )


def _paint_categorical(lab: NDArray, cell_ids: NDArray[np.int64], codes: NDArray[np.int64], n_cats: int):
    lut_size = int(max(lab.max(initial=0), cell_ids.max(initial=0), 0)) + 1
    lut = np.zeros(lut_size, dtype=np.int32)
    valid = (cell_ids >= 0) & (cell_ids < lut_size)
    lut[cell_ids[valid]] = (codes[valid] + 1).astype(np.int32)
    mapped = np.zeros_like(lab, dtype=np.int32)
    ok = (lab >= 0) & (lab < lut_size)
    mapped[ok] = lut[lab[ok]]
    palette = mpl_colormaps["tab20"].resampled(max(n_cats, 1))(np.linspace(0, 1, max(n_cats, 1)))
    rgba = np.zeros((*lab.shape, 4), dtype=np.float32)
    for i in range(n_cats):
        rgba[mapped == i + 1] = (*palette[i, :3], 1.0)
    return rgba, palette


def _paint_continuous(
    lab: NDArray,
    cell_ids: NDArray[np.int64],
    values: NDArray[np.float64],
    percentile: tuple[float, float],
    *,
    vmin: float | None = None,
    vmax: float | None = None,
):
    lut_size = int(max(lab.max(initial=0), cell_ids.max(initial=0), 0)) + 1
    lut = np.full(lut_size, np.nan, dtype=np.float32)
    valid = (cell_ids >= 0) & (cell_ids < lut_size) & np.isfinite(values)
    lut[cell_ids[valid]] = values[valid].astype(np.float32)
    img = np.full(lab.shape, np.nan, dtype=np.float32)
    ok = (lab >= 0) & (lab < lut_size)
    img[ok] = lut[lab[ok]]
    finite = lut[np.isfinite(lut)]
    if vmin is None or vmax is None:
        if finite.size == 0:
            return img, np.nan, np.nan
        lo = float(np.percentile(finite, percentile[0]))
        hi = float(np.percentile(finite, percentile[1]))
        if vmin is None:
            vmin = lo
        if vmax is None:
            vmax = hi
    vmin_f = float(vmin)
    vmax_f = float(vmax)
    if vmax_f <= vmin_f:
        vmax_f = vmin_f + 1e-6
    return img, vmin_f, vmax_f


def _draw_cells_on_ax(
    ax: Axes,
    rgb: NDArray,
    lab: NDArray,
    cell_ids: NDArray[np.int64],
    values: pd.Series,
    color: str,
    *,
    cmap: str,
    alpha: float,
    outlines: bool,
    outline_width: float,
    metric_percentile: tuple[float, float],
    title: str,
    vmin: float | None = None,
    vmax: float | None = None,
    categories: Sequence[str] | None = None,
    show_legend: bool = True,
    show_cbar: bool = True,
):
    ax.imshow(rgb)
    if _is_categorical(values) or not is_numeric_dtype(values):
        if categories is not None:
            cats = pd.Categorical(values.astype(str), categories=[str(c) for c in categories])
            codes = np.asarray(cats.codes)
            cat_labels = [str(c) for c in cats.categories]
        else:
            cats = values.astype("category")
            codes = np.asarray(cats.cat.codes)
            cat_labels = [str(c) for c in cats.cat.categories]
        overlay, palette = _paint_categorical(lab, cell_ids, codes, len(cat_labels))
        overlay[..., 3] *= alpha
        overlay[lab == 0, 3] = 0
        ax.imshow(overlay)
        handles = [Patch(facecolor=palette[i], label=cat_labels[i]) for i in range(len(cat_labels))]
        if show_legend and handles:
            ax.legend(handles=handles, loc="upper right", frameon=True, fontsize=7, title=color)
    else:
        numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=np.float64)
        metric_img, vmin, vmax = _paint_continuous(
            lab, cell_ids, numeric, metric_percentile, vmin=vmin, vmax=vmax
        )
        im = ax.imshow(
            np.ma.masked_invalid(metric_img),
            cmap=mpl_colormaps[cmap] if isinstance(cmap, str) else cmap,
            norm=Normalize(vmin=vmin, vmax=vmax),
            interpolation="nearest",
            alpha=alpha,
        )
        if show_cbar and np.isfinite(vmin) and np.isfinite(vmax):
            cbar = ax.figure.colorbar(im, ax=ax, fraction=0.035, pad=0.03, shrink=0.8)
            cbar.set_label(color)
    if outlines:
        _draw_outlines(ax, lab, width=outline_width)
    ax.set_title(title)
    ax.axis("off")
