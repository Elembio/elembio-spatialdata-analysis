"""Well-image plotters."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from math import ceil

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Patch
from numpy.typing import NDArray

from elembio_spatialdata_analysis.vis.spatial._access import (
    _match_hw,
    _name_list,
    cell_shape,
    find_cell_well,
    full_to_level,
    has_image,
    image_level_info,
    label_image,
    list_image_wells,
)
from elembio_spatialdata_analysis.vis.spatial._layout import resolve_flowcell_grid
from elembio_spatialdata_analysis.vis.spatial._overlay import (
    _color_series,
    _draw_cells_on_ax,
    _draw_outlines,
    _draw_shape_outline,
    _is_categorical,
)
from elembio_spatialdata_analysis.vis.spatial._raster import (
    _composite,
    _resolve_channel_names,
    _resolve_tints,
    image_planes,
    image_planes_window,
)

_DEFAULT_FIGSIZE = (8.0, 8.0)
_DEFAULT_PCT = (2.0, 99.8)


def _figure_axes(
    n: int,
    *,
    ncols: int | None,
    figsize: tuple[float, float],
    ax: Axes | Sequence[Axes] | None,
    nrows: int | None = None,
):
    import matplotlib.pyplot as plt

    if n < 1:
        raise ValueError("Need at least one panel.")
    if ax is not None:
        axes = [ax] if isinstance(ax, Axes) else list(ax)
        if len(axes) != n:
            raise ValueError(f"Expected {n} axes, got {len(axes)}.")
        return axes[0].figure, axes
    cols = ncols if ncols is not None else min(n, 3)
    cols = max(1, cols)
    rows = nrows if nrows is not None else ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(figsize[0] * cols, figsize[1] * rows))
    flat = [axes] if isinstance(axes, Axes) else list(np.asarray(axes, dtype=object).ravel())
    for extra in flat[n:]:
        extra.set_visible(False)
    return fig, flat[:n]


def _panel_figsize(figsize: tuple[float, float], n_slots: int) -> tuple[float, float]:
    if n_slots > 1 and tuple(figsize) == _DEFAULT_FIGSIZE:
        return (3.2, 3.2)
    return figsize


def _channel_legend(ax: Axes, names: Sequence[str], tints: Sequence[NDArray]) -> None:
    handles = [Patch(facecolor=tint, label=name) for name, tint in zip(names, tints, strict=True)]
    if handles:
        ax.legend(handles=handles, loc="lower right", frameon=True, fontsize=7, title="channels")


def _load_labels(sdata, well: str, *, compartment: str, scale: int | str | None, rgb: NDArray):
    lab, _key = label_image(sdata, well, compartment=compartment, scale=scale)
    return _match_hw(rgb, lab)


def _blank_rgb(lab: NDArray) -> NDArray[np.float32]:
    """Black canvas with the label mask's spatial extent."""
    return np.zeros((*lab.shape, 3), dtype=np.float32)


def _validate_image_visibility(
    *,
    show_image: bool,
    channels,
    channel_colors,
    pct: tuple[float, float],
    separate: bool,
    color,
    outlines: bool,
) -> None:
    if show_image:
        return
    incompatible: list[str] = []
    if channels is not None:
        incompatible.append("channels=")
    if channel_colors is not None:
        incompatible.append("channel_colors=")
    if pct != _DEFAULT_PCT:
        incompatible.append("pct=")
    if separate:
        incompatible.append("separate=True")
    if incompatible:
        options = ", ".join(incompatible)
        raise ValueError(f"show_image=False cannot be combined with {options}.")
    if color is None and not outlines:
        raise ValueError("show_image=False requires color= and/or outlines=True.")


def _subset_adata(adata, well: str):
    if "Well" not in adata.obs.columns:
        return adata
    return adata[adata.obs["Well"].astype(str) == str(well)]


def _shared_continuous_limits(adata, wells: Sequence[str], color: str, percentile: tuple[float, float]):
    chunks: list[np.ndarray] = []
    for well in wells:
        sub = _subset_adata(adata, well)
        if getattr(sub, "n_obs", 0) == 0:
            continue
        numeric = pd.to_numeric(_color_series(sub, color), errors="coerce").to_numpy(dtype=np.float64)
        chunks.append(numeric)
    if not chunks:
        return None, None
    finite = np.concatenate(chunks)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return None, None
    vmin = float(np.percentile(finite, percentile[0]))
    vmax = float(np.percentile(finite, percentile[1]))
    if vmax <= vmin:
        vmax = vmin + 1e-6
    return vmin, vmax


def _shared_categories(adata, wells: Sequence[str], color: str) -> list[str]:
    seen: list[str] = []
    for well in wells:
        sub = _subset_adata(adata, well)
        if getattr(sub, "n_obs", 0) == 0:
            continue
        values = _color_series(sub, color).astype(str)
        for v in pd.unique(values):
            if v not in seen:
                seen.append(str(v))
    return seen


def _load_well_view(
    sdata,
    well: str,
    names: Sequence[str],
    tints: Sequence[NDArray],
    *,
    scale: int | str | None,
    pct: tuple[float, float],
    need_labels: bool,
    compartment: str,
    show_image: bool,
):
    if not show_image:
        lab, label_key = label_image(sdata, well, compartment=compartment, scale=scale)
        return _blank_rgb(lab), lab, label_key
    _resolve_channel_names(sdata, well, names)
    planes, image_key = image_planes(sdata, well, names, scale=scale, pct=pct)
    rgb = _composite(planes, tints)
    lab = None
    if need_labels:
        rgb, lab = _load_labels(sdata, well, compartment=compartment, scale=scale, rgb=rgb)
        h, w = rgb.shape[:2]
        planes = [p[:h, :w] for p in planes]
        rgb = _composite(planes, tints)
    return rgb, lab, image_key


def _plot_single_well(
    sdata,
    well: str,
    *,
    channels,
    channel_colors,
    scale,
    pct,
    outlines,
    outline_width,
    compartment,
    adata,
    color,
    cmap,
    alpha,
    metric_percentile,
    separate,
    show_image,
    ax,
    figsize,
    title,
    ncols,
):
    if show_image:
        names = _resolve_channel_names(sdata, well, channels)
        planes, image_key = image_planes(sdata, well, names, scale=scale, pct=pct)
        tints = _resolve_tints(names, channel_colors)
        lab = None
        if outlines or color is not None:
            rgb0 = _composite(planes, tints)
            rgb0, lab = _load_labels(sdata, well, compartment=compartment, scale=scale, rgb=rgb0)
            h, w = rgb0.shape[:2]
            planes = [p[:h, :w] for p in planes]
        rgb = _composite(planes, tints)
    else:
        lab, image_key = label_image(sdata, well, compartment=compartment, scale=scale)
        rgb = _blank_rgb(lab)
        names, planes, tints = [], [], []

    if color is not None:
        if adata is None:
            raise ValueError("color= requires adata (the cell table from Loader.load_tables).")
        if "Cell" not in adata.obs.columns:
            raise KeyError("Filling cells requires obs['Cell'] (matches well-keyed label pixels).")
        color_names = _name_list(color, what="color")
        sub = _subset_adata(adata, well)
        if sub.n_obs == 0:
            raise ValueError(f"No cells with Well == {well!r}.")
        cell_ids = pd.to_numeric(sub.obs["Cell"], errors="coerce").to_numpy()
        finite = np.isfinite(cell_ids)
        cell_ids = np.where(finite, cell_ids, -1).astype(np.int64)
        fig, axes = _figure_axes(len(color_names), ncols=ncols, figsize=figsize, ax=ax)
        for ax_i, name in zip(axes, color_names, strict=True):
            values = _color_series(sub, name)
            panel_title = f"{name}  |  {image_key}" if title is None else f"{title}  |  {name}"
            _draw_cells_on_ax(
                ax_i,
                rgb,
                lab,
                cell_ids,
                values,
                name,
                cmap=cmap,
                alpha=alpha,
                outlines=outlines,
                outline_width=outline_width,
                metric_percentile=metric_percentile,
                title=panel_title,
            )
        fig.tight_layout()
        return axes[0] if isinstance(color, str) else axes

    if separate:
        fig, axes = _figure_axes(len(names), ncols=ncols, figsize=figsize, ax=ax)
        for ax_i, plane, name, tint in zip(axes, planes, names, tints, strict=True):
            ax_i.imshow(_composite([plane], [tint]))
            if outlines:
                _draw_outlines(ax_i, lab, width=outline_width)
            ax_i.set_title(name if title is None else f"{title}  |  {name}")
            ax_i.axis("off")
        fig.tight_layout()
        return axes[0] if isinstance(channels, str) or len(names) == 1 else axes

    _, axes = _figure_axes(1, ncols=1, figsize=figsize, ax=ax)
    ax0 = axes[0]
    ax0.imshow(rgb)
    if outlines:
        _draw_outlines(ax0, lab, width=outline_width)
    if show_image:
        _channel_legend(ax0, names, tints)
    ax0.set_title(title or image_key)
    ax0.axis("off")
    return ax0


def _centered_window(center: float, size: int, extent: int) -> tuple[int, int, int]:
    """``(start, stop, offset)`` of a fixed-size window centered on ``center``.

    ``start``/``stop`` are clipped to ``extent``; ``offset`` is where the read
    region belongs inside the padded canvas, so the cell stays centered even
    at the edge of a well.
    """
    start = round(center) - size // 2
    read_start = max(start, 0)
    read_stop = min(start + size, extent)
    if read_stop <= read_start:
        raise ValueError("Crop window falls outside the image.")
    return read_start, read_stop, read_start - start


def plot_cell(
    sdata,
    cell_id: int,
    *,
    well: str | None = None,
    channels: str | Sequence[str] | None = None,
    channel_colors: str | Sequence[str] | Mapping[str, str] | None = None,
    scale: int | str | None = 0,
    crop_size: int = 256,
    pct: tuple[float, float] = _DEFAULT_PCT,
    outline: bool = True,
    outline_color: str = "white",
    outline_width: float = 1.5,
    smooth: float = 6.0,
    compartment: str = "Cell",
    ax: Axes | None = None,
    figsize: tuple[float, float] = (5.0, 5.0),
    title: str | None = None,
) -> Axes:
    """Crop the well image around one cell and outline that cell.

    Parameters
    ----------
    sdata
        Opened ``SpatialData`` (``spatialdata.read_zarr`` on a directory store).
    cell_id
        Value from ``obs["Cell"]``. Ids are unique across the store, so the
        well is looked up from a table on ``sdata`` unless ``well=`` is given.
    well
        Skip the table lookup and read this well directly.
    channels, channel_colors
        Same meaning as ``plot_image``.
    scale
        Pyramid level. Default ``0`` is full resolution. Only the crop window
        is read, so this is cheap. Smaller numbers are sharper; larger numbers
        are coarser. The field of view does not change with ``scale``.
    crop_size
        Side of the square field in **full-resolution** pixels. Default
        ``256``. The same value covers the same area at every ``scale``.
        Crops at the edge of a well are padded, keeping the cell centered.
    outline
        Draw ``{well}-shapes-{compartment}`` for this cell as a vector
        contour. Neighbors stay visible as image context but are not drawn.
    smooth
        Contour smoothing, in full-resolution pixels. Stored shapes trace the
        segmentation mask pixel by pixel, so the raw boundary is a staircase;
        the default (``6``) rounds it off for figures. ``0`` draws the exact
        stored boundary, which is what to use when the outline will be measured.

    Returns
    -------
    Axes
        The axes the crop was drawn on.
    """
    if crop_size < 1:
        raise ValueError("crop_size must be at least 1 pixel.")
    if well is None:
        well = find_cell_well(sdata, cell_id)

    geom, _shapes_key = cell_shape(sdata, well, cell_id, compartment=compartment)
    image_key, (h, w), y_ratio, x_ratio = image_level_info(sdata, well, scale)
    centroid = geom.centroid
    cy = full_to_level(centroid.y, y_ratio)
    cx = full_to_level(centroid.x, x_ratio)
    level_h = max(1, round(crop_size / y_ratio))
    level_w = max(1, round(crop_size / x_ratio))
    y0, y1, y_off = _centered_window(cy, level_h, h)
    x0, x1, x_off = _centered_window(cx, level_w, w)

    names = _resolve_channel_names(sdata, well, channels)
    tints = _resolve_tints(names, channel_colors)
    planes, _ = image_planes_window(sdata, well, names, window=(y0, y1, x0, x1), scale=scale, pct=pct)
    read = _composite(planes, tints)
    rgb = np.zeros((level_h, level_w, 3), dtype=read.dtype)
    rgb[y_off : y_off + read.shape[0], x_off : x_off + read.shape[1]] = read

    _, axes = _figure_axes(1, ncols=1, figsize=figsize, ax=ax)
    ax0 = axes[0]
    ax0.imshow(rgb)
    if outline:

        def _to_display(xy: NDArray) -> NDArray:
            out = np.empty_like(xy)
            out[:, 0] = full_to_level(xy[:, 0], x_ratio) - x0 + x_off
            out[:, 1] = full_to_level(xy[:, 1], y_ratio) - y0 + y_off
            return out

        _draw_shape_outline(
            ax0,
            geom,
            to_display=_to_display,
            color=outline_color,
            width=outline_width,
            smooth=smooth,
        )
    _channel_legend(ax0, names, tints)
    ax0.set_title(title if title is not None else f"Cell {cell_id}  |  {image_key}")
    ax0.axis("off")
    return ax0


def plot_image(
    sdata,
    well: str | None = None,
    *,
    wells: Sequence[str] | None = None,
    channels: str | Sequence[str] | None = None,
    channel_colors: str | Sequence[str] | Mapping[str, str] | None = None,
    scale: int | str | None = None,
    pct: tuple[float, float] = _DEFAULT_PCT,
    outlines: bool = False,
    outline_width: float = 1.0,
    compartment: str = "Cell",
    adata=None,
    color: str | Sequence[str] | None = None,
    cmap: str = "viridis",
    alpha: float = 0.55,
    metric_percentile: tuple[float, float] = (2.0, 98.0),
    separate: bool = False,
    show_image: bool = True,
    ax: Axes | Sequence[Axes] | None = None,
    figsize: tuple[float, float] = _DEFAULT_FIGSIZE,
    title: str | None = None,
    ncols: int | None = None,
) -> Axes | list[Axes]:
    """Plot fluorescence well images.

    Omit ``well`` to lay out every well that has ``{well}-images-cell_paint``
    in a plate grid (letter columns × number rows). Pass ``wells=`` for a
    subset. Pass ``well=`` for a single well (one panel, same as before).

    Parameters
    ----------
    sdata
        Opened ``SpatialData`` (``spatialdata.read_zarr`` on a directory store).
    well
        One well id (``"A1"``). Looks up ``{well}-images-cell_paint``.
        ``None`` plots the plate.
    wells
        Optional subset of well ids for the plate view. Cannot be combined
        with ``well=``.
    channels
        One name or a list. ``None`` is Nucleus / Actin / Cell-Membrane
        (whichever exist). Intensities are tinted and overlaid. Shared
        across the plate.
    channel_colors
        Tint for each channel: one color, a list in channel order, or a
        ``{channel: color}`` dict. Defaults: Nucleus red, Actin green,
        Cell-Membrane blue, then a fluorescence palette.
    outlines
        Draw ``{well}-labels-Cell`` borders as vector strokes. No table needed.
    outline_width
        Line weight (matplotlib points). ``0.5`` is thinner than ``1``.
        ``0`` draws none. Default ``1``.
    adata, color
        Optional. Fill each segmented cell by an ``obs`` / ``var`` /
        ``qc_metrics`` value. On the plate, each panel subsets
        ``Well == that well``. Continuous fills share one color scale.
        A list of colors is single-well only.
    alpha
        Table-fill opacity. Default ``0.55`` lets fluorescence remain visible.
        With ``show_image=False``, pass ``alpha=1`` for a solid fill.
    separate
        If True, one panel per channel instead of a single overlay.
        Single-well only.
    show_image
        Whether to draw fluorescence beneath segmentation. Default ``True``.
        Set ``False`` to draw only a table fill and/or label outlines on a
        black background; fluorescence pixels are not loaded. In this mode,
        ``channels``, ``channel_colors``, ``pct``, and ``separate`` do not
        apply and cannot be passed.
    scale
        Pyramid level. ``None`` is the coarsest level; ``0`` is full resolution.
        ``scale=0`` is rejected when more than one well is drawn.

    Returns
    -------
    Axes or list of Axes
        A single ``Axes`` when one panel is drawn; otherwise a list with one
        ``Axes`` per panel (well, channel, or color column).
    """
    _validate_image_visibility(
        show_image=show_image,
        channels=channels,
        channel_colors=channel_colors,
        pct=pct,
        separate=separate,
        color=color,
        outlines=outlines,
    )
    if well is not None and wells is not None:
        raise ValueError("Pass well= or wells=, not both.")
    if well is not None:
        return _plot_single_well(
            sdata,
            well,
            channels=channels,
            channel_colors=channel_colors,
            scale=scale,
            pct=pct,
            outlines=outlines,
            outline_width=outline_width,
            compartment=compartment,
            adata=adata,
            color=color,
            cmap=cmap,
            alpha=alpha,
            metric_percentile=metric_percentile,
            separate=separate,
            show_image=show_image,
            ax=ax,
            figsize=figsize,
            title=title,
            ncols=ncols,
        )

    discovered = list_image_wells(sdata)
    if wells is not None:
        requested = [str(w) for w in wells]
        if not requested:
            raise ValueError("wells= is empty.")
        present_for_layout = requested
    else:
        present_for_layout = discovered
        if not present_for_layout:
            raise ValueError("No wells with images.")

    slot_order, n_rows, n_cols = resolve_flowcell_grid(discovered, wells=present_for_layout)
    allowed = set(present_for_layout)
    present_slots = [w for w in slot_order if has_image(sdata, w) and w in allowed]
    if scale == 0 and len(present_slots) > 1:
        raise ValueError("scale=0 is full resolution; refuse it for more than one well. Use scale=None.")
    if separate and len(slot_order) > 1:
        raise ValueError("separate=True is single-well only.")
    if color is not None and not isinstance(color, str) and len(slot_order) > 1:
        raise ValueError("A list of color= columns is single-well only.")

    if not present_slots:
        raise ValueError("No wells with images in the requested set.")

    if show_image:
        names = _resolve_channel_names(sdata, present_slots[0], channels)
        tints = _resolve_tints(names, channel_colors)
    else:
        names, tints = [], []
    need_labels = bool(outlines or color is not None or not show_image)

    fill_vmin = fill_vmax = None
    categories = None
    if color is not None:
        if adata is None:
            raise ValueError("color= requires adata (the cell table from Loader.load_tables).")
        if "Cell" not in adata.obs.columns:
            raise KeyError("Filling cells requires obs['Cell'] (matches well-keyed label pixels).")
        sample = None
        for slot in present_slots:
            sub = _subset_adata(adata, slot)
            if getattr(sub, "n_obs", 0) > 0:
                sample = _color_series(sub, color)
                break
        if sample is not None and (_is_categorical(sample) or not pd.api.types.is_numeric_dtype(sample)):
            categories = _shared_categories(adata, present_slots, color)
        elif sample is not None:
            fill_vmin, fill_vmax = _shared_continuous_limits(adata, present_slots, color, metric_percentile)

    n_slots = n_rows * n_cols
    fig, axes = _figure_axes(
        n_slots,
        ncols=n_cols,
        nrows=n_rows,
        figsize=_panel_figsize(figsize, n_slots),
        ax=ax,
    )
    last_fill_i = max(
        (i for i, w in enumerate(slot_order) if has_image(sdata, w) and w in allowed), default=-1
    )
    for i, (ax_i, slot) in enumerate(zip(axes, slot_order, strict=False)):
        if not has_image(sdata, slot) or slot not in allowed:
            ax_i.set_title(slot)
            ax_i.axis("off")
            continue
        rgb, lab, image_key = _load_well_view(
            sdata,
            slot,
            names,
            tints,
            scale=scale,
            pct=pct,
            need_labels=need_labels,
            compartment=compartment,
            show_image=show_image,
        )
        panel_title = title or image_key
        if color is not None:
            sub = _subset_adata(adata, slot)
            if sub.n_obs == 0:
                ax_i.imshow(rgb)
                if outlines:
                    _draw_outlines(ax_i, lab, width=outline_width)
                ax_i.set_title(panel_title)
                ax_i.axis("off")
                continue
            cell_ids = pd.to_numeric(sub.obs["Cell"], errors="coerce").to_numpy()
            finite = np.isfinite(cell_ids)
            cell_ids = np.where(finite, cell_ids, -1).astype(np.int64)
            values = _color_series(sub, color)
            _draw_cells_on_ax(
                ax_i,
                rgb,
                lab,
                cell_ids,
                values,
                color,
                cmap=cmap,
                alpha=alpha,
                outlines=outlines,
                outline_width=outline_width,
                metric_percentile=metric_percentile,
                title=panel_title,
                vmin=fill_vmin,
                vmax=fill_vmax,
                categories=categories,
                show_legend=i == last_fill_i,
                show_cbar=i == last_fill_i,
            )
        else:
            ax_i.imshow(rgb)
            if outlines:
                _draw_outlines(ax_i, lab, width=outline_width)
            ax_i.set_title(panel_title)
            ax_i.axis("off")
            if show_image and i == last_fill_i:
                _channel_legend(ax_i, names, tints)
    for extra in axes[len(slot_order) :]:
        extra.set_visible(False)
    fig.tight_layout()
    return axes
