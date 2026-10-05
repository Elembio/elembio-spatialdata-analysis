"""Public plotting functions for filtering QC diagnostics.

All six functions (plus the ``plot_cell_paint_*`` aliases) are re-exported at
``ebsd.vis.*`` and ``ebsd.vis.filter.*``.

Prerequisite
------------
Call ``ebsd.pp.compute_qc_metrics(adata, target_type)`` before calling any of
these functions (except ``plot_total_counts_by_batch``), or the metric lookup
will raise a clear ``KeyError``.
"""

from __future__ import annotations

import warnings
from collections.abc import Iterable, Mapping, Sequence

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp

from elembio_spatialdata_analysis.pp.filter._qc_metrics import _obsm_to_df
from elembio_spatialdata_analysis.pp.filter._quantile import QuantileSpec, _resolve_well_batch_quantiles
from elembio_spatialdata_analysis.vis.filter._engine import (
    _QC_METRIC_META,
    _obs_col_to_series,
    _plot_per_well_histogram,
    _prepare_well_grid,
    _resolve_qc_values,
    _resolve_well_layout,
    _row_sums,
    _sorted_well_labels,
)
from elembio_spatialdata_analysis.vis.filter._violin import (
    _normalize_quantile_spec,
    _normalize_quantile_value,
    _render_violin_grid,
    _resolve_default,
)

# ---------------------------------------------------------------------------
# plot_qc_per_well
# ---------------------------------------------------------------------------


def plot_qc_per_well(
    adata: ad.AnnData,
    metric: str = "total_counts",
    *,
    layer: str | None = None,
    well_colname: str = "Well",
    wells: Iterable[str] | None = None,
    orient_flowcell: str | None = None,
    bins: int | str = 40,
    cols: int = 4,
    figsize: tuple[float, float] | None = None,
    sharex: bool = True,
    sharey: bool = True,
    title: str | None = None,
    xlabel: str | None = None,
    facecolor: str = "#4C78A8",
    alpha: float = 0.45,
    density: bool = False,
    logx: bool = False,
    show_median: bool = False,
    show_n: bool = False,
    grid_y: bool = True,
    show_threshold: bool = False,
    threshold: float | None = None,
    threshold_color: str = "red",
    threshold_alpha: float = 0.7,
    threshold_ls: str = "--",
    threshold_label_fs: int = 9,
    show_percentiles: bool = False,
    lower_pct: float | None = 0.025,
    upper_pct: float | None = 0.975,
    percentile_color: str = "red",
    percentile_alpha: float = 0.7,
    percentile_ls: str = "--",
    percentile_label_fs: int = 9,
    show_percentile: bool = False,
    max_pct: float = 0.99,
) -> plt.Figure:
    """Per-well histograms for any QC metric.

    This is the single entry point for all per-well QC histogram plots.
    Pass ``metric`` to select what to plot — any column in
    ``adata.obsm['qc_metrics']``, any column in ``adata.obs``, or the
    special value ``"total_counts"`` to sum ``adata.X`` (or ``layer``).

    Requires ``ebsd.pp.compute_qc_metrics(adata, target_type)`` to have
    been called first for metrics sourced from ``obsm['qc_metrics']``.

    Parameters
    ----------
    adata
        Annotated data matrix.
    metric
        Which QC metric to plot.  Built-in recognised names:
        ``"total_counts"``, ``"total_assignment_rate"``,
        ``"total_nsb_counts"``, ``"total_gapdh_counts"``,
        ``"nuclear_frac"``, ``"n_genes_by_counts"``,
        ``"pct_counts_mt"``, ``"Area"`` (cell area from ``adata.obs``).
        Any other string is looked up in ``adata.obsm['qc_metrics']``
        then ``adata.obs``.
    layer
        Expression layer to use when ``metric="total_counts"``.
    well_colname
        ``adata.obs`` column holding well labels.
    wells
        Explicit well order; defaults to alphanumeric sort.
    orient_flowcell
        Named flowcell layout (e.g. ``"12_well"``).
    bins
        Number of histogram bins or ``"fd"`` for Freedman-Diaconis.
    cols
        Maximum subplot columns.
    figsize
        Figure size ``(width, height)``; auto-computed if ``None``.
    sharex, sharey
        Share axis limits across wells.
    title
        Figure-level super-title.
    xlabel
        Override the x-axis label (auto-derived from ``metric`` if ``None``).
    facecolor, alpha
        Bar fill colour and opacity.
    density
        Plot probability density instead of counts.
    logx
        Use log-scale x axis.
    show_median
        Draw a vertical line at the per-well median.
    show_n
        Annotate each panel with the cell count.
    grid_y
        Show horizontal grid lines.
    show_threshold, threshold
        Draw a single absolute threshold line.
    threshold_color, threshold_alpha, threshold_ls, threshold_label_fs
        Style for the threshold line.
    show_percentiles, lower_pct, upper_pct
        Draw symmetric lower/upper percentile bounds.
    percentile_color, percentile_alpha, percentile_ls, percentile_label_fs
        Style for the percentile lines.
    show_percentile, max_pct
        Draw a single upper-tail percentile line.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if well_colname not in adata.obs:
        raise KeyError(f"obs['{well_colname}'] not found in AnnData.")
    if show_percentiles:
        for pct, label in ((lower_pct, "lower_pct"), (upper_pct, "upper_pct")):
            if pct is not None and not (0 <= pct <= 1):
                raise ValueError(f"{label} must be in [0,1] or None.")
    if show_percentile and not (0 <= max_pct <= 1):
        raise ValueError("max_pct must be in [0, 1].")

    values_all, auto_xlabel = _resolve_qc_values(adata, metric, layer)

    meta = _QC_METRIC_META.get(metric)
    clip_0_1 = meta[1] if meta else False
    effective_logx = logx or (meta[2] if meta else False)

    well_labels = _obs_col_to_series(adata, well_colname)
    well_order, n_rows, n_cols, figsize = _prepare_well_grid(
        well_labels, wells=wells, orient_flowcell=orient_flowcell, max_cols=cols, figsize=figsize
    )
    return _plot_per_well_histogram(
        values_all,
        well_labels,
        well_order,
        n_rows,
        n_cols,
        figsize,
        xlabel=xlabel or auto_xlabel,
        bins=bins,
        sharex=sharex,
        sharey=sharey,
        clip_0_1=clip_0_1,
        title=title,
        facecolor=facecolor,
        alpha=alpha,
        density=density,
        logx=effective_logx,
        show_median=show_median,
        show_n=show_n,
        grid_y=grid_y,
        show_threshold=show_threshold,
        threshold=threshold,
        threshold_color=threshold_color,
        threshold_alpha=threshold_alpha,
        threshold_ls=threshold_ls,
        threshold_label_fs=threshold_label_fs,
        show_percentiles=show_percentiles,
        lower_pct=lower_pct,
        upper_pct=upper_pct,
        percentile_color=percentile_color,
        percentile_alpha=percentile_alpha,
        percentile_ls=percentile_ls,
        percentile_label_fs=percentile_label_fs,
        show_percentile=show_percentile,
        max_pct=max_pct,
    )


# ---------------------------------------------------------------------------
# plot_violin_per_well
# ---------------------------------------------------------------------------


def plot_violin_per_well(
    adata: ad.AnnData,
    metric: str,
    *,
    well_colname: str = "Well",
    wells: Iterable[str] | None = None,
    scale: float = 1.0,
    ylabel: str | None = None,
    title: str | None = None,
    figsize: tuple[float, float] = (5.5, 3.4),
    color: str = "#4C78A8",
    alpha: float = 0.55,
    show_median: bool = True,
    hline: float | None = None,
    hline_label: str | None = None,
    hline_color: str = "#e53e3e",
    dpi: int = 150,
) -> plt.Figure:
    """Per-well violin plot for any QC metric.

    Reads *metric* from ``adata.obsm['qc_metrics']``, then ``adata.obs``,
    then ``adata.obsm['cell_geometry']`` — the same resolution order as
    :func:`plot_qc_per_well`.

    Parameters
    ----------
    adata
        Annotated data matrix.
    metric
        QC metric column name (e.g. ``"pct_counts_mt"``).
    well_colname
        ``adata.obs`` column holding well labels.
    wells
        Explicit well order; defaults to alphanumeric sort.
    scale
        Multiply resolved values by this factor before plotting — e.g. a
        px²→µm² conversion (``PixelSizeUm**2``).  Set *ylabel* to match.
    ylabel
        Y-axis label; defaults to *metric*.
    title
        Figure title; defaults to ``"{ylabel} per Well"``.
    figsize
        Figure size in inches ``(width, height)``.
    color
        Violin fill colour.
    alpha
        Violin fill opacity.
    show_median
        Draw a horizontal line at each per-well median.
    hline
        Optional absolute horizontal reference line (e.g. filter threshold).
    hline_label
        Legend label for the reference line.
    hline_color
        Colour for the reference line.
    dpi
        Figure resolution (used by callers via ``save_or_encode_figure``).

    Returns
    -------
    matplotlib.figure.Figure
    """
    if well_colname not in adata.obs:
        raise KeyError(f"obs['{well_colname}'] not found in AnnData.")

    values_all, auto_ylabel = _resolve_qc_values(adata, metric, None)
    if scale != 1.0:
        values_all = np.asarray(values_all, dtype=np.float64) * scale
    effective_ylabel = ylabel or auto_ylabel
    effective_title = title or f"{effective_ylabel} per Well"

    well_labels = _obs_col_to_series(adata, well_colname)
    well_order = list(wells) if wells is not None else sorted(well_labels.unique().tolist())

    df_v = pd.DataFrame(
        {
            well_colname: well_labels.values,
            metric: pd.to_numeric(values_all, errors="coerce"),
        }
    ).dropna()

    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("#f8f9fa")

    n = len(well_order)
    data_by_well = [df_v.loc[df_v[well_colname] == w, metric].values for w in well_order]

    if any(len(d) > 0 for d in data_by_well):
        parts_v = ax.violinplot(
            data_by_well,
            positions=range(n),
            widths=0.65,
            showmedians=show_median,
            showextrema=False,
        )
        for body in parts_v["bodies"]:
            body.set_facecolor(color)
            body.set_alpha(alpha)
            body.set_edgecolor(_darken(color))
            body.set_linewidth(0.8)
        if show_median and "cmedians" in parts_v:
            parts_v["cmedians"].set_color(_darken(color))
            parts_v["cmedians"].set_linewidth(1.8)

    if hline is not None:
        ax.axhline(
            hline,
            color=hline_color,
            linestyle="--",
            linewidth=1.2,
            label=hline_label or f"{hline:.3g}",
            zorder=4,
        )
        ax.legend(fontsize=7, framealpha=0.8)

    ax.set_xticks(range(n))
    ax.set_xticklabels(well_order, rotation=0, ha="center", fontsize=9, fontweight="500")
    ax.set_ylabel(effective_ylabel, fontsize=9, color="#555")
    ax.set_title(effective_title, fontsize=10.5, fontweight="600", color="#222", pad=8)
    ax.tick_params(axis="y", labelsize=8, colors="#666")
    ax.tick_params(axis="x", length=0)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.6, color="#ddd", zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
    fig.tight_layout()
    return fig


def _darken(hex_color: str, factor: float = 0.65) -> str:
    """Return a darkened version of a hex colour for edges/medians."""
    hex_color = hex_color.lstrip("#")
    r, g, b = (int(hex_color[i : i + 2], 16) for i in (0, 2, 4))
    return f"#{int(r * factor):02x}{int(g * factor):02x}{int(b * factor):02x}"


# ---------------------------------------------------------------------------
# plot_filtered_histograms_per_well
# ---------------------------------------------------------------------------


def plot_filtered_histograms_per_well(
    adata: ad.AnnData,
    flags: pd.DataFrame | pd.Series,
    *,
    metric: str = "total_counts",
    obs_col: str | None = None,
    layer: str | None = None,
    flag_col: str = "passes_all_filters",
    well_col: str = "Well",
    wells: Iterable[str] | None = None,
    bins: int = 40,
    cols: int = 6,
    figsize: tuple[float, float] = (18.0, 6.0),
    colors: Mapping[str, str] | None = None,
    sharex: bool = True,
    sharey: bool = True,
    show_raw: bool = True,
    x_label: str | None = None,
    title: str | None = None,
) -> plt.Figure:
    """Per-well histograms overlaying raw vs. filtered cell distributions.

    Parameters
    ----------
    adata
        Input data.
    flags
        Filter flags aligned to ``adata.obs_names``; the column identified by
        ``flag_col`` must be boolean-castable.
    metric
        What to plot.  ``"total_counts"`` reads from the expression matrix
        (``adata.X`` or ``layer``).  Any other value reads
        ``adata.obsm['qc_metrics']``, then ``adata.obs[obs_col or metric]``.
    obs_col
        Explicit ``adata.obs`` column override when ``metric`` is not found
        in ``obsm['qc_metrics']``.
    layer
        Expression layer for ``metric="total_counts"``.
    flag_col
        Column in ``flags`` indicating passing cells.
        Default ``"passes_all_filters"``.
    well_col
        Observation column for well assignment. Default ``"Well"``.
    wells
        Explicit well order; defaults to alphanumeric sort.
    bins, cols, figsize, colors, sharex, sharey, show_raw, x_label, title
        Standard plotting options.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if isinstance(flags, pd.DataFrame):
        if flag_col not in flags.columns:
            raise KeyError(f"Column '{flag_col}' not found in flags DataFrame.")
        flag_series = flags[flag_col]
    else:
        flag_series = flags

    flag_series = flag_series.reindex(adata.obs_names)
    if flag_series.isna().any():
        raise ValueError("`flags` index does not match `adata.obs_names`.")
    mask_filtered = flag_series.astype(bool).to_numpy()

    # Use _resolve_qc_values for all metrics (including "total_counts") so that
    # pre-computed obsm['qc_metrics'] is always preferred over summing adata.X.
    # Only fall through to a raw X sum when (a) qc_metrics is absent AND
    # (b) metric == "total_counts" AND (c) an explicit layer override was given.
    try:
        values_all, auto_xlabel = _resolve_qc_values(adata, metric, layer)
        xlabel = x_label or auto_xlabel
    except KeyError:
        if metric == "total_counts":
            values_all = _row_sums(adata.layers[layer] if layer is not None else adata.X)
            xlabel = x_label or "Total counts per cell"
        else:
            col = obs_col if obs_col is not None else metric
            if col not in adata.obs:
                raise KeyError(f"obs['{col}'] not found in AnnData.") from None
            values_all = adata.obs[col].to_numpy(dtype=float, na_value=np.nan)
            xlabel = x_label or col

    well_labels = _obs_col_to_series(adata, well_col)
    well_order = list(wells) if wells is not None else _sorted_well_labels(pd.unique(well_labels))
    color_map = colors or {"original": "#4C78A8", "filtered": "#2ca02c"}

    n_wells = len(well_order)
    n_cols = min(cols, max(1, n_wells))
    n_rows = int(np.ceil(n_wells / n_cols))
    well_labels_arr = well_labels.to_numpy()

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(figsize[0], figsize[1] * n_rows),
        sharex=sharex,
        sharey=sharey,
        squeeze=False,
    )

    legend_axis = None

    for idx, well in enumerate(well_order):
        row, col = divmod(idx, n_cols)
        ax = axes[row][col]
        well_mask = well_labels_arr == str(well)

        if not well_mask.any():
            ax.axis("off")
            continue

        filtered_vals = values_all[well_mask & mask_filtered]
        filtered_vals = filtered_vals[np.isfinite(filtered_vals)]
        raw_vals = values_all[well_mask] if show_raw else None
        if raw_vals is not None:
            raw_vals = raw_vals[np.isfinite(raw_vals)]

        components = [filtered_vals]
        if show_raw and raw_vals is not None and raw_vals.size:
            components.append(raw_vals)

        if not any(a.size for a in components):
            ax.axis("off")
            continue

        finite_combined = np.concatenate(components)
        nbins = max(1, min(bins, finite_combined.size))
        vmin, vmax = finite_combined.min(), finite_combined.max()
        bin_edges = (
            np.array([vmin - 0.5, vmax + 0.5])
            if np.isclose(vmin, vmax)
            else np.linspace(vmin, vmax, nbins + 1)
        )

        if show_raw and raw_vals is not None and raw_vals.size:
            ax.hist(
                raw_vals,
                bins=bin_edges,
                color=color_map.get("original", "#4C78A8"),
                alpha=0.25,
                edgecolor="none",
                label="All cells (raw)",
            )

        if filtered_vals.size:
            counts, edges = np.histogram(filtered_vals, bins=bin_edges)
            if show_raw:
                ax.step(
                    edges[:-1],
                    counts,
                    where="post",
                    color=color_map.get("filtered", "#2ca02c"),
                    linewidth=1.2,
                    label="Cells passing filters",
                )
            else:
                ax.hist(
                    filtered_vals,
                    bins=bin_edges,
                    color=color_map.get("filtered", "#2ca02c"),
                    alpha=0.6,
                    edgecolor=color_map.get("filtered", "#2ca02c"),
                    linewidth=0.8,
                    label="Cells passing filters",
                )

        if legend_axis is None:
            legend_axis = ax

        ax.set_title(str(well))
        if col == 0:
            ax.set_ylabel("Cell count", fontsize=12)
        if row == n_rows - 1:
            ax.set_xlabel(xlabel, fontsize=12)
        ax.tick_params(axis="both", labelsize=10)

    for idx in range(n_wells, n_rows * n_cols):
        row, col = divmod(idx, n_cols)
        axes[row][col].axis("off")

    if legend_axis is not None:
        legend_axis.legend(loc="upper right", fontsize=11)

    if title is not None:
        fig.suptitle(title, fontsize=16)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
    else:
        fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# plot_total_counts_by_batch
# ---------------------------------------------------------------------------


def plot_total_counts_by_batch(
    adata: ad.AnnData,
    wells: Iterable[str] | None = None,
    batch_colname: str = "batch_name",
    well_colname: str = "Well",
    cmap: str = "Set2",
    ylim_percentile: float | None = 99.5,
    ncols: int = 4,
    title: str | None = None,
    show_percentiles: bool = False,
    lower_pct: QuantileSpec = 0.025,
    upper_pct: QuantileSpec = 0.975,
    line_color: str = "black",
    line_alpha: float = 0.7,
    line_width: float = 1.5,
    label_fontsize: int = 8,
    label_offset: float = 0.1,
    orient_flowcell: str | None = None,
    logy: bool = False,
    max_violin_points: int = 2000,
) -> plt.Figure:
    """Violin plots of per-cell total counts grouped by batch across wells.

    Percentile overlays respect well/batch overrides from
    ``QuantileSpec`` mappings — mirrors the filtering threshold resolution.

    Parameters
    ----------
    adata
        Annotated data matrix. Requires ``adata.var[batch_colname]``.
    wells
        Explicit well order; defaults to alphanumeric sort.
    batch_colname
        ``adata.var`` column identifying feature batches.
    well_colname
        ``adata.obs`` column identifying wells.
    cmap
        Matplotlib colormap name.
    ylim_percentile
        Upper y-axis percentile cap (``None`` = no cap).
    ncols
        Maximum number of well columns in the grid.
    title
        Figure super-title.
    show_percentiles
        Draw quantile threshold lines on violins.
    lower_pct, upper_pct
        Lower/upper threshold quantiles as scalars or ``QuantileSpec`` mappings.
    line_color, line_alpha, line_width
        Style for threshold lines.
    label_fontsize, label_offset
        Style for threshold labels.
    orient_flowcell
        Optional flowcell orientation key (e.g. ``"12_well"``).
    logy
        Use log scale on the y-axis.
    max_violin_points
        Downsample each per-well/per-batch array for KDE speed.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if well_colname not in adata.obs:
        raise KeyError(f"obs['{well_colname}'] not found")
    if batch_colname not in adata.var:
        raise KeyError(f"var['{batch_colname}'] not found")

    lower_spec = _normalize_quantile_spec(lower_pct, label="lower_pct")
    upper_spec = _normalize_quantile_spec(upper_pct, label="upper_pct")
    lower_default = _resolve_default(lower_spec, fallback=0.025)
    upper_default = _resolve_default(upper_spec, fallback=None)

    if (
        isinstance(lower_spec, (int, float))
        and isinstance(upper_spec, (int, float))
        and upper_spec < 1.0
        and not (lower_spec < upper_spec)
    ):
        raise ValueError("Require lower_pct < upper_pct when both are finite probabilities.")

    obs_wells = _obs_col_to_series(adata, well_colname)
    well_order, grid_rows, grid_cols = _resolve_well_layout(
        obs_wells, wells=wells, orient_flowcell=orient_flowcell, max_cols=ncols
    )
    nrows, ncols = grid_rows, grid_cols

    batches = adata.var[batch_colname].astype(str)
    unique_batches = batches.unique().tolist()
    if not unique_batches:
        raise ValueError("No batches found.")

    count_cols = [f"total_counts_{b}" for b in unique_batches]
    lower_levels = _resolve_well_batch_quantiles(
        lower_spec,
        count_cols=count_cols,
        well_labels=obs_wells,
        default_value=lower_default,
        spec_name="lower_pct",
    )
    upper_levels = _resolve_well_batch_quantiles(
        upper_spec,
        count_cols=count_cols,
        well_labels=obs_wells,
        default_value=upper_default,
        spec_name="upper_pct",
    )
    batch_to_col = {b: f"total_counts_{b}" for b in unique_batches}

    obs_wells_arr = obs_wells.to_numpy()

    batch_counts_cache: dict[str, np.ndarray] = {}
    if "total_counts" in adata.obsm:
        tc_df = _obsm_to_df(adata.obsm["total_counts"], adata.obs_names)
        for batch in unique_batches:
            if batch in tc_df.columns:
                batch_counts_cache[batch] = tc_df[batch].to_numpy(dtype=float)
    else:
        X = adata.X
        is_sparse = sp.issparse(X)
        batches_arr = batches.to_numpy()
        for batch in unique_batches:
            var_mask = batches_arr == batch
            if not var_mask.any():
                continue
            subX = X[:, var_mask]
            batch_counts_cache[batch] = (
                np.asarray(subX.sum(axis=1)).ravel()
                if is_sparse
                else np.asarray(subX.sum(axis=1), dtype=float).ravel()
            )

    _rng = np.random.default_rng(0)
    per_well_data: list[tuple[str, list[np.ndarray], list[str]]] = []
    all_vals: list[np.ndarray] = []

    for well in well_order:
        cell_mask = obs_wells_arr == well
        if not cell_mask.any():
            per_well_data.append((well, [], []))
            continue
        vdata: list[np.ndarray] = []
        labels: list[str] = []
        for batch in unique_batches:
            all_counts = batch_counts_cache.get(batch)
            if all_counts is None:
                continue
            vals = all_counts[cell_mask]
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                continue
            if max_violin_points > 0 and vals.size > max_violin_points:
                vals = _rng.choice(vals, max_violin_points, replace=False)
            vdata.append(vals)
            labels.append(batch)
            all_vals.append(vals)
        per_well_data.append((well, vdata, labels))

    if not all_vals:
        raise ValueError("No values found to plot after filtering.")

    all_vals_concat = np.concatenate(all_vals)
    if logy:
        positive_vals = all_vals_concat[all_vals_concat > 0]
        y_min: float | None = positive_vals.min() if positive_vals.size else 1e-3
        y_max = positive_vals.max() if positive_vals.size else 1.0
    else:
        y_max = float(
            np.max(all_vals_concat)
            if ylim_percentile is None
            else np.percentile(all_vals_concat, ylim_percentile)
        )
        y_min = None

    return _render_violin_grid(
        per_well_data,
        nrows=nrows,
        ncols=ncols,
        y_min=y_min,
        y_max=y_max,
        logy=logy,
        cmap=cmap,
        ylabel="Total counts per cell",
        title=title,
        show_percentiles=show_percentiles,
        lower_levels=lower_levels,
        upper_levels=upper_levels,
        col_key_fn=lambda lbl: batch_to_col.get(lbl, lbl),
        line_color=line_color,
        line_alpha=line_alpha,
        line_width=line_width,
        label_fontsize=label_fontsize,
        label_offset=label_offset,
    )


# ---------------------------------------------------------------------------
# plot_cell_paint_intensities
# ---------------------------------------------------------------------------


def plot_imaging_intensities(
    adata: ad.AnnData,
    wells: Iterable[str] | None = None,
    target_type_colname: str = "target_type",
    target_type_value: str | None = None,
    well_colname: str = "Well",
    cmap: str = "Set2",
    ylim_percentile: float | None = 99.5,
    ncols: int = 4,
    title: str | None = None,
    show_percentiles: bool = False,
    lower_pct: QuantileSpec = 0.025,
    upper_pct: QuantileSpec = 0.975,
    line_color: str = "black",
    line_alpha: float = 0.7,
    line_width: float = 1.5,
    label_fontsize: int = 8,
    label_offset: float = 0.1,
    orient_flowcell: str | None = None,
    logy: bool = False,
    max_violin_points: int = 2000,
    channels: Sequence[str] | None = None,
) -> plt.Figure:
    """Violin plots of the structural imaging (CellPaint) channel intensities per well.

    Each violin is one morphology channel. Channels are selected by
    :func:`~elembio_spatialdata_analysis.tl.morphology.select_morphology_channels`
    (the CellPaint ``CP*`` batch), so the protein/antibody panel, CellProfiler
    features, compartment variants, and backgrounds are excluded. Pass
    ``channels`` to plot an explicit set instead.

    Parameters
    ----------
    adata
        Imaging AnnData object.
    wells
        Explicit well order; defaults to alphanumeric sort.
    target_type_colname
        Deprecated; ignored. Retained for signature compatibility.
    target_type_value
        Deprecated; ignored (channel selection is batch-anchored). Passing a
        non-``None`` value emits a ``DeprecationWarning``; use ``channels`` to
        override.
    well_colname
        ``adata.obs`` column identifying wells.
    cmap
        Matplotlib colormap name.
    ylim_percentile
        Upper y-axis percentile cap (``None`` = no cap).
    ncols
        Maximum number of well columns in the grid.
    title
        Figure super-title.
    show_percentiles
        Draw quantile threshold lines on violins.
    lower_pct, upper_pct
        Lower/upper threshold quantiles as scalars or ``QuantileSpec`` mappings.
    line_color, line_alpha, line_width
        Style for threshold lines.
    label_fontsize, label_offset
        Style for threshold labels.
    orient_flowcell
        Optional flowcell orientation key (e.g. ``"12_well"``).
    logy
        Use log scale on the y-axis.
    max_violin_points
        Downsample each per-well/per-channel array for KDE speed.
    channels
        Explicit channel (``var_names``) list to plot; overrides the
        batch-anchored morphology-channel selection.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if well_colname not in adata.obs:
        raise KeyError(f"obs['{well_colname}'] not found")
    if target_type_value is not None:
        warnings.warn(
            "`target_type_value` is deprecated and ignored; imaging morphology channels "
            "are selected from the CP* batch via select_morphology_channels. "
            "Pass `channels=[...]` to plot an explicit set.",
            DeprecationWarning,
            stacklevel=2,
        )

    if channels is not None:
        primary_var_names: list[str] = list(channels)
    else:
        from elembio_spatialdata_analysis.tl.morphology import select_morphology_channels

        primary_var_names = select_morphology_channels(adata.var)
    if not primary_var_names:
        raise ValueError(
            "No CellPaint morphology channels found (need a CP* batch). "
            "Pass `channels=[...]` to plot an explicit set."
        )

    obs_wells = _obs_col_to_series(adata, well_colname)
    well_order, grid_rows, grid_cols = _resolve_well_layout(
        obs_wells, wells=wells, orient_flowcell=orient_flowcell, max_cols=ncols
    )
    nrows, ncols = grid_rows, grid_cols

    lower_default = _normalize_quantile_value(
        lower_pct if not isinstance(lower_pct, dict) else 0.025, label="lower_pct"
    )
    upper_default = _normalize_quantile_value(
        upper_pct if not isinstance(upper_pct, dict) else 0.975, label="upper_pct"
    )

    if (
        lower_default is not None
        and upper_default is not None
        and 0 < lower_default < 1
        and 0 < upper_default < 1
    ) and lower_default >= upper_default:
        raise ValueError("Require lower_pct < upper_pct when both are finite probabilities.")

    lower_levels = _resolve_well_batch_quantiles(
        lower_pct,
        count_cols=primary_var_names,
        well_labels=obs_wells,
        default_value=lower_default,
        spec_name="lower_pct",
    )
    upper_levels = _resolve_well_batch_quantiles(
        upper_pct,
        count_cols=primary_var_names,
        well_labels=obs_wells,
        default_value=upper_default,
        spec_name="upper_pct",
    )

    col_indices = [adata.var_names.get_loc(ch) for ch in primary_var_names]
    X_primary = adata.X[:, col_indices]
    X_primary = X_primary.toarray() if sp.issparse(X_primary) else np.asarray(X_primary, dtype=float)

    obs_wells_arr = obs_wells.to_numpy()
    _rng = np.random.default_rng(0)
    per_well_data: list[tuple[str, list[np.ndarray], list[str]]] = []
    all_vals: list[np.ndarray] = []

    for well in well_order:
        cell_mask = obs_wells_arr == well
        if not cell_mask.any():
            per_well_data.append((well, [], []))
            continue
        X_well = X_primary[cell_mask]
        vdata: list[np.ndarray] = []
        labels: list[str] = []
        for j, ch in enumerate(primary_var_names):
            vals = X_well[:, j].ravel().astype(float)
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                continue
            if max_violin_points > 0 and vals.size > max_violin_points:
                vals = _rng.choice(vals, max_violin_points, replace=False)
            vdata.append(vals)
            labels.append(ch)
            all_vals.append(vals)
        per_well_data.append((well, vdata, labels))

    if not all_vals:
        raise ValueError("No values found to plot after filtering.")

    all_vals_concat = np.concatenate(all_vals)
    if logy:
        positive_vals = all_vals_concat[all_vals_concat > 0]
        y_min: float | None = positive_vals.min() if positive_vals.size else 1e-3
        y_max = positive_vals.max() if positive_vals.size else 1.0
    else:
        y_max = float(
            np.max(all_vals_concat)
            if ylim_percentile is None
            else np.percentile(all_vals_concat, ylim_percentile)
        )
        y_min = None

    return _render_violin_grid(
        per_well_data,
        nrows=nrows,
        ncols=ncols,
        y_min=y_min,
        y_max=y_max,
        logy=logy,
        cmap=cmap,
        ylabel="Intensity",
        title=title,
        show_percentiles=show_percentiles,
        lower_levels=lower_levels,
        upper_levels=upper_levels,
        col_key_fn=lambda ch: ch,
        line_color=line_color,
        line_alpha=line_alpha,
        line_width=line_width,
        label_fontsize=label_fontsize,
        label_offset=label_offset,
    )


# ---------------------------------------------------------------------------
# plot_imaging_morphology
# ---------------------------------------------------------------------------


def plot_imaging_morphology(
    adata: ad.AnnData,
    metric: str,
    wells: Iterable[str] | None = None,
    channels: Sequence[str] | None = None,
    well_colname: str = "Well",
    cmap: str = "Set2",
    ylim_percentile: float | None = 99.5,
    ncols: int = 4,
    title: str | None = None,
    orient_flowcell: str | None = None,
    logy: bool = False,
    max_violin_points: int = 2000,
) -> plt.Figure:
    """Violin plots of a CellPaint morphology metric per channel per well.

    Morphology features are derived CellProfiler measurements (e.g. Texture,
    Granularity, Intensity statistics, RadialDistribution) whose ``var_names``
    follow the pattern ``<MetricName>_<Channel>`` or contain the channel name
    as a suffix after the last underscore.

    Parameters
    ----------
    adata
        CellPaint AnnData object.
    metric
        The morphology metric prefix to plot.  Examples:

        * ``"MeanIntensity"``  → matches ``Intensity_MeanIntensity_Actin``, etc.
        * ``"Granularity_3"``  → matches ``Granularity_3_Actin``, etc.
        * ``"Texture_Entropy_3_00_256"`` → full texture feature prefix.

        The match is a case-insensitive substring search on ``var_names``; one
        violin is drawn per channel that has a matching feature.
    channels
        Explicit list of channel suffixes to include, e.g.
        ``["Actin", "Nucleus", "Cell-Membrane"]``.  When ``None`` (default)
        all channels found in matching features are used.
    well_colname
        Column in ``adata.obs`` that identifies the well.
    cmap
        Matplotlib colormap name for channel colours.
    ylim_percentile
        Upper y-axis percentile cap (``None`` = no cap).
    ncols
        Maximum number of well columns in the grid.
    title
        Figure super-title.
    orient_flowcell
        Optional flowcell orientation key (e.g. ``"12_well"``).
    logy
        Use log scale on the y-axis.
    max_violin_points
        Downsample each per-well/per-channel array to this many points for
        KDE speed.  Set to 0 to disable.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if well_colname not in adata.obs:
        raise KeyError(f"obs['{well_colname}'] not found")

    metric_lower = metric.lower()
    matching = [v for v in adata.var_names if metric_lower in v.lower()]
    if not matching:
        raise ValueError(
            f"No features matching '{metric}' found in adata.var_names. "
            "Check the metric string (case-insensitive substring match)."
        )

    channel_to_var: dict[str, str] = {}
    for var_name in matching:
        ch = var_name.rsplit("_", 1)[-1]
        if ch not in channel_to_var:
            channel_to_var[ch] = var_name

    if channels is not None:
        channel_to_var = {ch: v for ch, v in channel_to_var.items() if ch in channels}
    if not channel_to_var:
        raise ValueError(f"No channels remain after filtering. Available: {list(channel_to_var.keys())}")

    channel_list = list(channel_to_var.keys())
    var_list = [channel_to_var[ch] for ch in channel_list]

    obs_wells = _obs_col_to_series(adata, well_colname)
    well_order, grid_rows, grid_cols = _resolve_well_layout(
        obs_wells, wells=wells, orient_flowcell=orient_flowcell, max_cols=ncols
    )
    nrows, ncols = grid_rows, grid_cols

    col_indices = [adata.var_names.get_loc(v) for v in var_list]
    X_morph = adata.X[:, col_indices]
    X_morph = X_morph.toarray() if sp.issparse(X_morph) else np.asarray(X_morph, dtype=float)

    obs_wells_arr = obs_wells.to_numpy()
    _rng = np.random.default_rng(0)
    per_well_data: list[tuple[str, list[np.ndarray], list[str]]] = []
    all_vals: list[np.ndarray] = []

    for well in well_order:
        cell_mask = obs_wells_arr == well
        if not cell_mask.any():
            per_well_data.append((well, [], []))
            continue
        X_well = X_morph[cell_mask]
        vdata: list[np.ndarray] = []
        labels: list[str] = []
        for j, ch in enumerate(channel_list):
            vals = X_well[:, j].ravel().astype(float)
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                continue
            if max_violin_points > 0 and vals.size > max_violin_points:
                vals = _rng.choice(vals, max_violin_points, replace=False)
            vdata.append(vals)
            labels.append(ch)
            all_vals.append(vals)
        per_well_data.append((well, vdata, labels))

    if not all_vals:
        raise ValueError("No finite values found to plot.")

    all_vals_concat = np.concatenate(all_vals)
    if logy:
        positive_vals = all_vals_concat[all_vals_concat > 0]
        y_min: float | None = positive_vals.min() if positive_vals.size else 1e-3
        y_max = positive_vals.max() if positive_vals.size else 1.0
    else:
        y_max = float(
            np.max(all_vals_concat)
            if ylim_percentile is None
            else np.percentile(all_vals_concat, ylim_percentile)
        )
        y_min = None

    return _render_violin_grid(
        per_well_data,
        nrows=nrows,
        ncols=ncols,
        y_min=y_min,
        y_max=y_max,
        logy=logy,
        cmap=cmap,
        ylabel=metric,
        title=title or f"CellPaint Morphology: {metric}",
        show_percentiles=False,
        lower_levels={},
        upper_levels={},
        col_key_fn=lambda ch: ch,
        line_color="black",
        line_alpha=0.7,
        line_width=1.5,
        label_fontsize=8,
        label_offset=0.1,
    )


# Backward-compat aliases (CellPaint -> Imaging rename).
plot_cell_paint_intensities = plot_imaging_intensities
plot_cell_paint_morphology = plot_imaging_morphology
