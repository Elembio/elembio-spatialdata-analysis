"""Private violin-plot engine shared by plot_total_counts_by_batch and plot_cell_paint_intensities.

Nothing in this module is part of the public API.
"""

from __future__ import annotations

from collections.abc import Mapping

import matplotlib.pyplot as plt
import numpy as np

from elembio_spatialdata_analysis.vis.filter._engine import _display_pct, _fs

# ---------------------------------------------------------------------------
# Quantile normalisation helpers (used by violin public functions)
# ---------------------------------------------------------------------------


def _normalize_quantile_value(value: object, *, label: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, str):
        s = value.strip().lower().replace("%", "").replace("th", "")
        if not s:
            raise ValueError(f"Quantile '{label}' is empty after stripping.")
        value = float(s)
    try:
        q = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise TypeError(f"Quantile '{label}' must be numeric or None.") from exc
    if q < 0:
        raise ValueError(f"Quantile '{label}' cannot be negative (got {q}).")
    if q > 1.0:
        if q <= 100.0:
            q /= 100.0
        else:
            raise ValueError(f"Quantile '{label}' percent must be ≤ 100 (got {q}).")
    return q


def _normalize_quantile_spec(value: object, *, label: str) -> object:
    if value is None or isinstance(value, (int, float, str)):
        return _normalize_quantile_value(value, label=label)
    if not isinstance(value, Mapping):
        raise TypeError(f"Quantile specification '{label}' must be numeric, mapping, or None.")
    normalized: dict[str, object] = {}
    for key, val in value.items():
        key_label = f"{label}[{key}]"
        normalized[str(key)] = (
            _normalize_quantile_spec(val, label=key_label)
            if isinstance(val, Mapping)
            else _normalize_quantile_value(val, label=key_label)
        )
    return normalized


def _resolve_default(spec: object, *, fallback: float | None) -> float | None:
    if spec is None:
        return None
    if isinstance(spec, Mapping):
        for k in ("__default__", "default"):
            if k in spec and spec[k] is not None:
                return float(spec[k])  # type: ignore[arg-type]
        return fallback
    return float(spec)  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# Shared violin rendering engine
# ---------------------------------------------------------------------------


def _render_violin_grid(
    per_well_data: list[tuple[str, list[np.ndarray], list[str]]],
    *,
    nrows: int,
    ncols: int,
    y_min: float | None,
    y_max: float,
    logy: bool,
    cmap: str,
    ylabel: str,
    title: str | None,
    show_percentiles: bool,
    lower_levels: dict,
    upper_levels: dict,
    col_key_fn: object,  # callable(label) -> col_key for threshold lookup
    line_color: str,
    line_alpha: float,
    line_width: float,
    label_fontsize: int,
    label_offset: float,
) -> plt.Figure:
    """Render a per-well violin grid from pre-computed per_well_data.

    Parameters
    ----------
    per_well_data
        List of ``(well_label, [array_per_group], [group_label])`` tuples.
    nrows, ncols
        Grid dimensions.
    y_min, y_max
        Y-axis bounds (``y_min`` only used in log mode).
    logy
        Use log y-axis.
    cmap
        Matplotlib colormap name.
    ylabel
        Y-axis label for the leftmost column.
    title
        Figure super-title.
    show_percentiles
        Draw threshold lines on violins.
    lower_levels, upper_levels
        Dicts ``{well: {col_key: quantile_value}}`` from
        ``_resolve_well_batch_quantiles``.
    col_key_fn
        Callable mapping a group label to the key used in
        ``lower_levels``/``upper_levels``.
    line_color, line_alpha, line_width
        Style for threshold lines.
    label_fontsize, label_offset
        Style for threshold labels.
    """
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(max(8, 4 * ncols), 4 * nrows),
        squeeze=False,
        sharey=True,
    )
    cmap_obj = plt.get_cmap(cmap)
    label_fontsize_scaled = _fs(label_fontsize)
    tick_fs = _fs(10)
    title_fs = _fs(12)
    ylabel_fs = _fs(11)

    for idx, (well, vdata, labels) in enumerate(per_well_data):
        ax = axes[idx // ncols, idx % ncols]
        if not vdata:
            ax.set_visible(False)
            continue

        parts = ax.violinplot(vdata, showmeans=False, showmedians=True, showextrema=False)
        colors = [cmap_obj(j / max(1, len(labels))) for j in range(len(labels))]
        for body, color in zip(parts["bodies"], colors, strict=False):
            body.set_facecolor(color)
            body.set_edgecolor("black")
            body.set_alpha(0.8)
        if "cmedians" in parts:
            parts["cmedians"].set_color("black")
            parts["cmedians"].set_linewidth(1.2)

        if logy:
            lower = max((y_min or 1e-3) * 0.9, 1e-6)
            upper = max(y_max * 1.1, lower * 10)
            ax.set_yscale("log")
            ax.set_ylim(lower, upper)
        else:
            ax.set_ylim(-0.05 * y_max, y_max)

        if show_percentiles:
            lower_targets = lower_levels.get(well, {}) or {}
            upper_targets = upper_levels.get(well, {}) or {}
            for j, (vals, lbl) in enumerate(zip(vdata, labels, strict=False)):
                x_left = j + 0.9
                x_right = j + 1.1
                col_key = col_key_fn(lbl)  # type: ignore[operator]
                min_q = lower_targets.get(col_key)
                max_q = upper_targets.get(col_key)
                for q_val in [min_q, max_q]:
                    if q_val is None or (q_val == max_q and q_val >= 1.0):
                        continue
                    y_line = np.quantile(vals, q_val)
                    ax.hlines(
                        y_line, x_left, x_right, color=line_color, linewidth=line_width, alpha=line_alpha
                    )
                    ax.text(
                        x_right + label_offset,
                        y_line,
                        _display_pct(q_val),
                        va="center",
                        ha="left",
                        fontsize=label_fontsize_scaled,
                        alpha=line_alpha,
                    )

        ax.set_xticks(np.arange(1, len(labels) + 1))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=tick_fs)
        ax.set_title(f"Well {well}", fontsize=title_fs)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if idx % ncols == 0:
            ax.set_ylabel(ylabel, fontsize=ylabel_fs)

    for idx in range(len(per_well_data), nrows * ncols):
        axes[idx // ncols, idx % ncols].set_visible(False)

    if title:
        fig.suptitle(title, y=0.995, fontsize=_fs(12))
    fig.tight_layout()
    return fig
