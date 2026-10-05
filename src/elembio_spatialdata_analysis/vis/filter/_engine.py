"""Private histogram engine and layout helpers for vis.filter.

Nothing in this module is part of the public API.  All symbols are prefixed
with ``_`` or are module-level constants consumed only by ``plots.py``.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp

from elembio_spatialdata_analysis.pp.filter._constants import OBSM_CELL_GEOMETRY, QC_METRICS_OBSM_KEY
from elembio_spatialdata_analysis.pp.filter._qc_metrics import _obsm_to_df

# ---------------------------------------------------------------------------
# Flowcell orientation registry
# ---------------------------------------------------------------------------

# Columns are alphabetic (A..F), rows are numeric (1, 2)
FLOWCELL_ORIENTATIONS: dict[str, dict[str, Sequence[str]]] = {
    "12_well": {
        "cols": list("ABCDEF"),
        "rows": ["1", "2"],
    },
    "6_well": {
        "cols": list("ABC"),
        "rows": ["1", "2"],
    },
}

# Map from expected well set (frozenset) → orientation key, for auto-detection.
_FLOWCELL_WELL_SETS: dict[frozenset[str], str] = {
    frozenset(f"{col}{row}" for row in layout["rows"] for col in layout["cols"]): key
    for key, layout in FLOWCELL_ORIENTATIONS.items()
}


def _auto_detect_flowcell(well_labels: Sequence[object]) -> str | None:
    """Infer flowcell orientation from the set of well labels.

    Returns the orientation key (e.g. ``"12_well"``) when the observed wells
    are a subset of a known layout's wells, ``None`` otherwise.
    """
    observed = frozenset(str(lab).upper() for lab in well_labels if str(lab).strip())
    for known_wells, key in _FLOWCELL_WELL_SETS.items():
        if observed <= known_wells:
            return key
    return None


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------


def _fs(sz: float, min_sz: int = 10) -> int:
    """Scale a font size to at least ``min_sz``."""
    return max(int(sz), min_sz)


def _display_pct(q: float) -> str:
    """Format a quantile fraction as a percentage string."""
    return f"{q * 100:g}%"


# ---------------------------------------------------------------------------
# Obs materialization (lazy/xarray-backed obs)
# ---------------------------------------------------------------------------


def _obs_col_to_series(adata: ad.AnnData, col: str) -> pd.Series:
    """Materialize an obs column to a pandas Series.

    With ``read_lazy``, ``adata.obs`` is xarray-backed; iterating over a column
    triggers per-element loads and is extremely slow. This helper forces a
    single materialization up front.
    """
    obs_col = adata.obs[col]
    if isinstance(obs_col, pd.Series):
        return obs_col.astype(str)
    # xarray or other lazy: materialize via to_numpy/values, then wrap
    if hasattr(obs_col, "to_numpy"):
        arr = obs_col.to_numpy()
    elif hasattr(obs_col, "values"):
        arr = obs_col.values
    else:
        arr = np.asarray(obs_col)
    # Avoid numpy astype(str) on pandas StringDType; use pandas for conversion
    out = pd.Series(arr, index=adata.obs_names)
    return out.astype(str)


# ---------------------------------------------------------------------------
# Well layout helpers
# ---------------------------------------------------------------------------


def _sorted_well_labels(labels: Sequence[object]) -> list[str]:
    """Sort wells by numeric suffix then alphabetic prefix (A1, B1, …, A2, …)."""

    def _key(label: object) -> tuple[int, str]:
        text = str(label)
        prefix = "".join(filter(str.isalpha, text)).upper()
        digits = "".join(filter(str.isdigit, text))
        suffix = int(digits) if digits else 0
        return suffix, prefix

    return sorted((str(lab) for lab in labels), key=_key)


def _resolve_well_layout(
    well_labels: Sequence[object],
    *,
    wells: Iterable[str] | None,
    orient_flowcell: str | None,
    max_cols: int,
) -> tuple[list[str], int, int]:
    """Determine well ordering and subplot grid size.  Returns (well_order, n_rows, n_cols).

    When ``orient_flowcell`` is ``None`` (default) the layout is auto-detected
    from the observed well labels — e.g. A1–F2 → ``"12_well"``.  Pass an
    explicit key to override, or ``"none"`` to suppress auto-detection.
    """
    normalized_labels = [str(label) for label in well_labels]
    unique_labels = list(np.unique(normalized_labels))

    if wells is not None:
        well_order = [str(well) for well in wells]
        n_wells = len(well_order)
        n_cols = min(max_cols, max(1, n_wells))
        n_rows = int(np.ceil(n_wells / n_cols)) if n_cols else 1
        return well_order, n_rows, n_cols

    # Resolve orientation key: explicit > auto-detect > None
    if orient_flowcell is not None:
        key = orient_flowcell.strip().lower().replace("-", "_")
        resolved_key = None if key == "none" else key
    else:
        resolved_key = _auto_detect_flowcell(unique_labels)

    if resolved_key:
        layout = FLOWCELL_ORIENTATIONS.get(resolved_key)
        if layout is None:
            known = ", ".join(sorted(FLOWCELL_ORIENTATIONS))
            raise KeyError(
                f"Unknown flowcell orientation '{resolved_key}'. Available options: {known or 'none'}."
            )
        cols_order = [str(col).upper() for col in layout["cols"]]
        rows_order = [str(row) for row in layout["rows"]]
        if not rows_order or not cols_order:
            raise ValueError(f"Flowcell orientation '{resolved_key}' is malformed.")
        base_order = [f"{col}{row}" for row in rows_order for col in cols_order]
        extra_labels = [lab for lab in _sorted_well_labels(unique_labels) if lab not in base_order]
        well_order = base_order + extra_labels
        n_cols = len(cols_order)
        n_rows = int(np.ceil(len(well_order) / n_cols))
        return well_order, n_rows, n_cols

    well_order = _sorted_well_labels(unique_labels)
    n_wells = len(well_order)
    n_cols = min(max_cols, max(1, n_wells))
    n_rows = int(np.ceil(n_wells / n_cols)) if n_cols else 1
    return well_order, n_rows, n_cols


def _default_figsize(n_rows: int, n_cols: int) -> tuple[float, float]:
    """Figure size for an ``n_rows × n_cols`` per-well grid.

    A lone subplot — a single-well ("mega well") run — gets a comfortable
    standalone size instead of the tiny per-cell size that suits a dense grid.
    Every multi-cell grid (6-/12-well plates, etc.) keeps the 3.0 × 2.6
    per-cell sizing unchanged.
    """
    if n_rows == 1 and n_cols == 1:
        return 6.0, 4.0
    return 3.0 * n_cols, 2.6 * n_rows


def _prepare_well_grid(
    well_labels: pd.Series,
    *,
    wells: Iterable[str] | None,
    orient_flowcell: str | None,
    max_cols: int,
    figsize: tuple[float, float] | None,
) -> tuple[list[str], int, int, tuple[float, float]]:
    well_order, n_rows, n_cols = _resolve_well_layout(
        well_labels, wells=wells, orient_flowcell=orient_flowcell, max_cols=max_cols
    )
    return well_order, n_rows, n_cols, figsize or _default_figsize(n_rows, n_cols)


# ---------------------------------------------------------------------------
# Binning helpers
# ---------------------------------------------------------------------------


def _freedman_diaconis_bins(values: np.ndarray) -> int:
    """Freedman-Diaconis rule on linear scale; guards against degenerate inputs."""
    filtered = np.asarray(values, dtype=np.float64)
    filtered = filtered[np.isfinite(filtered)]
    if filtered.size < 2:
        return 1
    q75, q25 = np.percentile(filtered, [75, 25])
    iqr = q75 - q25
    if iqr <= 0:
        return max(1, min(40, filtered.size))
    bin_width = 2 * iqr / np.cbrt(filtered.size)
    if bin_width <= 0:
        return max(1, min(40, filtered.size))
    span = filtered.max() - filtered.min()
    if np.isclose(span, 0):
        return 1
    return int(np.clip(np.ceil(span / bin_width), 1, 512))


def _freedman_diaconis_bins_log10(values: np.ndarray) -> int:
    """Freedman-Diaconis rule on log10 scale."""
    filtered = np.asarray(values, dtype=np.float64)
    filtered = filtered[np.isfinite(filtered) & (filtered > 0)]
    if filtered.size < 2:
        return 1
    log_vals = np.log10(filtered)
    q75, q25 = np.percentile(log_vals, [75, 25])
    iqr = q75 - q25
    if iqr <= 0:
        return max(1, min(40, log_vals.size))
    bin_width = 2 * iqr / np.cbrt(log_vals.size)
    if bin_width <= 0:
        return max(1, min(40, log_vals.size))
    span = log_vals.max() - log_vals.min()
    if np.isclose(span, 0):
        return 1
    return int(np.clip(np.ceil(span / bin_width), 1, 512))


def _log_bin_edges(values_list: list[np.ndarray], nbins: int | str) -> np.ndarray:
    """Compute log-spaced bin edges from a list of per-well arrays."""
    positives = []
    for raw_arr in values_list:
        arr = np.asarray(raw_arr, dtype=float)
        arr = arr[np.isfinite(arr) & (arr > 0)]
        if arr.size:
            positives.append(arr)
    if not positives:
        return np.array([1.0, 10.0])
    allpos = np.concatenate(positives)
    lo, hi = allpos.min(), allpos.max()
    if np.isclose(lo, hi):
        return np.array([lo / 1.2, hi * 1.2])
    k = (
        _freedman_diaconis_bins_log10(allpos)
        if (isinstance(nbins, str) and nbins.lower() == "fd")
        else max(1, int(nbins))
    )
    return np.logspace(np.log10(lo), np.log10(hi), k + 1)


def _linear_bin_edges(
    values_all: np.ndarray,
    values_per_well: list[np.ndarray],
    nbins: int | str,
) -> np.ndarray | None:
    """Compute linear shared bin edges from per-well value lists."""
    present = [(v.min(), v.max()) for v in values_per_well if v.size]
    if not present:
        return None
    gmin = min(a for a, _ in present)
    gmax = max(b for _, b in present)
    if np.isclose(gmin, gmax):
        return np.array([gmin - 0.5, gmax + 0.5])
    nb = (
        _freedman_diaconis_bins(values_all)
        if (isinstance(nbins, str) and nbins.lower() == "fd")
        else max(1, int(nbins))
    )
    return np.linspace(gmin, gmax, nb + 1)


def _per_panel_bin_edges(
    vals: np.ndarray,
    nbins: int | str,
    logx: bool,
) -> np.ndarray:
    """Compute bin edges for a single panel when sharex=False."""
    if logx:
        return _log_bin_edges([vals], nbins)
    vmin, vmax = vals.min(), vals.max()
    if np.isclose(vmin, vmax):
        return np.array([vmin - 0.5, vmax + 0.5])
    nb = (
        _freedman_diaconis_bins(vals)
        if (isinstance(nbins, str) and nbins.lower() == "fd")
        else max(1, int(nbins))
    )
    return np.linspace(vmin, vmax, nb + 1)


# ---------------------------------------------------------------------------
# Matrix helpers
# ---------------------------------------------------------------------------


def _row_sums(mat: Any) -> np.ndarray:
    """Sparse-safe row sums."""
    if sp.issparse(mat):
        return np.asarray(mat.sum(axis=1)).ravel().astype(float)
    return np.nansum(np.asarray(mat, dtype=float), axis=1)


def _compute_totals(adata: ad.AnnData, layer: str | None) -> np.ndarray:
    matrix = adata.layers[layer] if layer is not None else adata.X
    return _row_sums(matrix)


# ---------------------------------------------------------------------------
# QC metric registry and resolver
# ---------------------------------------------------------------------------

# Each entry: (xlabel, clip_0_1, default_logx)
_QC_METRIC_META: dict[str, tuple[str, bool, bool]] = {
    "total_counts": ("Total counts per cell", False, False),
    "total_assignment_rate": ("Assignment rate", True, False),
    "total_nsb_counts": ("NSB counts per cell", False, False),
    "total_gapdh_counts": ("GAPDH counts per cell", False, False),
    "nuclear_frac": ("Nuclear assignment rate", True, False),
    "n_genes_by_counts": ("Number of genes detected per cell", False, False),
    "pct_counts_mt": ("Mitochondrial gene percentage", True, False),
    # Read from obsm['cell_geometry'] (or obs), not qc_metrics
    "Area": ("Cell area (px²)", False, False),
}


def _require_qc_metric(adata: ad.AnnData, col: str) -> np.ndarray:
    """Read ``col`` from ``obsm['qc_metrics']`` or raise with a clear message."""
    if QC_METRICS_OBSM_KEY in adata.obsm:
        qc = _obsm_to_df(adata.obsm[QC_METRICS_OBSM_KEY], adata.obs_names)
        if col in qc.columns:
            return qc[col].reindex(adata.obs_names).to_numpy()
        raise KeyError(
            f"Column '{col}' not found in adata.obsm['{QC_METRICS_OBSM_KEY}']. "
            f"Available: {list(qc.columns)}. "
            "Re-run ebsd.pp.compute_qc_metrics(adata, target_type) to recompute."
        )
    raise KeyError(
        f"adata.obsm['{QC_METRICS_OBSM_KEY}'] not found. "
        "Run ebsd.pp.compute_qc_metrics(adata, target_type) before plotting."
    )


def _resolve_qc_values(adata: ad.AnnData, metric: str, layer: str | None) -> tuple[np.ndarray, str]:
    """Return ``(values, xlabel)`` for any supported QC metric.

    Resolution order
    ----------------
    1. ``"total_counts"`` with an explicit ``layer`` override — always sums that
       layer directly (the user wants a specific matrix, not the cached metric).
    2. ``"total_counts"`` without a layer override — prefer
       ``obsm['qc_metrics']['total_counts']`` if present (avoids touching
       lazy/zarr-backed ``X``); fall back to row sums of ``adata.X``.
    3. Any other metric found in ``adata.obsm['qc_metrics']``.
    4. Any other metric found in ``adata.obs``.
    5. Any other metric found in ``adata.obsm['cell_geometry']`` (e.g. ``"Area"``).
    6. Raise ``KeyError``.
    """
    meta = _QC_METRIC_META.get(metric)
    default_xlabel = meta[0] if meta else metric

    if metric == "total_counts":
        # Explicit layer override — user wants this specific matrix.
        if layer is not None:
            return _compute_totals(adata, layer), default_xlabel
        # Prefer pre-computed total_counts to avoid touching lazy/zarr-backed X.
        if QC_METRICS_OBSM_KEY in adata.obsm:
            qc = _obsm_to_df(adata.obsm[QC_METRICS_OBSM_KEY], adata.obs_names)
            if "total_counts" in qc.columns:
                return qc["total_counts"].to_numpy(), default_xlabel
        return _compute_totals(adata, None), default_xlabel

    if QC_METRICS_OBSM_KEY in adata.obsm:
        qc = _obsm_to_df(adata.obsm[QC_METRICS_OBSM_KEY], adata.obs_names)
        if metric in qc.columns:
            return qc[metric].to_numpy(), default_xlabel

    if metric in adata.obs.columns:
        return adata.obs[metric].to_numpy(dtype=float, na_value=np.nan), default_xlabel

    if OBSM_CELL_GEOMETRY in adata.obsm:
        cg = _obsm_to_df(adata.obsm[OBSM_CELL_GEOMETRY], adata.obs_names)
        if metric in cg.columns:
            return cg[metric].to_numpy(dtype=float, na_value=np.nan), default_xlabel

    raise KeyError(
        f"Metric '{metric}' not found in adata.obsm['{QC_METRICS_OBSM_KEY}'], adata.obs, "
        f"or adata.obsm['{OBSM_CELL_GEOMETRY}']. "
        "Run ebsd.pp.compute_qc_metrics(adata, target_type) first."
    )


# ---------------------------------------------------------------------------
# Core histogram engine
# ---------------------------------------------------------------------------


def _plot_per_well_histogram(
    values_all: np.ndarray,
    well_labels: pd.Series,
    well_order: list[str],
    n_rows: int,
    n_cols: int,
    figsize: tuple[float, float],
    *,
    xlabel: str,
    bins: int | str = 40,
    sharex: bool = True,
    sharey: bool = True,
    title: str | None = None,
    facecolor: str = "#4C78A8",
    alpha: float = 0.45,
    density: bool = False,
    logx: bool = False,
    clip_0_1: bool = False,
    show_median: bool = False,
    show_n: bool = False,
    grid_y: bool = True,
    show_percentiles: bool = False,
    lower_pct: float | None = 0.025,
    upper_pct: float | None = 0.975,
    percentile_color: str = "red",
    percentile_alpha: float = 0.7,
    percentile_ls: str = "--",
    percentile_label_fs: int = 9,
    show_threshold: bool = False,
    threshold: float | None = None,
    threshold_color: str = "red",
    threshold_alpha: float = 0.7,
    threshold_ls: str = "--",
    threshold_label_fs: int = 9,
    show_percentile: bool = False,
    max_pct: float = 0.99,
) -> plt.Figure:
    """Core histogram engine: draws a per-well grid from a flat value array."""
    well_labels_arr = well_labels.to_numpy()
    finite_all = values_all[np.isfinite(values_all)]
    if finite_all.size == 0:
        raise ValueError(f"No finite values to plot for '{xlabel}'.")

    if sharex:
        vals_per_well = []
        for w in well_order:
            v = values_all[well_labels_arr == w]
            v = v[np.isfinite(v)]
            if v.size:
                vals_per_well.append(v)
        if logx:
            bin_edges_global: np.ndarray | None = _log_bin_edges(vals_per_well, bins)
        else:
            bin_edges_global = _linear_bin_edges(finite_all, vals_per_well, bins)
    else:
        bin_edges_global = None

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=figsize,
        sharex=sharex,
        sharey=sharey,
        squeeze=False,
    )

    for idx, well in enumerate(well_order):
        r, c = divmod(idx, n_cols)
        ax = axes[r][c]

        mask = well_labels_arr == str(well)
        if not mask.any():
            ax.axis("off")
            continue

        vals = values_all[mask]
        vals = vals[np.isfinite(vals)]

        if clip_0_1:
            vals = np.clip(vals, 0.0, 1.0)

        if logx:
            pos = vals[vals > 0]
            if pos.size == 0:
                ax.axis("off")
                continue
            eps = float(pos.min()) / 10.0
            vals = np.where(vals <= 0, eps, vals)

        if vals.size == 0:
            ax.axis("off")
            continue

        bin_edges = (
            bin_edges_global if bin_edges_global is not None else _per_panel_bin_edges(vals, bins, logx)
        )

        ax.hist(
            vals,
            bins=bin_edges,
            color=facecolor,
            alpha=alpha,
            edgecolor="black",
            linewidth=0.6,
            density=density,
        )

        if show_median:
            med = float(np.median(vals))
            ax.axvline(med, color="black", linewidth=1.0, alpha=0.6)
            ytop = ax.get_ylim()[1]
            ax.text(
                med,
                ytop * 0.97,
                "median",
                color="black",
                fontsize=_fs(7),
                alpha=0.7,
                rotation=90,
                va="top",
                ha="right",
            )

        if show_percentiles and (lower_pct is not None or upper_pct is not None):
            ytop = ax.get_ylim()[1]

            def _pct_line(x: float, label: str, ha: str, *, ax=ax, ytop=ytop) -> None:
                ax.axvline(
                    x, color=percentile_color, linestyle=percentile_ls, linewidth=1.2, alpha=percentile_alpha
                )
                ax.text(
                    x,
                    ytop * 0.9,
                    label,
                    color=percentile_color,
                    fontsize=_fs(percentile_label_fs),
                    rotation=90,
                    va="top",
                    ha=ha,
                    alpha=percentile_alpha,
                )

            if lower_pct is not None:
                _pct_line(float(np.quantile(vals, lower_pct)), f"{lower_pct * 100:g}%", ha="right")
            if upper_pct is not None:
                _pct_line(float(np.quantile(vals, upper_pct)), f"{upper_pct * 100:g}%", ha="left")

        if show_threshold and threshold is not None:
            ax.axvline(
                threshold, color=threshold_color, linestyle=threshold_ls, linewidth=1.2, alpha=threshold_alpha
            )
            ytop = ax.get_ylim()[1]
            ax.text(
                threshold,
                ytop * 0.9,
                f"{threshold:.2f}",
                color=threshold_color,
                fontsize=_fs(threshold_label_fs),
                rotation=90,
                va="top",
                ha="right",
                alpha=threshold_alpha,
            )

        if show_percentile:
            thr = float(np.quantile(vals, max_pct))
            ax.axvline(
                thr, color=threshold_color, linestyle=threshold_ls, linewidth=1.2, alpha=threshold_alpha
            )
            ytop = ax.get_ylim()[1]
            ax.text(
                thr,
                ytop * 0.9,
                f"{max_pct * 100:g}%",
                color=threshold_color,
                fontsize=_fs(threshold_label_fs),
                rotation=90,
                va="top",
                ha="left",
                alpha=threshold_alpha,
            )

        if grid_y:
            ax.grid(axis="y", alpha=0.15, linestyle="-", linewidth=0.6)
        panel_title = f"Well {well}" + (f" (n={vals.size:,})" if show_n else "")
        ax.set_title(panel_title, fontsize=_fs(9), fontweight="500", color="#333")
        ax.tick_params(axis="both", labelsize=_fs(8), colors="#666")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#ddd")
        ax.spines["bottom"].set_color("#ddd")

        if logx:
            ax.set_xscale("log")

        if c == 0:
            ax.set_ylabel("Density" if density else "Cell count", fontsize=_fs(9), color="#555")

    if hasattr(fig, "supxlabel"):
        fig.supxlabel(xlabel, fontsize=_fs(9), color="#555")
    else:
        fig.text(0.5, 0.02, xlabel, ha="center", fontsize=_fs(9), color="#555")

    for idx in range(len(well_order), n_rows * n_cols):
        r, c = divmod(idx, n_cols)
        axes[r][c].axis("off")

    if title:
        fig.suptitle(title, fontsize=_fs(11), fontweight="600", color="#222")

    fig.tight_layout()
    return fig
