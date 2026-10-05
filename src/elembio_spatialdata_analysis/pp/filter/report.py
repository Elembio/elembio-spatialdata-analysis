"""Filter report generation: generate_filter_report."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import anndata as ad
import numpy as np
import pandas as pd

from elembio_spatialdata_analysis.pp.filter._constants import (
    FILTER_FLAGS_FULL_UNS_KEY,
    FILTER_FLAGS_OBSM_KEY,
    FILTER_PARAMS_UNS_KEY,
)
from elembio_spatialdata_analysis.pp.filter.filter_cells import FilterDiagnostics

if TYPE_CHECKING:
    pass


def generate_filter_report(
    filter_results: ad.AnnData | FilterDiagnostics | pd.DataFrame,
    *,
    filter_order: list[str] | None = None,
    include_total: bool = True,
    exclude_filters: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Summarize filtering outcomes, highlighting unique vs shared failures per filter.

    Parameters
    ----------
    filter_results : AnnData, FilterDiagnostics, or pandas.DataFrame
        - ``AnnData``: reads ``uns['filter_flags_full']`` (written by ``filter_cells``
          before inplace subsetting) with fallback to ``obsm['filter_flags']``, plus
          ``uns['filter_params']`` for threshold metadata.
        - ``FilterDiagnostics``: uses the ``flags`` and ``thresholds`` fields directly.
        - ``pandas.DataFrame``: raw boolean filter-flags dataframe (no threshold info).
    filter_order : list[str], optional
        Preferred ordering of filter columns.
    include_total : bool, default True
        If True, append a ``TOTAL`` row summarizing failures across all filters.
    exclude_filters : sequence of str, optional
        Filter columns to skip. Defaults to excluding ``passes_all_modalities`` and
        ``passes_any_modalities`` if present.

    Returns
    -------
    pandas.DataFrame
        Table with per-filter pass/fail counts and threshold metadata.
    """
    if isinstance(filter_results, ad.AnnData):
        adata = filter_results
        # Prefer the full pre-subset flags stored in uns (written by filter_cells before
        # _inplace_subset_obs), which contains all cells including those that were removed.
        # Fall back to obsm['filter_flags'] for backward compat or compute_filter_diagnostics usage.
        if FILTER_FLAGS_FULL_UNS_KEY in adata.uns:
            filter_flags = adata.uns[FILTER_FLAGS_FULL_UNS_KEY]
        elif FILTER_FLAGS_OBSM_KEY in adata.obsm:
            filter_flags = adata.obsm[FILTER_FLAGS_OBSM_KEY]
        else:
            raise KeyError("No filter flags found in adata. Run ebsd.pp.filter_cells() first.")
        thresholds = adata.uns.get(FILTER_PARAMS_UNS_KEY)
    elif isinstance(filter_results, FilterDiagnostics):
        filter_flags = filter_results.flags
        thresholds = filter_results.thresholds
    else:
        filter_flags = filter_results
        thresholds = None

    default_order = [
        "passes_nan_filter",
        "passes_empty_filter",
        "passes_counts",
        "passes_counts_per_batch",  # legacy compat (0.11.0)
        "passes_imaging_intensity",
        "passes_imaging_counts",  # legacy compat (0.10.0)
        "passes_cell_paint_counts",  # legacy compat
        "passes_assignment_rate",
        "passes_area",
        "passes_nsb_outlier",
        "passes_gapdh_outlier",
        "passes_mt_outlier",
        "passes_spatial_artifacts",
    ]
    preferred_order = filter_order or default_order

    exclusion_set = set(exclude_filters or ("passes_all_modalities", "passes_any_modalities"))
    exclusion_set.add("passes_all_filters")

    # Automatically hide filters whose threshold(s) are all None — they were
    # disabled and would show 0 failures, cluttering the report.
    _disabled_threshold_map: dict[str, tuple[str, ...]] = {
        "passes_assignment_rate": ("min_assignment_rate",),
        "passes_nsb_outlier": ("max_quant_nsb",),
        "passes_gapdh_outlier": ("max_quant_gapdh",),
        "passes_mt_outlier": ("max_pct_mt",),
        "passes_counts": ("min_quant_count", "max_quant_count"),
        "passes_counts_per_batch": ("min_quant_count", "max_quant_count"),  # legacy compat (0.11.0)
        "passes_area": ("min_quant_area", "max_quant_area"),
        "passes_spatial_artifacts": ("filter_spatial_artifacts",),
    }
    if thresholds is not None:
        for flag_col, param_keys in _disabled_threshold_map.items():
            # None disables quantile/rate gates; False disables the boolean
            # spatial-artifact gate (default off).
            if all(thresholds.get(k) is None or thresholds.get(k) is False for k in param_keys):
                exclusion_set.add(flag_col)

    filter_cols = [col for col in preferred_order if col not in exclusion_set and col in filter_flags.columns]
    extra_filters = [
        col
        for col in filter_flags.columns
        if col.startswith("passes_") and col not in filter_cols and col not in exclusion_set
    ]
    filter_cols.extend(extra_filters)

    if not filter_cols:
        raise ValueError("No filter columns found in `filter_flags`.")

    filter_matrix = filter_flags[filter_cols].fillna(False).astype(bool)
    fail_matrix = ~filter_matrix

    all_filters_mask = filter_flags.get("passes_all_filters")
    expected_all = filter_matrix.all(axis=1)
    if all_filters_mask is not None and not all_filters_mask.fillna(False).equals(expected_all):
        raise ValueError("`passes_all_filters` column is inconsistent with filter flags.")

    total_cells = len(filter_flags)

    filter_threshold_map: dict[str, tuple[str, ...]] = {
        "passes_counts": ("min_quant_count", "max_quant_count"),
        "passes_counts_per_batch": ("min_quant_count", "max_quant_count"),  # legacy compat (0.11.0)
        "passes_assignment_rate": ("min_assignment_rate",),
        "passes_area": ("min_quant_area", "max_quant_area"),
        "passes_nsb_outlier": ("max_quant_nsb",),
        "passes_gapdh_outlier": ("max_quant_gapdh",),
        "passes_mt_outlier": ("max_pct_mt",),
        "passes_nan_filter": (),
        "passes_empty_filter": ("min_total_count",),
        "passes_spatial_artifacts": ("filter_spatial_artifacts",),
    }

    def _format_threshold_summary(filter_name: str) -> str | None:
        if thresholds is None:
            return None
        keys = filter_threshold_map.get(filter_name)
        if not keys:
            return None
        pairs = []
        for key in keys:
            if key in thresholds:
                pairs.append(f"{key}={thresholds[key]!r}")
        if not pairs:
            return None
        return ", ".join(pairs)

    report_rows: list[dict[str, object]] = []

    for col in filter_cols:
        failures = fail_matrix[col]
        fail_count = int(failures.sum())

        if len(filter_cols) == 1:
            other_failures = pd.Series(False, index=failures.index)
        else:
            other_failures = fail_matrix.drop(columns=col).any(axis=1)

        unique_failures = failures & (~other_failures)
        overlap_failures = failures & other_failures

        unique_count = int(unique_failures.sum())
        overlap_count = int(overlap_failures.sum())
        pass_count = total_cells - fail_count

        report_rows.append(
            {
                "filter": col,
                "cells_passed": pass_count,
                "cells_failed": fail_count,
                "cells_failed_unique": unique_count,
                "cells_failed_overlap": overlap_count,
                "pct_failed_total": fail_count / total_cells if total_cells else np.nan,
                "pct_failed_unique": unique_count / total_cells if total_cells else np.nan,
                "pct_failed_overlap": overlap_count / total_cells if total_cells else np.nan,
                "pct_passed_total": pass_count / total_cells if total_cells else np.nan,
                "thresholds": _format_threshold_summary(col),
            }
        )

    report = pd.DataFrame(report_rows)

    if include_total:
        fail_counts_per_cell = fail_matrix.sum(axis=1)
        total_fail = int((fail_counts_per_cell >= 1).sum())
        total_unique = int((fail_counts_per_cell == 1).sum())
        total_overlap = int((fail_counts_per_cell >= 2).sum())
        total_pass = total_cells - total_fail

        total_row = {
            "filter": "TOTAL",
            "cells_passed": total_pass,
            "cells_failed": total_fail,
            "cells_failed_unique": total_unique,
            "cells_failed_overlap": total_overlap,
            "pct_failed_total": total_fail / total_cells if total_cells else np.nan,
            "pct_failed_unique": total_unique / total_cells if total_cells else np.nan,
            "pct_failed_overlap": total_overlap / total_cells if total_cells else np.nan,
            "pct_passed_total": total_pass / total_cells if total_cells else np.nan,
            "thresholds": None,
        }
        report = pd.concat([report, pd.DataFrame([total_row])], ignore_index=True)

    return report
