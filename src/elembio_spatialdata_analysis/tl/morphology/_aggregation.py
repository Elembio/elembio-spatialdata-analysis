"""CellPaint morphology feature aggregation.

Reduces dimensionality of CellPaint features by collapsing groups of related
measurements into a single representative value.  The default operation is to
aggregate **texture** features only (4 orientations → 1 rotation-invariant
value).  Radial and granularity features are kept separate by default because
they encode spatial gradients and scale-specific information respectively.
"""

from __future__ import annotations

from typing import Literal

import anndata as ad
import numpy as np
import pandas as pd

from elembio_spatialdata_analysis.tl.morphology._feature_groups import get_morphology_feature_groups


def aggregate_morphology_features(
    adata: ad.AnnData,
    aggregation_method: Literal["mean", "median", "max", "min", "std"] = "mean",
    feature_types: list[str] | None = None,
) -> ad.AnnData:
    """Aggregate related morphology features to reduce dimensionality.

    By default only aggregates **texture** features (4 orientations → 1 value)
    to create rotation-invariant features.  Radial and granularity features are
    kept in their original form.

    Parameters
    ----------
    adata:
        AnnData object with CellPaint morphology features in ``X`` / ``var``.
    aggregation_method:
        Reduction applied within each feature group.
    feature_types:
        Feature types to aggregate.  Defaults to ``["texture"]``.
        Options: ``"texture"``, ``"radial"``, ``"granularity"``.

        .. warning::
            Aggregating ``"radial"`` or ``"granularity"`` features may destroy
            spatial / scale information.  Use with caution.

    Returns
    -------
    New AnnData with aggregated features.  ``obs`` is preserved unchanged.
    """
    if feature_types is None:
        feature_types = ["texture"]

    feature_groups = get_morphology_feature_groups(adata)
    texture_groups = feature_groups["texture_groups"]
    radial_groups = feature_groups["radial_groups"]
    granularity_groups = feature_groups["granularity_groups"]

    X_orig = adata.X
    if hasattr(X_orig, "toarray"):
        X_orig = X_orig.toarray()

    feature_names_list = adata.var_names.tolist()
    feature_idx_map = {name: idx for idx, name in enumerate(feature_names_list)}

    agg_funcs = {
        "mean": np.mean,
        "median": np.median,
        "max": np.max,
        "min": np.min,
        "std": np.std,
    }
    if aggregation_method not in agg_funcs:
        raise ValueError(
            f"Unknown aggregation method: {aggregation_method!r}. Choose from: {list(agg_funcs)}"
        )
    agg_func = agg_funcs[aggregation_method]

    aggregated_cols: list[np.ndarray] = []
    aggregated_names: list[str] = []
    processed: set[str] = set()

    def _add_group(key: str, feat_list: list[str]) -> None:
        processed.update(feat_list)
        if len(feat_list) > 1:
            cols = [feature_idx_map[f] for f in feat_list]
            aggregated_cols.append(agg_func(X_orig[:, cols], axis=1))
            aggregated_names.append(key)
        else:
            aggregated_cols.append(X_orig[:, feature_idx_map[feat_list[0]]])
            aggregated_names.append(feat_list[0])

    if "texture" in feature_types:
        for key, feat_list in sorted(texture_groups.items()):
            _add_group(key, feat_list)

    if "radial" in feature_types:
        zone_based = ("FracAtD", "MeanFrac", "RadialCV")
        for key, feat_list in sorted(radial_groups.items()):
            processed.update(feat_list)
            if any(z in key for z in zone_based):
                for feat in feat_list:
                    aggregated_cols.append(X_orig[:, feature_idx_map[feat]])
                    aggregated_names.append(feat)
            elif len(feat_list) > 1:
                cols = [feature_idx_map[f] for f in feat_list]
                aggregated_cols.append(agg_func(X_orig[:, cols], axis=1))
                aggregated_names.append(key)
            else:
                aggregated_cols.append(X_orig[:, feature_idx_map[feat_list[0]]])
                aggregated_names.append(feat_list[0])

    if "granularity" in feature_types:
        for key, feat_list in sorted(granularity_groups.items()):
            processed.update(feat_list)
            scale_1 = [f for f in feat_list if "_1_" in f or f.endswith("_1")]
            other_scales = [f for f in feat_list if f not in scale_1]
            for feat in scale_1:
                aggregated_cols.append(X_orig[:, feature_idx_map[feat]])
                aggregated_names.append(feat)
            if len(other_scales) > 1:
                cols = [feature_idx_map[f] for f in other_scales]
                aggregated_cols.append(agg_func(X_orig[:, cols], axis=1))
                aggregated_names.append(f"{key}_scales2to8")
            elif len(other_scales) == 1:
                aggregated_cols.append(X_orig[:, feature_idx_map[other_scales[0]]])
                aggregated_names.append(other_scales[0])

    for feat in feature_names_list:
        if feat not in processed:
            aggregated_cols.append(X_orig[:, feature_idx_map[feat]])
            aggregated_names.append(feat)

    if not aggregated_cols:
        raise ValueError("No features remaining after aggregation. Check feature_types and input data.")

    X_agg = np.column_stack([c if c.ndim == 2 else c.reshape(-1, 1) for c in aggregated_cols])

    return ad.AnnData(
        X=X_agg,
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=aggregated_names),
    )
