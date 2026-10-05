"""Multi-modality filtering: filter_modalities."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import anndata as ad
import pandas as pd

from elembio_spatialdata_analysis.pp.filter._spatialdata import _update_spatialdata_attrs
from elembio_spatialdata_analysis.pp.filter.filter_cells import (
    ModalityFilterResult,
    filter_cells,
)
from elembio_spatialdata_analysis.tl.constants import normalize_target_type

# Modalities that own cell-segmentation area QC in a multimodal run.
# All other modalities share the same segmented cells and should not
# re-apply area filtering — Transcript already enforces area quality.
_AREA_OWNER_MODALITIES = frozenset({"transcript", "threeprimeuntargeted"})


def filter_modalities(
    adata_by_modality: Mapping[str, ad.AnnData],
    *,
    modality_order: Sequence[str] | None = None,
    modality_kwargs: Mapping[str, Mapping[str, Any]] | None = None,
    return_pass_summary: bool = True,
    disable_area_for_non_owners: bool = True,
    **common_kwargs: Any,
) -> tuple[dict[str, ModalityFilterResult], pd.DataFrame | None]:
    """Filter multiple modalities and collect their diagnostics.

    Parameters
    ----------
    adata_by_modality : Mapping[str, AnnData]
        Mapping from modality label to the corresponding AnnData object.
    modality_order : sequence of str, optional
        Override the iteration order of modalities.
    modality_kwargs : Mapping[str, Mapping[str, Any]], optional
        Per-modality overrides merged with ``common_kwargs`` before invoking
        ``filter_cells``.
    return_pass_summary : bool, default True
        If True, return a pass matrix describing how each cell fares across modalities.
    disable_area_for_non_owners : bool, default True
        When True, area filtering (``min_quant_area`` / ``max_quant_area``) is
        automatically disabled for modalities that do not own cell-segmentation
        quality (i.e. everything except Transcript / ThreePrimeUntargeted).
        Set to False to apply uniform area defaults across all modalities.
    **common_kwargs : Any
        Additional keyword arguments forwarded to each ``filter_cells`` call.

    Returns
    -------
    (results_by_modality, pass_summary)
        ``results_by_modality`` maps the modality label to ``ModalityFilterResult``.
        ``pass_summary`` is a DataFrame with columns ``passes_<modality>`` plus summary
        statistics; ``None`` when ``return_pass_summary=False``.
    """
    if not adata_by_modality:
        raise ValueError("`adata_by_modality` must contain at least one modality.")

    if modality_order is None:
        modality_order = list(adata_by_modality.keys())
    else:
        missing = [m for m in modality_order if m not in adata_by_modality]
        if missing:
            raise KeyError(f"Modalities {', '.join(missing)} are missing from `adata_by_modality`.")

    results: dict[str, ModalityFilterResult] = {}
    pass_columns: dict[str, pd.Series] = {}

    for modality in modality_order:
        adata_obj = adata_by_modality[modality]
        kwargs = dict(common_kwargs)

        # In a multimodal run, only the area-owning modalities should enforce
        # cell-segmentation area bounds.  Non-owning modalities share the same
        # segmented cells; their area quality is guaranteed by the owner's filter.
        if disable_area_for_non_owners:
            normalized_mod = normalize_target_type(modality)
            if normalized_mod not in _AREA_OWNER_MODALITIES:
                kwargs.setdefault("min_quant_area", None)
                kwargs.setdefault("max_quant_area", None)

        if modality_kwargs and modality in modality_kwargs:
            kwargs.update(modality_kwargs[modality])

        filtered, diag = filter_cells(adata_obj, modality, diagnostics=True, **kwargs)
        results[modality] = ModalityFilterResult(adata=filtered, diagnostics=diag)

        normalized = diag.thresholds.get("normalized_target_type", normalize_target_type(modality))
        pass_col = f"passes_{normalized}"
        passes = diag.flags["passes_all_filters"].rename(pass_col)
        pass_columns[pass_col] = passes

    pass_summary = pd.DataFrame(pass_columns).fillna(False).infer_objects(copy=False)
    pass_cols = list(pass_columns.keys())

    pass_summary["n_modalities_passed"] = pass_summary[pass_cols].sum(axis=1)
    pass_summary["passes_all_modalities"] = pass_summary[pass_cols].all(axis=1)
    pass_summary["passes_any_modalities"] = pass_summary[pass_cols].any(axis=1)

    passes_all_mask = pass_summary["passes_all_modalities"]
    cells_intersection = pass_summary.index[passes_all_mask]

    for result in results.values():
        adata_mod = result.adata
        mask = adata_mod.obs_names.isin(cells_intersection)
        result.adata = adata_mod[mask].copy()
        _update_spatialdata_attrs(result.adata)

    if not return_pass_summary:
        return results, None

    return results, pass_summary
