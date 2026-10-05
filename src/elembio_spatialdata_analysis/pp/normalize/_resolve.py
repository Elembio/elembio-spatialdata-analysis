"""Well-group resolution helpers used by normalization primitives."""

from __future__ import annotations

import numpy as np
import pandas as pd

# Candidate column names tried in priority order when run_colname is None.
# "run_name" is written by the Loader (load_tables_concatenated).
# "run_id"   is a legacy column name used in older single-run tables.
_RUN_COLNAME_CANDIDATES: tuple[str, ...] = ("run_name", "run_id")


def _obs_col_to_array(obs_df, col: str) -> np.ndarray:
    """Materialize an obs column to a numpy object array.

    With lazy loading, ``adata.obs`` is xarray-backed; column access returns
    a DataArray whose ``to_numpy()`` does not accept ``dtype`` or ``copy``.
    """
    col_data = obs_df[col]
    if isinstance(col_data, pd.Series):
        return col_data.to_numpy(dtype=object, copy=False)
    arr = col_data.to_numpy() if hasattr(col_data, "to_numpy") else np.asarray(col_data)
    return np.asarray(arr, dtype=object)


def _resolve_run_colname(obs_df: pd.DataFrame, run_colname: str | None) -> str | None:
    """Return the effective run column name to use for well-group disambiguation.

    Parameters
    ----------
    obs_df:
        ``adata.obs`` DataFrame.
    run_colname:
        Explicit column name, or ``None`` to auto-detect.

    Returns
    -------
    The first matching column name found in *obs_df*, or ``None`` when no
    run column is present (single-run data — no disambiguation needed).

    Behaviour
    ---------
    * ``run_colname=None`` (default): auto-detect — tries ``"run_name"``
      then ``"run_id"`` in order.  If neither exists the function returns
      ``None`` and well groups are formed from the well label alone.
    * ``run_colname="some_col"``: use that column exactly.  Raises
      ``KeyError`` if the column is absent, since the caller explicitly
      requested it.
    """
    if run_colname is not None:
        if run_colname not in obs_df.columns:
            raise KeyError(
                f"run_colname='{run_colname}' not found in AnnData.obs. "
                f"Available columns: {list(obs_df.columns)}"
            )
        return run_colname

    # Auto-detect: return the first candidate present
    for candidate in _RUN_COLNAME_CANDIDATES:
        if candidate in obs_df.columns:
            return candidate
    return None


def _resolve_well_groups(
    obs_df: pd.DataFrame,
    well_colname: str,
    run_colname: str | None = None,
) -> tuple[np.ndarray, list[np.ndarray]]:
    """Build well-group labels and row-index arrays.

    When a run column is present (or auto-detected), well groups are
    partitioned by the ``(run, well)`` combination so that repeated well
    identifiers across runs are treated independently.

    Parameters
    ----------
    obs_df:
        ``adata.obs`` DataFrame.
    well_colname:
        Column that identifies individual wells.
    run_colname:
        Column that identifies the run/experiment.

        * ``None`` (default): auto-detect — tries ``"run_name"`` (written
          by ``Loader.load(..., concatenated=True)``) then ``"run_id"``
          (legacy).  If neither exists, wells are grouped by label alone
          (correct for single-run data).
        * Any string: use that column exactly; raises ``KeyError`` if absent.

    Returns
    -------
    wells:
        Unique group labels — plain well labels for single-run data, or
        ``(run, well)`` tuples for multi-run data.
    well_indices:
        List of integer row-index arrays, one per group.
    """
    if well_colname not in obs_df.columns:
        raise KeyError(f"Column '{well_colname}' not found in AnnData.obs")

    effective_run_col = _resolve_run_colname(obs_df, run_colname)

    well_values = _obs_col_to_array(obs_df, well_colname)
    well_values = np.array([np.nan if pd.isna(val) else val for val in well_values], dtype=object)

    if effective_run_col is not None:
        run_values = _obs_col_to_array(obs_df, effective_run_col)
        run_values = np.array([None if pd.isna(val) else val for val in run_values], dtype=object)
        labels = np.empty(well_values.shape[0], dtype=object)
        for idx, (run_val, well_val) in enumerate(zip(run_values, well_values, strict=False)):
            labels[idx] = (run_val, well_val)
    else:
        labels = well_values

    label_series = pd.Series(labels, copy=False)
    wells = label_series.unique()
    well_indices = [label_series.index[label_series == well].to_numpy(dtype=int) for well in wells]
    return wells, well_indices
