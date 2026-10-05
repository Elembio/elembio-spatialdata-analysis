"""Merge and concatenation utilities for AnnData tables and Polars frames.

This module is backend-agnostic and reusable by both :mod:`tables` and
:mod:`points` modules.
"""

from __future__ import annotations

import re
import warnings
from collections.abc import MutableMapping

import anndata as ad
import numpy as np
import pandas as pd
import polars as pl

# ---------------------------------------------------------------------------
# Column mismatch diagnostics
# ---------------------------------------------------------------------------


def _check_var_mismatch(tables: dict[str, ad.AnnData]) -> bool:
    """Check whether ``.var_names`` are identical across all tables.

    Returns ``True`` if all tables share the same var_names, ``False``
    otherwise.  On mismatch, emits a warning listing the differences
    and noting that ``.var`` metadata will be cleared.
    """
    if len(tables) < 2:
        return True

    keys = list(tables.keys())
    ref_key = keys[0]
    ref_vars = set(tables[ref_key].var_names)

    for key in keys[1:]:
        other_vars = set(tables[key].var_names)
        only_ref = ref_vars - other_vars
        only_other = other_vars - ref_vars
        if only_ref or only_other:
            parts: list[str] = []
            if only_ref:
                sample = sorted(only_ref)[:5]
                parts.append(f"only in {ref_key!r}: {sample}" + (" ..." if len(only_ref) > 5 else ""))
            if only_other:
                sample = sorted(only_other)[:5]
                parts.append(f"only in {key!r}: {sample}" + (" ..." if len(only_other) > 5 else ""))
            warnings.warn(
                f"var_names differ across runs ({'; '.join(parts)}). "
                ".var metadata (target_id, batch_name, etc.) will be cleared "
                "because it no longer matches across runs.",
                UserWarning,
                stacklevel=3,
            )
            return False

    return True


# ---------------------------------------------------------------------------
# obsm / uns merge diagnostics
# ---------------------------------------------------------------------------


def _check_obsm_mismatch(
    tables: dict[str, ad.AnnData],
) -> list[str]:
    """Warn about differing ``.obsm`` keys or column schemas across tables."""
    if len(tables) < 2:
        return []

    keys = list(tables.keys())
    ref_key = keys[0]
    ref_obsm_keys = set(tables[ref_key].obsm.keys())
    issues: list[str] = []

    for key in keys[1:]:
        other_keys = set(tables[key].obsm.keys())
        if ref_obsm_keys != other_keys:
            issues.append(
                f"obsm keys differ between {ref_key!r} ({sorted(ref_obsm_keys)}) "
                f"and {key!r} ({sorted(other_keys)})"
            )
            continue

        for obsm_key in ref_obsm_keys:
            ref_val = tables[ref_key].obsm[obsm_key]
            other_val = tables[key].obsm[obsm_key]
            if (
                isinstance(ref_val, pd.DataFrame)
                and isinstance(other_val, pd.DataFrame)
                and list(ref_val.columns) != list(other_val.columns)
            ):
                issues.append(f"obsm[{obsm_key!r}] columns differ between {ref_key!r} and {key!r}")

    for msg in issues:
        warnings.warn(msg, stacklevel=3)
    return issues


# ---------------------------------------------------------------------------
# Row-wise (axis=0) concatenation of AnnData tables
# ---------------------------------------------------------------------------


def _downcast_stringdtype(df: pd.DataFrame) -> None:
    """Convert ``StringDType`` columns to ``object`` in-place.

    Numpy's ``StringDType`` cannot be promoted with float/int dtypes,
    which causes ``ad.concat`` to fail when the same column has different
    dtypes across runs.  Casting to ``object`` lets pandas handle the
    mixed-type concatenation gracefully.
    """
    for col in df.columns:
        if isinstance(df[col].dtype, np.dtypes.StringDType):
            df[col] = df[col].astype(object)


def _concat_tables_rowwise(
    tables: dict[str, ad.AnnData],
    *,
    source_stores: dict[str, str] | None = None,
) -> ad.AnnData:
    """Concatenate AnnData objects across runs (row-wise, axis=0).

    Merge strategy:

    - ``.var_names`` are outer-joined (union of all targets).
    - ``.var`` columns are kept only when var_names are identical across
      runs; otherwise cleared (values like target_id, batch_name no longer
      match).
    - ``.varm`` DataFrames are column-stacked with ``{run_name}_{col}``
      prefixes when var_names are identical; otherwise dropped.
    - ``.obsm`` DataFrames are outer-merged (missing columns filled with
      NaN).
    - ``.uns`` entries are kept only if identical across all runs
      (``"same"`` strategy); differing entries are dropped.

    Adds provenance columns to ``.obs``:

    - ``source_key``: the key in *tables*
    - ``run_name``: *source_key* with any ``__N`` dedup suffix stripped
    - ``source_store``: the original path/URI (from *source_stores* dict)

    Parameters
    ----------
    tables
        ``{source_key: adata}``.
    source_stores
        ``{source_key: path_or_uri}``.  Optional; used to populate
        ``source_store`` in ``.obs``.
    """
    if not tables:
        raise ValueError("No tables provided for concatenation.")

    var_identical = _check_var_mismatch(tables)
    _check_obsm_mismatch(tables)

    # Normalize StringDType columns to object to avoid numpy DType
    # promotion errors when concatenating across runs with mixed dtypes.
    for adata in tables.values():
        _downcast_stringdtype(adata.obs)
        _downcast_stringdtype(adata.var)

    # Pre-extract varm before ad.concat (which uses merge="first" and
    # silently drops all but the first run's varm).
    merged_varm = _merge_varm(tables) if var_identical else {}

    combined = ad.concat(
        tables,
        axis=0,
        join="outer",
        merge="first",
        uns_merge="same",
        label="source_key",
        index_unique="-",
    )

    if not var_identical:
        combined.var = combined.var.drop(columns=combined.var.columns)

    # Drop uns entries that differ across runs.
    _drop_differing_uns(combined.uns, tables)

    # Replace varm with the properly merged version.
    for key in list(combined.varm.keys()):
        del combined.varm[key]
    for key, df in merged_varm.items():
        combined.varm[key] = df

    source_keys = combined.obs["source_key"].to_numpy().astype(str)
    combined.obs["run_name"] = [re.sub(r"__\d+$", "", sk) for sk in source_keys]

    if source_stores:
        combined.obs["source_store"] = [source_stores.get(sk, sk) for sk in source_keys]

    return combined


def _drop_differing_uns(
    combined_uns: MutableMapping,
    tables: dict[str, ad.AnnData],
) -> None:
    """Remove ``.uns`` entries that are not identical across all runs."""
    if len(tables) < 2:
        return
    adatas = list(tables.values())
    for key in list(combined_uns.keys()):
        values = [a.uns[key] for a in adatas if key in a.uns]
        if len(values) < len(adatas):
            del combined_uns[key]
            continue
        ref = values[0]
        if not all(_uns_values_equal(ref, v) for v in values[1:]):
            del combined_uns[key]


def _uns_values_equal(a: object, b: object) -> bool:
    """Recursively compare two uns values for equality."""
    if type(a) is not type(b):
        return False
    if isinstance(a, dict):
        if a.keys() != b.keys():  # type: ignore[union-attr]
            return False
        return all(_uns_values_equal(a[k], b[k]) for k in a)  # type: ignore[index]
    if isinstance(a, pd.DataFrame):
        return a.shape == b.shape and (a == b).all().all()  # type: ignore[union-attr]
    if isinstance(a, np.ndarray):
        return np.array_equal(a, b)
    try:
        return bool(a == b)
    except (ValueError, TypeError):
        return False


def _to_dataframe(v: object) -> pd.DataFrame | None:
    """Convert a varm value to DataFrame, handling Dataset2D from lazy reads."""
    if isinstance(v, pd.DataFrame):
        return v
    if hasattr(v, "to_memory"):
        materialized = v.to_memory()
        if isinstance(materialized, pd.DataFrame):
            return materialized
    return None


def _merge_varm(tables: dict[str, ad.AnnData]) -> dict[str, pd.DataFrame]:
    """Merge ``.varm`` DataFrames across runs, prefixing columns with run name.

    Each column ``{well_id}`` becomes ``{run_name}_{well_id}``.
    Only DataFrame-compatible varm entries present in all tables are merged.
    """
    if not tables:
        return {}

    keys = list(tables.keys())
    run_names = {key: re.sub(r"__\d+$", "", key) for key in keys}

    # Find varm keys common to all tables.
    common_varm_keys: set[str] | None = None
    for adata in tables.values():
        varm_keys = {k for k in adata.varm if _to_dataframe(adata.varm[k]) is not None}
        common_varm_keys = varm_keys if common_varm_keys is None else common_varm_keys & varm_keys
    if not common_varm_keys:
        return {}

    result: dict[str, pd.DataFrame] = {}
    for varm_key in sorted(common_varm_keys):
        frames: list[pd.DataFrame] = []
        for source_key in keys:
            df = _to_dataframe(tables[source_key].varm[varm_key])
            assert df is not None
            prefix = run_names[source_key]
            renamed = df.rename(columns={col: f"{prefix}_{col}" for col in df.columns})
            frames.append(renamed)
        result[varm_key] = pd.concat(frames, axis=1)

    return result


# ---------------------------------------------------------------------------
# Polars frame concatenation (for points)
# ---------------------------------------------------------------------------


def _concat_frames(
    frames: dict[str, pl.LazyFrame],
) -> pl.LazyFrame:
    """Concatenate Polars LazyFrames with provenance tracking.

    Adds a ``source_key`` column before concatenating so each row can be
    traced back to its origin.

    Parameters
    ----------
    frames
        ``{source_key: LazyFrame}``.
    """
    if not frames:
        raise ValueError("No frames provided for concatenation.")

    tagged: list[pl.LazyFrame] = []
    for key, lf in frames.items():
        tagged.append(lf.with_columns(pl.lit(key).alias("source_key")))

    return pl.concat(tagged, how="diagonal_relaxed")
