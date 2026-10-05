"""Quantile resolution helpers for per-well and per-batch threshold configuration."""

from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

QuantileSpec = float | Mapping[str, "float | Mapping[str, float]"] | None


def _resolve_quantile_overrides(
    value: QuantileSpec,
    *,
    spec_name: str,
) -> dict[str, float | None]:
    """Coerce a scalar or mapping into a normalized override dictionary."""

    def _validate(val: float | None, key: str) -> float | None:
        if val is None:
            return None
        if not 0.0 <= val <= 1.0:
            raise ValueError(f"{spec_name} for '{key}' must be within [0, 1]; got {val}.")
        return float(val)

    if isinstance(value, Mapping):
        overrides: dict[str, float] = {}
        default = value.get("__default__", value.get("default"))
        if default is not None:
            overrides["__default__"] = _validate(float(default), f"{spec_name}.__default__")
        for key, val in value.items():
            if key in {"__default__", "default"}:
                continue
            if isinstance(val, Mapping):
                raise TypeError(
                    f"{spec_name} does not support nested mappings beyond one level (offending key '{key}')."
                )
            overrides[str(key)] = _validate(float(val), f"{spec_name}[{key}]")
        return overrides

    return {"__default__": _validate(float(value), spec_name)}


def _lookup_override(
    overrides: Mapping[str, float | None] | None,
    column: str,
) -> float | None:
    if not overrides:
        return None
    if column in overrides:
        return overrides[column]
    batch_key = column.removeprefix("total_counts_")
    if batch_key in overrides:
        return overrides[batch_key]
    return None


def _resolve_well_batch_quantiles(
    quantile_spec: QuantileSpec,
    *,
    count_cols: list[str],
    well_labels: pd.Series,
    default_value: float | None,
    spec_name: str,
) -> dict[str, dict[str, float | None]]:
    """
    Resolve quantile configuration into explicit well→batch thresholds (0-1 values).

    Precedence for a given (well, batch) pair:
        1. Well-specific overrides (exact column, batch name, then well default)
        2. Global batch overrides (exact column, batch name, then batch default)
        3. Top-level default (either from ``min_quant_count`` scalar or ``__default__``)
    """

    def _validate_default(val: float | None, key: str) -> float | None:
        if val is None:
            return None
        if not 0.0 <= val <= 1.0:
            raise ValueError(f"{spec_name} default for '{key}' must be within [0, 1]; got {val}.")
        return float(val)

    wells = pd.Index(well_labels.astype(str)).unique()
    default = _validate_default(default_value, "__default__")
    batch_overrides: dict[str, float | None] = {}
    well_overrides: dict[str, dict[str, float | None]] = {}

    if isinstance(quantile_spec, Mapping):
        for key, val in quantile_spec.items():
            key_str = str(key)
            if key_str in {"__default__", "default"}:
                default = _validate_default(float(val), key_str)
            elif key_str == "__batch__":
                batch_overrides.update(_resolve_quantile_overrides(val, spec_name=f"{spec_name}[__batch__]"))
            elif isinstance(val, Mapping):
                well_overrides[key_str] = _resolve_quantile_overrides(
                    val, spec_name=f"{spec_name}[{key_str}]"
                )
            else:
                target_is_well = key_str in wells
                overrides = _resolve_quantile_overrides(val, spec_name=f"{spec_name}[{key_str}]")
                if target_is_well:
                    well_overrides[key_str] = overrides
                else:
                    batch_overrides.update(overrides)
    elif quantile_spec is not None and not isinstance(quantile_spec, (int, float)):
        raise TypeError(f"{spec_name} must be a float or mapping; received type {type(quantile_spec)}.")

    well_quantiles: dict[str, dict[str, float | None]] = {}

    for well in wells:
        per_well: dict[str, float] = {}
        well_override = well_overrides.get(well, {})
        well_default = well_override.get("__default__")
        for col in count_cols:
            val = _lookup_override(well_override, col)
            if val is None and well_default is not None:
                val = well_default
            if val is None:
                val = _lookup_override(batch_overrides, col)
            if val is None and batch_overrides.get("__default__") is not None:
                val = batch_overrides["__default__"]
            if val is None:
                val = default
            per_well[col] = val
        well_quantiles[well] = per_well

    return well_quantiles


def _resolve_threshold_defaults(
    normalized_target_type: str,
    overrides: dict,
) -> dict:
    """Merge modality defaults with explicit overrides for filtering thresholds."""
    from elembio_spatialdata_analysis.pp.filter._constants import (
        _FILTER_PARAM_DEFAULTS,
        _UNSET,
        DEFAULT_BOUNDS_PER_MODALITY,
    )

    modality_defaults = DEFAULT_BOUNDS_PER_MODALITY.get(normalized_target_type, {})
    resolved: dict = {}

    for key, base_value in _FILTER_PARAM_DEFAULTS.items():
        override_value = overrides.get(key, _UNSET)
        if override_value is not _UNSET:
            resolved[key] = override_value
            continue
        if key in modality_defaults:
            resolved[key] = modality_defaults[key]
        else:
            resolved[key] = base_value

    return resolved
