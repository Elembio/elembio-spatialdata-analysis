"""YAML-driven traffic-light evaluation for QC metrics.

Each gate maps a metric name to threshold rules (``warn_below`` / ``fail_below``
/ ``warn_above`` / ``fail_above``) and produces a ``GateResult`` with a
severity (``pass`` / ``warn`` / ``fail`` / ``unknown``) for any observed value.

The default gates config ships with the package.  Additional or stricter
thresholds can be layered on top with :meth:`QcGates.overlay`.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import yaml

_logger = logging.getLogger(__name__)

Severity = Literal["pass", "warn", "fail", "unknown"]

_DEFAULT_GATES_PATH = Path(__file__).resolve().parent / "config" / "qc_gates.yml"


@dataclass
class GateResult:
    """Outcome of evaluating a single metric against its gate."""

    metric: str
    value: float | int | None
    severity: Severity
    threshold: dict[str, float]
    message: str | None = None


class QcGates:
    """In-memory store of metric → threshold rules."""

    def __init__(self, rules: dict[str, dict[str, Any]] | None = None) -> None:
        self._rules: dict[str, dict[str, Any]] = rules or {}

    @classmethod
    def from_yaml(cls, path: Path | None = None) -> QcGates:
        """Load gates from a YAML file. Falls back to the bundled default."""
        target = path if path is not None else _DEFAULT_GATES_PATH
        if not target.is_file():
            _logger.debug("QC gates file %s not found; returning empty rule set", target)
            return cls({})
        try:
            with target.open("r", encoding="utf-8") as f:
                data = yaml.safe_load(f) or {}
        except Exception as exc:
            _logger.warning("Failed to parse QC gates file %s: %s", target, exc)
            return cls({})
        if not isinstance(data, dict):
            return cls({})
        return cls({str(k): v for k, v in data.items() if isinstance(v, dict)})

    def overlay(self, other: QcGates) -> QcGates:
        """Return a new :class:`QcGates` where ``other`` overrides matching rules."""
        merged = dict(self._rules)
        merged.update(other._rules)
        return QcGates(merged)

    def evaluate(self, metric: str, value: float | int | None) -> GateResult:
        """Evaluate one metric value against its rule.

        Severity ladder: any ``fail_*`` violation → ``fail``; any ``warn_*``
        violation → ``warn``; otherwise ``pass`` if a rule exists.  Returns
        ``unknown`` when no rule is registered or the value is missing.
        """
        rule = self._rules.get(metric)
        if rule is None or value is None:
            return GateResult(metric=metric, value=value, severity="unknown", threshold={})

        try:
            v = float(value)
        except (TypeError, ValueError):
            return GateResult(metric=metric, value=value, severity="unknown", threshold=dict(rule))

        thresholds: dict[str, float] = {
            k: float(rule[k]) for k in ("warn_below", "fail_below", "warn_above", "fail_above") if k in rule
        }

        if "fail_below" in thresholds and v < thresholds["fail_below"]:
            return GateResult(
                metric, value, "fail", thresholds, _msg(metric, v, "below", thresholds["fail_below"])
            )
        if "fail_above" in thresholds and v > thresholds["fail_above"]:
            return GateResult(
                metric, value, "fail", thresholds, _msg(metric, v, "above", thresholds["fail_above"])
            )
        if "warn_below" in thresholds and v < thresholds["warn_below"]:
            return GateResult(
                metric, value, "warn", thresholds, _msg(metric, v, "below", thresholds["warn_below"])
            )
        if "warn_above" in thresholds and v > thresholds["warn_above"]:
            return GateResult(
                metric, value, "warn", thresholds, _msg(metric, v, "above", thresholds["warn_above"])
            )
        return GateResult(metric=metric, value=value, severity="pass", threshold=thresholds)

    def has_rule(self, metric: str) -> bool:
        return metric in self._rules

    def metrics(self) -> list[str]:
        return sorted(self._rules.keys())


def _msg(metric: str, observed: float, direction: str, limit: float) -> str:
    return f"{metric} {observed:g} {direction} threshold {limit:g}"


def load_default_gates() -> QcGates:
    """Convenience: load the bundled default ``qc/config/qc_gates.yml``."""
    return QcGates.from_yaml(_DEFAULT_GATES_PATH)


__all__ = [
    "GateResult",
    "QcGates",
    "Severity",
    "load_default_gates",
]
