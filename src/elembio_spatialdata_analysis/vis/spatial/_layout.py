"""Flowcell well grid (letter columns × number rows). Does not import vis.filter."""

from __future__ import annotations

from collections.abc import Sequence

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

_FLOWCELL_WELL_SETS: dict[frozenset[str], str] = {
    frozenset(f"{col}{row}" for row in layout["rows"] for col in layout["cols"]): key
    for key, layout in FLOWCELL_ORIENTATIONS.items()
}


def _sorted_well_ids(labels: Sequence[object]) -> list[str]:
    def _key(label: object) -> tuple[int, str]:
        text = str(label)
        prefix = "".join(filter(str.isalpha, text)).upper()
        digits = "".join(filter(str.isdigit, text))
        suffix = int(digits) if digits else 0
        return suffix, prefix

    return sorted((str(lab) for lab in labels), key=_key)


def _auto_detect_flowcell(well_ids: Sequence[object]) -> str | None:
    observed = frozenset(str(lab).upper() for lab in well_ids if str(lab).strip())
    for known_wells, key in _FLOWCELL_WELL_SETS.items():
        if observed <= known_wells:
            return key
    return None


def resolve_flowcell_grid(
    present: Sequence[str],
    *,
    wells: Sequence[str] | None = None,
) -> tuple[list[str], int, int]:
    """Return ``(slot_order, n_rows, n_cols)`` for a plate of well images.

    Slots follow vis.filter: alphabetic columns, numeric rows (12-well is
    2 × 6: A1–F1 then A2–F2). Missing wells stay in the grid as empty slots.
    ``wells`` restricts which ids are considered present for auto-detect.
    """
    present_ids = [str(w) for w in present]
    if wells is not None:
        present_ids = [str(w) for w in wells]
    if not present_ids:
        raise ValueError("No wells with images.")
    key = _auto_detect_flowcell(present_ids)
    if key is None:
        order = _sorted_well_ids(present_ids)
        n = len(order)
        n_cols = min(3, max(1, n))
        n_rows = (n + n_cols - 1) // n_cols
        return order, n_rows, n_cols
    layout = FLOWCELL_ORIENTATIONS[key]
    cols_order = [str(c).upper() for c in layout["cols"]]
    rows_order = [str(r) for r in layout["rows"]]
    slot_order = [f"{col}{row}" for row in rows_order for col in cols_order]
    extra = [w for w in _sorted_well_ids(present_ids) if w not in slot_order]
    slot_order = slot_order + extra
    n_cols = len(cols_order)
    n_rows = (len(slot_order) + n_cols - 1) // n_cols
    return slot_order, n_rows, n_cols
