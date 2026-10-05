"""Modality routing constants for the SpatialData analysis toolkit.

This module owns target-type aliases and table-name mappings.
- Morphology-channel selection lives in :mod:`elembio_spatialdata_analysis.tl.morphology`.

Imaging modality
~~~~~~~~~~~~~~~~
The pipeline emits a single merged table
``global-table-imaging`` (lowercase) for all imaging batches, with
``uns["metadata_attrs"]["target_type"] = "imaging"``.  Per-feature
sub-types (``"CellPaint"``, ``"ProteinIF"``) are stored in
``var["target_type"]``.

``"CellPaint"`` and ``"cell_paint"`` remain valid inputs throughout the
API and normalize to ``"imaging"`` for backward compatibility.
"""

# Valid target type inputs mapped to their internal normalized form.
# Accepts capitalized table-name forms (matching zarr store keys), the
# lowercase internal form used in routing logic, and backward-compat
# aliases.
VALID_TARGET_TYPES: dict[str, str] = {
    # Capitalized — matches zarr table keys (e.g. "global-table-Transcript")
    "Transcript": "transcript",
    "Protein": "protein",
    "ThreePrimeUntargeted": "threeprimeuntargeted",
    # Imaging — the primary name for the unified imaging modality.
    # The merged table is ``global-table-imaging`` (lowercase, unlike the
    # other modalities) and ``uns["metadata_attrs"]["target_type"]`` is
    # ``"imaging"``.  CellPaint and ProteinIF are sub-types distinguished
    # by ``var["target_type"]`` within the merged table.
    "Imaging": "imaging",
    "imaging": "imaging",
    # Backward-compat aliases — old code that passes ``"CellPaint"`` or
    # ``"cell_paint"`` continues to work; both normalize to ``"imaging"``.
    "CellPaint": "imaging",
    "cell_paint": "imaging",
    # ProteinIF — per-feature sub-type within the imaging table.
    # Accepted as an input alias so callers can pass the value read
    # directly from ``var["target_type"]``.
    "ProteinIF": "imaging",
    # OPS — written as ``target_type="OPS"`` by the pipeline in every
    # per-target-site OPS table's ``uns['metadata_attrs']``.
    # ``"SpecializedTargeted"`` is accepted as a symmetric alias because
    # it appears as the *batch_type* (pass-through from Panel.json
    # ``AssayType``) in global metadata — never as a target_type in the
    # AnnData itself, so accepting it is harmless cosmetics.
    "OPS": "ops",
    "SpecializedTargeted": "ops",
    # Protein / Transcript carry a ``*Barcoding`` *batch_type* (pass-through
    # from Panel.json ``AssayType``) distinct from their ``Protein`` /
    # ``Transcript`` target_type — accept both, symmetric with
    # ``SpecializedTargeted`` for OPS above.
    "ProteinBarcoding": "protein",
    "TranscriptBarcoding": "transcript",
    # Lowercase canonical form (also accepted as input)
    "transcript": "transcript",
    "protein": "protein",
    "threeprimeuntargeted": "threeprimeuntargeted",
    "ops": "ops",
}

# Supported normalized target types (lowercase).
SUPPORTED_TARGET_TYPES: tuple[str, ...] = (
    "transcript",
    "protein",
    "imaging",
    "threeprimeuntargeted",
    "ops",
)

# Mapping from normalized target types to table-name suffixes used in SpatialData.
# Used when constructing table names like ``f"global-table-{suffix}"``.
#
# NOTE: ``"imaging"`` maps to lowercase ``"imaging"`` — unlike the other
# modalities which are capitalized.  This matches the build repo's
# ``_IMAGING_TABLE_TYPE = "imaging"`` (``global-table-imaging``).
#
# NOTE: OPS is intentionally absent.  OPS tables are emitted by the
# pipeline with per-target-site keys
# (``global-table-{batch_name}-{target_site}``, e.g. ``global-table-B01-A549``);
# no aggregated ``global-table-OPS`` ever exists.  Adding ``"ops": "OPS"``
# here would silently produce broken table-key lookups in any caller
# that does ``f"global-table-{normalized_to_table_name(t)}"``.  OPS
# discovery uses metadata-driven enumeration instead — see
# :meth:`elembio_spatialdata_analysis.load.backends.SpatialBackend.table_suffixes_for_target_type`.
NORMALIZED_TO_TABLE_NAME: dict[str, str] = {
    "transcript": "Transcript",
    "protein": "Protein",
    "imaging": "imaging",
    "threeprimeuntargeted": "ThreePrimeUntargeted",
}

# Reverse mapping: table-name suffixes to normalized types.
# Includes ``"CellPaint"`` as a legacy alias so that old-format SDOs
# (``global-table-CellPaint``) are still recognised.
TABLE_NAME_TO_NORMALIZED: dict[str, str] = {v: k for k, v in NORMALIZED_TO_TABLE_NAME.items()}
TABLE_NAME_TO_NORMALIZED["CellPaint"] = "imaging"  # legacy SDO compat

# Table-name suffixes as used in SpatialData.
TABLE_NAMES: list[str] = list(NORMALIZED_TO_TABLE_NAME.values())

# Known suffixes for global table keys (used when discovering modalities)
_TABLE_KEY_SUFFIXES = ("-Filtered", "-corrected")


def parse_table_key_to_normalized(table_key: str) -> str | None:
    """Parse a SpatialData table key into a normalized target type, if it is a known modality.

    Handles keys like ``"global-table-Transcript"``,
    ``"global-table-imaging"``, ``"global-table-CellPaint"`` (legacy),
    ``"global-table-ThreePrimeUntargeted-Filtered"``, or
    ``"global-table-Protein-corrected"``.  Returns ``None`` for unknown
    or per-batch tables (e.g. ``"global-table-B06"``).

    Parameters
    ----------
    table_key : str
        Key from sdata.tables or folder name under tables/.

    Returns
    -------
    str or None
        Normalized target type (e.g. ``"transcript"``, ``"imaging"``),
        or ``None`` if unknown.
    """
    if not isinstance(table_key, str) or not table_key.startswith("global-table-"):
        return None
    base = table_key.replace("global-table-", "", 1)
    for suffix in _TABLE_KEY_SUFFIXES:
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break
    if base not in TABLE_NAME_TO_NORMALIZED:
        return None
    return TABLE_NAME_TO_NORMALIZED[base]


def normalize_target_type(target_type: str) -> str:
    """Normalize a target type string to its internal lowercase form.

    Accepts the 5 canonical modality names in either their capitalized
    table-name form, their lowercase internal form, or backward-compat
    aliases.  Any other value raises ``ValueError``.

    Parameters
    ----------
    target_type : str
        One of:

        - ``"Transcript"`` / ``"transcript"``
        - ``"Protein"`` / ``"protein"``
        - ``"Imaging"`` / ``"imaging"`` / ``"CellPaint"`` / ``"cell_paint"``
          / ``"ProteinIF"``
        - ``"ThreePrimeUntargeted"`` / ``"threeprimeuntargeted"``
        - ``"OPS"`` / ``"SpecializedTargeted"`` / ``"ops"``

    Returns
    -------
    str
        Lowercase internal form, e.g. ``"transcript"``, ``"imaging"``,
        ``"ops"``.

    Raises
    ------
    ValueError
        If *target_type* is not one of the accepted strings above.

    Examples
    --------
    >>> normalize_target_type("Transcript")
    'transcript'
    >>> normalize_target_type("Imaging")
    'imaging'
    >>> normalize_target_type("CellPaint")
    'imaging'
    >>> normalize_target_type("cell_paint")
    'imaging'
    >>> normalize_target_type("ProteinIF")
    'imaging'
    >>> normalize_target_type("OPS")
    'ops'
    """
    if not isinstance(target_type, str) or not target_type.strip():
        raise ValueError("`target_type` must be a non-empty string.")

    normalized = VALID_TARGET_TYPES.get(target_type.strip())
    if normalized is None:
        raise ValueError(
            f"Invalid target_type {target_type!r}. "
            'Valid values: "Transcript", "Protein", "Imaging" '
            '(or "CellPaint"), "ThreePrimeUntargeted", "OPS" '
            '(or "SpecializedTargeted"), '
            "and their lowercase internal forms."
        )
    return normalized


def validate_target_type(target_type: str) -> str:
    """Validate *target_type* and return it unchanged.

    Unlike :func:`normalize_target_type`, this preserves the original casing so
    the caller can use it directly as a zarr table-key suffix
    (``f"global-table-{target_type}"``).  Raises ``ValueError`` for any
    unrecognised value.

    Parameters
    ----------
    target_type : str
        Any value accepted by :func:`normalize_target_type`.

    Returns
    -------
    str
        The input string, unchanged.
    """
    normalize_target_type(target_type)  # raises if invalid
    return target_type


def table_name_to_normalized(table_name: str) -> str:
    """
    Convert a capitalized table name to normalized target type.

    Parameters
    ----------
    table_name : str
        Capitalized table name (e.g., "Transcript", "Protein").

    Returns
    -------
    str
        Normalized target type (e.g., "transcript", "protein").
    """
    return TABLE_NAME_TO_NORMALIZED.get(table_name, normalize_target_type(table_name))


def normalized_to_table_name(normalized_type: str) -> str:
    """
    Convert a normalized target type to capitalized table name.

    Parameters
    ----------
    normalized_type : str
        Normalized target type (e.g., "transcript", "protein").

    Returns
    -------
    str
        Capitalized table name (e.g., "Transcript", "Protein").
    """
    return NORMALIZED_TO_TABLE_NAME.get(normalized_type, normalized_type.title())
