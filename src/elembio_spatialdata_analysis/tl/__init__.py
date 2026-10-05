"""Analysis tools: imaging morphology helpers and modality-name utilities.

Includes the morphology functions from ``tl.morphology`` and helpers such as
``normalize_target_type`` that map modality names (``"Transcript"``,
``"imaging"``, ``"OPS"``, ...) to a canonical form.
"""

from elembio_spatialdata_analysis.tl.constants import (
    TABLE_NAMES,
    normalize_target_type,
    normalized_to_table_name,
    table_name_to_normalized,
    validate_target_type,
)
from elembio_spatialdata_analysis.tl.morphology import (
    aggregate_morphology_features,
    get_morphology_feature_groups,
    select_morphology_channels,
)

__all__ = [
    "TABLE_NAMES",
    "aggregate_morphology_features",
    "get_morphology_feature_groups",
    "normalize_target_type",
    "normalized_to_table_name",
    "select_morphology_channels",
    "table_name_to_normalized",
    "validate_target_type",
]
