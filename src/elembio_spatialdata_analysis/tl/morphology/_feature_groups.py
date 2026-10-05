"""CellPaint morphology feature grouping.

By default, cell geometry / AreaShape features live in ``adata.obsm['cell_geometry']``.
If ``Loader.load_tables(..., include_cell_geometry=True)`` was used, the same
columns are also appended to ``adata.var`` / ``X``. Otherwise, ``adata.X`` /
``adata.var`` are expected to hold only CellProfiler morphology features:
Intensity, Texture, Granularity, RadialDistribution, Location, and the six
primary channel raw values.

``get_morphology_feature_groups`` categorises whatever is present in
``adata.var_names`` — it does not read ``obsm['cell_geometry']``.
"""

from __future__ import annotations

import anndata as ad


def get_morphology_feature_groups(adata: ad.AnnData) -> dict:
    """Identify and group CellPaint morphology features by type.

    Groups features found in ``adata.var_names`` into:

    - ``basic_channels``: raw channel intensity values (no CellProfiler prefix,
      or ending in ``_Nuclear``)
    - ``granularity_groups``: ``Granularity_<Scale>_<Channel>`` grouped by channel
    - ``texture_groups``: ``Texture_<Type>_3_XX_256_<Channel>`` grouped by type × channel
    - ``intensity_features``: ``Intensity_*`` features (flat list)
    - ``location_features``: ``Location_*`` features (flat list)
    - ``radial_groups``: ``RadialDistribution_<Type>_<Zone>of4_<Channel>`` grouped
      by type × channel
    - ``other_features``: anything that does not match the above patterns
    - ``all_channels``: set of all unique channel names identified

    .. note::
        ``AreaShape_*``, ``AreaUm``, and ``NuclearAreaUm`` are usually only in
        ``adata.obsm['cell_geometry']``. They appear in ``adata.var_names`` when
        the table was loaded with ``include_cell_geometry=True`` or in legacy
        layouts. In those cases they are grouped under ``other_features``.

    Parameters
    ----------
    adata:
        AnnData object with CellPaint features in ``var_names``.

    Returns
    -------
    dict with keys listed above.
    """
    feature_names = adata.var_names.tolist()

    basic_channels: list[str] = []
    granularity_groups: dict[str, list[str]] = {}
    texture_groups: dict[str, list[str]] = {}
    intensity_features: list[str] = []
    location_features: list[str] = []
    radial_groups: dict[str, list[str]] = {}
    other_features: list[str] = []
    all_channels: set[str] = set()

    for feat in feature_names:
        if feat.startswith("Granularity_"):
            rest = feat[len("Granularity_") :]
            parts = rest.split("_", 1)
            if len(parts) == 2:
                channel = parts[1]
                key = f"Granularity_{channel}"
                granularity_groups.setdefault(key, []).append(feat)
                all_channels.add(channel)
            else:
                other_features.append(feat)

        elif feat.startswith("Texture_"):
            rest = feat[len("Texture_") :]
            parts = rest.split("_")
            type_end_idx = None
            for i, part in enumerate(parts):
                if (
                    part == "3"
                    and i + 2 < len(parts)
                    and parts[i + 1] in ("00", "01", "02", "03")
                    and parts[i + 2] == "256"
                ):
                    type_end_idx = i
                    break
            if type_end_idx is not None and type_end_idx > 0:
                texture_type = "_".join(parts[:type_end_idx])
                channel = "_".join(parts[type_end_idx + 3 :])
                key = f"Texture_{texture_type}_{channel}"
                texture_groups.setdefault(key, []).append(feat)
                all_channels.add(channel)
            else:
                texture_groups.setdefault(feat, []).append(feat)

        elif feat.startswith("Intensity_"):
            intensity_features.append(feat)

        elif feat.startswith("Location_"):
            location_features.append(feat)

        elif feat.startswith("RadialDistribution_"):
            rest = feat[len("RadialDistribution_") :]
            parts = rest.split("_")
            zone_idx = next(
                (i for i, p in enumerate(parts) if p in ("1of4", "2of4", "3of4", "4of4")),
                None,
            )
            if zone_idx is not None and zone_idx > 0:
                radial_type = "_".join(parts[:zone_idx])
                channel = "_".join(parts[zone_idx + 1 :])
                key = f"RadialDistribution_{radial_type}_{channel}"
                radial_groups.setdefault(key, []).append(feat)
                all_channels.add(channel)
            else:
                radial_groups.setdefault(feat, []).append(feat)

        elif not feat.startswith(("AreaShape_",)) and feat not in ("AreaUm", "NuclearAreaUm"):
            # Basic channel: no recognised prefix, or ends in _Nuclear
            if "_" not in feat or feat.endswith("_Nuclear"):
                basic_channels.append(feat)
                all_channels.add(feat)
            else:
                other_features.append(feat)

        else:
            # AreaShape_* / AreaUm / NuclearAreaUm — legacy layout, stash in other
            other_features.append(feat)

    return {
        "basic_channels": basic_channels,
        "granularity_groups": granularity_groups,
        "texture_groups": texture_groups,
        "intensity_features": intensity_features,
        "location_features": location_features,
        "radial_groups": radial_groups,
        "other_features": other_features,
        "all_channels": all_channels,
    }
