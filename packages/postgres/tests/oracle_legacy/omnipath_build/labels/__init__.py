"""Stored human-readable entity labels.

Every entity gets a best label by a per-entity-type cascade, precomputed at
build (never per request) and stored on ``entity.label`` with the producing rule
recorded in ``entity.label_rule``. This package holds the label builders run in
``derive``:

* gene-based entities → gene symbol (``populate_entity_labels``);
* chemicals → brevity-first name cascade (``populate_chemical_labels``);
* lipids → Goslin short names (omnipath-metabo postbuild, future).

A universal identifier fallback guarantees a non-empty label for every entity;
the richer per-type rules overwrite it as they land.
"""

from __future__ import annotations

from omnipath_build.labels.chemical_labels import (
    ChemicalLabelStats,
    populate_chemical_labels,
)
from omnipath_build.labels.entity_labels import (
    EntityLabelStats,
    populate_entity_labels,
)
from omnipath_build.labels.entity_name import (
    EntityNameStats,
    populate_entity_name,
)

__all__ = [
    'ChemicalLabelStats',
    'EntityLabelStats',
    'EntityNameStats',
    'populate_chemical_labels',
    'populate_entity_labels',
    'populate_entity_name',
]
