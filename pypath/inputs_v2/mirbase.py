"""
miRBase: pre-miRNA and mature miRNA as distinct entities + maturation edges.

miRBase assigns two different accession namespaces to the two maturation
stages of a microRNA: ``MI#`` to the precursor (pre-miRNA stem-loop) and
``MIMAT#`` to each mature product cleaved from it. This module models them as
two distinct MicroRNA entities, distinguished by their source accession namespaces.
Each mature product derives_from its precursor.

Data come from the legacy ``pypath.inputs.mirbase`` tables (already on the
dlmachine download stack); all organisms are emitted (the miRBase name itself
encodes the organism, e.g. ``hsa-mir-21`` / ``hsa-miR-21-5p``).

Data source: https://www.mirbase.org/
"""

from __future__ import annotations

from biolink_model.datamodel.model import MicroRNA, slots
from omnipath_core.naming import Namespace

import collections
from collections.abc import Generator
from typing import Any

from pypath.internals.cv_terms import (
    LicenseCV,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    AssociationBuilder,
    AssociationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.inputs_v2.base import Dataset, Resource, ResourceConfig
from pypath.inputs.mirbase import (
    mirbase_mirna,
    mirbase_mirna_mature,
    mirbase_mirna_pre_mature,
)


config = ResourceConfig(
    id=ResourceCv.MIRBASE,
    name='miRBase',
    url='https://www.mirbase.org/',
    license=LicenseCV.PUBLIC,
    update_category=UpdateCategoryCV.REGULAR,
    primary_category='mirna',
    short='miRBase',
    pubmed='30423142',
    description=(
        'miRBase is the primary public repository and online resource for '
        'microRNA sequences and annotation. This inputs_v2 module emits '
        'precursor (MI#) and mature (MIMAT#) miRNAs as distinct entities '
        'joined by maturation relations.'
    ),
)


# =============================================================================
# Raw parsers (data come from the legacy mirbase tables, not a single file)
# =============================================================================


def _precursor_to_matures() -> dict[str, list[str]]:
    """Map each precursor MI# to its mature MIMAT# products."""
    matures: dict[str, list[str]] = collections.defaultdict(list)
    for precursor_row, mature_row in mirbase_mirna_pre_mature(None):
        mi_accession = precursor_row[1]
        mimat_accession = mature_row[3]
        if mi_accession and mimat_accession:
            if mimat_accession not in matures[mi_accession]:
                matures[mi_accession].append(mimat_accession)
    return dict(matures)


def _precursors_raw(
    opener: Any = None,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield one row per pre-miRNA, carrying its mature products for the edge."""
    matures = _precursor_to_matures()
    for row in mirbase_mirna(None):
        mi_accession = row[1]
        if not mi_accession:
            continue
        yield {
            'mirbase_pre': mi_accession,
            'name': row[2] if len(row) > 2 else None,
            'synonym': row[3] if len(row) > 3 else None,
            'description': row[4] if len(row) > 4 else None,
            'matures': matures.get(mi_accession, []),
        }


def _matures_raw(
    opener: Any = None,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield one row per mature miRNA (MIMAT#).

    Only the mature's own name (``row[1]``) is kept as the entity name; the
    parent precursor name (``row[2]``) is deliberately not emitted as a
    synonym, so precursor names resolve to MI# and mature names to MIMAT#
    without cross-contamination.
    """
    precursors = collections.defaultdict(list)
    for precursor, matures in _precursor_to_matures().items():
        for mature in matures:
            precursors[mature].append(precursor)
    for row in mirbase_mirna_mature(None):
        mimat_accession = row[3] if len(row) > 3 else None
        if not mimat_accession:
            continue
        yield {
            'mirbase_mat': mimat_accession,
            'precursors': precursors.get(mimat_accession, []),
            'name': row[1] if len(row) > 1 else None,
        }


# =============================================================================
# Schemas
# =============================================================================

f = FieldConfig()


precursors_schema = EntityBuilder(
    entity_type=MicroRNA,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.MIRBASE_PRECURSOR, value=f('mirbase_pre')),
        CV(term=Namespace.NAME, value=f('name')),
        CV(term=Namespace.SYNONYM, value=f('synonym')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.description, value=f('description')),
    ),
)


matures_schema = EntityBuilder(
    entity_type=MicroRNA,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.MIRBASE_MATURE, value=f('mirbase_mat')),
        CV(term=Namespace.NAME, value=f('name')),
    ),
    annotations=AnnotationsBuilder(),
    associations=AssociationsBuilder(
        AssociationBuilder(
            object_entity_type=MicroRNA,
            object_identifier_type=Namespace.MIRBASE_PRECURSOR,
            object_identifier=f('precursors'),
            predicate=slots.derives_from,
        ),
    ),
)


# =============================================================================
# Resource definition
# =============================================================================

resource = Resource(
    config,
    precursors=Dataset(
        download=None,
        mapper=precursors_schema,
        raw_parser=_precursors_raw,
    ),
    matures=Dataset(
        download=None,
        mapper=matures_schema,
        raw_parser=_matures_raw,
    ),
)


__all__ = ['config', 'resource']
