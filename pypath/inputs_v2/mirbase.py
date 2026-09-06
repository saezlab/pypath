"""
miRBase: pre-miRNA and mature miRNA as distinct entities + maturation edges.

miRBase assigns two different accession namespaces to the two maturation
stages of a microRNA: ``MI#`` to the precursor (pre-miRNA stem-loop) and
``MIMAT#`` to each mature product cleaved from it. This module models them as
two distinct MicroRNA entities, distinguished by their source accession namespaces.
Each mature product derives_from its precursor.

Data come from the official release-22 EMBL archive, pinned explicitly because
release 23 does not currently expose its full annotation download. All organisms
and explicit precursor/product links are retained. Raw evidence records carry
the source release and URL independently of the OmniPath build version.

Data source: https://www.mirbase.org/
"""

from __future__ import annotations

from biolink_model.datamodel.model import MicroRNA, slots
from omnipath_core.naming import Namespace

import re
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
from pypath.inputs_v2.base import (
    Dataset, Download, Resource, ResourceConfig, _first_handle,
)


SOURCE_RELEASE = '22'
SOURCE_URL = f'https://www.mirbase.org/download_version_files/{SOURCE_RELEASE}/miRNA.dat'

download = Download(
    url=SOURCE_URL,
    filename='miRNA.dat',
    subfolder=f'mirbase/release-{SOURCE_RELEASE}',
    large=True,
    ext='.dat',
    default_mode='r',
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
        'joined by maturation relations. Source: official miRBase release 22 EMBL archive.'
    ),
)


# =============================================================================
# Raw parsers: release-pinned EMBL records
# =============================================================================


def _embl_records(opener):
    """Stream complete EMBL entries; reject error pages and truncated records."""
    handle = _first_handle(opener)
    if handle is None:
        raise ValueError('miRBase EMBL download is empty')
    record = []
    count = 0
    for line in handle:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.rstrip('\r\n')
        if not record:
            if not line.strip():
                continue
            if not line.startswith('ID   '):
                raise ValueError('Expected miRBase EMBL ID record')
        elif line.startswith('ID   '):
            raise ValueError('Unterminated miRBase EMBL entry')
        record.append(line)
        if line == '//':
            yield _parse_entry(record)
            count += 1
            record = []
    if record or not count:
        raise ValueError('Truncated or empty miRBase EMBL archive')


def _parse_entry(lines):
    name = lines[0][5:].split()[0]
    accessions = [
        accession.strip()
        for line in lines if line.startswith('AC   ')
        for accession in line[5:].split(';') if accession.strip()
    ]
    if len(accessions) != 1 or not re.fullmatch(r'MI\d+', accessions[0]):
        raise ValueError(f'Invalid precursor accession for {name}')
    features = []
    feature = None
    for line in lines:
        if not line.startswith('FT   '):
            continue
        key = line[5:21].strip()
        if key:
            feature = [] if key == 'miRNA' else None
            if feature is not None:
                features.append(feature)
        elif feature is not None:
            feature.append(line[21:].strip())
    products = {}
    for feature in features:
        qualifiers = ' '.join(feature)
        accession = re.search(r'/accession="(MIMAT\d+)"', qualifiers)
        product = re.search(r'/product="([^"\n]+)"', qualifiers)
        if not accession or not product:
            raise ValueError(f'Missing mature accession/product for {name}')
        accession, product = accession[1], product[1]
        if accession in products and products[accession] != product:
            raise ValueError(f'Conflicting mature names for {accession}')
        products[accession] = product
    return {
        'mirbase_pre': accessions[0],
        'name': name,
        'description': ' '.join(line[5:] for line in lines if line.startswith('DE   ')),
        'products': products,
        'source_release': SOURCE_RELEASE,
        'source_url': SOURCE_URL,
    }


def _precursors_raw(
    opener: Any = None,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    for entry in _embl_records(opener):
        products = entry.pop('products')
        yield {**entry, 'synonym': None, 'matures': list(products)}


def _matures_raw(
    opener: Any = None,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Aggregate all explicit parents; never use precursor names as aliases."""
    matures = {}
    for entry in _embl_records(opener):
        for accession, name in entry['products'].items():
            row = matures.setdefault(accession, {
                'mirbase_mat': accession,
                'name': name,
                'precursors': [],
                'source_release': SOURCE_RELEASE,
                'source_url': SOURCE_URL,
            })
            if row['name'] != name:
                raise ValueError(f'Conflicting mature names for {accession}')
            if entry['mirbase_pre'] not in row['precursors']:
                row['precursors'].append(entry['mirbase_pre'])
    yield from matures.values()


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
        download=download,
        mapper=precursors_schema,
        raw_parser=_precursors_raw,
    ),
    matures=Dataset(
        download=download,
        mapper=matures_schema,
        raw_parser=_matures_raw,
    ),
)


__all__ = ['config', 'resource']
