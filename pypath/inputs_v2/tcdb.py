"""
Parse TCDB (Transporter Classification Database) data and emit Entity records.

Two datasets are provided:

- transporters: proteins annotated with their TC number and transporter family
  name, derived from acc2tcid.py (UniProt → TCID) joined with families.py
  (family prefix → family name).

- substrates: transporter-substrate interactions linking UniProt accessions to
  CHEBI-identified small molecules, derived from getSubstrates.py (TCID →
  substrates) joined with acc2tcid.py (TCID → UniProt).
"""

from __future__ import annotations

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content
from omnipath_core.interaction_profiles import TRANSPORT_QUALIFIERS

import functools
from pathlib import Path
from pypath.inputs_v2.parsers import protein_ids
from pypath.inputs_v2.parsers.protein_ids import accession_namespace as _accession_namespace

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers import tcdb as _parsers
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import Entity, Identifier, Membership
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)


config = ResourceConfig(
    id=ResourceCv.TCDB,
    name='TCDB',
    url='https://www.tcdb.org/',
    license=LicenseCV.CC_BY_SA_3_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='33170213',
    content=(Content.MOLECULAR_TRANSPORT, Content.MOLECULAR_ANNOTATIONS),
    summary='Transporter protein classifications and transporter–substrate relationships.',
    description=(
        'The Transporter Classification Database (TCDB) provides a '
        'comprehensive classification system for membrane transport proteins '
        'based on both function and phylogeny. It assigns TC numbers to '
        'transporters and curates their substrates using CHEBI identifiers.'
    ),
)

f = FieldConfig(
    extract={
        'chebi': r'CHEBI:(\d+)',
    },
)

def preparation_inputs():
    return [Path(protein_ids.__file__)]


_transporters_schema = EntityBuilder(
    entity_type=model.Protein,
    identifiers=IdentifiersBuilder(
        CV(term=lambda row: _accession_namespace(row.get('uniprot')), value=f('uniprot'))
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.has_topic,
            value=f('tcid', transform=lambda v: 'TC:' + str(v)),
        )
    ),
)

# TC numbers classify proteins; they do not identify a transport event.
_transport_schema = RelationBuilder(
    subject=EntityBuilder(
        entity_type=model.Protein,
        identifiers=IdentifiersBuilder(
            CV(term=lambda row: _accession_namespace(row.get('transporter_uniprot')), value=f('transporter_uniprot'))
        ),
    ),
    predicate=slots.affects,
    object=EntityBuilder(
        entity_type=model.ChemicalEntity,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.CHEBI, value=f('substrate_id', extract='chebi')),
            CV(term=Namespace.NAME, value=f('substrate_name')),
        ),
        annotations=AnnotationsBuilder(),
    ),
    annotations=AnnotationsBuilder(
        *(CV(term=term, value=value) for term, value in TRANSPORT_QUALIFIERS),
        CV(
            term=slots.has_topic,
            value=f('tcid', transform=lambda v: 'TC:' + str(v)),
        )
    ),
)


_families_download = Download(
    url='http://www.tcdb.org/cgi-bin/projectv/public/families.py',
    filename='tcdb_families.html',
    subfolder='tcdb',
    ext='html',
)
FAMILY_NAMESPACE = 'tcdb_family'


def _iter_families(opener, max_records: int | None = None, force_refresh: bool = False, **_kwargs):
    """Each family of families.py (a TC prefix such as ``2.A.1`` and its
    name) with the proteins whose TC numbers are in it."""
    if not opener or not opener.result:
        return
    families_opener = _families_download.open(force_refresh=force_refresh)
    if not families_opener or not families_opener.result:
        return
    names = _parsers._parse_families(families_opener.result)
    members: dict[str, list[str]] = {}
    for tcid, uniprots in _parsers._build_tcid_to_uniprots(opener.result).items():
        prefix = _parsers._family_prefix(tcid)
        if prefix:
            members.setdefault(prefix, []).extend(uniprots)
    for count, (prefix, uniprots) in enumerate(sorted(members.items())):
        if max_records is not None and count >= max_records:
            return
        yield {'id': prefix, 'name': names.get(prefix, ''), 'members': sorted(set(uniprots))}


def map_family(row: dict) -> Entity:
    """A transporter family with its proteins as members."""
    identifiers = [Identifier(type=FAMILY_NAMESPACE, value=row['id'])]
    if row['name']:
        identifiers.append(Identifier(type=Namespace.NAME, value=row['name']))
    return Entity(
        type=model.ProteinFamily,
        identifiers=identifiers,
        membership=[
            Membership(
                member=Entity(
                    type=model.Protein,
                    identifiers=[Identifier(type=namespace, value=accession)],
                )
            )
            for accession in row['members']
            if (namespace := _accession_namespace(accession))
        ],
    )


resource = Resource(
    config,
    transporters=Dataset(
        download=(
            _acc2tc := Download(
                url='http://www.tcdb.org/cgi-bin/projectv/public/acc2tcid.py',
                filename='tcdb_acc2tc.tsv',
                subfolder='tcdb',
                ext='tsv',
            )
        ),
        mapper=_transporters_schema,
        raw_parser=functools.partial(
            _parsers.transporters,
            families_download=_families_download,
        ),
    ),
    families=Dataset(
        download=_acc2tc,
        mapper=map_family,
        raw_parser=_iter_families,
    ),
    transports=Dataset(
        download=Download(
            url='https://www.tcdb.org/cgi-bin/substrates/getSubstrates.py',
            filename='tcdb_substrates.tsv',
            subfolder='tcdb',
            ext='tsv',
        ),
        mapper=_transport_schema,
        raw_parser=functools.partial(
            _parsers.transport_substrates,
            acc2tc_download=_acc2tc,
        ),
    ),
)
