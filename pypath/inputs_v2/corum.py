"""
Parse CORUM data and emit Entity records.

This module converts CORUM protein complex data into Entity records using the
declarative schema pattern defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin

from pypath.internals.cv_terms import (
    LicenseCV,
    OntologyCv,
    UpdateCategoryCV,
    ResourceCv,
)
from pypath.internals.tabular_builder import (
    AssociationBuilder,
    AssociationsBuilder,
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembershipBuilder,
    MembersFromList,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.corum import _raw
from pypath.internals.silver_schema import Entity, Identifier


config = ResourceConfig(
    id=ResourceCv.CORUM,
    name='CORUM',
    url='https://mips.helmholtz-muenchen.de/corum/',
    license=LicenseCV.CC_BY_NC_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='30357367',
    annotation_ontologies=(OntologyCv.GENE_ONTOLOGY,),
    content=(Content.PROTEIN_COMPLEXES,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Curated mammalian protein complexes and their subunits.',
    description=(
        'CORUM is a manually curated repository of experimentally characterized '
        'protein complexes from mammalian organisms. Each complex includes '
        'information about subunit composition, functional annotations, and '
        'literature references.'
    ),
)

f = FieldConfig(
    delimiter=';',
    extract={
        'uniprot': r'((?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?(?:-PRO_\d+)?)',
    },
    map={
        'organism_taxid': {
            'Human': '9606',
            'Homo sapiens': '9606',
            'Mouse': '10090',
            'Mus musculus': '10090',
            'Rat': '10116',
            'Rattus norvegicus': '10116',
        },
    },
)

FUNCAT_NAMESPACE = 'funcat'


def _funcat(row: dict) -> list[tuple[str, str]]:
    """FunCat IDs with their descriptions, both ';'-separated in the same order."""
    ids = [v.strip() for v in str(row.get('FunCat ID') or '').split(';')]
    names = [v.strip() for v in str(row.get('FunCat description') or '').split(';')]
    if len(names) != len(ids):
        names = [''] * len(ids)
    return [(i, n) for i, n in zip(ids, names) if i and i != 'None']


def _iter_funcat(opener, max_records=None, **kwargs):
    seen = set()
    for row in _raw(opener, **kwargs):
        for term in _funcat(row):
            if term[0] not in seen:
                seen.add(term[0])
                yield {'id': term[0], 'name': term[1]}
                if max_records is not None and len(seen) >= max_records:
                    return


def map_funcat(row: dict) -> Entity:
    identifiers = [Identifier(type=FUNCAT_NAMESPACE, value=row['id'])]
    if row['name']:
        identifiers.append(Identifier(type=Namespace.NAME, value=row['name']))
    return Entity(type=model.OntologyClass, identifiers=identifiers)


complexes_schema = EntityBuilder(
    entity_type=model.MacromolecularComplex,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.CORUM, value=f('ComplexID')),
        CV(term=Namespace.NAME, value=f('ComplexName')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.publications,
            value=f(
                'PubMed ID',
                delimiter=';',
                transform=lambda v: 'PMID:' + str(v).removeprefix('PMID:')
                if str(v).removeprefix('PMID:').isdigit()
                and int(str(v).removeprefix('PMID:')) > 0
                else None,
            ),
        )
    ),
    associations=AssociationsBuilder(
        AssociationBuilder(
            object_entity_type=model.OntologyClass,
            object_identifier_type=Namespace.GO,
            object_identifier=f('GO ID', delimiter=';'),
        ),
        # The complex's MIPS FunCat categories ("10.01.09.05"); the terms
        # dataset names them from the descriptions beside the IDs.
        AssociationBuilder(
            predicate=slots.associated_with,
            object_entity_type=model.OntologyClass,
            object_identifier_type=FUNCAT_NAMESPACE,
            object_identifier=lambda row: [i for i, _ in _funcat(row)],
        ),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=model.Protein,
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.UNIPROT,
                    value=f(
                        'subunits(UniProt IDs)',
                        delimiter=';',
                        extract='uniprot',
                    ),
                )
            ),
            entity_annotations=AnnotationsBuilder(
                CV(
                    term=slots.in_taxon,
                    value=f(
                        'Organism',
                        map='organism_taxid',
                        transform=lambda v: 'NCBITaxon:'
                        + str(v).removeprefix('NCBITaxon:')
                        if str(v).removeprefix('NCBITaxon:').isdigit()
                        and int(str(v).removeprefix('NCBITaxon:')) > 0
                        else None,
                    ),
                )
            ),
        )
    ),
)

_download = Download(
    url='https://rescued.omnipathdb.org/CORUM_allComplexes.txt.zip',
    filename='CORUM_allComplexes.txt.zip',
    subfolder='corum',
    large=True,
    ext='zip',
)
resource = Resource(
    config,
    complexes=Dataset(
        download=_download,
        mapper=complexes_schema,
        raw_parser=_raw,
    ),
    terms=Dataset(
        download=_download,
        mapper=map_funcat,
        raw_parser=_iter_funcat,
    ),
)
