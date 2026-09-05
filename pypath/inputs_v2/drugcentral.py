"""DrugCentral mappings preserve action codes without guessing target complexes."""

from __future__ import annotations

import csv
import hashlib

from biolink_model.datamodel.model import (
    ChemicalEntity,
    DirectionQualifierEnum,
    NamedThing,
    Protein,
    slots,
)
from omnipath_core.naming import Namespace

from pypath.inputs_v2._measurements import measurement as _measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_tsv
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Entity,
    Identifier,
    Membership,
)
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.share.downloads import download_and_open

DRUGCENTRAL_STRUCTURES_URL = 'https://unmtid-shinyapps.net/download/DrugCentral/2021_09_01/structures.smiles.tsv'
config = ResourceConfig(
    id=ResourceCv.DRUGCENTRAL,
    name='DrugCentral',
    url='https://drugcentral.org/',
    license=LicenseCV.CC_BY_SA_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='36484092',
    primary_category='interactions',
    description='DrugCentral is an online drug information resource with active ingredients, chemical entities, pharmacologic action, indications, mechanism of action and drug-target interaction data.',
)
AFFINITY_TERMS = {
    'IC50': 'BAO:0000190',
    'EC50': 'BAO:0000188',
    'KI': 'BAO:0000192',
    'KD': 'BAO:0000034',
    'KON': 'BAO:0000480',
    'KOFF': 'BAO:0000479',
}
f = FieldConfig()
ACTION_DIRECTION = {
    'AGONIST': DirectionQualifierEnum.increased,
    'ACTIVATOR': DirectionQualifierEnum.increased,
    'ANTAGONIST': DirectionQualifierEnum.decreased,
    'INHIBITOR': DirectionQualifierEnum.decreased,
    'BLOCKER': DirectionQualifierEnum.decreased,
}


def _clean(value):
    return str(value or '').strip().strip('"')


def _values(value):
    # These lists are positionally aligned; dropping blanks changes identity.
    return [v.strip() for v in _clean(value).split('|')]


def _taxon(row):
    organism = _clean(row.get('ORGANISM'))
    taxon = {
        'Homo sapiens': '9606',
        'Mus musculus': '10090',
        'Rattus norvegicus': '10116',
    }.get(organism)
    return f'NCBITaxon:{taxon}' if taxon else None


chemical_builder = EntityBuilder(
    entity_type=ChemicalEntity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.DRUGCENTRAL, value=f('STRUCT_ID')),
        CV(term=Namespace.INCHIKEY, value=f('structure_inchikey')),
        CV(term=Namespace.NAME, value=f('DRUG_NAME')),
    ),
    annotations=AnnotationsBuilder(),
)
protein_builder = EntityBuilder(
    entity_type=Protein,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.UNIPROT, value=f('accession')),
        CV(term=Namespace.GENESYMBOL, value=f('gene')),
        CV(term=Namespace.UNIPROT_ENTRY, value=f('swissprot')),
        CV(term=Namespace.NAME, value=f('name')),
    ),
    annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=_taxon)),
)


def _target_entity(row):
    accessions = _values(row.get('ACCESSION'))
    genes = _values(row.get('GENE'))
    entries = _values(row.get('SWISSPROT'))
    proteins = [
        protein_builder.build(
            {
                **row,
                'accession': accession,
                'gene': genes[i] if i < len(genes) else None,
                'swissprot': entries[i] if i < len(entries) else None,
                'name': row.get('TARGET_NAME')
                if len(accessions) == 1
                else None,
            }
        )
        for i, accession in enumerate(accessions)
    ]
    proteins = [p for p in proteins if p is not None]
    if len(proteins) == 1:
        return proteins[0]
    if not proteins:
        return None
    return Entity(
        type=NamedThing,
        identifiers=[
            Identifier(type='drugcentral.target_group', value=hashlib.sha256('|'.join(sorted(accessions)).encode()).hexdigest()),
            Identifier(
                type=Namespace.NAME, value=_clean(row.get('TARGET_NAME'))
            )
        ],
        annotations=protein_builder.annotations.build(row),
        membership=[
            Membership(member=p, predicate=slots.has_member) for p in proteins
        ],
    )


def _predicate(row):
    action = _clean(row.get('ACTION_TYPE')).upper()
    if action in ACTION_DIRECTION:
        return slots.affects
    if action == 'BINDER':
        return slots.physically_interacts_with
    return slots.associated_with


interactions_schema = RelationBuilder(
    subject=chemical_builder,
    predicate=_predicate,
    object=_target_entity,
    annotations=AnnotationsBuilder(
        CV(
            term=slots.object_direction_qualifier,
            value=lambda row: ACTION_DIRECTION.get(
                _clean(row.get('ACTION_TYPE')).upper()
            ),
        ),
        CV(
            term=lambda row: AFFINITY_TERMS.get(
                str(row.get('ACT_TYPE') or '').upper()
            ),
            value=lambda row: _measurement(
                row.get('ACT_VALUE'),
                row.get('ACT_UNIT'),
                row.get('ACT_TYPE'),
                row.get('RELATION'),
            ),
        ),
        CV(term=slots.description, value=f('ACT_COMMENT')),
        CV(term=slots.supporting_data_source, value=f('ACT_SOURCE')),
        CV(term=slots.supporting_data_source, value=f('MOA_SOURCE')),
        CV(term=slots.source_record_urls, value=f('ACT_SOURCE_URL')),
        CV(term=slots.source_record_urls, value=f('MOA_SOURCE_URL')),
    ),
    identifiers=IdentifiersBuilder(),
)
map_drugcentral_interaction = interactions_schema


def _raw(opener, **kwargs):
    structures = download_and_open(
        url=DRUGCENTRAL_STRUCTURES_URL,
        filename='structures.smiles.tsv',
        subfolder='drugcentral',
        large=True,
        default_mode='r',
    )
    try:
        inchikeys = {
            _clean(r.get('ID')): _clean(r.get('InChIKey'))
            for r in csv.DictReader(structures.result, delimiter='\t')
        }
    finally:
        structures.close()
    for row in iter_tsv(opener, **kwargs):
        yield {
            **row,
            'structure_inchikey': inchikeys.get(_clean(row.get('STRUCT_ID'))),
        }


resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            url='https://unmtid-dbs.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz',
            filename='drug.target.interaction.tsv.gz',
            subfolder='drugcentral',
            ext='gz',
        ),
        mapper=interactions_schema,
        raw_parser=_raw,
    ),
)
