"""BRENDA: enzyme classes, molecular forms, and the compounds that regulate enzymes.

Identifiers use UniProt and EC; classification is an association with an explicit
source concept, not an invented catalytic reaction.

Regulation: each activator, inhibitor and cofactor BRENDA states for a protein
record is a relation from the compound to the protein's activity. The flat file
names compounds only; BRENDA's own ligand table (its knowledge graph's names,
synonyms and structures) gives a name its ligand ID, InChIKey and InChI where it
names exactly one ligand. Inhibitors carry their Ki values (mM).
"""

from __future__ import annotations
import functools
from pathlib import Path
from urllib.parse import urlencode

from biolink_model.datamodel import model
from biolink_model.datamodel.model import (
    CausalMechanismQualifierEnum,
    DirectionQualifierEnum,
    GeneOrGeneProductOrChemicalEntityAspectEnum,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, Topic

from omnipath_resources.inputs_v2._measurements import measurement
from omnipath_resources.inputs_v2._molecular_forms import protein_variants
from omnipath_core.molecular_forms import normalize_molecular_form

from omnipath_resources.inputs_v2.base import (
    ResourceConfig,
    Download,
    Resource,
    Dataset,
    ontology_entity_mapper,
)
from omnipath_resources.inputs_v2.parsers import brenda as _parsers
from omnipath_resources.inputs_v2.parsers import protein_ids
from omnipath_resources.internals.tabular_builder import (
    AssociationBuilder,
    AnnotationsBuilder,
    AssociationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from omnipath_resources.internals.cv_terms import (
    LicenseCV,
    OntologyCv,
    UpdateCategoryCV,
    ResourceCv,
)

# =================================== SET-UP ===================================

URL = 'https://www.brenda-enzymes.org/download.php'
_ENZYME_CLASSIFICATION_ONTOLOGY_ID = 'enzyme_classification'

config = ResourceConfig(
    id=ResourceCv.BRENDA,
    name='BRENDA',
    url='https://www.brenda-enzymes.org/index.php',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='41206471',
    annotation_ontologies=(OntologyCv.ENZYME_CLASSIFICATION,),
    topics=(Topic.METABOLISM,),
    content=(Content.COMPOUND_TARGET_INTERACTIONS, Content.MOLECULAR_ANNOTATIONS, Content.ONTOLOGIES),
    summary=(
        'Enzyme activators, inhibitors and cofactors, protein enzyme-class annotations '
        'and the enzyme classification hierarchy.'
    ),
    description=(
        'BRENDA is the main collection of enzyme functional data available to '
        'the scientific community.'
    ),
)


download = Download(
    url=URL,
    filename='brenda.txt.tar.gz',
    subfolder='brenda',
    large=True,
    ext='tar.gz',
    default_mode='rb',
    download_kwargs={
        'post': True,
        'query': {'dlfile': 'dl-textfile', 'accept-license': '1'},
    },
)

# BRENDA's ligand table, from its knowledge graph (https://sparql.dsmz.de/brenda).
_SPARQL = 'https://sparql.dsmz.de/api/brenda'
_PREFIXES = (
    'PREFIX d3o: <https://purl.dsmz.de/schema/> '
    'PREFIX dct: <http://purl.org/dc/terms/> '
    'PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> '
)
_COMPOUND_NAMES_QUERY = _PREFIXES + (
    'SELECT ?id ?kind ?name WHERE { ?c a d3o:ChemicalCompound ; dct:hasIdentifier ?id . '
    '{ ?c rdfs:label ?name BIND("label" AS ?kind) } UNION '
    '{ ?c d3o:hasSynonym ?name BIND("synonym" AS ?kind) } }'
)
_COMPOUND_STRUCTURES_QUERY = _PREFIXES + (
    'SELECT ?id ?inchikey ?inchi WHERE { ?c a d3o:ChemicalCompound ; dct:hasIdentifier ?id ; '
    'd3o:hasStructure ?s . ?s d3o:hasInChIKey ?inchikey ; d3o:hasInChI ?inchi }'
)


def _sparql_download(query: str, filename: str) -> Download:
    return Download(
        url=f'{_SPARQL}?{urlencode({"query": query, "action": "tsv_export"})}',
        filename=filename,
        subfolder='brenda',
        large=True,
        ext='tsv',
        default_mode='r',
    )


compound_names_download = _sparql_download(_COMPOUND_NAMES_QUERY, 'brenda_compound_names.tsv')
compound_structures_download = _sparql_download(_COMPOUND_STRUCTURES_QUERY, 'brenda_compound_structures.tsv')

# =================================== SCHEMA ===================================

f = FieldConfig(
    extract={},
    map={},
    transform={},
)

enzyme_ontology_schema = ontology_entity_mapper(
    _parsers.term_record_to_term,
    ontology_id=_ENZYME_CLASSIFICATION_ONTOLOGY_ID,
    identifier_type=Namespace.EC,
)

schema = EntityBuilder(
    entity_type=model.Protein,
    identifiers=IdentifiersBuilder(
        CV(term=lambda row: protein_ids.accession_namespace(row.get('UniProt')), value=f('UniProt'))
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.publications,
            value=f(
                'Refs',
                transform=lambda v: 'PMID:' + str(v).removeprefix('PMID:')
                if str(v).removeprefix('PMID:').isdigit()
                else None,
            ),
        )
    ),
    associations=AssociationsBuilder(
        AssociationBuilder(
            object_entity_type=model.OntologyClass,
            object_identifier_type=Namespace.EC,
            object_identifier=f(_parsers.ec_term_id),
        )
    ),
)


def _molecular_form(row):
    descriptor = row['descriptor']
    if row['section'] == 'EN':
        if descriptor.lower() in {
            'wild-type',
            'wild type',
            'more',
            'additional information',
            'no data',
        }:
            return None
        return normalize_molecular_form(
            {
                'variants': protein_variants(
                    descriptor,
                    coordinate_reference={
                        'identifier': None,
                        'coordinate_system': 'protein',
                        'position_base': 1,
                    },
                )
            },
            allow_resolved=False,
        )
    # Negated, uncertain and generic PM statements retain source context only.
    # We don't translate arbitrary narrative words to controlled PTM terms.
    if descriptor not in {
        'glycoprotein',
        'phosphoprotein',
        'lipoprotein',
        'side-chain modification',
        'proteolytic modification',
    }:
        return None
    return normalize_molecular_form(
        {
            'modifications': [
                {
                    'term': descriptor,
                    'description': '; '.join(row['notes']) or None,
                }
            ]
        },
        allow_resolved=False,
    )


molecular_forms_schema = EntityBuilder(
    entity_type=model.Protein,
    molecular_form=_molecular_form,
    identifiers=IdentifiersBuilder(
        CV(term=lambda row: protein_ids.accession_namespace(row.get('UniProt')), value=f('UniProt')),
        CV(term='brenda_protein_record', value=lambda row: row['EC'] + '#' + row['protein_record_id']),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term='brenda:molecular_observation',
            value=lambda row: _parsers.observation_texts(
                row['descriptor'], row['notes']
            ),
        ),
        CV(term=slots.in_taxon_label, value=f('organism')),
        CV(term='brenda:protein_note', value=f('protein_notes')),
        CV(term=slots.publications, value=lambda row: ['PMID:' + p for p in row.get('Refs', [])]),
    ),
)

_REGULATION_QUALIFIERS = {
    'activator': (DirectionQualifierEnum.increased, CausalMechanismQualifierEnum.activation),
    'inhibitor': (DirectionQualifierEnum.decreased, CausalMechanismQualifierEnum.inhibition),
    'cofactor': (None, CausalMechanismQualifierEnum.cofactor),
}

regulation_schema = RelationBuilder(
    subject=EntityBuilder(
        entity_type=model.ChemicalEntity,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.INCHIKEY, value=f('inchikey')),
            CV(term=Namespace.INCHI, value=f('inchi')),
            CV(term='brenda_ligand', value=f('brenda_ligand')),
            CV(term=Namespace.NAME, value=f('compound_name')),
        ),
        annotations=AnnotationsBuilder(),
    ),
    predicate=slots.affects,
    object=EntityBuilder(
        entity_type=model.Protein,
        identifiers=IdentifiersBuilder(
            CV(term=lambda row: protein_ids.accession_namespace(row.get('UniProt')), value=f('UniProt')),
            CV(term='brenda_protein_record', value=lambda row: row['EC'] + '#' + row['protein_record_id']),
        ),
        annotations=AnnotationsBuilder(CV(term=slots.in_taxon_label, value=f('organism'))),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.object_aspect_qualifier, value=GeneOrGeneProductOrChemicalEntityAspectEnum.activity),
        CV(term=slots.object_direction_qualifier, value=lambda row: _REGULATION_QUALIFIERS[row['role']][0]),
        CV(term=slots.causal_mechanism_qualifier, value=lambda row: _REGULATION_QUALIFIERS[row['role']][1]),
        CV(term='brenda:role', value=f('role')),
        CV(term='brenda:note', value=lambda row: row.get('notes') or None),
        CV(term=slots.publications, value=lambda row: ['PMID:' + p for p in row.get('Refs', [])]),
        CV(term='BAO:0000192', value=lambda row: [
            m for v in row.get('ki') or () if (m := measurement(v, 'mM', 'Ki'))
        ] or None),
    ),
)

# ================================= RESOURCE ===================================

def preparation_inputs():
    return [Path(protein_ids.__file__)]


resource = Resource(
    config=config,
    enzyme_ontology=Dataset(
        download=download,
        mapper=enzyme_ontology_schema,
        raw_parser=_parsers.iter_enzyme_ontology,
    ),
    protein_enzyme_class_annotations=Dataset(
        download=download,
        mapper=schema,
        raw_parser=_parsers.iter_protein_enzyme_class_annotations,
    ),
    molecular_forms=Dataset(download=download, mapper=molecular_forms_schema, raw_parser=_parsers.iter_molecular_forms),
    regulation=Dataset(
        download=download,
        mapper=regulation_schema,
        raw_parser=functools.partial(
            _parsers.iter_regulation,
            names_download=compound_names_download,
            structures_download=compound_structures_download,
        ),
    ),
)
