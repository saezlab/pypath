"""WikiPathways RDF pathways and diagram interactions in native Biolink.

Only exact interaction classes supply effect direction. Unsigned conversion
diagram edges remain broad interactions until a reaction-node model is supplied.
"""

from __future__ import annotations

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace

from pypath.internals.cv_terms import (
    LicenseCV,
    OntologyCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    AssociationBuilder,
    AssociationsBuilder,
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.wikipathways import _raw, current_rdf_url


config = ResourceConfig(
    id=ResourceCv.WIKIPATHWAYS,
    name='WikiPathways',
    url='https://www.wikipathways.org/',
    license=LicenseCV.CC0_1_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='37941138',
    primary_category='pathways',
    annotation_ontologies=(OntologyCv.WIKIPATHWAYS,),
    description=(
        'WikiPathways is a community-curated pathway database. '
        'This inputs_v2 module parses the current RDF pathway export and '
        'emits pathway and directed interaction entities.'
    ),
)


entity_type_map = {
    'protein': model.Protein,
    'gene': model.Gene,
    'rna': model.RNAProduct,
    'dna': model.NucleicAcidEntity,
    'complex': model.MacromolecularComplex,
    'chemical': model.ChemicalEntity,
    'physical_entity': model.PhysicalEntity,
    'protein_family': model.ProteinFamily,
    'reaction': model.MolecularActivity,
    'degradation': model.BiologicalProcess,
    'catalysis': model.MolecularActivity,
    'control': model.BiologicalProcess,
    'interaction': model.BiologicalProcess,
    'cv_term': model.OntologyClass,
    'pathway': model.Pathway,
}
_TAXON_SCOPED_ENTITY_TYPES = {
    model.Protein,
    model.Gene,
    model.RNAProduct,
    model.NucleicAcidEntity,
}

f = FieldConfig(
    extract={
        'ensembl_id': '^(ENS[A-Z0-9]*\\d+(?:\\.\\d+)?)$',
        'uniprot_id': '^((?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\\d+)?)$',
        'entrez_id': '^(\\d+)$',
        'chebi': '^(?:CHEBI:)?(\\d+)$',
    },
    map={
        'entity_type': lambda value: entity_type_map.get(
            value, model.PhysicalEntity
        )
    },
)


def _member_taxon_id(prefix: str):
    def _value(row):
        entity_type = entity_type_map.get(
            row.get(f'{prefix}_entity_type'), model.PhysicalEntity
        )
        return (
            'NCBITaxon:' + str(row.get('taxon_id'))
            if entity_type in _TAXON_SCOPED_ENTITY_TYPES
            and str(row.get('taxon_id') or '').isdigit()
            and int(row['taxon_id']) > 0
            else None
        )

    return _value


pathways_schema = EntityBuilder(
    entity_type=model.Pathway,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.WIKIPATHWAYS, value=f('pathway_id')),
        CV(term=Namespace.WIKIPATHWAYS_VERSION, value=f('pathway_version_id')),
        CV(term=Namespace.NAME, value=f('title')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.in_taxon,
            value=f(
                'taxon_id',
                transform=lambda v: 'NCBITaxon:'
                + str(v).removeprefix('NCBITaxon:')
                if str(v).removeprefix('NCBITaxon:').isdigit()
                and int(str(v).removeprefix('NCBITaxon:')) > 0
                else None,
            ),
        ),
        CV(
            term=slots.publications,
            value=f(
                'pubmed_ids',
                delimiter=';',
                transform=lambda v: 'PMID:' + str(v).removeprefix('PMID:')
                if str(v).removeprefix('PMID:').isdigit()
                and int(str(v).removeprefix('PMID:')) > 0
                else None,
            ),
        ),
        CV(term=slots.description, value=f('description')),
    ),
)


def _participant_builder(prefix: str, role) -> EntityBuilder:
    return EntityBuilder(
        entity_type=f(f'{prefix}_entity_type', map='entity_type'),
        identifiers=IdentifiersBuilder(
            CV(
                term=Namespace.UNIPROT,
                value=f(
                    f'{prefix}_uniprot', delimiter=';', extract='uniprot_id'
                ),
            ),
            CV(
                term=Namespace.ENTREZ,
                value=f(f'{prefix}_entrez', delimiter=';', extract='entrez_id'),
            ),
            CV(
                term=Namespace.ENSEMBL,
                value=f(
                    f'{prefix}_ensembl', delimiter=';', extract='ensembl_id'
                ),
            ),
            CV(
                term=Namespace.CHEBI,
                value=f(f'{prefix}_chebi', delimiter=';', extract='chebi'),
            ),
            CV(term=Namespace.HMDB, value=f(f'{prefix}_hmdb', delimiter=';')),
            CV(
                term=Namespace.KEGG,
                value=f(f'{prefix}_kegg_compound', delimiter=';'),
            ),
            CV(
                term=Namespace.PUBCHEM,
                value=f(f'{prefix}_pubchem_compound', delimiter=';'),
            ),
            CV(
                term=Namespace.GENESYMBOL,
                value=f(f'{prefix}_hgnc', delimiter=';'),
            ),
            CV(term=Namespace.NAME, value=f(f'{prefix}_label')),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=_member_taxon_id(prefix))
        ),
        associations=AssociationsBuilder(
            AssociationBuilder(
                object_entity_type=model.Pathway,
                object_identifier_type=Namespace.WIKIPATHWAYS,
                object_identifier=f('pathway_term_accession'),
            )
        ),
    )


def _interaction_types(row):
    values = row.get('interaction_types') or []
    return set(values.split(';') if isinstance(values, str) else values)


def _interaction_direction(row):
    values = _interaction_types(row)
    if 'Stimulation' in values and 'Inhibition' not in values:
        return model.DirectionQualifierEnum.increased
    if 'Inhibition' in values and 'Stimulation' not in values:
        return model.DirectionQualifierEnum.decreased
    return None


def wikipathways_predicate(row):
    values = _interaction_types(row)
    if values & {'Stimulation', 'Inhibition', 'TranscriptionTranslation'}:
        return slots.affects
    if values & {'Binding', 'ComplexBinding'}:
        return slots.physically_interacts_with
    # Conversion and catalysis diagram edges need a reaction node to be more specific.
    return slots.interacts_with


interactions_schema = RelationBuilder(
    subject=_participant_builder('source', 'source'),
    predicate=wikipathways_predicate,
    object=_participant_builder('target', 'target'),
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.NAME, value=f('interaction_local_id'))
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.object_direction_qualifier, value=_interaction_direction),
        CV(term=slots.has_topic, value=f('pathway_version_id')),
        CV(
            term=slots.in_taxon,
            value=f(
                'taxon_id',
                transform=lambda v: 'NCBITaxon:'
                + str(v).removeprefix('NCBITaxon:')
                if str(v).removeprefix('NCBITaxon:').isdigit()
                and int(str(v).removeprefix('NCBITaxon:')) > 0
                else None,
            ),
        ),
    ),
)


download = Download(
    url=current_rdf_url,
    filename='wikipathways_rdf_wp.zip',
    subfolder='wikipathways',
    large=True,
    ext='zip',
    default_mode='rb',
)


resource = Resource(
    config,
    pathways=Dataset(
        download=download,
        mapper=pathways_schema,
        raw_parser=lambda opener, **kwargs: _raw(
            opener, data_type='pathways', **kwargs
        ),
    ),
    interactions=Dataset(
        download=download,
        mapper=interactions_schema,
        raw_parser=lambda opener, **kwargs: _raw(
            opener, data_type='interactions', **kwargs
        ),
    ),
)
