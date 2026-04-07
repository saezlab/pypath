"""
Parse Reactome BioPAX data and emit Entity records.

This module converts Reactome BioPAX data into Entity records using the
 declarative schema pattern.
"""

from __future__ import annotations

from pypath.internals.cv_terms import (
    BiologicalRoleCv,
    EntityTypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    LicenseCV,
    MoleculeAnnotationsCv,
    OntologyCv,
    ParticipantMetadataCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    Member,
    MembersFromList,
    MembershipBuilder,
)
from pypath.internals.ontology_builder import OntologyBuilder, RelationshipBuilder
from pypath.internals.ontology_schema import OntologyDocument, OntologyTypedef
from pypath.inputs_v2.base import Dataset, Download, OntologyDataset, Resource, ResourceConfig
from pypath.inputs_v2.parsers.reactome import _raw


config = ResourceConfig(
    id=ResourceCv.REACTOME,
    name='Reactome',
    url='https://reactome.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='34788843',
    primary_category='pathways',
    annotation_ontologies=(OntologyCv.GENE_ONTOLOGY, OntologyCv.REACTOME_PATHWAYS),
    description=(
        'Reactome is a free, open-source, curated and peer-reviewed '
        'pathway database.'
    ),
)


_MISSING_VALUE = '__MISSING__'

entity_type_map = {v.value: v for v in EntityTypeCv}
role_map = {
    'reactant': BiologicalRoleCv.REACTANT,
    'product': BiologicalRoleCv.PRODUCT,
    'template': BiologicalRoleCv.TEMPLATE,
    'controller': BiologicalRoleCv.CONTROLLER,
    'controlled': BiologicalRoleCv.CONTROLLED,
}

f = FieldConfig(
    delimiter=';',
    map={
        'entity_type': lambda value: (
            EntityTypeCv.PHYSICAL_ENTITY
            if not value or value == _MISSING_VALUE
            else entity_type_map.get(value, EntityTypeCv.PHYSICAL_ENTITY)
        ),
        'role': lambda value: (
            None
            if not value or value == _MISSING_VALUE
            else role_map.get(value)
        ),
        'missing': lambda value: '' if not value or value == _MISSING_VALUE else value,
        'split': lambda value: [] if not value or value == _MISSING_VALUE else [
            item for item in value.split(';') if item
        ],
    },
)


reactions_schema = EntityBuilder(
    entity_type=f('entity_type', map=entity_type_map),
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('display_name')),
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('synonyms')),
        CV(term=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=f('reactome_stable_id')),
        CV(term=IdentifierNamespaceCv.REACTOME_ID, value=f('reactome_id')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('pubmed')),
        CV(term=MoleculeAnnotationsCv.EC_NUMBER, value=f('ec_number')),
        CV(term=InteractionMetadataCv.CONVERSION_DIRECTION, value=f('direction')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('pathway_term_accession', map='split')),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=f('participant_entity_type', delimiter='||', map='entity_type'),
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.NAME, value=f('participant_display_name', delimiter='||', map='missing')),
                CV(term=IdentifierNamespaceCv.SYNONYM, value=f('participant_synonyms', delimiter='||', map='split')),
                CV(
                    term=IdentifierNamespaceCv.REACTOME_STABLE_ID,
                    value=f('participant_reactome_stable_id', delimiter='||', map='split'),
                ),
                CV(term=IdentifierNamespaceCv.UNIPROT, value=f('participant_uniprot', delimiter='||', map='split')),
                CV(term=IdentifierNamespaceCv.CHEBI, value=f('participant_chebi', delimiter='||', map='split')),
                CV(
                    term=IdentifierNamespaceCv.PUBCHEM_COMPOUND,
                    value=f('participant_pubchem_compound', delimiter='||', map='split'),
                ),
                CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('participant_kegg', delimiter='||', map='split')),
                CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('participant_go', delimiter='||', map='split')),
            ),
            annotations=AnnotationsBuilder(
                CV(term=f('participant_role', delimiter='||', map='role')),
                CV(
                    term=ParticipantMetadataCv.STOICHIOMETRY,
                    value=f('participant_stoichiometry', delimiter='||', map='missing'),
                ),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(
                    term=IdentifierNamespaceCv.NCBI_TAX_ID,
                    value=f('participant_ncbi_tax_id', delimiter='||', map='missing'),
                ),
                CV(
                    term=IdentifierNamespaceCv.CV_TERM_ACCESSION,
                    value=f('participant_pathway_term_accession', delimiter='||', map='split'),
                ),
            ),
        ),
    ),
)


controls_schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('display_name')),
        CV(term=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=f('reactome_stable_id')),
        CV(term=IdentifierNamespaceCv.REACTOME_ID, value=f('reactome_id')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('go')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=InteractionMetadataCv.CONTROL_TYPE, value=f('control_type')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('pathway_term_accession', map='split')),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=f('controller_entity_type', map='entity_type'),
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.NAME, value=f('controller_display_name', map='missing')),
                    CV(term=IdentifierNamespaceCv.SYNONYM, value=f('controller_synonyms', map='split')),
                    CV(term=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=f('controller_reactome_stable_id', map='split')),
                    CV(term=IdentifierNamespaceCv.UNIPROT, value=f('controller_uniprot', map='split')),
                    CV(term=IdentifierNamespaceCv.CHEBI, value=f('controller_chebi', map='split')),
                    CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=f('controller_pubchem_compound', map='split')),
                    CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('controller_kegg', map='split')),
                    CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('controller_go', map='split')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('controller_ncbi_tax_id', map='missing')),
                    CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('controller_pathway_term_accession', map='split')),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=BiologicalRoleCv.CONTROLLER),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=f('controlled_entity_type', map='entity_type'),
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.NAME, value=f('controlled_display_name', map='missing')),
                    CV(term=IdentifierNamespaceCv.SYNONYM, value=f('controlled_synonyms', map='split')),
                    CV(term=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=f('controlled_reactome_stable_id', map='split')),
                    CV(term=IdentifierNamespaceCv.UNIPROT, value=f('controlled_uniprot', map='split')),
                    CV(term=IdentifierNamespaceCv.CHEBI, value=f('controlled_chebi', map='split')),
                    CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=f('controlled_pubchem_compound', map='split')),
                    CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('controlled_kegg', map='split')),
                    CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('controlled_go', map='split')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('controlled_ncbi_tax_id', map='missing')),
                    CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('controlled_pathway_term_accession', map='split')),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=BiologicalRoleCv.CONTROLLED),
            ),
        ),
    ),
)


control_groups_schema = EntityBuilder(
    entity_type=f('controller_entity_type', map='entity_type'),
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('controller_display_name', map='missing')),
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('controller_synonyms', map='split')),
        CV(term=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=f('controller_reactome_stable_id', map='split')),
        CV(term=IdentifierNamespaceCv.UNIPROT, value=f('controller_uniprot', map='split')),
        CV(term=IdentifierNamespaceCv.CHEBI, value=f('controller_chebi', map='split')),
        CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=f('controller_pubchem_compound', map='split')),
        CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('controller_kegg', map='split')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('controller_go', map='split')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('controller_ncbi_tax_id', map='missing')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('controller_pathway_term_accession', map='split')),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=f('controller_member_entity_type', delimiter='||', map='entity_type'),
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.NAME, value=f('controller_member_display_name', delimiter='||', map='missing')),
                CV(term=IdentifierNamespaceCv.SYNONYM, value=f('controller_member_synonyms', delimiter='||', map='split')),
                CV(term=IdentifierNamespaceCv.REACTOME_STABLE_ID, value=f('controller_member_reactome_stable_id', delimiter='||', map='split')),
                CV(term=IdentifierNamespaceCv.UNIPROT, value=f('controller_member_uniprot', delimiter='||', map='split')),
                CV(term=IdentifierNamespaceCv.CHEBI, value=f('controller_member_chebi', delimiter='||', map='split')),
                CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=f('controller_member_pubchem_compound', delimiter='||', map='split')),
                CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('controller_member_kegg', delimiter='||', map='split')),
                CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('controller_member_go', delimiter='||', map='split')),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('controller_member_ncbi_tax_id', delimiter='||', map='missing')),
                CV(
                    term=IdentifierNamespaceCv.CV_TERM_ACCESSION,
                    value=f('controller_member_pathway_term_accession', delimiter='||', map='split'),
                ),
            ),
        ),
    ),
)


pathway_ontology_schema = OntologyBuilder(
    id='id',
    name='name',
    definition='definition',
    synonyms=f('synonyms', delimiter=';'),
    comments=f('comments', delimiter=';'),
    xrefs=f('xrefs', delimiter=';'),
    relationships=[
        RelationshipBuilder(
            type='part_of',
            target=f('parent_ids', delimiter=';'),
            target_name=f('parent_names', delimiter=';'),
        ),
    ],
)

pathway_ontology_document = OntologyDocument(
    ontology='reactome_pathways',
    default_namespace='reactome_pathways',
    remark='Reactome pathway ontology exported from Reactome BioPAX via pypath.',
    typedefs=[OntologyTypedef(id='part_of', name='part_of')],
)


download = Download(
    url='https://reactome.org/download/current/biopax.zip',
    filename='reactome_biopax.zip',
    subfolder='reactome',
    large=True,
    ext='zip',
    default_mode='rb',
)

resource = Resource(
    config,
    reactions=Dataset(
        download=download,
        mapper=reactions_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='reactions', **kwargs),
    ),
    controls=Dataset(
        download=download,
        mapper=controls_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='controls', **kwargs),
    ),
    control_groups=Dataset(
        download=download,
        mapper=control_groups_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='control_groups', **kwargs),
    ),
    pathway_ontology=OntologyDataset(
        download=download,
        mapper=pathway_ontology_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='pathway_terms', **kwargs),
        document=pathway_ontology_document,
        extension='obo',
        file_stem='reactome_pathways',
    ),
)
