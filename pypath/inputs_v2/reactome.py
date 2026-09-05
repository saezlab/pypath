"""Reactome BioPAX pathways, activities and explicit regulatory controls.

Parsers preserve source type tokens; native classes and predicates are chosen
here. Degradation controls target the process rather than an inferred substrate
quantity change.
"""

from __future__ import annotations
from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from functools import partial
import re
from pypath.internals.cv_terms import (
    LicenseCV,
    OntologyCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.silver_schema import EntityRef, OntologyRelation
from pypath.internals.tabular_builder import (
    AssociationBuilder,
    AssociationsBuilder,
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembersFromList,
    MembershipBuilder,
    RelationBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.reactome import _raw

config = ResourceConfig(
    id=ResourceCv.REACTOME,
    name='Reactome',
    url='https://reactome.org/',
    license=LicenseCV.CC0_1_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='34788843',
    primary_category='pathways',
    annotation_ontologies=(
        OntologyCv.GENE_ONTOLOGY,
        OntologyCv.REACTOME_PATHWAYS,
    ),
    description='Reactome is a free, open-source, curated and peer-reviewed pathway database.',
)
_MISSING_VALUE = '__MISSING__'
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
role_map = {
    'reactant': slots.has_input,
    'product': slots.has_output,
    'template': slots.has_input,
}


def _control_effect(row):
    return {
        'ACTIVATION': model.DirectionQualifierEnum.increased,
        'INHIBITION': model.DirectionQualifierEnum.decreased,
    }.get(row.get('control_type'))


_TAXON_SCOPED_ENTITY_TYPES = {
    model.Protein,
    model.Gene,
    model.RNAProduct,
    model.NucleicAcidEntity,
}
_REACTOME_PATHWAY_ONTOLOGY_ID = 'reactome_pathways'
f = FieldConfig(
    delimiter=';',
    map={
        'entity_type': lambda value: model.PhysicalEntity
        if not value or value == _MISSING_VALUE
        else entity_type_map.get(value, model.PhysicalEntity),
        'role': lambda value: None
        if not value or value == _MISSING_VALUE
        else role_map.get(value),
        'missing': lambda value: ''
        if not value or value == _MISSING_VALUE
        else value,
        'split': lambda value: []
        if not value or value == _MISSING_VALUE
        else [item for item in value.split(';') if item],
        'split_chebi': lambda value: []
        if not value or value == _MISSING_VALUE
        else [
            match.group(1)
            for item in value.split(';')
            if (match := re.fullmatch('(?:CHEBI:)?(\\d+)', item.strip()))
        ],
    },
)


def _participant_taxon_ids(row):
    entity_types = f(
        'participant_entity_type',
        delimiter='||',
        map='entity_type',
        preserve_indices=True,
    ).extract(row)
    taxon_ids = f(
        'participant_ncbi_tax_id',
        delimiter='||',
        map='missing',
        preserve_indices=True,
    ).extract(row)
    size = max(len(entity_types), len(taxon_ids))
    values = []
    for idx in range(size):
        entity_type = entity_types[idx] if idx < len(entity_types) else None
        taxon_id = taxon_ids[idx] if idx < len(taxon_ids) else None
        values.append(
            next(iter(_taxon_id_for_entity_type(entity_type, taxon_id)), None)
        )
    return values


def _taxon_id_for_entity_type(entity_type, taxon_id):
    mapped = (
        entity_type
        if isinstance(entity_type, type)
        else entity_type_map.get(entity_type, model.PhysicalEntity)
    )
    values = taxon_id if isinstance(taxon_id, list) else [taxon_id]
    return [
        'NCBITaxon:' + str(v).removeprefix('NCBITaxon:')
        for v in values
        if mapped in _TAXON_SCOPED_ENTITY_TYPES
        and v
        and str(v).removeprefix('NCBITaxon:').isdigit()
    ]


def _entity_taxon_id(entity_type_field, taxon_field):
    def extractor(row):
        entity_type = row.get(entity_type_field)
        taxon_id = f(taxon_field, map='missing').extract(row)
        return _taxon_id_for_entity_type(entity_type, taxon_id)

    return extractor


def _controller_member_taxon_ids(row):
    entity_types = f(
        'controller_member_entity_type',
        delimiter='||',
        map='entity_type',
        preserve_indices=True,
    ).extract(row)
    taxon_ids = f(
        'controller_member_ncbi_tax_id',
        delimiter='||',
        map='missing',
        preserve_indices=True,
    ).extract(row)
    size = max(len(entity_types), len(taxon_ids))
    values = []
    for idx in range(size):
        entity_type = entity_types[idx] if idx < len(entity_types) else None
        taxon_id = taxon_ids[idx] if idx < len(taxon_ids) else None
        values.append(
            next(iter(_taxon_id_for_entity_type(entity_type, taxon_id)), None)
        )
    return values


def _split_values(value: object, *, delimiter: str = ';') -> list[str]:
    if not value or value == _MISSING_VALUE:
        return []
    return [
        item
        for item in (part.strip() for part in str(value).split(delimiter))
        if item and item != _MISSING_VALUE
    ]


def _pathway_ref(stable_id: str) -> EntityRef:
    return EntityRef(
        type=model.Pathway,
        identifier_type=Namespace.REACTOME,
        identifier=stable_id,
    )


def _pathway_association(object_identifier) -> AssociationBuilder:
    return AssociationBuilder(
        object_entity_type=model.Pathway,
        object_identifier_type=Namespace.REACTOME,
        object_identifier=object_identifier,
        predicate=slots.part_of,
    )


def _pathway_associations(object_identifier) -> AssociationsBuilder:
    return AssociationsBuilder(_pathway_association(object_identifier))


def _cv_term_association(object_identifier) -> AssociationBuilder:
    return AssociationBuilder(
        object_entity_type=model.OntologyClass,
        object_identifier_type=Namespace.GO,
        object_identifier=object_identifier,
    )


def _cv_term_associations(object_identifier) -> AssociationsBuilder:
    return AssociationsBuilder(_cv_term_association(object_identifier))


def _combined_associations(
    *associations: AssociationBuilder,
) -> AssociationsBuilder:
    return AssociationsBuilder(*associations)


def _pathway_ontology_relations(row: dict) -> list[OntologyRelation]:
    relations: list[OntologyRelation] = []
    seen: set[tuple[str, str]] = set()
    for parent_group in _split_values(
        row.get('parent_pathway_reactome_stable_id'), delimiter='||'
    ):
        for parent_id in _split_values(parent_group):
            key = ('part_of', parent_id)
            if key in seen:
                continue
            seen.add(key)
            relations.append(
                OntologyRelation(
                    predicate='part_of',
                    object=_pathway_ref(parent_id),
                    ontology_id=_REACTOME_PATHWAY_ONTOLOGY_ID,
                )
            )
    return relations


reactions_schema = EntityBuilder(
    entity_type=f('entity_type', map=entity_type_map),
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.REACTOME, value=f('reactome_stable_id')),
        CV(term=Namespace.REACTOME_ID, value=f('reactome_id')),
        CV(term=Namespace.NAME, value=f('display_name')),
        CV(term=Namespace.SYNONYM, value=f('synonyms')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.publications,
            value=f(
                'pubmed',
                transform=lambda v: 'PMID:' + str(v).removeprefix('PMID:')
                if str(v).removeprefix('PMID:').isdigit()
                and int(str(v).removeprefix('PMID:')) > 0
                else None,
            ),
        ),
        CV(
            term=slots.has_topic,
            value=lambda row, source=f('ec_number'): [
                'EC:' + str(v).removeprefix('EC:') for v in source.extract(row)
            ],
        ),
    ),
    associations=_pathway_associations(
        f('pathway_term_accession', map='split')
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=f(
                'participant_entity_type', delimiter='||', map='entity_type'
            ),
            predicate=f(
                'participant_role',
                delimiter='||',
                map=lambda v: role_map.get(v),
                preserve_indices=True,
            ),
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.REACTOME,
                    value=f(
                        'participant_reactome_stable_id',
                        delimiter='||',
                        map='split',
                    ),
                ),
                CV(
                    term=Namespace.UNIPROT,
                    value=f('participant_uniprot', delimiter='||', map='split'),
                ),
                CV(
                    term=Namespace.CHEBI,
                    value=f(
                        'participant_chebi', delimiter='||', map='split_chebi'
                    ),
                ),
                CV(
                    term=Namespace.PUBCHEM,
                    value=f(
                        'participant_pubchem_compound',
                        delimiter='||',
                        map='split',
                    ),
                ),
                CV(
                    term=Namespace.KEGG,
                    value=f('participant_kegg', delimiter='||', map='split'),
                ),
                CV(
                    term=Namespace.NAME,
                    value=f(
                        'participant_display_name',
                        delimiter='||',
                        map='missing',
                    ),
                ),
                CV(
                    term=Namespace.SYNONYM,
                    value=f(
                        'participant_synonyms', delimiter='||', map='split'
                    ),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(
                    term=slots.stoichiometry,
                    value=f(
                        'participant_stoichiometry',
                        delimiter='||',
                        map='missing',
                    ),
                )
            ),
            entity_annotations=AnnotationsBuilder(
                CV(term=slots.in_taxon, value=_participant_taxon_ids)
            ),
            entity_associations=_combined_associations(
                _pathway_association(
                    f(
                        'participant_pathway_term_accession',
                        delimiter='||',
                        map='split',
                    )
                ),
                _cv_term_association(
                    f('participant_go', delimiter='||', map='split')
                ),
            ),
        )
    ),
)


def reactome_predicate(row):
    return (
        slots.catalyzes
        if row.get('control_class') == 'Catalysis'
        else slots.regulates
    )


controller_builder = EntityBuilder(
    entity_type=f('controller_entity_type', map='entity_type'),
    identifiers=IdentifiersBuilder(
        CV(
            term=Namespace.REACTOME,
            value=f('controller_reactome_stable_id', map='split'),
        ),
        CV(term=Namespace.UNIPROT, value=f('controller_uniprot', map='split')),
        CV(
            term=Namespace.CHEBI, value=f('controller_chebi', map='split_chebi')
        ),
        CV(
            term=Namespace.PUBCHEM,
            value=f('controller_pubchem_compound', map='split'),
        ),
        CV(term=Namespace.KEGG, value=f('controller_kegg', map='split')),
        CV(
            term=Namespace.NAME,
            value=f('controller_display_name', map='missing'),
        ),
        CV(term=Namespace.SYNONYM, value=f('controller_synonyms', map='split')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.in_taxon,
            value=_entity_taxon_id(
                'controller_entity_type', 'controller_ncbi_tax_id'
            ),
        )
    ),
    associations=_combined_associations(
        _pathway_association(
            f('controller_pathway_term_accession', map='split')
        ),
        _cv_term_association(f('controller_go', map='split')),
    ),
)
controlled_builder = EntityBuilder(
    entity_type=f('controlled_entity_type', map='entity_type'),
    identifiers=IdentifiersBuilder(
        CV(
            term=Namespace.REACTOME,
            value=f('controlled_reactome_stable_id', map='split'),
        ),
        CV(term=Namespace.UNIPROT, value=f('controlled_uniprot', map='split')),
        CV(
            term=Namespace.CHEBI, value=f('controlled_chebi', map='split_chebi')
        ),
        CV(
            term=Namespace.PUBCHEM,
            value=f('controlled_pubchem_compound', map='split'),
        ),
        CV(term=Namespace.KEGG, value=f('controlled_kegg', map='split')),
        CV(
            term=Namespace.NAME,
            value=f('controlled_display_name', map='missing'),
        ),
        CV(term=Namespace.SYNONYM, value=f('controlled_synonyms', map='split')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.in_taxon,
            value=_entity_taxon_id(
                'controlled_entity_type', 'controlled_ncbi_tax_id'
            ),
        )
    ),
    associations=_combined_associations(
        _pathway_association(
            f('controlled_pathway_term_accession', map='split')
        ),
        _cv_term_association(f('controlled_go', map='split')),
    ),
)
controls_schema = RelationBuilder(
    subject=controller_builder,
    predicate=reactome_predicate,
    object=controlled_builder,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.REACTOME, value=f('reactome_stable_id')),
        CV(term=Namespace.REACTOME_ID, value=f('reactome_id')),
        CV(term=Namespace.NAME, value=f('display_name')),
    ),
    annotations=AnnotationsBuilder(),
)
pathways_schema = EntityBuilder(
    entity_type=model.Pathway,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.REACTOME, value=f('reactome_stable_id')),
        CV(term=Namespace.REACTOME_ID, value=f('reactome_id')),
        CV(term=Namespace.NAME, value=f('display_name')),
        CV(term=Namespace.SYNONYM, value=f('synonyms')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.publications,
            value=f(
                'pubmed',
                transform=lambda v: 'PMID:' + str(v).removeprefix('PMID:')
                if str(v).removeprefix('PMID:').isdigit()
                and int(str(v).removeprefix('PMID:')) > 0
                else None,
            ),
        ),
        CV(
            term=slots.in_taxon,
            value=f(
                'ncbi_tax_id',
                transform=lambda v: 'NCBITaxon:'
                + str(v).removeprefix('NCBITaxon:')
                if str(v).removeprefix('NCBITaxon:').isdigit()
                and int(str(v).removeprefix('NCBITaxon:')) > 0
                else None,
            ),
        ),
        CV(term=slots.description, value=f('definition')),
        CV(term=slots.description, value=f('comments')),
    ),
    associations=_cv_term_associations(f('go')),
    ontology_relations=_pathway_ontology_relations,
)
control_groups_schema = EntityBuilder(
    entity_type=f('controller_entity_type', map='entity_type'),
    identifiers=IdentifiersBuilder(
        CV(
            term=Namespace.REACTOME,
            value=f('controller_reactome_stable_id', map='split'),
        ),
        CV(term=Namespace.UNIPROT, value=f('controller_uniprot', map='split')),
        CV(
            term=Namespace.CHEBI, value=f('controller_chebi', map='split_chebi')
        ),
        CV(
            term=Namespace.PUBCHEM,
            value=f('controller_pubchem_compound', map='split'),
        ),
        CV(term=Namespace.KEGG, value=f('controller_kegg', map='split')),
        CV(
            term=Namespace.NAME,
            value=f('controller_display_name', map='missing'),
        ),
        CV(term=Namespace.SYNONYM, value=f('controller_synonyms', map='split')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.in_taxon,
            value=_entity_taxon_id(
                'controller_entity_type', 'controller_ncbi_tax_id'
            ),
        )
    ),
    associations=_combined_associations(
        _pathway_association(
            f('controller_pathway_term_accession', map='split')
        ),
        _cv_term_association(f('controller_go', map='split')),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=f(
                'controller_member_entity_type',
                delimiter='||',
                map='entity_type',
            ),
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.REACTOME,
                    value=f(
                        'controller_member_reactome_stable_id',
                        delimiter='||',
                        map='split',
                    ),
                ),
                CV(
                    term=Namespace.UNIPROT,
                    value=f(
                        'controller_member_uniprot', delimiter='||', map='split'
                    ),
                ),
                CV(
                    term=Namespace.CHEBI,
                    value=f(
                        'controller_member_chebi',
                        delimiter='||',
                        map='split_chebi',
                    ),
                ),
                CV(
                    term=Namespace.PUBCHEM,
                    value=f(
                        'controller_member_pubchem_compound',
                        delimiter='||',
                        map='split',
                    ),
                ),
                CV(
                    term=Namespace.KEGG,
                    value=f(
                        'controller_member_kegg', delimiter='||', map='split'
                    ),
                ),
                CV(
                    term=Namespace.NAME,
                    value=f(
                        'controller_member_display_name',
                        delimiter='||',
                        map='missing',
                    ),
                ),
                CV(
                    term=Namespace.SYNONYM,
                    value=f(
                        'controller_member_synonyms',
                        delimiter='||',
                        map='split',
                    ),
                ),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(term=slots.in_taxon, value=_controller_member_taxon_ids)
            ),
            entity_associations=_combined_associations(
                _pathway_association(
                    f(
                        'controller_member_pathway_term_accession',
                        delimiter='||',
                        map='split',
                    )
                ),
                _cv_term_association(
                    f('controller_member_go', delimiter='||', map='split')
                ),
            ),
        )
    ),
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
        raw_parser=partial(_raw, data_type='reactions'),
    ),
    pathways=Dataset(
        download=download,
        mapper=pathways_schema,
        raw_parser=partial(_raw, data_type='pathways'),
    ),
    controls=Dataset(
        download=download,
        mapper=controls_schema,
        raw_parser=partial(_raw, data_type='controls'),
    ),
    control_groups=Dataset(
        download=download,
        mapper=control_groups_schema,
        raw_parser=partial(_raw, data_type='control_groups'),
    ),
)
