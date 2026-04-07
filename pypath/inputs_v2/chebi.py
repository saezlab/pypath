"""Parse ChEBI OBO data into molecule entities and a standalone ontology export."""

from __future__ import annotations

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    MoleculeAnnotationsCv,
    OntologyCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.ontology_schema import OntologyDocument, OntologyRelationship, OntologyTerm
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, OntologyDataset, Resource, ResourceConfig
from pypath.inputs_v2.parsers.chebi import _raw


config = ResourceConfig(
    id=ResourceCv.CHEBI,
    name='ChEBI',
    url='https://www.ebi.ac.uk/chebi/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='17932068',
    primary_category='small_molecules',
    annotation_ontologies=(OntologyCv.CHEBI,),
    description=(
        'ChEBI is a manually curated database and ontology of small molecular '
        'entities. This inputs_v2 module emits searchable small-molecule '
        'entities plus a standalone ChEBI ontology export derived from the '
        'official OBO release.'
    ),
)

f = FieldConfig()


molecules_schema = EntityBuilder(
    entity_type=EntityTypeCv.SMALL_MOLECULE,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.CHEBI, value=f('chebi_id')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('name')),
        CV(term=IdentifierNamespaceCv.SYNONYM, value=lambda row: row.get('synonyms', [])),
        CV(term=IdentifierNamespaceCv.SMILES, value=f('smiles')),
        CV(term=IdentifierNamespaceCv.STANDARD_INCHI, value=f('inchi')),
        CV(term=IdentifierNamespaceCv.STANDARD_INCHI_KEY, value=f('inchikey')),
        CV(term=IdentifierNamespaceCv.MOLECULAR_FORMULA, value=f('formula')),
        CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=lambda row: row.get('kegg_compound', [])),
        CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=lambda row: row.get('pubchem_compound', [])),
        CV(term=IdentifierNamespaceCv.HMDB, value=lambda row: row.get('hmdb', [])),
        CV(term=IdentifierNamespaceCv.LIPIDMAPS, value=lambda row: row.get('lipidmaps', [])),
        CV(term=IdentifierNamespaceCv.CAS, value=lambda row: row.get('cas', [])),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.DESCRIPTION, value=f('definition')),
        CV(term=MoleculeAnnotationsCv.MASS_DALTON, value=f('mass')),
        CV(term=MoleculeAnnotationsCv.MW_MONOISOTOPIC, value=f('monoisotopic_mass')),
        CV(term=MoleculeAnnotationsCv.MOLECULAR_CHARGE, value=f('charge')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=lambda row: row.get('ancestor_terms', [])),
    ),
)


pathway_ontology_document = OntologyDocument(
    ontology='chebi',
    default_namespace='chebi',
    remark='ChEBI ontology exported from chebi.obo.gz via pypath.',
)


def _id_translation_mapper(row: dict) -> dict | None:
    chebi_id = row.get('chebi_id')
    standard_inchi = row.get('inchi')
    if not chebi_id or not standard_inchi:
        return None
    return {
        'source': 'chebi',
        'key_type': 'MI:0474:Chebi',
        'key_value': chebi_id,
        'standard_inchi': standard_inchi,
    }


def _ontology_mapper(row: dict) -> OntologyTerm | None:
    term_id = row.get('id')
    name = row.get('name')

    if not term_id or not name:
        return None

    relationships = [
        OntologyRelationship(
            type=relationship['type'],
            target=relationship['target'],
            target_name=relationship.get('target_name') or None,
        )
        for relationship in row.get('relationships', [])
        if relationship.get('type') and relationship.get('target')
    ]

    return OntologyTerm(
        id=term_id,
        name=name,
        definition=row.get('definition') or None,
        synonyms=row.get('synonyms') or None,
        comments=row.get('comments') or None,
        xrefs=row.get('xrefs') or None,
        is_a=row.get('is_a') or None,
        relationships=relationships or None,
        is_obsolete=bool(row.get('is_obsolete')),
    )


download = Download(
    url='https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.gz',
    filename='chebi.obo.gz',
    subfolder='chebi',
    large=True,
    ext='gz',
)


resource = Resource(
    config,
    molecules=Dataset(
        download=download,
        mapper=molecules_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='molecules', **kwargs),
    ),
    id_translation=Dataset(
        download=download,
        mapper=lambda row: row,
        raw_parser=lambda opener, **kwargs: (
            row
            for raw_row in _raw(opener, data_type='molecules', **kwargs)
            if (row := _id_translation_mapper(raw_row)) is not None
        ),
        kind='id_translation',
    ),
    ontology=OntologyDataset(
        download=download,
        mapper=_ontology_mapper,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='ontology_terms', **kwargs),
        document=pathway_ontology_document,
        extension='obo',
        file_stem='chebi',
    ),
)
