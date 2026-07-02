"""Parse ChEBI OBO data into molecule entities with ontology hierarchy edges."""

from __future__ import annotations

from collections.abc import Iterable
import re

from pypath.inputs_v2.base import (
    Dataset,
    Download,
    Resource,
    ResourceConfig,
)
from pypath.inputs_v2.parsers.chebi import _raw
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    MoleculeAnnotationsCv,
    OntologyAnnotationCv,
    OntologyCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.silver_schema import EntityRef, OntologyRelation
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)

config = ResourceConfig(
    id=ResourceCv.CHEBI,
    name='ChEBI',
    url='https://www.ebi.ac.uk/chebi/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='41312627',
    primary_category='small_molecules',
    annotation_ontologies=(OntologyCv.CHEBI,),
    description=(
        'ChEBI is a manually curated database and ontology of small molecular '
        'entities. This inputs_v2 module emits searchable small-molecule '
        'entities with ontology hierarchy edges derived from the official OBO '
        'release.'
    ),
)

f = FieldConfig(
    extract={
        'chebi': r'^(?:CHEBI:)?(\d+)$',
    },
)


molecules_schema = EntityBuilder(
    entity_type=EntityTypeCv.CHEMICAL,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.CHEBI, value=f('chebi_id', extract='chebi')),
        CV(
            term=IdentifierNamespaceCv.CHEBI,
            value=lambda row: [
                value
                for value in (
                    _chebi_identifier(alt_id)
                    for alt_id in row.get('alt_ids', [])
                )
                if value
            ],
        ),
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
        CV(term=OntologyAnnotationCv.DEFINITION, value=f('definition')),
        CV(term=MoleculeAnnotationsCv.MASS_DALTON, value=f('mass')),
        CV(term=MoleculeAnnotationsCv.MW_MONOISOTOPIC, value=f('monoisotopic_mass')),
        CV(term=MoleculeAnnotationsCv.MOLECULAR_CHARGE, value=f('charge')),
    ),
    ontology_relations=lambda row: _ontology_relations(row),
)


def _chebi_identifier(value: object) -> str | None:
    match = re.fullmatch(r'(?:CHEBI:)?(\d+)', str(value or '').strip())
    return match.group(1) if match else None


def _ontology_relations(row: dict) -> list[OntologyRelation]:
    relations: list[OntologyRelation] = []
    seen: set[tuple[str, str]] = set()
    for parent in row.get('is_a') or []:
        parent_id = _chebi_identifier(parent)
        if not parent_id:
            continue
        key = ('is_a', parent_id)
        if key in seen:
            continue
        seen.add(key)
        relations.append(
            OntologyRelation(
                predicate='is_a',
                object=EntityRef(
                    type=EntityTypeCv.CHEMICAL,
                    identifier_type=IdentifierNamespaceCv.CHEBI,
                    identifier=parent_id,
                ),
                ontology_id='chebi',
            )
        )
    for relationship in row.get('relationships') or []:
        predicate = relationship.get('type')
        target_id = _chebi_identifier(relationship.get('target'))
        if not predicate or not target_id:
            continue
        key = (predicate, target_id)
        if key in seen:
            continue
        seen.add(key)
        relations.append(
            OntologyRelation(
                predicate=predicate,
                object=EntityRef(
                    type=EntityTypeCv.CHEMICAL,
                    identifier_type=IdentifierNamespaceCv.CHEBI,
                    identifier=target_id,
                ),
                ontology_id='chebi',
            )
        )
    return relations


def _id_translation_rows(row: dict) -> Iterable[dict]:
    chebi_match = re.fullmatch(
        r'(?:CHEBI:)?(\d+)',
        str(row.get('chebi_id') or '').strip(),
    )
    chebi_id = chebi_match.group(1) if chebi_match else None
    standard_inchi = row.get('inchi')
    if not chebi_id or not standard_inchi:
        return

    seen = set()
    for raw_id in (row.get('chebi_id'), *(row.get('alt_ids') or [])):
        match = re.fullmatch(r'(?:CHEBI:)?(\d+)', str(raw_id or '').strip())
        if not match:
            continue
        key_value = match.group(1)
        if key_value in seen:
            continue
        seen.add(key_value)
        yield {
            'source': 'chebi',
            'key_type': 'MI:0474:Chebi',
            'key_value': key_value,
            'standard_inchi': standard_inchi,
        }


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
            id_row
            for raw_row in _raw(opener, data_type='molecules', **kwargs)
            for id_row in _id_translation_rows(raw_row)
        ),
        kind='id_translation',
    ),
)
