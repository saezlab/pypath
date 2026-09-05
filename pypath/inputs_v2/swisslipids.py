"""
Parse SwissLipids data and emit Entity records.

This module converts SwissLipids lipid data into Entity records using
the schema defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from biolink_model.datamodel.model import ChemicalEntity, slots
from omnipath_core.naming import Namespace

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_tsv
from pypath.internals.cv_terms import (
    LicenseCV,
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
    id=ResourceCv.SWISSLIPIDS,
    name='SwissLipids',
    url='https://www.swisslipids.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='25943471',
    primary_category='lipids',
    annotation_ontologies=(OntologyCv.SWISSLIPIDS,),
    description='SwissLipids is a curated resource providing a framework for the annotation of mass spectrometry data. It provides over 750,000 lipid structures with expert curation of lipid classes and nomenclature, hierarchical organization, cross-references to other databases (ChEBI, LIPID MAPS, HMDB), and integration with mass spectrometry tools. The database covers all major lipid categories including fatty acyls, glycerolipids, glycerophospholipids, sphingolipids, sterol lipids, prenol lipids, saccharolipids, and polyketides.',
)
f = FieldConfig(
    extract={'chebi': '^(?:CHEBI:)?(\\d+)$'},
    transform={
        'inchi': lambda v: None
        if not v or str(v).strip().lower() in {'none', 'inchi=none'}
        else str(v).strip(),
        'inchikey': lambda v: None
        if not v or str(v).strip().lower() in {'none', 'inchikey=none'}
        else str(v).strip().removeprefix('InChIKey='),
    },
)
lipids_schema = EntityBuilder(
    entity_type=ChemicalEntity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.SWISSLIPIDS, value=f('Lipid ID')),
        CV(
            term=Namespace.INCHIKEY,
            value=f('InChI key (pH7.3)', transform='inchikey'),
        ),
        CV(term=Namespace.INCHI, value=f('InChI (pH7.3)', transform='inchi')),
        CV(term=Namespace.SMILES, value=f('SMILES (pH7.3)')),
        CV(term=Namespace.CHEBI, value=f('CHEBI', delimiter='|', extract='chebi')),
        CV(term=Namespace.LIPIDMAPS, value=f('LIPID MAPS', delimiter='|')),
        CV(term=Namespace.HMDB, value=f('HMDB', delimiter='|')),
        CV(term=Namespace.METANETX, value=f('MetaNetX', delimiter='|')),
        CV(term=Namespace.NAME, value=f('Name')),
        CV(term=Namespace.SYNONYM, value=f('Synonyms*', delimiter=';')),
        CV(term=Namespace.SYNONYM, value=f('Abbreviation*')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.has_chemical_formula, value=f('Formula (pH7.3)')),
        CV(term='chemrof:charge', value=f('Charge (pH7.3)')),
        CV(
            term='chemrof:monoisotopic_mass',
            value=f('Exact Mass (neutral form)'),
        ),
        CV(
            term=slots.publications,
            value=f(
                'PMID',
                delimiter='|',
                transform=lambda v: 'PMID:' + str(v).removeprefix('PMID:'),
            ),
        ),
    ),
    ontology_relations=lambda row: _ontology_relations(row),
)


def _swisslipids_identifier(value: object) -> str | None:
    identifier = str(value or '').strip()
    if not identifier or identifier.lower() in {'none', 'nan'}:
        return None
    return identifier


def _ontology_relations(row: dict) -> list[OntologyRelation]:
    parent_id = _swisslipids_identifier(row.get('Parent'))
    lipid_id = _swisslipids_identifier(row.get('Lipid ID'))
    if not parent_id or parent_id == lipid_id:
        return []
    return [
        OntologyRelation(
            predicate=slots.subclass_of,
            object=EntityRef(
                type=ChemicalEntity,
                identifier_type=Namespace.SWISSLIPIDS,
                identifier=parent_id,
            ),
            ontology_id='swisslipids',
        )
    ]


download = Download(
    url='https://swisslipids.org/api/file.php?cas=download_files&file=lipids.tsv',
    filename='lipids.tsv.gz',
    subfolder='swisslipids',
    encoding='latin-1',
)


def _id_translation_row(row: dict) -> dict | None:
    swisslipids_id = row.get('Lipid ID')
    standard_inchi = f.transform['inchi'](row.get('InChI (pH7.3)'))
    if not swisslipids_id or not standard_inchi:
        return None
    return {
        'source': 'swisslipids',
        'key_type': Namespace.SWISSLIPIDS,
        'key_value': swisslipids_id,
        'standard_inchi': standard_inchi,
    }


def _id_translation_raw(opener, max_records: int | None = None, **kwargs):
    emitted = 0
    for raw_row in iter_tsv(opener, **kwargs):
        row = _id_translation_row(raw_row)
        if row is None:
            continue
        yield row
        emitted += 1
        if max_records is not None and emitted >= max_records:
            break


resource = Resource(
    config,
    lipids=Dataset(
        download=download, mapper=lipids_schema, raw_parser=iter_tsv
    ),
    id_translation=Dataset(
        download=download,
        mapper=lambda row: row,
        raw_parser=_id_translation_raw,
        kind='id_translation',
    ),
)
