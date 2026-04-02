"""
Parse SwissLipids data and emit Entity records.

This module converts SwissLipids lipid data into Entity records using
the schema defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    MoleculeAnnotationsCv,
    ResourceCv,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_tsv


config = ResourceConfig(
    id=ResourceCv.SWISSLIPIDS,
    name='SwissLipids',
    url='https://www.swisslipids.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='25943471',
    description=(
        'SwissLipids is a curated resource providing a framework for the '
        'annotation of mass spectrometry data. It provides over 750,000 lipid '
        'structures with expert curation of lipid classes and nomenclature, '
        'hierarchical organization, cross-references to other databases (ChEBI, '
        'LIPID MAPS, HMDB), and integration with mass spectrometry tools. '
        'The database covers all major lipid categories including fatty acyls, '
        'glycerolipids, glycerophospholipids, sphingolipids, sterol lipids, '
        'prenol lipids, saccharolipids, and polyketides.'
    ),
)

f = FieldConfig(
    extract={
        'chebi': r'^(?:CHEBI:)?(\d+)$',
    },
    transform={
        'inchi': lambda v: None if not v or str(v).strip().lower() in {'none', 'inchi=none'} else str(v).strip(),
        'inchikey': lambda v: None if not v or str(v).strip().lower() in {'none', 'inchikey=none'} else str(v).strip().removeprefix('InChIKey='),
        'chebi': lambda v: f'CHEBI:{v}' if v else None,
    },
)

lipids_schema = EntityBuilder(
    entity_type=EntityTypeCv.LIPID,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.SWISSLIPIDS, value=f('Lipid ID')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('Name')),
        CV(
            term=IdentifierNamespaceCv.STANDARD_INCHI_KEY,
            value=f('InChI key (pH7.3)', transform='inchikey'),
        ),
        CV(
            term=IdentifierNamespaceCv.STANDARD_INCHI,
            value=f('InChI (pH7.3)', transform='inchi'),
        ),
        CV(term=IdentifierNamespaceCv.SMILES, value=f('SMILES (pH7.3)')),
        CV(
            term=IdentifierNamespaceCv.CHEBI,
            value=f('CHEBI', extract='chebi', transform='chebi'),
        ),
        CV(term=IdentifierNamespaceCv.LIPIDMAPS, value=f('LIPID MAPS')),
        CV(term=IdentifierNamespaceCv.HMDB, value=f('HMDB')),
        CV(term=IdentifierNamespaceCv.METANETX, value=f('MetaNetX')),
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('Synonyms*', delimiter=';')),
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('Abbreviation*')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.LIPID_HIERARCHY_LEVEL, value=f('Level')),
        CV(term=MoleculeAnnotationsCv.LIPID_MAIN_CLASS, value=f('Lipid class*')),
        CV(term=IdentifierNamespaceCv.SWISSLIPIDS, value=f('Parent')),
        CV(term=MoleculeAnnotationsCv.LIPID_STRUCTURAL_COMPONENTS, value=f('Components*')),
        CV(term=MoleculeAnnotationsCv.MOLECULAR_CHARGE, value=f('Charge (pH7.3)')),
        CV(term=MoleculeAnnotationsCv.MASS_DALTON, value=f('Exact Mass (neutral form)')),
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('PMID', delimiter='|')),
    ),
)

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
        'key_type': 'OM:0009:Swisslipids',
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
        download=download,
        mapper=lipids_schema,
        raw_parser=iter_tsv,
    ),
    id_translation=Dataset(
        download=download,
        mapper=lambda row: row,
        raw_parser=_id_translation_raw,
        kind='id_translation',
    ),
)
