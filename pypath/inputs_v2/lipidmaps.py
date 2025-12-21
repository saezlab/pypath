"""
Parse LIPID MAPS Structure Database (LMSD) data and emit Entity records.

This module downloads the LIPID MAPS Structure Database in SDF format and
converts lipid structures into Entity records using the schema defined in
pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    MoleculeAnnotationsCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.formats.sdf import SdfReader
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig


def _parse_lipids_sdf(opener, **_kwargs: object) -> Generator[dict[str, str], None, None]:
    if not opener or not opener.result:
        return

    for file_handle in opener.result.values():
        sdf_reader = SdfReader(
            file_handle,
            names={
                'HMDB_ID': 'HMDB_ID',
                'PUBCHEM_CID': 'PUBCHEM_CID',
                'SWISSLIPIDS_ID': 'SWISSLIPIDS_ID',
                'LM_ID': 'LM_ID',
                'ABBREVIATION': 'ABBREVIATION',
                'CHEBI_ID': 'CHEBI_ID',
                'SYNONYMS': 'SYNONYMS',
                'INCHI': 'INCHI',
                'INCHI_KEY': 'INCHI_KEY',
                'COMMON_NAME': 'COMMON_NAME',
                'SYSTEMATIC_NAME': 'SYSTEMATIC_NAME',
                'SMILES': 'SMILES',
                'FORMULA': 'FORMULA',
            },
            fields={
                'EXACT_MASS',
                'CATEGORY',
                'MAIN_CLASS',
                'SUB_CLASS',
                'NAME',
            },
        )

        for record in sdf_reader:
            flat_record = {}
            if 'name' in record:
                flat_record.update(record['name'])
            if 'annot' in record:
                flat_record.update(record['annot'])
            yield flat_record
        break


config = ResourceConfig(
    id=ResourceCv.LIPIDMAPS,
    name='LIPID MAPS Structure Database',
    url='https://lipidmaps.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='37855672',
    description=(
        'The LIPID MAPS Structure Database (LMSD) is a comprehensive database '
        'of lipid structures, annotations, and cross-references. It contains '
        'over 47,000 unique lipid structures, classified according to a '
        'comprehensive lipid classification system. The database includes '
        'structure-based identifiers (InChI, InChIKey, SMILES), chemical '
        'properties (formula, exact mass), lipid classification (category, '
        'main class, sub class), and cross-references to other databases '
        '(ChEBI, PubChem, SwissLipids, HMDB).'
    ),
)

f = FieldConfig(
    extract={
        'chebi': r'^(?:CHEBI:)?(\d+)$',
    },
    transform={
        'chebi': lambda v: f'CHEBI:{v}' if v else None,
    },
)

lipids_schema = EntityBuilder(
    entity_type=EntityTypeCv.LIPID,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.LIPIDMAPS, value=f('LM_ID')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('COMMON_NAME')),
        CV(term=IdentifierNamespaceCv.IUPAC_NAME, value=f('SYSTEMATIC_NAME')),
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('ABBREVIATION')),
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('SYNONYMS', delimiter=';')),
        CV(term=IdentifierNamespaceCv.STANDARD_INCHI_KEY, value=f('INCHI_KEY')),
        CV(term=IdentifierNamespaceCv.STANDARD_INCHI, value=f('INCHI')),
        CV(term=IdentifierNamespaceCv.SMILES, value=f('SMILES')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('FORMULA')),
        CV(
            term=IdentifierNamespaceCv.CHEBI,
            value=f('CHEBI_ID', extract='chebi'),
        ),
        CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=f('PUBCHEM_CID')),
        CV(term=IdentifierNamespaceCv.HMDB, value=f('HMDB_ID')),
        CV(term=IdentifierNamespaceCv.SWISSLIPIDS, value=f('SWISSLIPIDS_ID')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.MASS_DALTON, value=f('EXACT_MASS')),
        CV(term=MoleculeAnnotationsCv.LIPID_CATEGORY, value=f('CATEGORY')),
        CV(term=MoleculeAnnotationsCv.LIPID_MAIN_CLASS, value=f('MAIN_CLASS')),
        CV(term=MoleculeAnnotationsCv.LIPID_SUB_CLASS, value=f('SUB_CLASS')),
    ),
)

resource = Resource(
    config,
    lipids=Dataset(
        download=Download(
            url='https://lipidmaps.org/files/?file=LMSD&ext=sdf.zip',
            filename='structures.zip',
            subfolder='lipidmaps',
            large=True,
            ext='zip',
            default_mode='rb',
        ),
        mapper=lipids_schema,
        raw_parser=_parse_lipids_sdf,
    ),
)

lipids = resource.lipids
