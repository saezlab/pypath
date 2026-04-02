"""
Parse HMDB (Human Metabolome Database) data and emit Entity records.

This module downloads HMDB metabolite XML data and converts it into Entity
records using the declarative schema pattern from tabular_builder.
"""

from __future__ import annotations

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
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.hmdb import _raw


config = ResourceConfig(
    id=ResourceCv.HMDB,
    name='Human Metabolome Database',
    url='https://hmdb.ca/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='37953221',
    description=(
        'The Human Metabolome Database (HMDB) is a comprehensive database '
        'containing detailed information about small molecule metabolites '
        'found in the human body. It includes chemical, clinical, and '
        'biochemical/molecular biology data, with over 220,000 metabolite '
        'entries including both water-soluble and lipid-soluble metabolites.'
    ),
)

f = FieldConfig(
    extract={
        'chebi': r'^(?:CHEBI:)?(\d+)$',
        'drugbank': r'^(DB\d+)$',
        'kegg_compound': r'^([CDGcdg])(\d{4,5})$',
    },
    transform={
        'chebi': lambda v: f'CHEBI:{v}' if v else None,
        'kegg_compound': lambda v: f'{v[0].upper()}{v[1].zfill(5)}' if v and len(v) == 2 else None,
    },
)

metabolites_schema = EntityBuilder(
    entity_type=EntityTypeCv.SMALL_MOLECULE,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.HMDB, value=f('accession')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('name')),
        CV(term=IdentifierNamespaceCv.IUPAC_TRADITIONAL_NAME, value=f('traditional_iupac')),
        CV(term=IdentifierNamespaceCv.IUPAC_NAME, value=f('iupac_name')),
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('synonyms', delimiter=';')),
        CV(term=IdentifierNamespaceCv.STANDARD_INCHI_KEY, value=f('inchikey')),
        CV(term=IdentifierNamespaceCv.STANDARD_INCHI, value=f('inchi')),
        CV(term=IdentifierNamespaceCv.SMILES, value=f('smiles')),
        CV(
            term=IdentifierNamespaceCv.CHEBI,
            value=f('chebi_id', extract='chebi'),
        ),
        CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=f('pubchem_compound_id')),
        CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('kegg_id', extract='kegg_compound', transform='kegg_compound')),
        CV(term=IdentifierNamespaceCv.DRUGBANK, value=f('drugbank_id', extract='drugbank')),
        CV(term=IdentifierNamespaceCv.CAS, value=f('cas_registry_number')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.DESCRIPTION, value=f('description')),
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('pubmed_ids', delimiter=';')),
    ),
)

download = Download(
    url='https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip',
    filename='hmdb_metabolites.zip',
    subfolder='hmdb',
    large=True,
    ext='zip',
    default_mode='rb',
)


def _id_translation_row(row: dict) -> dict | None:
    hmdb_id = row.get('accession')
    standard_inchi = row.get('inchi')
    if not hmdb_id or not standard_inchi:
        return None
    return {
        'source': 'hmdb',
        'key_type': 'OM:0004:Hmdb',
        'key_value': hmdb_id,
        'standard_inchi': standard_inchi,
    }


resource = Resource(
    config,
    metabolites=Dataset(
        download=download,
        mapper=metabolites_schema,
        raw_parser=_raw,
    ),
    id_translation=Dataset(
        download=download,
        mapper=lambda row: row,
        raw_parser=lambda opener, **kwargs: (
            row
            for raw_row in _raw(opener, **kwargs)
            if (row := _id_translation_row(raw_row)) is not None
        ),
    ),
)
