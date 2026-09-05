"""
Parse HMDB (Human Metabolome Database) data and emit Entity records.

This module downloads HMDB metabolite XML data and converts it into Entity
records using the declarative schema pattern from tabular_builder.
"""

from __future__ import annotations

from biolink_model.datamodel.model import ChemicalEntity, OntologyClass, slots
from omnipath_core.naming import Namespace

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.hmdb import _raw, chemont_name_map
from pypath.internals.cv_terms import (
    LicenseCV,
    OntologyCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    AssociationBuilder,
    AssociationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)

config = ResourceConfig(
    id=ResourceCv.HMDB,
    name='Human Metabolome Database',
    short='HMDB',
    url='https://hmdb.ca/',
    license=LicenseCV.CC_BY_NC_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='34986597',
    primary_category='small_molecules',
    annotation_ontologies=(OntologyCv.CHEMONT,),
    description='The Human Metabolome Database (HMDB) is a comprehensive database containing detailed information about small molecule metabolites found in the human body. It includes chemical, clinical, and biochemical/molecular biology data, with over 220,000 metabolite entries including both water-soluble and lipid-soluble metabolites.',
)
f = FieldConfig(
    extract={
        'chebi': '^(?:CHEBI:)?(\\d+)$',
        'drugbank': '^(DB\\d+)$',
        'kegg_compound': '^([CDGcdg])(\\d{4,5})$',
    },
    transform={
        'kegg_compound': lambda v: f'{v[0].upper()}{v[1].zfill(5)}'
        if v and len(v) == 2
        else None
    },
)
metabolites_schema = EntityBuilder(
    entity_type=ChemicalEntity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.HMDB, value=f('accession')),
        CV(term=Namespace.INCHIKEY, value=f('inchikey')),
        CV(term=Namespace.INCHI, value=f('inchi')),
        CV(term=Namespace.SMILES, value=f('smiles')),
        CV(term=Namespace.CHEBI, value=f('chebi_id', extract='chebi')),
        CV(term=Namespace.PUBCHEM, value=f('pubchem_compound_id')),
        CV(
            term=Namespace.KEGG,
            value=f(
                'kegg_id', extract='kegg_compound', transform='kegg_compound'
            ),
        ),
        CV(term=Namespace.DRUGBANK, value=f('drugbank_id', extract='drugbank')),
        CV(term=Namespace.CAS, value=f('cas_registry_number')),
        CV(term=Namespace.NAME, value=f('name')),
        CV(term=Namespace.SYNONYM, value=f('traditional_iupac')),
        CV(term=Namespace.SYNONYM, value=f('iupac_name')),
        CV(term=Namespace.SYNONYM, value=f('synonyms', delimiter=';')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.description, value=f('description')),
        CV(
            term=slots.publications,
            value=f(
                'pubmed_ids',
                delimiter=';',
                transform=lambda v: 'PMID:' + str(v).removeprefix('PMID:'),
            ),
        ),
    ),
    associations=AssociationsBuilder(
        AssociationBuilder(
            object_entity_type=OntologyClass,
            object_identifier_type=Namespace.CHEMONT,
            object_identifier=f('chemont_ids', delimiter=';'),
        )
    ),
)
download = Download(
    url='https://rescued.omnipathdb.org/hmdb_metabolites.zip',
    filename='hmdb_metabolites.zip',
    subfolder='hmdb',
    large=True,
    ext='.zip',
    default_mode='rb',
)
chemont_download = Download(
    url='http://classyfire.wishartlab.com/system/downloads/1_0/chemont/ChemOnt_2_1.obo.zip',
    filename='ChemOnt_2_1.obo.zip',
    subfolder='chemont',
    large=True,
    ext='zip',
    needed=['ChemOnt_2_1.obo'],
)


def _raw_with_ontology_maps(opener, force_refresh: bool = False, **kwargs):
    chemont_map = chemont_name_map(
        chemont_download.open(force_refresh=force_refresh)
    )
    yield from _raw(
        opener,
        force_refresh=force_refresh,
        chemont_name_map=chemont_map,
        **kwargs,
    )


def _id_translation_row(row: dict) -> dict | None:
    hmdb_id = row.get('accession')
    standard_inchi = row.get('inchi')
    if not hmdb_id or not standard_inchi:
        return None
    return {
        'source': 'hmdb',
        'key_type': Namespace.HMDB,
        'key_value': hmdb_id,
        'standard_inchi': standard_inchi,
    }


resource = Resource(
    config,
    metabolites=Dataset(
        download=download,
        mapper=metabolites_schema,
        raw_parser=_raw_with_ontology_maps,
    ),
    id_translation=Dataset(
        download=download,
        mapper=lambda row: row,
        raw_parser=lambda opener, **kwargs: (
            row
            for raw_row in _raw(opener, **kwargs)
            if (row := _id_translation_row(raw_row)) is not None
        ),
        kind='id_translation',
    ),
)
