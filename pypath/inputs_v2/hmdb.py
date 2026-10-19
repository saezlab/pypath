"""
Parse HMDB (Human Metabolome Database) data and emit Entity records.

This module downloads HMDB metabolite XML data and converts it into Entity
records using the declarative schema pattern from tabular_builder.

A metabolite is ``located_in`` the cellular locations, tissues and
biospecimens HMDB lists for it, and ``associated_with`` its diseases. Biolink
maps HMDB's ``at_cellular_location``, ``at_tissue`` and ``in_biospecimen`` to
``located_in``. A location is the GO, UBERON or CL term of the crosswalk below
where the two name the same thing, otherwise a term of HMDB's own; the
``terms`` dataset names HMDB's terms.
"""

from __future__ import annotations

from biolink_model.datamodel.model import ChemicalEntity, OntologyClass, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.hmdb import _raw, chemont_name_map, iter_terms
from pypath.internals.silver_schema import Entity, Identifier
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
    name='HMDB',
    full='Human Metabolome Database',
    short='HMDB',
    url='https://hmdb.ca/',
    license=LicenseCV.CC_BY_NC_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='34986597',
    annotation_ontologies=(OntologyCv.CHEMONT,),
    topics=(Topic.METABOLISM,),
    content=(Content.MOLECULAR_ANNOTATIONS, Content.IDENTIFIER_MAPPINGS),
    summary='Human metabolite identities, chemical properties and biological annotations.',
    description='The Human Metabolome Database (HMDB) is a comprehensive database containing detailed information about small molecule metabolites found in the human body. It includes chemical, clinical, and biochemical/molecular biology data, with over 220,000 metabolite entries including both water-soluble and lipid-soluble metabolites.',
)
# The vocabularies of HMDB's biological properties: 11 cellular locations, 42
# tissues and 18 biospecimens in the release of the cache (2026-10-19), each
# reviewed against the term it is given. A value not here keeps HMDB's term:
# "All Tissues", "Microsomes" (GO's term is obsolete), "Breath", and the
# pathological "Pericardial Effusion" and "Ascites Fluid".
CELLULAR_LOCATION_TERMS = {
    'Cytoplasm': 'GO:0005737',
    'Endoplasmic reticulum': 'GO:0005783',
    'Extracellular': 'GO:0005576',  # extracellular region
    'Golgi apparatus': 'GO:0005794',
    'Inner mitochondrial membrane': 'GO:0005743',
    'Lysosome': 'GO:0005764',
    'Membrane': 'GO:0016020',
    'Mitochondria': 'GO:0005739',
    'Nucleus': 'GO:0005634',
    'Peroxisome': 'GO:0005777',
}
# Tissues HMDB names by a cell type take the Cell Ontology term.
TISSUE_TERMS = {
    'Adipose Tissue': 'UBERON:0001013',
    'Adrenal Cortex': 'UBERON:0001235',
    'Adrenal Gland': 'UBERON:0002369',
    'Adrenal Medulla': 'UBERON:0001236',
    'Basal Ganglia': 'UBERON:0010011',
    'Bile': 'UBERON:0001970',
    # HMDB lists the gall bladder apart, so this is the urinary bladder.
    'Bladder': 'UBERON:0001255',
    'Blood': 'UBERON:0000178',
    'Bone Marrow': 'UBERON:0002371',
    'Brain': 'UBERON:0000955',
    'Cartilage': 'UBERON:0002418',
    'Epidermis': 'UBERON:0001003',
    'Erythrocyte': 'CL:0000232',
    'Eye Lens': 'UBERON:0000965',
    'Fibroblasts': 'CL:0000057',
    'Gall Bladder': 'UBERON:0002110',
    'Hair': 'UBERON:0001037',
    'Heart': 'UBERON:0000948',
    'Intestine': 'UBERON:0000160',
    'Kidney': 'UBERON:0002113',
    'Leukocyte': 'CL:0000738',
    'Liver': 'UBERON:0002107',
    'Lung': 'UBERON:0002048',
    'Neuron': 'CL:0000540',
    'Ovary': 'UBERON:0000992',
    'Pancreas': 'UBERON:0001264',
    'Parathyroid Gland': 'UBERON:0001132',
    'Pineal Gland': 'UBERON:0001905',
    'Placenta': 'UBERON:0001987',
    'Platelet': 'CL:0000233',
    'Prostate': 'UBERON:0002367',
    'Retina': 'UBERON:0000966',
    'Semen': 'UBERON:0001968',
    'Skeletal Muscle': 'UBERON:0001134',
    'Smooth Muscle': 'UBERON:0001135',
    'Spleen': 'UBERON:0002106',
    'Testis': 'UBERON:0000473',
    'Thyroid Gland': 'UBERON:0002046',
    'Umbilical cord': 'UBERON:0002331',
    'Urine': 'UBERON:0001088',
    'Vitreous humor': 'UBERON:0001797',
}
BIOSPECIMEN_TERMS = {
    'Amniotic Fluid': 'UBERON:0000173',
    'Aqueous Humour': 'UBERON:0001796',
    'Bile': 'UBERON:0001970',
    'Blood': 'UBERON:0000178',
    'Breast Milk': 'UBERON:0001913',  # milk; HMDB is human
    'Cellular Cytoplasm': 'GO:0005737',
    'Cerebrospinal Fluid (CSF)': 'UBERON:0001359',
    'Feces': 'UBERON:0001988',
    'Lymph': 'UBERON:0002391',
    'Prostate Tissue': 'UBERON:0002367',
    'Saliva': 'UBERON:0001836',
    'Semen': 'UBERON:0001968',
    'Sweat': 'UBERON:0001089',
    'Tears': 'UBERON:0001827',
    'Urine': 'UBERON:0001088',
}
# Field of the parser: namespace of HMDB's own terms, crosswalk.
LOCATIONS = {
    'cellular_locations': ('hmdb_cellular_location', CELLULAR_LOCATION_TERMS),
    'tissue_locations': ('hmdb_tissue', TISSUE_TERMS),
    'biospecimen_locations': ('hmdb_biospecimen', BIOSPECIMEN_TERMS),
}
DISEASE_NAMESPACE = 'hmdb_disease'
_PREFIX_NAMESPACE = {'GO': Namespace.GO, 'UBERON': Namespace.UBERON, 'CL': Namespace.CL}


def _names(row: dict, field: str) -> list[str]:
    """The parser's list of names; a single name is a list of one."""
    value = row.get(field) or []
    return [value] if isinstance(value, str) else list(value)


def _location_associations(field: str) -> list[AssociationBuilder]:
    """Statements to the mapped terms, one builder per ontology, and to HMDB's."""
    namespace, crosswalk = LOCATIONS[field]
    builders = [
        AssociationBuilder(
            predicate=slots.located_in,
            object_entity_type=OntologyClass,
            object_identifier_type=_PREFIX_NAMESPACE[prefix],
            object_identifier=lambda row, prefix=prefix: [
                term
                for name in _names(row, field)
                if (term := crosswalk.get(name)) and term.startswith(prefix + ':')
            ],
        )
        for prefix in sorted({term.split(':')[0] for term in crosswalk.values()})
    ]
    builders.append(
        AssociationBuilder(
            predicate=slots.located_in,
            object_entity_type=OntologyClass,
            object_identifier_type=namespace,
            object_identifier=lambda row: [
                name for name in _names(row, field) if name not in crosswalk
            ],
        )
    )
    return builders


def map_term(row: dict) -> Entity | None:
    """A location, biospecimen or disease of HMDB's own vocabulary."""
    if row['field'] == 'diseases':
        identifiers = [Identifier(type=DISEASE_NAMESPACE, value=row['name'])]
        if row.get('omim'):
            identifiers.append(Identifier(type='omim', value=row['omim']))
    else:
        namespace, crosswalk = LOCATIONS[row['field']]
        if row['name'] in crosswalk:
            return None
        identifiers = [Identifier(type=namespace, value=row['name'])]
    identifiers.append(Identifier(type=Namespace.NAME, value=row['name']))
    return Entity(type=OntologyClass, identifiers=identifiers)


f = FieldConfig(
    extract={
        'chebi': '^(?:CHEBI:)?(\\d+)$',
        'drugbank': '^(DB\\d+)$',
        # One group: the builder takes the first group of a pattern.
        'kegg_compound': '^([CDGcdg]\\d{4,5})$',
    },
    transform={
        'kegg_compound': lambda v: f'{v[0].upper()}{v[1:].zfill(5)}' if v else None,
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
        ),
        *_location_associations('cellular_locations'),
        *_location_associations('tissue_locations'),
        *_location_associations('biospecimen_locations'),
        AssociationBuilder(
            predicate=slots.associated_with,
            object_entity_type=OntologyClass,
            object_identifier_type=DISEASE_NAMESPACE,
            object_identifier=lambda row: _names(row, 'diseases'),
        ),
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
    terms=Dataset(
        download=download,
        mapper=map_term,
        raw_parser=iter_terms,
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
