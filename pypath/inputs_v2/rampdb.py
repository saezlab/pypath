"""
Parse RaMP-DB data and emit Entity records.

This module converts annotations of lipids and metabolites into Entity records
using the declarative schema pattern.
"""

from functools import partial
import gzip
import os
from pathlib import Path
import re
import shutil

from biolink_model.datamodel.model import ChemicalEntity, Gene, Pathway, slots
from bs4 import BeautifulSoup
from omnipath_core.naming import Namespace
import requests

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_sqlite
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

BASE_URL = 'https://github.com/ncats/RaMP-DB/raw/refs/heads/%s/db/'


def get_ramp_latest_ver(branch='main'):
    """
    Retrieves the URL for the latest version of the database SQL file
    """
    res = requests.get(BASE_URL % branch)
    soup = BeautifulSoup(res.text, 'html.parser')
    files = sorted(
        {f.text for f in soup.find_all(title=re.compile('\\.sqlite.gz$'))}
    )
    ver = re.search('(\\d+\\.\\d+\\.\\d+)', files[-1]).group(1)
    return (ver, BASE_URL % branch + files[-1])


config = ResourceConfig(
    id=ResourceCv.RAMPDB,
    name='RaMP-DB',
    url='https://rampdb.nih.gov/',
    license=LicenseCV.GPL_2_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='36373969',
    primary_category='pathways',
    description='RaMP-DB (Relational database of Metabolomic Pathways) is a multi-sourced integrated database with comprehensive annotations on biological pathways, structure/chemistry, disease and ontology annotations for genes, proteins, and metabolites.',
)


class RampDownload:
    """Resolve the remote release only when a download is requested."""

    def open(self, *, force_refresh=False, **kwargs):
        version, url = get_ramp_latest_ver()
        return Download(
            url=url,
            filename=os.path.basename(url),
            subfolder='ramp',
            large=True,
            ext='gz',
            default_mode='rb',
        ).open(force_refresh=force_refresh, **kwargs)


download = RampDownload()
table_names = [
    'analyte',
    'metabolite_class',
    'analytehaspathway',
    'pathway',
    'source',
    'chem_props',
]


def parser(opener, table=None, **_kwargs):
    path = Path(opener.path.replace('.gz', ''))
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(opener.path, 'rb') as src, path.open('wb') as dst:
            shutil.copyfileobj(src, dst, length=1024 * 1024)
    yield from iter_sqlite(opener, table_name=table, sqlite_path=path)


MOLECULE_TYPE_TO_ENTITY_TYPE = {'gene': Gene, 'compound': ChemicalEntity}
SOURCE_TO_ENTITY_TYPE = {'hmdb': ChemicalEntity, 'lipidmaps': ChemicalEntity}
SOURCE_TO_TERM = {
    'hmdb': Namespace.HMDB,
    'lipidmaps': Namespace.LIPIDMAPS,
    'cas': Namespace.CAS,
    'en': Namespace.EC,
    'brenda': Namespace.EC,
    'chebi': Namespace.CHEBI,
    'chemspider': Namespace.CHEMSPIDER,
    'ensembl': Namespace.ENSG,
    'entrez': Namespace.ENTREZ,
    'gene_symbol': Namespace.GENESYMBOL,
    'kegg': Namespace.KEGG,
    'reactome': Namespace.REACTOME,
    'wiki': Namespace.WIKIPATHWAYS,
    'pfocr': Namespace.PFOCR,
}
ID_TO_ENTITY = {'G': Gene, 'C': ChemicalEntity}
f = FieldConfig(
    extract={'rampID': '^RAMP_([A-Z]+)_\\d+$'},
    map={
        'type_to_entity': MOLECULE_TYPE_TO_ENTITY_TYPE,
        'source_to_entity': SOURCE_TO_ENTITY_TYPE,
        'source_to_term': SOURCE_TO_TERM,
        'id_to_entity': ID_TO_ENTITY,
    },
    transform={
        'lower': lambda x: x.lower(),
        'upper': lambda x: x.upper(),
        'postcolon': lambda x: x.split(':')[-1],
    },
)
analyte_schema = EntityBuilder(
    entity_type=f('type', map='type_to_entity'),
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.RAMP, value=f('rampId')),
        CV(term=Namespace.NAME, value=f('common_name')),
    ),
    annotations=AnnotationsBuilder(),
)
metabolite_class_schema = EntityBuilder(
    entity_type=f('source', map='source_to_entity'),
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.RAMP, value=f('ramp_id')),
        CV(
            term=f('source', map='source_to_term', transform='lower'),
            value=f('class_source_id', transform='postcolon'),
        ),
    ),
    annotations=AnnotationsBuilder(),
)
source_schema = EntityBuilder(
    entity_type=f('geneOrCompound', map='type_to_entity'),
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.RAMP, value=f('rampId')),
        CV(
            term=f('IDtype', map='source_to_term', transform='lower'),
            value=f('sourceId', transform='postcolon'),
        ),
        CV(term=Namespace.NAME, value=f('commonName')),
    ),
    annotations=AnnotationsBuilder(),
)
chem_props_schema = EntityBuilder(
    entity_type=ChemicalEntity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.RAMP, value=f('ramp_id')),
        CV(
            term=f('chem_data_source', map='source_to_term', transform='lower'),
            value=f('chem_source_id', transform='postcolon'),
        ),
        CV(term=Namespace.SMILES, value=f('iso_smiles')),
        CV(term=Namespace.INCHIKEY, value=f('inchi_key')),
        CV(term=Namespace.INCHI, value=f('inchi')),
        CV(term=Namespace.NAME, value=f('common_name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.has_chemical_formula, value=f('mol_formula')),
        CV(term='chemrof:mass', value=f('mw')),
        CV(term='chemrof:monoisotopic_mass', value=f('monoisotop_mass')),
    ),
)
analytehaspathway_schema = RelationBuilder(
    subject=EntityBuilder(
        entity_type=f('rampId', map='id_to_entity', extract='rampID'),
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.RAMP, value=f('rampId'))
        ),
    ),
    predicate=slots.participates_in,
    object=EntityBuilder(
        entity_type=Pathway,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.RAMP, value=f('pathwayRampId'))
        ),
    ),
    identifiers=IdentifiersBuilder(),
)
pathway_schema = EntityBuilder(
    entity_type=Pathway,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.RAMP, value=f('pathwayRampId'))
    ),
    annotations=AnnotationsBuilder(),
)
schemas = {
    'analyte': analyte_schema,
    'metabolite_class': metabolite_class_schema,
    'analytehaspathway': analytehaspathway_schema,
    'pathway': pathway_schema,
    'source': source_schema,
    'chem_props': chem_props_schema,
}
kwargs = {
    t: Dataset(
        download=download,
        mapper=schemas[t],
        raw_parser=partial(parser, table=t),
    )
    for t in table_names
}
resource = Resource(config=config, **kwargs)
