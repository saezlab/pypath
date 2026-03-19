"""
Parse ChEMBL data and emit Entity records.

This module converts ChEMBL ligand-target interaction data
into Entity records using the schema defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from functools import partial
from pathlib import Path

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_sqlite
from pypath.share import cache

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
    MoleculeSubtypeCv,
    ProteinFunctionalClassCv,
    MoleculeAnnotationsCv,
)

from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)


VERSION = 36
# Relative path to the SQLite database file within the downloaded archive
DB_REL_PATH = f'chembl_{VERSION}/chembl_{VERSION}_sqlite/chembl_{VERSION}.db'
# Local path to the cached SQLite database file
SQLITE_PATH = Path(cache.get_cachedir()) / f'ChEMBL_SQLite_{VERSION}.sqlite'


def _chembl_url(version: int = VERSION, **_kwargs: object) -> str:
    return ('https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/'
            f'chembl_{version:02d}/chembl_{version:02d}_sqlite.tar.gz')


def _files_needed(version: int = VERSION, **_kwargs: object) -> list[str]:
    return [f'chembl_{version:02d}/chembl_{version:02d}_sqlite/chembl_{version:02d}.db']


config = ResourceConfig(
    id=ResourceCv.CHEMBL,
    name="ChEMBL",
    url="https://www.ebi.ac.uk/chembl/",
    license=LicenseCV.CC_BY_SA_3_0,
    update_category=UpdateCategoryCV.REGULAR,
    description=(
        "ChEMBL is a manually curated chemical database of bioactive molecules "
        "with drug-like properties."
    ),
)

download = Download(
    url=_chembl_url,
    filename=f"chembl_{VERSION}_sqlite.tar.gz",
    subfolder="chembl",
    large=True,
    ext=".tar.gz",
    needed=_files_needed()
)

# Mapping of ChEMBL molecule_type strings to EntityTypeCv
MOLECULE_TYPE_TO_ENTITY_TYPE = {
    'Small molecule': EntityTypeCv.SMALL_MOLECULE,
    'Protein': EntityTypeCv.PROTEIN,
    'Antibody': EntityTypeCv.PROTEIN,
    'Enzyme': EntityTypeCv.PROTEIN,
    'Oligosaccharide': EntityTypeCv.SMALL_MOLECULE,
    'Oligonucleotide': EntityTypeCv.SMALL_MOLECULE,
    'Gene': EntityTypeCv.GENE,
    'Cell': EntityTypeCv.PHYSICAL_ENTITY,
    'Unknown': EntityTypeCv.PHYSICAL_ENTITY,
    'Unclassified': EntityTypeCv.PHYSICAL_ENTITY,
}

# Mapping of ChEMBL molecule_type strings to more specific subclasses
MOLECULE_TYPE_TO_SUBTYPE = {
    'Antibody': MoleculeSubtypeCv.ANTIBODY,
    'Enzyme': ProteinFunctionalClassCv.ENZYME,
    'Small molecule': MoleculeSubtypeCv.SYNTHETIC_ORGANIC,
}

f = FieldConfig(
    extract={
        'chebi': r'^(?:CHEBI:)?(\d+)$',
    },
    transform={
        'chebi': lambda v: f'CHEBI:{v}' if v else None,
        'bool_to_cv': lambda v, cv: cv if str(v) == '1' else None,
    },
    map={
        'entity_type': MOLECULE_TYPE_TO_ENTITY_TYPE,
        'subtype': MOLECULE_TYPE_TO_SUBTYPE,
    },
)

molecules_schema = EntityBuilder(
    entity_type=f('molecule_type', map='entity_type'),
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.CHEMBL, value=f('chembl_id')),
        CV(term=IdentifierNamespaceCv.CHEMBL_INTERNAL_ID, value=f('molregno')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('pref_name')),
        CV(term=IdentifierNamespaceCv.CHEBI, value=f('chebi_id', extract='chebi', transform='chebi')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.CLINICAL_PHASE, value=f('max_phase')),
        CV(term=f('molecule_type', map='subtype')),
        CV(term=f('therapeutic_flag', transform=lambda v: f.transform['bool_to_cv'](v, MoleculeAnnotationsCv.APPROVED))),
        CV(term=f('withdrawn_flag', transform=lambda v: f.transform['bool_to_cv'](v, MoleculeAnnotationsCv.WITHDRAWN))),
        CV(term=f('polymer_flag', transform=lambda v: f.transform['bool_to_cv'](v, MoleculeAnnotationsCv.POLYMER))),
        CV(term=f('inorganic_flag', transform=lambda v: f.transform['bool_to_cv'](v, MoleculeAnnotationsCv.INORGANIC))),
        CV(term=f('natural_product', transform=lambda v: f.transform['bool_to_cv'](v, MoleculeAnnotationsCv.NATURAL_PRODUCT))),
    ),
)

resource = Resource(
    config=config,
    molecules=Dataset(
        download=download,
        mapper=molecules_schema,
        raw_parser=partial(
            iter_sqlite,
            table_name='molecule_dictionary',
            sqlite_path=SQLITE_PATH,
            db_rel_path=DB_REL_PATH,
        ),
    ),
)
