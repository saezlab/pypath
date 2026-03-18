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

molecules_schema = EntityBuilder(
    entity_type=EntityTypeCv.SMALL_MOLECULE,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.CHEMBL, value='molecule_chembl_id'),
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
