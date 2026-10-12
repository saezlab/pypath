"""
Parse LRdb ligand-receptor pairs and emit Relation records.

LRdb is the ligand-receptor table of SingleCellSignalR. Each row is a human
ligand-receptor pair named by gene symbols, with protein names, the databases
the pair was taken from and the PubMed IDs supporting it.
"""

from __future__ import annotations
from collections.abc import Callable
import re
from typing import Any
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from biolink_model.datamodel.model import (
    Protein,
    slots,
)
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic
from omnipath_core.source_attributes import PARTICIPANT_ROLE
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.lrdb import iter_pairs

config = ResourceConfig(
    id=ResourceCv.LRDB,
    name='LRdb',
    url='https://github.com/SCA-IRCM/LRdb',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.GPL_3_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='32196115',
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(Content.LIGAND_RECEPTOR_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='Human ligand–receptor pairs with their supporting databases and publications.',
    description='LRdb is the ligand-receptor database of SingleCellSignalR. It combines human ligand-receptor pairs from FANTOM5, HPRD, HPMR, the Guide to Pharmacology, UniProt, Reactome and cellsignal.com with pairs the authors curated from the literature, and records for each pair the databases and the PubMed IDs that support it.',
)
download_interactions = Download(
    url='https://raw.githubusercontent.com/SCA-IRCM/LRdb/master/LRdb_122019.txt',
    filename='LRdb_122019.txt',
    subfolder='lrdb',
)
# The table has no species column. LRdb is a human database.
TAXON = 'NCBITaxon:9606'
# Every value of the ``source`` column that names a database, in the source's
# spelling. The remaining value, ``literature``, marks a pair the LRdb authors
# curated themselves; its evidence is the row's PubMed IDs. The reader raises
# on any other value, so none is dropped here.
DATABASES = {
    'fantom5': 'fantom5',
    'HPRD': 'HPRD',
    'HPMR': 'HPMR',
    'IUPHAR': 'IUPHAR',
    'reactome': 'reactome',
    'uniprot': 'uniprot',
    'cellsignal.com': 'cellsignal.com',
}
# UniProt lists a precursor's chains or domains after its name; the source
# copied that list, often cut short, into the name column.
_name_sections = re.compile(r' \[(?:Cleaved into|Includes): .*$')


def _protein_name(value: str) -> str | None:
    return _name_sections.sub('', value).strip() or None


def _synonym(value: str) -> str | None:
    # A few synonyms end in a stray semicolon.
    return value.rstrip(';').strip() or None


def _remark(row: dict[str, Any]) -> str | None:
    # One remark went through CSV quoting several times; undo the doubling.
    return re.sub(r'"{2,}', '"', row.get('remarks') or '').strip() or None


f = FieldConfig(
    extract={
        # Stricter than the shared PubMed identifier rule (1 to 8 digits), because
        # the file has cut numbers. Every intact PMID in the file has seven or
        # eight digits. Many cells were cut after 10 or 11 characters, leaving a
        # fragment of the next PMID (``9914169,952``): 273 of the 3,251 rows have
        # such a token, 268 of their cells are exactly 10 or 11 characters long,
        # and in 270 the token is the last of the cell. Of the 100 distinct
        # fragments, a random 30 were looked up at NCBI; all 30 are 1975 or 1976
        # papers on other subjects (``1`` is a formate assay in methanol
        # poisoning, ``952`` growth of Staphylococcus on frankfurters). They are
        # left out, not repaired. Three PMIDs carry a stray closing parenthesis.
        'pmid': r'^([1-9]\d{6,7})\)?$',
    },
    map={'database': DATABASES},
    transform={
        'pmid_curie': lambda v: f'PMID:{v}',
        'protein_name': _protein_name,
        'synonym': _synonym,
    },
)


def _participant(role: str) -> Callable[[dict[str, Any]], Any]:
    builder = EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f(role)),
            CV(
                term=Namespace.NAME,
                value=f(f'{role}.name', transform='protein_name'),
            ),
            CV(
                term=Namespace.SYNONYM,
                value=f(f'{role}.synonyms', delimiter=',', transform='synonym'),
            ),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=TAXON),
            # The columns are named for the role (`ligand`, `receptor`).
            CV(term=PARTICIPANT_ROLE, value=role),
        ),
    )

    def build(row: dict[str, Any]) -> Any:
        # The gene symbol identifies the participant; a name alone does not.
        return builder.build(row) if (row.get(role) or '').strip() else None

    return build


# The ligand is the subject and the receptor the object.
interactions_schema = RelationBuilder(
    subject=_participant('ligand'),
    predicate=slots.interacts_with,
    object=_participant('receptor'),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.supporting_data_source,
            value=f('source', delimiter=',', map='database'),
        ),
        CV(
            term=slots.publications,
            value=f(
                'PMIDs',
                delimiter=',',
                extract='pmid',
                transform='pmid_curie',
            ),
        ),
        CV(term=slots.description, value=_remark),
    ),
)
resource = Resource(
    config,
    interactions=Dataset(
        download=download_interactions,
        mapper=interactions_schema,
        raw_parser=iter_pairs,
    ),
)
