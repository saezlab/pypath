"""
Parse PMI-DB protein-metabolite interactions and emit Relation records.

PMI-DB collects protein-metabolite interactions and non-interactions measured
by chemoproteomics (LiP-SMap, CETSA/TPP, SPROX, lipid, sterol and flavonoid
probes) in E. coli, yeast, human and mouse. Each row is one pair tested in one
study. Only the pairs a study found to interact are imported; the tested
non-interactions are not statements about the pair.
"""

from __future__ import annotations

from typing import Any, Iterator

from biolink_model.datamodel.model import ChemicalEntity, Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_tsv
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.PMIDB,
    name='PMI-DB',
    url='http://easybioai.com/PMIDB/',
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='33554247',
    topics=(Topic.METABOLISM,),
    content=(Content.PHYSICAL_INTERACTIONS, Content.COMPOUND_TARGET_INTERACTIONS),
    knowledge_origin=(KnowledgeOrigin.EXPERIMENTAL,),
    summary='Protein–metabolite interactions measured by chemoproteomics in E. coli, yeast, human and mouse.',
    description=(
        'PMI-DB provides both interactions and non-interactions between '
        'proteins and metabolites, collected from chemoproteomics studies, '
        'for experimental biologists and for training machine learning models.'
    ),
)

BASE_URL = 'http://easybioai.com/PMIDB/static/%s.txt'
# The two information files have no header row; the names are the website's.
METABOLITE_FIELDS = ('kegg', 'hmdb', 'name', 'link', 'smiles')
PROTEIN_FIELDS = ('uniprot', 'gene', 'go', 'link')
# `Sample source` as the file writes it. "yeast" is Saccharomyces cerevisiae
# (the SPROX study it comes from, PMID:24741112).
TAXA = {
    'Escherichia coli': '562',
    'yeast': '4932',
    'human': '9606',
    'mouse': '10090',
}
INTERACTS = 'TRUE'


def _download(name: str) -> Download:
    return Download(
        url=BASE_URL % name,
        filename=f'{name}.txt',
        subfolder='pmidb',
        ext='txt',
    )


def _value(value: Any) -> str | None:
    """A cell, with the file's NA and empty cells as None."""
    text = str(value or '').strip()
    return None if text in ('', 'NA') else text


def _info(name: str, fields: tuple[str, ...], key: str, force_refresh: bool) -> dict[str, dict]:
    opener = _download(name).open(force_refresh=force_refresh)
    return {
        row[key]: row
        for row in iter_tsv(opener, fieldnames=fields, quotechar=None)
        if _value(row.get(key))
    }


def iter_pmidb(opener, force_refresh: bool = False, **kwargs: Any) -> Iterator[dict[str, Any]]:
    """The interacting pairs of `interaction.txt`, each with its metabolite's
    identifiers from `metabolites_infor.txt` and its protein's gene names from
    `protein_infor.txt`."""
    metabolites = _info('metabolites_infor', METABOLITE_FIELDS, 'kegg', force_refresh)
    proteins = _info('protein_infor', PROTEIN_FIELDS, 'uniprot', force_refresh)
    for row in iter_tsv(opener, quotechar=None, required=('KEGG', 'Uniprot_KB_id', 'interaction'), **kwargs):
        if row['interaction'] != INTERACTS:
            continue
        source = row.get('Sample source')
        if source not in TAXA:
            raise ValueError(f'PMI-DB: unknown sample source {source!r}')
        metabolite = metabolites.get(row['KEGG'], {})
        protein = proteins.get(row['Uniprot_KB_id'], {})
        yield {
            **row,
            'taxon_id': TAXA[source],
            'hmdb': _value(metabolite.get('hmdb')),
            'metabolite_name': _value(metabolite.get('name')),
            'smiles': _value(metabolite.get('smiles')),
            'gene': _value(protein.get('gene')),
        }


f = FieldConfig()

metabolite_builder = EntityBuilder(
    entity_type=ChemicalEntity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.KEGG, value=f('KEGG')),
        CV(term=Namespace.HMDB, value=f('hmdb')),
        CV(term=Namespace.SMILES, value=f('smiles')),
        CV(term=Namespace.NAME, value=f('metabolite_name')),
    ),
)
protein_builder = EntityBuilder(
    entity_type=Protein,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.UNIPROT, value=f('Uniprot_KB_id')),
        CV(term=Namespace.GENESYMBOL, value=f('gene', delimiter=',')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=lambda row: f"NCBITaxon:{row['taxon_id']}"),
    ),
)
interactions_schema = RelationBuilder(
    subject=metabolite_builder,
    predicate=slots.interacts_with,
    object=protein_builder,
    annotations=AnnotationsBuilder(
        CV(term=slots.publications, value=lambda row: pmid_curie(str(row.get('literature') or '').removeprefix('PMID:'))),
        # The study's method and thresholds, as written.
        CV(term='pmidb:technology', value=lambda row: _value(row.get('Technology'))),
        CV(term='pmidb:probe', value=lambda row: _value(row.get('probe'))),
        CV(term='pmidb:mean', value=lambda row: _value(row.get('mean'))),
        CV(term='pmidb:ms_quality_control_method', value=lambda row: _value(row.get('MS Quality control method'))),
        CV(term='pmidb:candidate_selection_cutoff', value=lambda row: _value(row.get('Candidate selection cutoff'))),
    ),
)

resource = Resource(
    config,
    interactions=Dataset(
        download=_download('interaction'),
        mapper=interactions_schema,
        raw_parser=iter_pmidb,
    ),
)
