"""
Parse the MIPS mammalian protein-protein interaction database (MPPI) and
emit Relation records.

MPPI is one PSI-MI XML level 1 file from 2004. Each ``<interaction>`` is one
experiment of one publication on two proteins: a PubMed ID, a PSI-MI
detection method and, for each protein, a sequence database accession, a name
and a species. The proteins are from several mammals, and the two proteins of
an interaction often from different ones. The file lists most interactions
in both orders; the second of such rows is not mapped.
"""

from __future__ import annotations

from collections.abc import Callable
import re
from typing import Any

from biolink_model.datamodel.model import Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin

from pypath.internals.cv_terms import (
    LicenseCV,
    OntologyCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.inputs_v2._publications import PMID_EXTRACT
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.mppi import iter_interactions


config = ResourceConfig(
    id=ResourceCv.MPPI,
    name='MPPI',
    full='MIPS Mammalian Protein-Protein Interaction Database',
    url='http://mips.helmholtz-muenchen.de/proj/ppi/',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.MPPI,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='15531608',
    annotation_ontologies=(OntologyCv.PSI_MI,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Manually curated mammalian protein-protein interactions from '
        'individual experiments, with their detection methods and '
        'publications.'
    ),
    description=(
        'The MIPS Mammalian Protein-Protein Interaction Database (MPPI) is a '
        'collection of physical protein interactions in mammals that expert '
        'curators extracted from the scientific literature. Only evidence '
        'from individually performed experiments was included, not '
        'large-scale surveys. Each record names the publication, the '
        'detection method and the two proteins with their species. The '
        'database was released in 2004 and is no longer maintained or served '
        'by MIPS.'
    ),
)

# The MIPS site no longer serves the file. This is OmniPath's copy of it,
# byte for byte the file that was at
# http://mips.helmholtz-muenchen.de/proj/ppi/data/mppi.gz from 2011 to 2022.
download_interactions = Download(
    url='https://rescued.omnipathdb.org/mppi.gz',
    filename='mppi.gz',
    subfolder='mppi',
    large=True,
    ext='gz',
    default_mode='rb',
    encoding=None,
)

# The accessions the source writes for each of its sequence databases
# (``primaryRef@db``). Swiss-Prot and TrEMBL are the two sections of UniProt;
# ``GB`` accessions are GenBank nucleotide records. PIR has no namespace yet.
_UNIPROT_RE = re.compile(
    r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}'
)
ACCESSIONS = {
    'SP': (Namespace.UNIPROT, _UNIPROT_RE),
    'TREMBL': (Namespace.UNIPROT, _UNIPROT_RE),
    'GB': (Namespace.GENBANK, re.compile(r'[A-Z]\d{5}|[A-Z]{2}\d{6}')),
    'PIR': ('pir', re.compile(r'[A-Z]\d{5}')),
}
# Gene names in front of a protein name: ``ORC2L;ORC2`` or ``CRK or CRKO``.
_GENE_NAMES_RE = re.compile(r'[\w./-]+(?:(?:;| or | OR )[\w./-]+)*')
_GENE_NAME_SEPARATOR_RE = re.compile(r';| or | OR ')
# A few cells list two or three PubMed IDs, separated by ``;`` or ``,``.
_PMID_SEPARATOR_RE = re.compile(r'[;,]')


# The forms of a UniProt accession that the source writes for a protein it
# names by more than the accession: ``P29354_2`` is the isoform ``P29354-2``,
# and ``DLK_MOUSE`` is an entry name (up to five characters for the protein,
# then the species code, with a letter in it).
_ISOFORM_FORM_RE = re.compile(rf'({_UNIPROT_RE.pattern})_(\d+)')
_ENTRY_NAME_RE = re.compile(r'[A-Z0-9]{1,10}_(?=[A-Z0-9]*[A-Z])[A-Z0-9]{1,5}')
_UNIPROT_DATABASES = ('SP', 'TREMBL')


def accession(row: dict[str, Any], side: str) -> tuple[Any, str] | None:
    """The participant's accession, when the source gives exactly one.

    Of the 3,628 references, 220 are not an accession of their database. Of
    these, two forms say exactly which UniProt record they mean and are read:
    ``P29354_2`` is the isoform ``P29354-2`` (``Namespace.UNIPROT``) and
    ``DLK_MOUSE`` is an entry name (``Namespace.UNIPROT_ENTRY``). The others
    (``-``, a missing ``id``, cells that name alternative proteins and were
    cut in the source such as ``P16092 (mouse) or sp`` and ``P27361 (ERK1) or
    sp``, a lower-case accession, numbers) do not identify one protein, so
    they give no identifier and are not repaired.
    """
    database = row.get(f'{side}_db')
    value = row.get(f'{side}_id') or ''
    namespace, pattern = ACCESSIONS.get(database, (None, None))
    if namespace is None:
        return None
    if pattern.fullmatch(value):
        return namespace, value
    if database in _UNIPROT_DATABASES:
        isoform = _ISOFORM_FORM_RE.fullmatch(value)
        if isoform:
            return Namespace.UNIPROT, f'{isoform.group(1)}-{isoform.group(2)}'
        if _ENTRY_NAME_RE.fullmatch(value):
            return Namespace.UNIPROT_ENTRY, value
    return None


def _name_parts(row: dict[str, Any], side: str) -> tuple[str, str]:
    """A name cell is ``gene names; protein name``; the gene names may be empty."""
    genes, _, name = (row.get(f'{side}_name') or '').partition('; ')
    return genes.strip(), name.strip()


def _gene_names(side: str) -> Callable[[dict[str, Any]], list[str]]:
    def _value(row: dict[str, Any]) -> list[str]:
        genes = _name_parts(row, side)[0]
        if not _GENE_NAMES_RE.fullmatch(genes):
            return []
        return _GENE_NAME_SEPARATOR_RE.split(genes)

    return _value


f = FieldConfig(
    extract={
        'pmid': PMID_EXTRACT,
        # Taxon 0 stands for an unknown species.
        'taxon': r'^([1-9]\d*)$',
    },
    transform={
        'pmid': lambda v: f'PMID:{v}',
        'taxon': lambda v: f'NCBITaxon:{v}',
    },
)


def _participant(side: str) -> EntityBuilder:
    def _pairs(row: dict[str, Any]) -> list[tuple[Any, str]]:
        found = accession(row, side)
        return [found] if found else []

    return EntityBuilder(
        entity_type=Protein,
        identifiers=IdentifiersBuilder(
            CV.from_pairs(_pairs),
            CV(
                term=Namespace.NAME,
                value=lambda row: _name_parts(row, side)[1] or None,
            ),
            CV(term=Namespace.SYNONYM, value=_gene_names(side)),
        ),
        annotations=AnnotationsBuilder(
            # The source's own taxon, which for 93 participants is Mammalia.
            CV(
                term=slots.in_taxon,
                value=f(f'{side}_taxon', extract='taxon', transform='taxon'),
            ),
        ),
    )


# MPPI records physical interactions and does not say which are direct; most
# were found by coimmunoprecipitation or copurification.
interactions_schema = RelationBuilder(
    subject=_participant('a'),
    predicate=slots.interacts_with,
    object=_participant('b'),
    annotations=AnnotationsBuilder(
        CV(term=slots.has_evidence_of_type, value=f('method_mi')),
        CV(
            term=slots.publications,
            value=f(
                lambda row: _PMID_SEPARATOR_RE.split(row.get('pubmed') or ''),
                extract='pmid',
                transform='pmid',
            ),
        ),
    ),
)


def map_interaction(row: dict[str, Any]):
    """A relation for an interaction of two proteins that both have an accession.

    The second row of a pair listed in both orders is not mapped: it repeats
    the first.
    """
    if row.get('participants') != 2 or row.get('mirror'):
        return None
    # A name alone does not identify a participant.
    if not (accession(row, 'a') and accession(row, 'b')):
        return None
    return interactions_schema(row)


resource = Resource(
    config,
    interactions=Dataset(
        download=download_interactions,
        mapper=map_interaction,
        raw_parser=iter_interactions,
    ),
)
