"""
Parse PepCyber:P~Pep interactions and emit Relation records.

PepCyber:P~Pep curated interactions between proteins with a phosphoprotein-
binding domain (PPBD: SH2, 14-3-3, PTB, FHA, BRCT, WW, ...) and their
substrates, the proteins that carry the phosphorylated peptide the domain
binds. Each interaction names the binding-domain family, the two proteins, the
binding peptide of the substrate with its phosphorylated residue and that
residue's position, the types of evidence and the source's category of the
interaction. The publications of the interaction's records are on the detail
pages, only some of which survive.

The module emits one relation per interaction: the PPBD protein is the
subject, and the substrate the object, whose molecular form carries the
phosphorylated residue where the peptide marks one. The website is gone; this
reads the copies that OmniPath keeps in its archive of rescued files.

The binding sequence is a relation attribute; the number of records of an
interaction and the descriptions and synonyms on the detail pages stay in the
parsed record.
"""

from __future__ import annotations

from collections.abc import Callable
import functools
import re
from typing import Any

from biolink_model.datamodel.model import Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.pepcyber import iter_pepcyber
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.PEPCYBER,
    name='PepCyber',
    full='PepCyber:P~Pep',
    url='http://www.pepcyber.org/PPEP/',
    # Gong et al. 2008 (PMC2238930): "available, free of charge, to
    # researchers" at academic and non-profit institutions. The site, which
    # said "All rights reserved" (2006-2010), is offline and states no other terms.
    license=LicenseCV.NONCOMMERCIAL,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='18160410',
    topics=(Topic.SIGNALING,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Interactions of phosphoprotein-binding domain proteins (SH2, 14-3-3, '
        'PTB, WW...) with phosphorylated substrates, with the phosphosite.'
    ),
    description=(
        'PepCyber:P~Pep was a database of interactions mediated by '
        'phosphoprotein-binding domains (PPBDs). It listed 11,269 records of '
        'interactions between 387 PPBD proteins and 1,471 substrate proteins '
        'in its last release (version 1.2, May 2010), curated from the '
        'literature, with the binding domain family, the phosphopeptide and '
        'the phosphorylated residue, the types of evidence and a category that '
        'says whether the domain is known to mediate the interaction. The '
        'website is no longer online.'
    ),
)

# OmniPath's copies of the search result page for all interactions, the
# directory of the detail pages that were rescued (704 of 9,530) and the pages.
download = Download(
    url='https://rescued.omnipathdb.org/pepcyber_search_result.html',
    filename='pepcyber_search_result.html',
    subfolder='pepcyber',
    large=True,
    encoding=None,
    default_mode='rb',
)
download_detail_listing = Download(
    url='https://rescued.omnipathdb.org/pepcyber_details/',
    filename='pepcyber_details.html',
    subfolder='pepcyber',
    large=True,
    encoding=None,
    default_mode='rb',
)
download_detail = Download(
    url=lambda number, **_: (
        f'https://rescued.omnipathdb.org/pepcyber_details/pepcyber_{number}.html'
    ),
    filename=lambda number, **_: f'pepcyber_{number}.html',
    subfolder='pepcyber/details',
    large=True,
    encoding=None,
    default_mode='rb',
)

_ACCESSION_RE = re.compile(
    r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
    r'(?:-\d+)?$'
)
# A RefSeq protein accession. The detail pages also give GenBank accessions
# (``AAA92686``) under the same heading, which are not RefSeq.
_REFSEQ_RE = re.compile(r'^[NXY]P_\d+(?:\.\d+)?$')
# PepCyber's numbers of its interactions, an identifier system without a
# registered namespace.
_PEPCYBER_INTERACTION = 'pepcyber'
# The binding domain family the source assigns to the PPBD protein ("PPBD
# class"), as an attribute of the relation; Biolink has no slot for it.
_BINDING_DOMAIN = 'pepcyber:binding_domain'
# The source's category of the interaction (Gong et al. 2008): A, the interaction
# is known to be mediated by a PPBD and the substrate peptide is identified; B,
# mediated by a PPBD, peptide not identified; C, not known whether a PPBD mediates it.
_CATEGORY = 'pepcyber:category'
_EVIDENCE_TYPE = 'pepcyber:evidence_type'
_BINDING_SEQUENCE = 'pepcyber:binding_sequence'

# The evidence types of the "Evidence" column, a comma separated list, by the
# lower case of the type as written. The source's own spelling of each is kept,
# including "Refernce" (four rows) and "NMR andITC".
_EVIDENCE_TYPES = {
    'inference': 'Inference',
    'mutagenesis': 'Mutagenesis',
    'peptide library': 'Peptide library',
    'peptide affinity isolation': 'Peptide affinity isolation',
    'structure': 'Structure',
    'pull down': 'Pull down',
    'mappit': 'MAPPIT',
    'immunoprecipitation': 'Immunoprecipitation',
    'other': 'Other',
    'refernce': 'Refernce',
    # One row each.
    'competitive inhibition assay': 'Competitive inhibition assay',
    'predicted': 'Predicted',
    'fluorescence polarization': 'Fluorescence Polarization',
    'nmr spectroscopy': 'NMR spectroscopy',
    'nmr anditc': 'NMR andITC',
    'dot blot assay': 'Dot blot assay',
}
# The categories of the source. A and B say that the interaction is mediated by a
# phosphoprotein-binding domain; C is not known to be. None of
# the categories says that the two proteins touch, and the evidence types
# include immunoprecipitation, pull down and inference, so the category gives
# ``interacts_with``. An experimental structure among the evidence types of
# the interaction shows the domain and the phosphopeptide in contact, which is
# a direct physical interaction.
_CATEGORIES = frozenset({'A', 'B', 'C'})
_STRUCTURE = 'Structure'


def _evidence_types(row: dict[str, Any]) -> list[str]:
    types = []
    for part in row.get('evidence') or ():
        try:
            types.append(_EVIDENCE_TYPES[part.strip().lower()])
        except KeyError:
            raise ValueError(f'Unknown PepCyber evidence type {part!r}') from None
    return list(dict.fromkeys(types))


def _predicate(row: dict[str, Any]) -> Any:
    if row.get('category') not in _CATEGORIES:
        raise ValueError(f'Unknown PepCyber category {row.get("category")!r}')
    if _STRUCTURE in _evidence_types(row):
        return slots.directly_physically_interacts_with
    return slots.interacts_with


def _publications(row: dict[str, Any]) -> list[str]:
    records = (row.get('detail') or {}).get('records') or ()
    return list(dict.fromkeys(c for r in records if (c := pmid_curie(r['pmid']))))


def _detail(row: dict[str, Any], role: str) -> dict[str, Any]:
    return (row.get('detail') or {}).get(role) or {}


def _accession(row: dict[str, Any], role: str) -> str | None:
    value = _detail(row, role).get('swissprot')
    return value if value and _ACCESSION_RE.match(value) else None


def _refseq(row: dict[str, Any], role: str) -> str | None:
    value = _detail(row, role).get('refseq')
    return value if value and _REFSEQ_RE.match(value) else None


def _phosphorylated_form(row: dict[str, Any]) -> dict | None:
    """The phosphorylated residue of the substrate, where the peptide marks it.

    The position is the row's index site. Where the detail page gives the
    substrate's RefSeq sequence, the position counts that sequence (the
    residue at it is the marked one in 646 of the 653 rows that can be
    checked); otherwise the sequence is not named.
    """
    residue, position = row.get('phospho_residue'), row.get('index_site')
    if residue is None or position is None:
        return None
    refseq = _refseq(row, 'substrate')
    return {
        'modifications': [
            {
                'term': 'phosphorylated residue',
                'residue': residue,
                'position': position,
                'end_position': position,
                'coordinate_reference': {
                    'identifier': {'ns': Namespace.REFSEQ_PROTEIN, 'id': refseq}
                    if refseq
                    else None,
                    'coordinate_system': 'protein',
                    'position_base': 1,
                },
                'description': None,
            }
        ],
    }


def _ppbd_identifiers(row: dict[str, Any]) -> list[tuple[str, str]]:
    pairs = [
        (Namespace.UNIPROT, _accession(row, 'ppbd')),
        (Namespace.REFSEQ_PROTEIN, _refseq(row, 'ppbd')),
        (Namespace.GENESYMBOL, row.get('ppbd_symbol')),
    ]
    return [(ns, v) for ns, v in pairs if v]


def _substrate_identifiers(row: dict[str, Any]) -> list[tuple[str, str]]:
    pairs = [
        (Namespace.UNIPROT, _accession(row, 'substrate')),
        (Namespace.ENTREZ, row.get('substrate_entrez')),
        (Namespace.REFSEQ_PROTEIN, _refseq(row, 'substrate')),
        (Namespace.GENESYMBOL, row.get('substrate_symbol')),
    ]
    return [(ns, v) for ns, v in pairs if v]


def _participant(
    identifiers: Callable[[dict[str, Any]], list[tuple[str, str]]],
    form: Callable[[dict[str, Any]], Any] | None,
    has_identity: Callable[[dict[str, Any]], bool],
) -> Callable[[dict[str, Any]], Any]:
    builder = EntityBuilder(
        entity_type=Protein,
        molecular_form=form,
        identifiers=IdentifiersBuilder(CV.from_pairs(identifiers)),
    )

    def build(row: dict[str, Any]) -> Any:
        return builder.build(row) if has_identity(row) else None

    return build


# The substrate is known by the Gene ID of its link, or by the UniProt accession
# of the detail page. Without either, the source's name is free text: a viral
# protein, a fusion gene, "Phosphatidic Acid", or a character-set artifact. The
# PPBD protein is named by its Entrez gene name (one row has none).
_ppbd = _participant(_ppbd_identifiers, None, lambda row: bool(row.get('ppbd_symbol')))
_substrate = _participant(
    _substrate_identifiers,
    _phosphorylated_form,
    lambda row: bool(row.get('substrate_entrez') or _accession(row, 'substrate')),
)


def _interaction_number(row: dict[str, Any]) -> str:
    return str(row['number'])


interactions_schema = RelationBuilder(
    subject=_ppbd,
    predicate=_predicate,
    object=_substrate,
    identifiers=IdentifiersBuilder(
        CV(term=_PEPCYBER_INTERACTION, value=_interaction_number),
    ),
    annotations=AnnotationsBuilder(
        CV(term=_BINDING_DOMAIN, value=lambda row: row.get('ppbd_class')),
        CV(term=_CATEGORY, value=lambda row: row.get('category')),
        CV(term=_EVIDENCE_TYPE, value=_evidence_types),
        CV(term=_BINDING_SEQUENCE, value=lambda row: row.get('binding_sequence')),
        CV(term=slots.publications, value=_publications),
    ),
)


resource = Resource(
    config,
    interactions=Dataset(
        download=download,
        mapper=interactions_schema,
        raw_parser=functools.partial(
            iter_pepcyber,
            details=download_detail,
            detail_numbers=download_detail_listing,
        ),
    ),
)
