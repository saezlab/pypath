"""
Parse DOMINO domain-peptide interactions and emit Relation records.

DOMINO is a database of interactions mediated by protein interaction
domains, curated by the MINT group. One row of its file is one experiment
that showed an interaction between two participants, with the detection
method, the publication, the binding site of each participant (residue
ranges and InterPro domains) and the modified residues the binding needs.

About half of the rows are a protein domain binding a synthetic or phage
displayed peptide that the file does not identify, by accession or by
sequence. This module emits the rows in which both participants are
identified, as an interaction between the two, with the binding sites and
modified residues as features of the participants. The residue of a
modification (its one-letter code) is read from the window of sequence that
the note of the modification gives, where the position and the window agree;
see ``modified_residue``.

The DOMINO website is no longer public. This module reads the last release
(22 October 2009) from OmniPath's archive of rescued files.
"""

from __future__ import annotations

import re
from typing import Any

from biolink_model.datamodel.model import Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

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
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.inputs_v2 import _mitab
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.domino import (
    iter_domino,
    parse_binding_site,
    parse_ptms,
)


config = ResourceConfig(
    id=ResourceCv.DOMINO,
    name='DOMINO',
    url='http://mint.bio.uniroma2.it/domino/search/searchWelcome.do',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.CC_BY_2_5,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='17135199',
    annotation_ontologies=(OntologyCv.PSI_MI,),
    topics=(Topic.SIGNALING,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.EXPERIMENTAL,
    ),
    summary=(
        'Protein-protein interactions mediated by interaction domains, with '
        'binding sites, detection methods and publications.'
    ),
    description=(
        'DOMINO is a database of interactions mediated by protein '
        'interaction domains such as SH3, SH2, PDZ, 14-3-3, PTB and WW, '
        'curated from the literature by the MINT group at the University of '
        'Rome Tor Vergata. It follows the PSI-MI standard and records, for '
        'each experiment, the shortest fragments of the two partners shown '
        'to be sufficient for binding, the domain they contain, and the '
        'modified residues and mutations that affect the interaction. Many '
        'of its entries are domains binding synthetic peptides. The last '
        'release is from 2009; its content was merged into MINT.'
    ),
)

# CYGD (the MIPS yeast genome database) names of yeast open reading frames
# (``YPR197C``), an identifier system without a registered namespace.
_CYGD = 'cygd'
_NAMESPACES = {**_mitab.MITAB_NAMESPACES, 'cygd': _CYGD}

# Every kind of alias in the file. A gene name is the gene symbol of UniProt;
# the synonyms mix gene synonyms with protein names.
_ALIAS_NAMESPACES = {
    'gene name': Namespace.GENESYMBOL,
    'gene name synonym': Namespace.SYNONYM,
    'orf name': Namespace.SYNONYM,
    'locus name': Namespace.SYNONYM,
    'isoform synonym': Namespace.SYNONYM,
}

# The interactor type columns hold plain words. A peptide is a fragment of a
# protein or a synthetic peptide; none in the file has an identifier.
_INTERACTOR_TYPES = {'protein': Protein, 'peptide': Protein}

# Interaction types of the file outside the direct interaction and
# association types that the shared MITAB rules map.
_BROAD_INTERACTIONS = {
    # physical interaction: obsolete, "can be direct or indirect"; PSI-MI
    # split it into association and physical association.
    'MI:0218',
}
_NOT_INTERACTIONS = {
    # colocalization: "a physical interaction among those molecules cannot
    # be inferred".
    'MI:0403',
}

_MINT_RECORD_RE = re.compile(r'MINT-\d+')
_ISOFORM_RE = re.compile(
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
    r'-\d+'
)
# What DOMINO calls the feature its ``binding site`` columns describe.
_BINDING_SITE = 'binding site'
# The residues that each modification, as DOMINO names it, can be on. The
# generic "phosphorylated residue" is on serine, threonine or tyrosine in
# the file. A modification not listed here gets no residue.
_RESIDUES = {
    'phosphorylated residue': 'STY',
    'o-phospho-serine': 'S',
    'o-phospho-threonine': 'T',
    "o4'-phospho-tyrosine": 'Y',
    'n6-methyl-lysine': 'K',
    'n6,n6-dimethyl-lysine': 'K',
    'n6,n6,n6-trimethyl-lysine': 'K',
}
# A note's window of the sequence is 20 residues: the modified one, the ten
# before it and the nine after it. It is shorter where the protein ends, and
# starts at the first residue when the modification is among the first ten.
_WINDOW_BEFORE = 10
_WINDOW_LENGTH = 20


def modified_residue(ptm: dict[str, Any]) -> str | None:
    """The one-letter code of the modified residue, from the note's window.

    The window states the sequence around the modified residue but not where
    in it the residue is; the layout of the 20-residue window (see above)
    does. The residue is read from there only when everything agrees: the
    modification has one position, the note gives that same single position,
    the window is no longer than the layout allows and reaches the residue,
    and the residue is one the modification can be on (``_RESIDUES``). In any
    other case, and when there is no window, the residue is None.
    """
    window = ptm.get('window')
    position = ptm['start']
    allowed = _RESIDUES.get(ptm['name'] or '')
    if not window or position is None or allowed is None:
        return None
    if not (ptm['end'] == window['start'] == window['end'] == position):
        return None
    offset = _WINDOW_BEFORE if position > _WINDOW_BEFORE else position - 1
    sequence = window['sequence']
    if offset >= len(sequence) or len(sequence) > min(
        _WINDOW_LENGTH, offset + 10
    ):
        return None
    return sequence[offset] if sequence[offset] in allowed else None


def _identifier_pairs(side: str):
    """The participant's identifier, then its gene name and synonyms."""

    def _pairs(row: dict[str, Any]) -> list[tuple[Any, str]]:
        pairs = _mitab.identifier_pairs(row.get(f'id_{side}'), _NAMESPACES)
        if not pairs:
            return []
        for item in _mitab.parse_mitab_cell(row.get(f'aliases_{side}')):
            namespace = _ALIAS_NAMESPACES.get(item.description or '')
            if namespace is not None and item.database == 'uniprotkb':
                pairs.append((namespace, item.accession))
        return list(dict.fromkeys(pairs))

    return _pairs


def _entity_type(side: str):
    column = f'interactor_type_{side}'
    return lambda row: _INTERACTOR_TYPES.get(str(row.get(column) or '').strip())


def _coordinate(row: dict[str, Any], side: str) -> dict[str, Any]:
    """Positions count residues of the participant, from 1.

    The sequence is named only when the participant is an isoform: an entry
    accession does not say which of its sequences the curator numbered.
    """
    isoform = next(
        (
            item.accession
            for item in _mitab.parse_mitab_cell(row.get(f'id_{side}'))
            if item.database == 'uniprotkb'
            and _ISOFORM_RE.fullmatch(item.accession)
        ),
        None,
    )
    return {
        'identifier': {'ns': 'uniprot', 'id': isoform} if isoform else None,
        'coordinate_system': 'protein',
        'position_base': 1,
    }


def _participant_form(side: str):
    """Binding site regions and modified residues of the participant.

    A binding site with one InterPro entry gives a region per residue range,
    each named by the entry. With several entries the file does not say
    which range is which domain: the ranges stay unnamed and every entry
    gives a region without a position. So does the entry of a site whose
    position is not known.
    """

    def _form(row: dict[str, Any]) -> dict[str, Any] | None:
        coordinate = _coordinate(row, side)
        regions = []
        site = parse_binding_site(row.get(f'binding_site_{side}'))
        if site:
            entries = [
                {'ns': 'interpro', 'id': accession}
                for accession, _ in site['domains']
            ]
            entry = entries[0] if len(entries) == 1 else None
            regions += [
                {
                    'type': _BINDING_SITE,
                    'identifier': entry,
                    'position': start,
                    'end_position': end,
                    'coordinate_reference': coordinate,
                }
                for start, end in site['ranges']
            ]
            if not site['ranges'] or entry is None:
                regions += [
                    {'type': _BINDING_SITE, 'identifier': identifier}
                    for identifier in entries
                ]
        modifications = [
            {
                # The modification as DOMINO names it; its PSI-MI accession
                # is an obsolete term.
                'term': ptm['name'],
                'residue': modified_residue(ptm),
                'position': ptm['start'],
                'end_position': ptm['end'],
                'coordinate_reference': (
                    coordinate if ptm['start'] is not None else None
                ),
            }
            for ptm in parse_ptms(row.get(f'ptm_{side}'))
            if ptm['name']
        ]
        if not regions and not modifications:
            return None
        return {'regions': regions, 'modifications': modifications}

    return _form


def _interaction_records(row: dict[str, Any]) -> list[tuple[Any, str]]:
    """The MINT and IMEx records of the row.

    The peptide rows have numbers of DOMINO itself under the same database
    name (``mint:DOMINO-1234``); those are not MINT accessions.
    """
    return [
        (namespace, accession)
        for namespace, accession in _mitab.identifier_pairs(
            row.get('interaction_ids')
        )
        if namespace != Namespace.MINT or _MINT_RECORD_RE.fullmatch(accession)
    ]


def _publications(row: dict[str, Any]) -> list[str]:
    """PubMed identifiers; a DOI for the articles the file cites by DOI."""
    cell = row.get('publications')
    return _mitab.pmid_curies(cell) + [
        f'doi:{item.accession}'
        for item in _mitab.parse_mitab_cell(cell)
        if item.database == 'doi' and item.accession.startswith('10.')
    ]


def _predicate(row: dict[str, Any]):
    """
    Predicate for the PSI-MI interaction type of the row.

    A direct interaction or one of its kinds (the enzymatic reactions) is a
    direct physical interaction; an association, a physical association or
    the obsolete "physical interaction" says only that the participants
    interact. Colocalization and any other type give None.
    """
    types = _mitab.mi_terms(row.get('interaction_types'))
    if not types or set(types) & _NOT_INTERACTIONS:
        return None
    predicate = _mitab.physical_predicate(types)
    if predicate is None and set(types) <= _BROAD_INTERACTIONS:
        return slots.interacts_with
    return predicate


def _participant(side: str) -> EntityBuilder:
    return EntityBuilder(
        entity_type=_entity_type(side),
        molecular_form=_participant_form(side),
        identifiers=IdentifiersBuilder(
            CV.from_pairs(_identifier_pairs(side)),
        ),
        annotations=AnnotationsBuilder(
            *_mitab.participant_annotations(
                taxon=f'taxid_{side}',
                biological_role=f'biological_role_{side}',
                experimental_role=f'experimental_role_{side}',
            ),
        ),
    )


interactions_schema = RelationBuilder(
    subject=_participant('a'),
    predicate=_predicate,
    object=_participant('b'),
    identifiers=IdentifiersBuilder(CV.from_pairs(_interaction_records)),
    annotations=AnnotationsBuilder(
        *_mitab.interaction_annotations(
            detection_method='detection_methods',
            interaction_type='interaction_types',
            source_database='source_databases',
            host_organism='host_organism',
        ),
        CV(term=slots.publications, value=_publications),
    ),
)


def exclusion_reason(row: dict[str, Any]) -> str | None:
    """Why a row gives no relation; None for a row that is mapped."""
    # A negative result must never become an interaction. The column holds
    # ``-`` or the sentence "the interactions are not observed under these
    # experimental conditions"; the caution column repeats it.
    if str(row.get('negative') or '').strip() not in ('', '-') or (
        'negative:' in str(row.get('caution') or '')
    ):
        return 'negative result'
    identified = {side: bool(_identifier_pairs(side)(row)) for side in 'ab'}
    for side in 'ab':
        if _entity_type(side)(row) is None:
            return 'unknown interactor type'
        peptide = str(row.get(f'interactor_type_{side}')).strip() == 'peptide'
        if peptide and not identified[side]:
            return 'peptide without an identifier'
    if not all(identified.values()):
        return 'participant without a known identifier'
    types = _mitab.mi_terms(row.get('interaction_types'))
    if set(types) & _NOT_INTERACTIONS:
        return 'colocalization'
    if _predicate(row) is None:
        return 'unknown interaction type'
    return None


def map_interaction(row: dict[str, Any]):
    return None if exclusion_reason(row) else interactions_schema(row)


resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            url=(
                'https://rescued.omnipathdb.org/'
                'domino_2009-10-22-full-binary.mitab26.gz'
            ),
            filename='domino_2009-10-22-full-binary.mitab26.gz',
            subfolder='domino',
            large=True,
            ext='gz',
        ),
        mapper=map_interaction,
        raw_parser=iter_domino,
    ),
)
