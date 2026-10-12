"""
Parse switches.ELM and emit Relation records of its motif-mediated interactions.

switches.ELM curates experimentally validated molecular switches in which a
short linear motif mediates an interaction between two proteins and the
interaction is turned on or off by a modification of one of the proteins
(``PTM-dependent``), by the binding of an effector, or by the inclusion or
removal of the motif in alternative isoforms. The table of the switches
(see ``parsers/switches_elm.py``) has one line for each interaction of a
switch: the two participants with the binding sites that mediate the
interaction (the motif class and its residues in A, the domain and its
residues in B), the type, subtype and mechanism of the switch, its outcome (the
interaction is induced, enhanced, inhibited or abrogated), the modification
and its sites, and the publications.

This module emits the interaction between the two proteins, with the binding
sites as regions of the participants, as ``elm.py`` does for ELM. The
interaction exists in one state of the switch. Where a modification induces or
enhances the interaction (the line says the outcome is ``Induction`` or
``Enhancement`` and the switch ``PTM-dependent``), the form of the participant
that the modification is added to carries the modified residues, as the
interaction needs them. Where the modification inhibits or abrogates the
interaction, or the line has an effector, a splice form or no modification,
the interaction is that of the other state and the form carries no
modification. The other columns of the table describe the switch, not the
interaction, and have no slot: they are not mapped.

Interactions of a participant that is not a protein (a small molecule, a
nucleic acid, a membrane) or a complex formed by another interaction of the
table, and the lines the table marks inactive, give no relation.
"""

from __future__ import annotations

import re
from typing import Any

from biolink_model.datamodel.model import Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import is_pmid, pmid_curie
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.switches_elm import iter_switches
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.SWITCHES_ELM,
    name='switches.ELM',
    full='switches.ELM: a compendium of conditional regulatory interaction interfaces',
    url='http://switches.elm.eu.org/',
    # "Data downloaded from switches.ELM is only free for academic use, not for
    # commercial use." (usage policy, http://switches.elm.eu.org/index.php?page=help)
    license=LicenseCV.NONCOMMERCIAL,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='23550212',
    topics=(Topic.SIGNALING,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Motif-mediated interactions between two proteins that a modification, '
        'an effector or alternative splicing switches on or off, with the '
        'binding sites, the modified residues and publications.'
    ),
    description=(
        'switches.ELM, of the ELM consortium at EMBL, curates experimentally '
        'validated motif-based molecular switches from the literature: '
        'interactions of a short linear motif with its binding domain whose '
        'occurrence is regulated, by a post-translational modification, the '
        'binding of an effector, competition between binding partners or '
        'alternative splicing. Each switch has a type and a subtype, the '
        'interactions it consists of with the binding sites and their '
        'positions, the modified residues and the enzymes that modify '
        'them, localisation, cell cycle phase and pathway, and the '
        'publications. The enzymes, localisation, cell cycle phase and '
        'pathway are in part inferred from other databases.'
    ),
)

# The file is the one the legacy reader used; the site still offers it.
download = Download(
    url='http://switches.elm.eu.org/downloads/switches.ELM-v1.txt',
    filename='switches.ELM-v1.txt',
    subfolder='switches_elm',
)

# switches.ELM's own numbers (``INTI000001``, ``SWTI000001``).
INTERACTION_ID = 'switches_elm_interaction'
SWITCH_ID = 'switches_elm_switch'
# ``Status`` -> whether the line is in use. The four lines of the table that
# are not (two switches) give no relation.
STATUSES = {'Active': True, 'Inactive': False}
# ``Switch mechanism``. Only ``PTM-dependent`` says that a modification
# switches the interaction.
MECHANISMS = frozenset(
    {
        'PTM-dependent',
        'Effector binding',
        'Alternative splicing',
        # The same mechanism, with another capitalisation.
        'Alternative Splicing',
        'Alternative promoter usage',
        'Alternative promoter usage and alternative splicing',
        'Alternative initiation',
        'RNA editing',
    }
)
# ``Switch outcome`` -> whether the interaction occurs in the modified (or
# effector bound) state. Induction and enhancement are the outcomes in which
# the modification makes or strengthens the interaction; with the others the
# interaction is that of the unmodified state.
OUTCOMES = {
    'Induction': True,
    'Enhancement': True,
    'Abrogation': False,
    'Inhibition': False,
}
# A modification of a residue is a PSI-MOD term (``MOD:00696``, phosphorylated
# residue). One line writes ``Proteolytic cleavage``, which names no residue
# modification and is not mapped.
_MOD_TERM = re.compile(r'MOD:\d+')
_BINDING_SITE = 'binding site'
# ELM class identifiers start with the type of the class.
_ELM_CLASS = re.compile(r'ELM:((?:CLV|DEG|DOC|LIG|MOD|TRG)_[\w-]+)')
_PFAM = re.compile(r'PFAM:(PF\d{5})')
_INTERPRO = re.compile(r'IPR:(IPR\d{6})')


def _check(row: dict[str, Any]) -> None:
    """Raise on a value of a categorical column that the crosswalks lack."""
    for name, value, known in (
        ('status', row.get('status'), STATUSES),
        ('switch mechanism', row.get('switch_mechanism'), MECHANISMS),
        ('switch outcome', row.get('outcome'), OUTCOMES),
    ):
        if (value is not None or name == 'status') and value not in known:
            raise ValueError(f'Unknown switches.ELM {name}: {value!r}')


def _site_identifier(text: str | None) -> dict[str, str] | None:
    """The ELM class, Pfam or InterPro accession of a binding site. Others
    (``coiled_coil``, ``OTHER:coiled coil``, a Pfam-B family) have none."""
    for pattern, namespace in (
        (_ELM_CLASS, 'elm'),
        (_PFAM, 'pfam'),
        (_INTERPRO, 'interpro'),
    ):
        match = pattern.fullmatch(text or '')
        if match:
            return {'ns': namespace, 'id': match[1]}
    return None


def _form(role: str):
    """The molecular form of participant ``role`` (``a`` or ``b``): its
    binding sites and, where the line makes the interaction depend on it, the
    modified residues."""

    def form(row: dict[str, Any]) -> dict[str, Any] | None:
        participant = row[role]
        accession = participant['id']
        # The sequence is named only when the participant is an isoform: an
        # entry accession does not say which of its sequences the positions
        # count.
        coordinate = {
            'identifier': {'ns': 'uniprot', 'id': accession} if '-' in accession else None,
            'coordinate_system': 'protein',
            'position_base': 1,
        }
        regions = []
        for site in participant['binding_sites']:
            identifier = _site_identifier(site['id'])
            located = site['start'] is not None
            regions.append(
                {
                    'type': _BINDING_SITE,
                    'identifier': identifier,
                    'position': site['start'],
                    'end_position': site['end'],
                    'coordinate_reference': coordinate if located else None,
                    'description': None if identifier else site['id'],
                }
            )
        modifications = _modifications(row, role, coordinate)
        if not (regions or modifications):
            return None
        return {'regions': regions or None, 'modifications': modifications or None}

    return form


def _modifications(
    row: dict[str, Any], role: str, coordinate: dict[str, Any]
) -> list[dict[str, Any]]:
    """The modified residues that the interaction of the line needs, on the
    participant that the modification is added to."""
    mechanism, outcome = row.get('switch_mechanism'), row.get('outcome')
    term = row.get('modification')
    if (
        mechanism != 'PTM-dependent'
        or not OUTCOMES.get(outcome)
        or row.get('affected') != role.upper()
        or not _MOD_TERM.fullmatch(term or '')
    ):
        return []
    return [
        {
            'term': term,
            'residue': site['residue'],
            'position': site['position'],
            'end_position': site['position'],
            'coordinate_reference': coordinate,
            'description': f'{site["residue"]}{site["position"]}',
        }
        for site in row['modification_sites']
    ]


def _protein(role: str):
    builder = EntityBuilder(
        entity_type=Protein,
        molecular_form=_form(role),
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.UNIPROT, value=lambda row: row[role]['id']),
        ),
    )

    def build(row: dict[str, Any]) -> Any:
        # Only a protein: another participant (a small molecule, a nucleic
        # acid, a membrane, a complex) is not an entity of this module.
        if not STATUSES[row['status']] or row[role]['kind'] != 'uniprot':
            return None
        return builder.build(row)

    return build


def _publications(row: dict[str, Any]) -> list[str]:
    return [pmid_curie(pmid) for pmid in row['pmids'] if is_pmid(pmid)]


# The protein with the motif is the subject and the protein with the domain
# the object, as in ``elm.py``. The table says that binding interfaces mediate
# the interaction, not that the contact is direct, so it is ``interacts_with``.
interactions_schema = RelationBuilder(
    subject=_protein('a'),
    predicate=slots.interacts_with,
    object=_protein('b'),
    identifiers=IdentifiersBuilder(
        CV(term=INTERACTION_ID, value=lambda row: row['interaction_id']),
        CV(term=SWITCH_ID, value=lambda row: row['switch_id']),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.publications, value=_publications),
    ),
)

def map_interaction(row: dict[str, Any]) -> Any:
    """The relation of one line; None for a line that is not of two proteins
    or is inactive. The builder drops an error of a participant, so the values
    are checked first."""
    _check(row)
    return interactions_schema.build(row)


resource = Resource(
    config,
    interactions=Dataset(
        download=download,
        mapper=map_interaction,
        raw_parser=iter_switches,
    ),
)
