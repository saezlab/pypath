"""
Reader for the SPIKE database XML export (``LatestSpikeDB.xml``).

The format is described in Appendix C of Elkon et al. 2008 (PMID 18289391).
A ``<SpikeDatabase>`` has four blocks:

- ``<SiteBlock>``: the laboratories that contribute data; not read.
- ``<BuildingBlock>``: the entities. ``<Gene>`` (a gene with its RNA and
  protein, cross-referenced to Entrez Gene), ``<SignallingMolecule>`` (a small
  molecule, cross-referenced to ChEBI), ``<Event>`` (a process, a stimulus or
  another entity outside the other types) and ``<Group>``, a complex or a
  family (``type``) with its ``<Member>`` entities.
- ``<RegulationBlock>``: ``<Regulation>`` elements, directed and signed. The
  ``<Target>`` is an entity or another regulation or interaction; the
  ``<PhysicalTarget>`` is the entity the biochemical mechanism acts on.
- ``<InteractionBlock>``: ``<Interaction>`` elements, undirected physical
  links between ``<ProteinA>`` and ``<ProteinB>``.

Elements refer to each other by the numeric ``id`` attribute. The readers
resolve these references and interpret nothing: codes (``effect``,
``integrity``, ``XRef/@db``) are passed on as written, and elements flagged
``isDeleted`` are yielded with ``is_deleted`` set.
"""

from __future__ import annotations

from collections.abc import Callable, Generator
from typing import Any
import xml.etree.ElementTree as ET

from pypath.inputs_v2.parsers.base import read_xml_root

__all__ = [
    'GROUP_TYPES',
    'QUALITY_LEVELS',
    'iter_groups',
    'iter_interactions',
    'iter_regulations',
    'iter_spike_groups',
    'iter_spike_interactions',
    'iter_spike_regulations',
]

# ``integrity``: the quality level of a relation, 1 the highest.
QUALITY_LEVELS = ('1', '2', '3', '4')
# ``Group/@type``: a complex or a family.
GROUP_TYPES = ('Complex', 'Family')

_ENTITY_TAGS = ('Gene', 'SignallingMolecule', 'Event', 'Group')


def _common(element: ET.Element) -> dict[str, Any]:
    """Attributes shared by entities, regulations and interactions."""
    return {
        'id': element.get('id'),
        'kind': element.tag,
        'accession': element.get('spikeAccession'),
        'data_source': element.get('dataSource'),
        'site': element.get('site'),
        'integrity': element.get('integrity'),
        'status': element.get('status'),
        'is_deleted': element.get('isDeleted') == '1',
    }


def _entity(element: ET.Element) -> dict[str, Any]:
    xref = element.find('XRef')
    return {
        **_common(element),
        'group_type': element.get('type') if element.tag == 'Group' else None,
        'name': (element.get('name') or '').strip() or None,
        'long_name': (element.get('longName') or '').strip() or None,
        'xref_db': xref.get('db') if xref is not None else None,
        'xref_id': xref.get('id') if xref is not None else None,
    }


def _description(element: ET.Element) -> str | None:
    return ' '.join((element.findtext('Description') or '').split()) or None


def _relation(element: ET.Element) -> dict[str, Any]:
    pubmed = []
    for reference in element.findall('Reference'):
        pmid = (reference.get('pmid') or '').strip()
        if pmid and pmid not in pubmed:
            pubmed.append(pmid)
    return {
        **_common(element),
        'description': _description(element),
        'pubmed': pubmed,
    }


class _Database:
    """The blocks of a SPIKE XML document, with the entities by id."""

    def __init__(self, opener) -> None:
        root = read_xml_root(opener)
        self.empty = root is None
        if self.empty:
            return
        if root.tag != 'SpikeDatabase':
            raise ValueError(
                f'Expected a SPIKE database XML document, found <{root.tag}>'
            )
        blocks = {
            name: root.find(name)
            for name in ('BuildingBlock', 'RegulationBlock', 'InteractionBlock')
        }
        missing = [name for name, block in blocks.items() if block is None]
        if missing:
            raise ValueError(f'SPIKE XML lacks {", ".join(missing)}')
        self.building, self.regulations, self.interactions = blocks.values()
        # A regulation can target a regulation or an interaction.
        self.relations = {
            element.get('id'): _common(element)
            for block in (self.regulations, self.interactions)
            for element in block
        }
        self.entities: dict[str, dict[str, Any]] = {}
        for element in self.building:
            if element.tag not in _ENTITY_TAGS:
                raise ValueError(f'Unknown SPIKE entity element: <{element.tag}>')
            self.entities[element.get('id')] = _entity(element)
        for element in self.building.iterfind('Group'):
            members = [
                self.reference(member.get('ref'))
                for member in element.findall('Member')
            ]
            self.entities[element.get('id')]['member_kinds'] = [
                {'kind': member.get('kind'), 'name': member.get('name')}
                for member in members
            ]

    def entity(self, element: ET.Element, tag: str) -> dict[str, Any]:
        """The entity that child ``tag`` of a relation refers to."""
        child = element.find(tag)
        if child is None or not child.get('ref'):
            raise ValueError(
                f'SPIKE {element.tag} {element.get("spikeAccession")!r} has no <{tag}>'
            )
        return self.reference(child.get('ref'))

    def reference(self, ref: str) -> dict[str, Any]:
        """The entity or relation with this id; ``kind`` is None if there is none."""
        if ref in self.entities:
            return self.entities[ref]
        if ref in self.relations:
            return self.relations[ref]
        return {'id': ref, 'kind': None}


def iter_spike_groups(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the complexes and families (``<Group>``) of a SPIKE XML document.

    Fields: ``id``, ``kind`` (``Group``), ``group_type`` (``Complex`` or
    ``Family``), ``accession``, ``name``, ``long_name``, ``description``,
    ``data_source``, ``site``, ``integrity``, ``status``, ``is_deleted``,
    ``xref_db``, ``xref_id`` and ``members``: the entities of the ``<Member>``
    elements, in document order, each with the same fields except
    ``description`` and ``members`` (``kind`` is the element name: ``Gene``,
    ``SignallingMolecule``, ``Event`` or ``Group``). A group, here and where
    a relation refers to it, also has ``member_kinds``: the ``kind`` and
    ``name`` of each of its members.
    """
    database = _Database(opener)
    if database.empty:
        return
    for element in database.building.iterfind('Group'):
        yield {
            **database.entities[element.get('id')],
            'description': _description(element),
            'members': [
                database.reference(member.get('ref'))
                for member in element.findall('Member')
            ],
        }


def iter_spike_regulations(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the ``<Regulation>`` elements of a SPIKE XML document.

    Fields: ``id``, ``kind`` (``Regulation``), ``accession``, ``data_source``,
    ``site``, ``integrity``, ``status``, ``is_deleted``, ``mechanism``,
    ``effect``, ``description``, ``pubmed`` (PubMed IDs of the references)
    and ``source``, ``target`` and ``physical_target``: the entities referred
    to, with the fields of a group member. A target that is a regulation or
    an interaction has ``kind`` ``Regulation`` or ``Interaction`` and no name.
    """
    database = _Database(opener)
    if database.empty:
        return
    for element in database.regulations.iterfind('Regulation'):
        yield {
            **_relation(element),
            'mechanism': element.get('mechanism'),
            'effect': element.get('effect'),
            'source': database.entity(element, 'Source'),
            'target': database.entity(element, 'Target'),
            'physical_target': database.entity(element, 'PhysicalTarget'),
        }


def iter_spike_interactions(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the ``<Interaction>`` elements of a SPIKE XML document.

    Fields: ``id``, ``kind`` (``Interaction``), ``accession``,
    ``data_source``, ``site``, ``integrity``, ``status``, ``is_deleted``,
    ``biological_assay``, ``description``, ``pubmed`` and ``protein_a`` and
    ``protein_b``: the entities referred to, with the fields of a group member.
    """
    database = _Database(opener)
    if database.empty:
        return
    for element in database.interactions.iterfind('Interaction'):
        yield {
            **_relation(element),
            'biological_assay': element.get('biologicalAssay'),
            'protein_a': database.entity(element, 'ProteinA'),
            'protein_b': database.entity(element, 'ProteinB'),
        }


def _of_level(rows, level: str) -> Generator[dict[str, Any], None, None]:
    """The rows of one quality level; a row of no level is an error."""
    if level not in QUALITY_LEVELS:
        raise ValueError(f'Not a SPIKE quality level: {level!r}')
    count = 0
    for row in rows:
        count += 1
        if row.get('integrity') not in QUALITY_LEVELS:
            raise ValueError(
                f'Unknown SPIKE quality level: {row.get("integrity")!r} '
                f'({row.get("accession")})'
            )
        if row['integrity'] == level:
            yield row
    if not count:
        raise ValueError('SPIKE XML has no relations of this kind')


def iter_regulations(
    opener,
    level: str,
    statements_of: Callable[[dict[str, Any]], list[str]],
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield one row per statement of the regulations of one quality level.

    ``statements_of`` takes a regulation (the fields of
    :func:`iter_spike_regulations`) and returns the statements it gives, as
    labels, none for a deleted regulation. ``STATEMENT`` is the label; a
    regulation without a statement is yielded once with ``STATEMENT`` None.
    """
    for row in _of_level(iter_spike_regulations(opener, **kwargs), level):
        for statement in statements_of(row) or [None]:
            yield {**row, 'STATEMENT': statement}


def iter_interactions(
    opener,
    level: str,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield the interactions of one quality level."""
    yield from _of_level(iter_spike_interactions(opener, **kwargs), level)


def iter_groups(opener, group_type: str, **kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the groups of one type (``Complex`` or ``Family``)."""
    count = 0
    for row in iter_spike_groups(opener, **kwargs):
        count += 1
        if row.get('group_type') not in GROUP_TYPES:
            raise ValueError(f'Unknown SPIKE group type: {row.get("group_type")!r}')
        if row['group_type'] == group_type:
            yield row
    if not count:
        raise ValueError('SPIKE XML has no complexes or families')
