"""Read the GPCR families of GPCRdb.

GPCRdb publishes its reference data in the ``protwis/gpcrdb_data`` repository.
``protein_data/proteins_and_families.txt`` is an indented outline of the
receptors:

    Class A (Rhodopsin) | gpcrdba
        Aminergic receptors
            5-Hydroxytryptamine receptors
                "gpcr","1","5-Hydroxytryptamine receptors","1","5-HT<sub>1A</sub> receptor",...

a receptor class (indented by 0), a ligand type (4), a receptor family (8) and
under it one CSV line (12) for each receptor. The line gives the IUPHAR family
and receptor numbers, the names, and for the human, rat and mouse gene the
HGNC, symbol, name, locus, RefSeq transcript and protein, UniProt accession and
Entrez Gene numbers, with the human gene first. The outline has no identifiers
for the class, the ligand type and the family.

GPCRdb's own identifiers of them are the slugs of its family tree
(``001`` for Class A, ``001_001`` for its aminergic receptors, ``001_001_001``
for the 5-hydroxytryptamine receptors), which the REST service gives for the
whole tree (``https://gpcrdb.org/services/proteinfamily/``: ``slug``, ``name``
and ``parent``). :func:`iter_families` reads both and joins the outline to the
tree by the names from the class down, which every level of the outline has in
the tree, and yields a row for each class, ligand type and family of the
outline.
"""

from __future__ import annotations

from collections.abc import Generator
import csv
import html
import json
import re
from typing import Any

from pypath.inputs_v2.base import Download
from pypath.inputs_v2.parsers.base import read_text

__all__ = ['SPECIES', 'iter_families']

# The columns of a receptor line that hold the UniProt accession of the human,
# the rat and the mouse receptor, and the organism of each.
SPECIES = (
    ('human', 15, '9606'),
    ('rat', 23, '10116'),
    ('mouse', 31, '10090'),
)
# A line of a receptor is the family number, the family name and the receptor
# number before the names; the shortest line is the one of a receptor that has
# only a human gene (up to its Entrez Gene number).
_MIN_FIELDS = 17
# A UniProt accession, with the isoform number some cells give (``P34998-2``).
_ACCESSION = re.compile(
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
    r'(?:-\d+)?'
)
_LEVELS = {0: 'class', 4: 'ligand_type', 8: 'family'}
_MARKUP = re.compile(r'<[^>]*>')


def _plain(name: str) -> str:
    """A name without the HTML markup GPCRdb writes in it: the subscripts
    (``GABA<sub>B</sub> receptors`` is "GABAB receptors") and the character
    entities (``&alpha;``)."""
    return html.unescape(_MARKUP.sub('', name))


def _tree_slugs(tree: Download, **kwargs: Any) -> dict[tuple[str, ...], str]:
    """The slug of each node of the family tree, by the names from its class down.

    The tree has a root, ``000``, "Parent family", which is not a part of the
    path; a name repeated at the same place of the tree is an error.
    """
    rows = json.loads(read_text(tree.open(**kwargs)))
    nodes = {row['slug']: row for row in rows}
    if not nodes:
        raise ValueError('GPCRdb family tree is empty')
    slugs: dict[tuple[str, ...], str] = {}
    for slug, node in nodes.items():
        path = []
        while node.get('parent'):
            path.append(node['name'])
            node = nodes[node['parent']['slug']]
        if not path:
            continue
        key = tuple(reversed(path))
        if key in slugs:
            raise ValueError(f'GPCRdb family tree has two nodes at {key}')
        slugs[key] = slug
    return slugs


def _outline(opener) -> Generator[tuple[str, str], None, None]:
    """The lines of the outline as ``(level, text)``, with ``receptor`` for the
    lines of receptors; blank lines are skipped, an indentation other than 0,
    4, 8 and 12 is an error."""
    for number, line in enumerate(read_text(opener).splitlines(), start=1):
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip(' '))
        if indent in _LEVELS:
            yield _LEVELS[indent], line.strip()
        elif indent == 12:
            yield 'receptor', line.strip()
        else:
            raise ValueError(f'GPCRdb outline line {number} is indented by {indent}')


def iter_families(
    opener, *, tree: Download, **kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Yield the classes, ligand types and families of the outline.

    ``opener`` is the outline, ``tree`` the download of the family tree. Every
    row has ``slug``, ``name`` and ``level`` (``class``, ``ligand_type`` or
    ``family``), then:

    * class and ligand type: ``children``, the slugs and names of the nodes of
      the next level that the outline puts below it, in outline order;
    * family: ``members``, the receptors the outline lists under it, in
      outline order, each with ``accession`` (UniProt), ``taxon`` (NCBI
      taxonomy number of the species of the column it is in) and ``species``.
      A receptor appears once for each accession of each species the line
      has (a cell can hold two, separated by a bar, and an isoform number);
      a column without an accession gives none.

    The class name of the outline ends with ``| gpcrdba`` (a letter that names
    the class the numbering scheme comes from); it is not part of the name. A
    node of the outline that the tree does not have, or that has no receptor
    below it, is an error.
    """
    forced = {key: kwargs[key] for key in ('force_refresh',) if key in kwargs}
    slugs = _tree_slugs(tree, **forced)
    path: list[str] = []
    nodes: dict[tuple[str, ...], dict[str, Any]] = {}
    order: list[tuple[str, ...]] = []

    for level, text in _outline(opener):
        if level == 'receptor':
            if len(path) != 3:
                raise ValueError(f'GPCRdb receptor line outside a family: {text[:60]!r}')
            fields = next(csv.reader([text]))
            if fields[0] != 'gpcr' or len(fields) < _MIN_FIELDS:
                raise ValueError(f'GPCRdb receptor line is not a gpcr record: {text[:60]!r}')
            members = nodes[tuple(path)]['members']
            for species, column, taxon in SPECIES:
                cell = fields[column] if len(fields) > column else ''
                # A receptor with two genes in a species has both accessions
                # in the cell, separated by a bar (the angiotensin AT1
                # receptor of rat and mouse).
                for accession in filter(None, cell.split('|')):
                    if not _ACCESSION.fullmatch(accession):
                        raise ValueError(
                            f'GPCRdb {species} accession {accession!r} in {fields[4]!r}'
                        )
                    member = {
                        'accession': accession, 'taxon': taxon, 'species': species,
                    }
                    if member not in members:
                        members.append(member)
            continue
        name = text.split(' | ')[0].strip() if level == 'class' else text
        depth = {'class': 1, 'ligand_type': 2, 'family': 3}[level]
        path = path[: depth - 1] + [name]
        key = tuple(path)
        if key not in slugs:
            raise ValueError(f'GPCRdb outline node {key} is not in the family tree')
        node = {'slug': slugs[key], 'name': _plain(name), 'level': level}
        if level == 'family':
            node['members'] = []
        else:
            node['children'] = []
        nodes[key] = node
        order.append(key)
        if depth > 1:
            parent = nodes[key[:-1]]
            parent['children'].append({'slug': slugs[key], 'name': _plain(name)})

    for key in order:
        node = nodes[key]
        if not node.get('members', node.get('children')):
            raise ValueError(f'GPCRdb node {key} has nothing below it')
        yield node
    if not order:
        raise ValueError('GPCRdb outline is empty')
