"""Read the classification of the human membrane proteome of Almen et al. 2009.

Additional file 1 of the paper is an ``.xls`` workbook; its ``Data`` sheet has one
row per human membrane protein, identified by an IPI accession:

    IPI Accession  IPI Description  Main Class  Classification  Protein coding
    Phobius TM  Phobius SP  Phobius Topology  SOSUI TM  SOSUI TOPOLOGY
    TMHMM TM  TMHMM Topology

``Main Class`` is one of five classes. ``Classification`` is the path of the
protein through the classes below it, as the paper's legend says: "EC number is
assigned for enzymes and semicolon separated classes and subclasses for all
other proteins" (``GPCR;Rhodopsin``, ``SLC;APC;SLC12``). The enzyme classes are
EC numbers of any depth (``2.7.11``), and a protein with two enzymatic
activities has two, joined by ``/``.

The in silico surfaceome table (``parsers/surfaceome.py``) carries the same
classification for its proteins, written in the same way and by UniProt
accession, so the functions that read a path are shared with it. This file is
read for the classes (:func:`iter_classes`); the proteins are read from that
table, because IPI accessions no longer exist and nothing resolves them.
"""

from __future__ import annotations

from collections.abc import Generator, Iterable
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_excel

COLUMNS = (
    'IPI Accession',
    'IPI Description',
    'Main Class',
    'Classification',
    'Protein coding',
    'Phobius TM',
    'Phobius SP',
    'Phobius Topology',
    'SOSUI TM',
    'SOSUI TOPOLOGY',
    'TMHMM TM',
    'TMHMM Topology',
)

MAIN_CLASSES = (
    'Enzymes',
    'Receptors',
    'Transporters',
    'Miscellaneous',
    'Unclassified',
)
ENZYMES = 'Enzymes'
UNCLASSIFIED = 'Unclassified'

# The two ways a classification is read: a path of the source's own classes, or
# the EC numbers of an enzyme.
CLASS = 'class'
EC = 'ec'

# Known defects of the workbook, repaired where the value is read. Both are
# spelling slips of a single class: the typed class has no other member, so the
# repaired value is the class the other rows of its family use (``Kinase``), and
# the EC class ``3.4.`` ends in a dot that no other EC number of the table has.
# The legacy OmniPath reader repaired the first one too.
CLASS_REPAIRS = {'KInase;PDGFR': 'Kinase;PDGFR'}
EC_REPAIRS = {'3.4.': '3.4'}

_EC_NUMBER = re.compile(r'^[1-7](\.\d+){0,3}$')


def class_path(main_class: str, classification: str) -> tuple[str, ...]:
    """The classes of a protein below its main class, from the top.

    A protein of the ``Unclassified`` main class has no classes below it
    (the source repeats the name in ``Classification``).
    """
    if main_class not in MAIN_CLASSES or main_class == ENZYMES:
        raise ValueError(f'Not a main class with a class path: {main_class!r}')
    if main_class == UNCLASSIFIED:
        if classification != UNCLASSIFIED:
            raise ValueError(
                f'Unclassified protein with classification {classification!r}'
            )
        return ()
    classification = CLASS_REPAIRS.get(classification, classification)
    path = tuple(classification.split(';'))
    if not all(path):
        raise ValueError(f'Empty class in classification {classification!r}')
    return path


def node_id(path: tuple[str, ...], depth: int) -> str:
    """Identifier of the class at ``depth`` (1 is the top one) of a path: the
    path to it as the source writes it."""
    return ';'.join(path[:depth])


def ec_numbers(classification: str) -> list[str]:
    """The EC numbers of an enzyme's classification, each without ``EC:``."""
    numbers = [
        EC_REPAIRS.get(number, number) for number in classification.split('/')
    ]
    for number in numbers:
        if not _EC_NUMBER.match(number):
            raise ValueError(f'Not an EC number: {number!r} in {classification!r}')
    return numbers


def statements(main_class: str, classification: str) -> list[tuple[str, str]]:
    """The classes a protein is placed in, as ``(kind, term)`` pairs.

    The most specific class of the source's own hierarchy is taken (the main
    class for an enzyme or an unclassified protein), and each EC number of an
    enzyme. The classes above the most specific one follow from the hierarchy
    of the ``terms``.
    """
    if main_class == ENZYMES:
        return [(CLASS, ENZYMES)] + [
            (EC, number) for number in ec_numbers(classification)
        ]
    path = class_path(main_class, classification)
    return [(CLASS, node_id(path, len(path)) if path else main_class)]


def class_nodes(
    pairs: Iterable[tuple[str, str]],
) -> Generator[dict[str, Any], None, None]:
    """The classes of the source's own hierarchy used by ``(main class,
    classification)`` pairs: every main class, and the nodes of every path.

    Yields ``id`` (the path as the source writes it, or the main class name),
    ``name`` (the path with `` > `` between levels), ``label`` (the last class)
    and ``parent`` (the id of the class above, None for a main class). A class
    path always has the same main class above it.
    """
    owner: dict[str, str] = {}
    nodes: dict[str, dict[str, Any]] = {
        name: {'id': name, 'name': name, 'label': name, 'parent': None}
        for name in MAIN_CLASSES
    }
    for main_class, classification in pairs:
        if main_class == ENZYMES:
            continue
        path = class_path(main_class, classification)
        if path and owner.setdefault(path[0], main_class) != main_class:
            raise ValueError(
                f'Class {path[0]!r} is under two main classes: '
                f'{owner[path[0]]!r} and {main_class!r}'
            )
        for depth in range(1, len(path) + 1):
            identifier = node_id(path, depth)
            nodes.setdefault(
                identifier,
                {
                    'id': identifier,
                    'name': ' > '.join(path[:depth]),
                    'label': path[depth - 1],
                    'parent': node_id(path, depth - 1) if depth > 1 else main_class,
                },
            )
    yield from nodes.values()


def iter_classes(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The classes of the hierarchy, see :func:`class_nodes`.

    An unknown main class, a classification that cannot be read or an IPI
    accession that occurs twice is an error.
    """
    seen: set[str] = set()
    pairs = []
    for row in iter_excel(opener, sheet='Data', required=COLUMNS):
        ipi = row['IPI Accession']
        if ipi in seen:
            raise ValueError(f'IPI accession {ipi!r} occurs twice')
        seen.add(ipi)
        main_class = row['Main Class']
        if main_class not in MAIN_CLASSES:
            raise ValueError(f'Unknown main class {main_class!r} of {ipi}')
        # Reading the classes of the row checks its classification.
        statements(main_class, row['Classification'])
        pairs.append((main_class, row['Classification']))
    yield from class_nodes(pairs)
