"""Read the in silico human surfaceome (SURFY) of Bausch-Fluck et al. 2018.

Table S3 of the paper is an ``.xlsx`` workbook with two sheets of 44 columns
(the first row is a title, the second the header). ``SurfaceomeMasterTable``
has every human UniProt protein and its features; ``in silico surfaceome
only`` has the 2,886 proteins that the random forest classifier SURFY
predicts to be at the cell surface, each with ``Surfaceome Label`` ``surface``.
The surfaceome is read from the second sheet (:func:`iter_proteins`). The
columns it reads, among those that describe the prediction:

    UniProt name, UniProt accession, UniProt gene, GeneID
    Surfaceome Label Source   how the protein got the label: the classifier
                              (``machine learning``), the positive training set
                              (``pos. trainingset``) or a UniProt GPI anchor
                              (``GPI (UniProt)``)
    MachineLearning trainingset, MachineLearning score
    MachineLearning FPR class (1=1%, 2=5%, 3=15%)
    TM domains, signalpeptide, topology, topology source
    Membranome Almen main-class, Membranome Almen sub-class
                              the class of the protein in Almen et al. 2009
                              (``parsers/almen2009.py``)

The first sheet gives the class of Almen et al. to 6,264 of its proteins, of
the surfaceome or not, by UniProt accession. Almen et al. identified their
proteins by IPI accessions, which no longer exist, so this sheet is where the
classification is read from for ``almen2009.py`` (:func:`iter_almen_classes`).

Excel has stored the EC class of a few enzymes as a number, and a number cannot
show a trailing zero: the 14 cells ``1.1`` of the first sheet are the 11
proteins of EC 1.1 and the 3 of EC 1.10 of Almen's own file, and nothing in
the sheet tells them apart. Such a protein is placed in the main class
``Enzymes`` only (see ``LOST_TRAILING_ZERO``).
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.almen2009 import (
    CLASS,
    ENZYMES,
    class_nodes,
    statements,
)
from pypath.inputs_v2.parsers.base import iter_excel

SHEET = 'in silico surfaceome only'
MASTER_SHEET = 'SurfaceomeMasterTable'
FPR_CLASS = 'MachineLearning FPR class (1=1%, 2=5%, 3=15%)'
COLUMNS = (
    'UniProt name',
    'UniProt accession',
    'UniProt gene',
    'GeneID',
    'Surfaceome Label',
    'Surfaceome Label Source',
    'MachineLearning trainingset',
    'MachineLearning score',
    FPR_CLASS,
    'TM domains',
    'signalpeptide',
    'topology',
    'topology source',
    'Membranome Almen main-class',
    'Membranome Almen sub-class',
)

# The columns read from the first sheet for the classes of Almen et al.
ALMEN_COLUMNS = (
    'UniProt name',
    'UniProt accession',
    'UniProt gene',
    'GeneID',
    'Membranome Almen main-class',
    'Membranome Almen sub-class',
)
# EC classes that a number in the sheet may stand for: the class of Almen's
# file whose last zero the number has lost, and the class it then looks like.
LOST_TRAILING_ZERO = {'1.1': ('1.1', '1.10')}

# The kind of the record that puts a protein in the list; the records of its
# classes have the kinds of ``parsers/almen2009.py`` (``class`` and ``ec``).
SURFACEOME = 'surfaceome'

LABEL = 'surface'
LABEL_SOURCES = ('machine learning', 'pos. trainingset', 'GPI (UniProt)')


def _sheet(opener) -> Generator[dict[str, Any], None, None]:
    # The first row of the sheet is a title.
    return iter_excel(opener, sheet=SHEET, header_row=2, required=COLUMNS)


def _text(value: Any) -> str | None:
    """A cell that holds a class as text or, where Excel has read an EC number
    like ``3.4`` as a number, as the number."""
    if value is None or isinstance(value, str):
        return value
    if isinstance(value, float) and value == int(value):
        return str(int(value))
    return str(value)


def _entrez(value: Any) -> list[str]:
    """Entrez Gene IDs of a protein: a number, or several joined by ``;``."""
    if value is None:
        return []
    ids = str(value).split(';')
    if not all(identifier.isdigit() for identifier in ids):
        raise ValueError(f'Not Entrez Gene IDs: {value!r}')
    return ids


def _classification(row: dict[str, Any]) -> tuple[str, str | None] | None:
    """The Almen class of a row, None where the source gives none.

    The sub-class of an enzyme is None where the cell is a number that stands
    for two EC classes (``LOST_TRAILING_ZERO``).
    """
    main_class = row['Membranome Almen main-class']
    sub_class = row['Membranome Almen sub-class']
    if main_class is None:
        if sub_class is not None:
            raise ValueError(f'Sub-class without main class: {sub_class!r}')
        return None
    # An EC number that Excel has read as a number: ``1.11`` and ``3.4`` are
    # EC classes, ``1`` the class of oxidoreductases. A number cannot show a
    # trailing zero: ``1.1`` may be EC 1.1 or EC 1.10.
    if main_class == ENZYMES:
        was_number = not isinstance(sub_class, str)
        sub_class = _text(sub_class)
        if was_number and sub_class in LOST_TRAILING_ZERO:
            return main_class, None
    if not isinstance(sub_class, str):
        raise ValueError(f'Class of {main_class!r} is not text: {sub_class!r}')
    return main_class, sub_class


def _class_statements(row: dict[str, Any]) -> list[tuple[str, str]]:
    """The classes of Almen et al. a row places its protein in, as ``(kind,
    term)`` pairs (see ``parsers/almen2009.py``); none where the row has no class."""
    classification = _classification(row)
    if classification is None:
        return []
    main_class, sub_class = classification
    if sub_class is None:
        return [(CLASS, main_class)]
    return statements(main_class, sub_class)


def iter_proteins(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """One record per statement about a protein of the surfaceome.

    Every protein gives a record of kind ``surfaceome`` (membership in the list,
    with the fields of the prediction) and, where the sheet has one, a record of
    kind ``class`` or ``ec`` for each class of Almen et al. (see
    ``parsers/almen2009.py``). Fields: ``uniprot``, ``entry_name``, ``gene``,
    ``entrez`` (a list; empty where the sheet has none), ``kind`` and ``term`` (for the
    classes), and for every record of the protein ``label_source``,
    ``trainingset`` (None unless ``surface``), ``score`` (the classifier's
    score), ``fpr_class`` (its false positive rate class: 1, 2, 3, or the same
    with ``nonTM``, None for a few proteins), ``tm_domains``, ``signal_peptide``,
    ``topology`` and ``topology_source``.
    """
    seen: set[str] = set()
    for row in _sheet(opener):
        accession = row['UniProt accession']
        if accession in seen:
            raise ValueError(f'UniProt accession {accession!r} occurs twice')
        seen.add(accession)
        if row['Surfaceome Label'] != LABEL:
            raise ValueError(
                f'{accession}: label {row["Surfaceome Label"]!r} in the surfaceome sheet'
            )
        if row['Surfaceome Label Source'] not in LABEL_SOURCES:
            raise ValueError(
                f'{accession}: unknown label source {row["Surfaceome Label Source"]!r}'
            )
        if row['signalpeptide'] not in (0, 1):
            raise ValueError(f'{accession}: signalpeptide {row["signalpeptide"]!r}')
        record = {
            'uniprot': accession,
            'entry_name': row['UniProt name'],
            'gene': row['UniProt gene'],
            'entrez': _entrez(row['GeneID']),
            'label_source': row['Surfaceome Label Source'],
            'trainingset': row['MachineLearning trainingset'],
            'score': float(row['MachineLearning score']),
            'fpr_class': None if row[FPR_CLASS] is None else str(row[FPR_CLASS]),
            'tm_domains': int(row['TM domains']),
            'signal_peptide': bool(row['signalpeptide']),
            'topology': row['topology'],
            'topology_source': row['topology source'],
        }
        yield {**record, 'kind': SURFACEOME, 'term': None}
        for kind, term in _class_statements(row):
            yield {**record, 'kind': kind, 'term': term}


def iter_almen_classes(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """One record per class of Almen et al. of a protein of the first sheet.

    Every protein the sheet gives a class to, whatever its surfaceome label,
    gives a record of kind ``class`` or ``ec`` for each of its classes (see
    ``parsers/almen2009.py``). Fields: ``uniprot``, ``entry_name``, ``gene``
    (None where the sheet has none), ``entrez`` (a list), ``kind`` and
    ``term``. A sheet without a classified protein is an error.
    """
    seen: set[str] = set()
    for row in iter_excel(opener, sheet=MASTER_SHEET, header_row=2, required=ALMEN_COLUMNS):
        classes = _class_statements(row)
        if not classes:
            continue
        accession = row['UniProt accession']
        if not accession or accession in seen:
            raise ValueError(f'UniProt accession {accession!r} is missing or occurs twice')
        seen.add(accession)
        record = {
            'uniprot': accession,
            'entry_name': row['UniProt name'],
            'gene': row['UniProt gene'],
            'entrez': _entrez(row['GeneID']),
        }
        for kind, term in classes:
            yield {**record, 'kind': kind, 'term': term}
    if not seen:
        raise ValueError('The surfaceome table gives no protein a class of Almen et al.')


def iter_terms(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The list of the surfaceome, and the classes of Almen et al. that the
    sheet's proteins are in (see ``class_nodes``). Classes have ``kind``
    ``class``; the list has ``kind`` ``surfaceome``."""
    yield {'kind': SURFACEOME}
    # An enzyme's sub-class takes no part in the hierarchy of the classes.
    pairs = [
        (pair[0], pair[1] or '')
        for row in _sheet(opener)
        if (pair := _classification(row)) is not None
    ]
    for node in class_nodes(pairs):
        yield {**node, 'kind': CLASS}
