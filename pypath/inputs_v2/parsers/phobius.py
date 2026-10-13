"""Read OmniPath's Phobius predictions for the human proteome.

OmniPath ran the Phobius predictor (Käll, Krogh and Sonnhammer 2004) on the
human Swiss-Prot proteome and published the output as a tab-separated file
(``phobius_human_proteome.tsv``, last modified 2020-03-30), with one row per
UniProt accession:

    Source  Uniprot_AC  Uniprot_ID  TM  SP  PREDICTION
    sp      P01920      DQB1_HUMAN  1   Y   n15-25c32/33o231-251i
    sp      P18440      ARY1_HUMAN  0   0   o

``TM`` is the number of predicted transmembrane helices, ``SP`` is ``Y`` for a
predicted signal peptide and ``0`` otherwise, and ``PREDICTION`` is the Phobius
topology: an optional signal peptide (``n15-25c32/33``: the n- and h-region
span, the c-region end and the cleavage site), then ``o`` or ``i`` for the side
of the membrane the sequence starts on (non-cytoplasmic or cytoplasmic) and the
helices (``231-251``) each followed by the side the sequence is on after it.
A topology of ``o`` or ``i`` alone is a prediction of no helix.

Source defect: in 317 rows a tab splits the entry name in two, and the column
is padded with spaces, as in ``CENL1_H<TAB>UMAN       0`` (295 rows, split after
``_H``) and ``TVA2_HU<TAB>MAN        0`` (22 rows, split after ``_HU``). The rows
still have six fields, so the shift is visible only in the third and fourth.
They are repaired by the rule in :func:`_repair`.

:func:`iter_proteins` yields one dictionary per row of the file, checked.
:func:`iter_predictions` yields one dictionary per statement: for each protein
the transmembrane prediction if it has at least one helix and the signal
peptide prediction if it has one. :func:`iter_terms` yields the two terms.
"""

from __future__ import annotations

from collections.abc import Generator
import logging
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

__all__ = [
    'SIGNAL_PEPTIDE',
    'TRANSMEMBRANE',
    'iter_predictions',
    'iter_proteins',
    'iter_terms',
]

logger = logging.getLogger(__name__)

TRANSMEMBRANE = 'transmembrane'
SIGNAL_PEPTIDE = 'signal_peptide'

COLUMNS = ('Source', 'Uniprot_AC', 'Uniprot_ID', 'TM', 'SP', 'PREDICTION')

_ACCESSION = re.compile(r'[A-Z0-9]{6,10}')
_ENTRY_NAME = re.compile(r'[A-Z0-9]+_HUMAN')
# The entry name as the tab leaves it in the third column, and the rest of it
# with the padding and the TM value in the fourth.
_SPLIT_NAME = re.compile(r'[A-Z0-9]+_H(?:U(?:M(?:A)?)?)?')
_SPLIT_TAIL = re.compile(r'(?P<tail>[A-Z]+)\s+(?P<tm>\d+)')
# A signal peptide, the side the sequence starts on, and helices each with the
# side it ends on.
_TOPOLOGY = re.compile(
    r'(?P<signal>n\d+-\d+c\d+/\d+)?(?P<start>[io])(?P<helices>(?:\d+-\d+[io])*)'
)
_HELIX = re.compile(r'\d+-\d+[io]')


def _repair(fields: list[str]) -> tuple[list[str], bool]:
    """The fields of a row with the entry name made whole.

    A row is broken when its fourth field is not a number. The rule: the third
    field is an entry name cut inside ``HUMAN`` (``CENL1_H``, ``TVA2_HU``), the
    fourth is the rest of the name, padding spaces and the TM value
    (``UMAN       0``); the entry name is the two parts joined and the TM value
    the number. A row that does not fit this rule raises.
    """
    if fields[3].isdigit():
        return fields, False
    tail = _SPLIT_TAIL.fullmatch(fields[3])
    if (
        tail is None
        or not _SPLIT_NAME.fullmatch(fields[2])
        or not _ENTRY_NAME.fullmatch(fields[2] + tail['tail'])
    ):
        raise ValueError(f'Phobius row is malformed and does not fit the repair rule: {fields!r}')
    return [*fields[:2], fields[2] + tail['tail'], tail['tm'], *fields[4:]], True


def _check(fields: list[str]) -> dict[str, Any]:
    source, accession, entry_name, tm, sp, prediction = fields
    if source != 'sp':
        raise ValueError(f'Phobius row is not from Swiss-Prot: {fields!r}')
    if not _ACCESSION.fullmatch(accession):
        raise ValueError(f'Phobius row has no UniProt accession: {fields!r}')
    # The file is a run on the human proteome, so every entry is a human one.
    if not _ENTRY_NAME.fullmatch(entry_name):
        raise ValueError(f'Phobius entry name does not end in _HUMAN: {fields!r}')
    if sp not in {'Y', '0'}:
        raise ValueError(f'Phobius row has an unknown SP value: {fields!r}')
    topology = _TOPOLOGY.fullmatch(prediction)
    if topology is None:
        raise ValueError(f'Phobius row has an unknown topology: {fields!r}')
    helices = len(_HELIX.findall(topology['helices']))
    if int(tm) != helices:
        raise ValueError(f'Phobius TM is not the number of helices of the topology: {fields!r}')
    signal = topology['signal'] is not None
    if (sp == 'Y') != signal:
        raise ValueError(f'Phobius SP does not agree with the topology: {fields!r}')
    return {
        'accession': accession,
        'entry_name': entry_name,
        'helices': helices,
        'signal_peptide': signal,
        'topology': prediction,
    }


def iter_proteins(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per row of the file: ``accession``, ``entry_name``,
    ``helices`` (the number of predicted transmembrane helices), ``signal_peptide``
    (a boolean), ``topology`` (the string as written) and ``repaired`` (True for
    a row whose entry name was split by a tab).

    A row that cannot be read, an entry name that does not end in ``_HUMAN``, a
    TM value that is not the number of helices in the topology, an SP value that
    does not agree with the topology (``Y`` when it starts with a signal
    peptide), an accession that is in the file twice and a file without rows are
    errors. The number of repaired rows is logged.
    """
    rows = iter_tsv(opener, quotechar=None, strict=True, required=COLUMNS)
    seen: set[str] = set()
    repaired = 0
    for row in rows:
        fields, fixed = _repair([row[column] for column in COLUMNS])
        protein = _check(fields)
        if protein['accession'] in seen:
            raise ValueError(f'Phobius accession {protein["accession"]} is in the file twice')
        seen.add(protein['accession'])
        repaired += fixed
        yield {**protein, 'repaired': fixed}
    logger.info('Phobius: %d rows, %d with the entry name repaired', len(seen), repaired)


def iter_predictions(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per statement: the row of :func:`iter_proteins` with
    ``term`` as ``transmembrane`` (for a protein with at least one helix) or
    ``signal_peptide``, a protein with both giving two. A protein with neither
    gives none: the file has no prediction to state for it.
    """
    for protein in iter_proteins(opener):
        if protein['helices']:
            yield {**protein, 'term': TRANSMEMBRANE}
        if protein['signal_peptide']:
            yield {**protein, 'term': SIGNAL_PEPTIDE}


def iter_terms(opener=None, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield the two terms of the module; the file does not list them."""
    for term in (TRANSMEMBRANE, SIGNAL_PEPTIDE):
        yield {'term': term}
