"""
Read the CollecTRI2 files and join them into one row per TF-target pair.

CollecTRI2 (Zenodo record 19773408) is a text-mined and curated collection of
transcriptional regulatory interactions between human transcription factors
and target genes, in two tab separated files, both exports of the ExTRI2
workflow (the first line of each is a comment that names the job):

``CollecTRI2.tsv.gz``, the aggregated interactions, with one line per pair and
56 columns: the pair (``TF:TG``, then the gene symbol of the transcription
factor and of the target) and, for each of the twelve sources, a ``present``
column that holds the source's name when it supports the pair, and the
source's own columns. Those are in every source the PubMed IDs of its
evidence (``PMID``, items joined by ``|``, and several IDs of one item by
``;``) and, for most sources, the sign of each item in the same order:

    ExTRI2 (text mining)   Sign                    UP, DOWN, UNKNOWN
    HTRI                   Technique, Confidence
    TRRUST                 Regulation              Activation, Repression, Unknown
    TFactS                 Sign, Species, Source   UP, DOWN
    GOA                    Sign                    UP, DOWN, Unknown
    IntAct                 Method ID
    SIGNOR                 Effect, Sign            UP, DOWN, Unknown
    CytReg                 Assay type, species, Activation/Repression, Year
    GEREDB                 Effect                  positive, negative, unknown
    NTNU Curated           Sign                    +, -
    Pavlidis2021           MoR                     activation, repression, not_applicable
    DoRothEA_A             Directed, Effect        Stimulate, Inhibit, Unknown

and the transcription factor's classification (Lambert, Lovering, GO terms,
TFClass).

``CollecTRI2_regulome.tsv.gz``, the signed regulome made from it for the
activity inference of transcription factors, one line per pair: ``ID``,
``source`` and ``target`` (the gene symbols), ``weight`` (1 or -1), ``TF Type``
(dbTF, coTF, coTF candidate or empty) and ``Auto-regulation``.

The weight is the sign that the ``regulome`` task of the SaezLab workflow
(https://github.com/Rbbt-Workflows/SaezLab) assigns to the pair. With the
parameters of this file, which reproduce all of its weights from the
aggregated file, the rule is:

1. Each item of evidence of each source is a PubMed ID (or another
   identifier) with a sign. The sign columns used are those listed above,
   except the ``MoR`` of Pavlidis2021 and the ``Directed`` of DoRothEA_A,
   which the workflow ignores. Sources without a sign column and items with an
   unknown sign give an evidence without sign. An item without an ID is an
   evidence of its own. Evidence IDs are counted once for each sign.
2. The pair is negative if it has evidence IDs of a negative sign and at
   least as many of them as of a positive sign.
3. A pair with no evidence of either sign is negative if most of the pairs of
   its transcription factor are negative (a factor's pair counts as negative
   if it has more negative than positive evidence IDs, and as positive
   otherwise), and positive if not.
4. Every other pair is positive.

So the weight of a pair is a statement of its evidence only where the signed
evidence has a strict majority. With a tie, or without a signed evidence, it is
a tie-break, a default or an inference from the other pairs of the factor.
:func:`iter_pairs` computes the counts of the signed evidence and says which
case a pair is.
"""

from __future__ import annotations

from collections.abc import Generator
import io
from types import SimpleNamespace
from typing import Any

from pypath.inputs_v2._publications import is_pmid
from pypath.inputs_v2.parsers.base import iter_tsv, read_text

MAIN_KEY = 'Transcription Factor (Associated Gene Name)'
MAIN_TARGET = 'Target Gene (Associated Gene Name)'
REGULOME_COLUMNS = ('ID', 'source', 'target', 'weight', 'TF Type', 'Auto-regulation')
# source -> (column of its sign, or None; column of its evidence IDs). The
# source is the name in the header of its columns (``[NTNU Curated] Sign``);
# its ``present`` column holds the label that the file writes in the pair.
SOURCES = {
    'ExTRI2': ('[ExTRI2] Sign', '[ExTRI2] PMID'),
    'HTRI': (None, '[HTRI] PMID'),
    'TRRUST': ('[TRRUST] Regulation', '[TRRUST] PMID'),
    'TFactS': ('[TFactS] Sign', '[TFactS] PMID'),
    'GOA': ('[GOA] Sign', '[GOA] PMID'),
    'IntAct': (None, '[IntAct] PMID'),
    # The columns ``Effect`` and ``Sign`` of SIGNOR both state the sign; the
    # workflow takes the last.
    'SIGNOR': ('[SIGNOR] Sign', '[SIGNOR] PMID'),
    'CytReg': ('[CytReg] Activation/Repression', '[CytReg] PMID'),
    'GEREDB': ('[GEREDB] Effect', '[GEREDB] PMID'),
    'NTNU Curated': ('[NTNU Curated] Sign', '[NTNU Curated] PMID'),
    'Pavlidis2021': (None, '[Pavlidis2021] PMID'),
    'DoRothEA_A': ('[DoRothEA_A] Effect', '[DoRothEA_A] PMID'),
}
# The values of the sign columns, lowercase -> the sign of the evidence; None
# is "no sign". A value that is not here is an error.
SIGNS = {
    'up': '+',
    'activation': '+',
    'positive': '+',
    '+': '+',
    'stimulate': '+',
    'down': '-',
    'repression': '-',
    'negative': '-',
    '-': '-',
    'inhibit': '-',
    'unknown': None,
    '?': None,
    '': None,
    # Both signs in one item (NTNU Curated); the workflow counts it as a
    # category of its own, which has no say in the sign of the pair.
    '+_-': None,
    '-_+': None,
}
# How the weight of a pair relates to its evidence (see the module docstring).
PMID_MAJORITY = 'pmid_majority'
PMID_TIE = 'pmid_tie'
NO_SIGNED_EVIDENCE = 'no_signed_evidence'
CONTRADICTED = 'contradicted'


def _read(opener, **kwargs: Any):
    """The lines of one file: the first line is a comment, the second the
    header, whose first name has a ``#``."""
    text = read_text(opener)
    if not text.startswith('#:'):
        raise ValueError('CollecTRI2 file does not start with the comment line of its export')
    return iter_tsv(
        SimpleNamespace(result=io.StringIO(text)),
        skip_rows=1,
        strict=True,
        quotechar=None,
        **kwargs,
    )


def _items(cell: object) -> list[str]:
    return str(cell or '').split('|') if cell else []


def _evidence(row: dict[str, Any], source: str) -> list[tuple[str, tuple[str, ...]]]:
    """The items of evidence of a source: the sign and the evidence IDs."""
    sign_column, id_column = SOURCES[source]
    ids = _items(row[id_column])
    signs = _items(row[sign_column]) if sign_column else []
    if not ids and not signs:
        ids, signs = [''], ['']
    elif not signs:
        signs = [''] * len(ids)
    elif not ids:
        ids = [''] * len(signs)
    if len(ids) != len(signs):
        raise ValueError(
            f'{row["TF:TG"]}: {source} has {len(ids)} evidence items and {len(signs)} signs'
        )
    evidence = []
    for number, (item, sign) in enumerate(zip(ids, signs)):
        if sign.lower() not in SIGNS:
            raise ValueError(f'Unknown CollecTRI2 sign of {source}: {sign!r}')
        # An item without an ID is an evidence of its own.
        found = tuple(
            token.strip() for token in item.split(';') if token.strip()
        ) or (f'{source}#{number}',)
        evidence.append((SIGNS[sign.lower()], found))
    return evidence


def _pair(
    row: dict[str, Any], regulome: dict[tuple[str, str], dict[str, Any]]
) -> tuple[dict[str, Any], dict[str, set[str]]]:
    """The row of a pair and the sets of evidence IDs by sign."""
    tf, target = row[MAIN_KEY], row[MAIN_TARGET]
    if (tf, target) not in regulome:
        raise ValueError(f'CollecTRI2 regulome lacks the pair {tf}:{target}')
    signed: dict[str, set[str]] = {'+': set(), '-': set()}
    sources, references = [], []
    for source in SOURCES:
        label = row[f'[{source}] present']
        if not label:
            continue
        sources.append(label)
        for sign, found in _evidence(row, source):
            if sign:
                signed[sign].update(found)
            references.extend(found)
    return {
        'tf': tf,
        'target': target,
        'sources': sources,
        # PubMed IDs, in the order of the file; other identifiers of the
        # evidence (PMC IDs, PAZAR project names) are not PubMed IDs.
        'pmids': list(dict.fromkeys(r for r in references if is_pmid(r))),
        'weight': regulome[(tf, target)]['weight'],
        'tf_type': regulome[(tf, target)]['TF Type'] or None,
        'auto_regulation': bool(regulome[(tf, target)]['Auto-regulation']),
        'up': len(signed['+']),
        'down': len(signed['-']),
    }, signed


def iter_pairs(
    opener,
    force_refresh: bool = False,
    regulome=None,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Yield the pairs of ``CollecTRI2.tsv.gz`` with their weight in the regulome.

    ``opener`` is the aggregated file and ``regulome`` the download of the
    regulome (anything with an ``open`` method that takes ``force_refresh``).
    A pair of the aggregated file that the regulome lacks, a pair twice in
    either file, a header without the expected columns, an unknown sign, a
    different number of evidence items and signs and a file without pairs
    raise ``ValueError``. Pairs of the regulome that the aggregated file lacks give
    no row.

    Keys: ``tf`` and ``target`` (gene symbols); ``sources`` (the labels of the
    sources that support the pair, as the file writes them); ``pmids`` (the
    PubMed IDs of their evidence); ``weight`` (``1`` or ``-1``, the
    regulome's consensus sign); ``tf_type`` (dbTF, coTF, coTF candidate or
    None); ``auto_regulation`` (bool); ``up`` and ``down`` (the numbers of
    distinct evidence IDs with a positive or a negative sign); and
    ``sign_decision``: ``pmid_majority`` if one sign has more evidence than
    the other and the weight is that sign, ``pmid_tie`` if as many of both,
    ``no_signed_evidence`` if none has a sign (the weight is then a default or
    an inference from the factor's other pairs), ``contradicted`` if the
    weight is the other sign than the majority of the evidence.
    """
    if regulome is None:
        raise ValueError('The CollecTRI2 regulome download is required')
    weights: dict[tuple[str, str], dict[str, Any]] = {}
    for line in _read(
        regulome.open(force_refresh=force_refresh),
        required=('#ID', 'source', 'target', 'weight', 'TF Type', 'Auto-regulation'),
    ):
        key = (line['source'], line['target'])
        if key in weights:
            raise ValueError(f'CollecTRI2 regulome lists {key} twice')
        weights[key] = line
    if not weights:
        raise ValueError('CollecTRI2 regulome has no pairs')
    required = ['#TF:TG', MAIN_KEY, MAIN_TARGET]
    for source, (sign_column, id_column) in SOURCES.items():
        required += [f'[{source}] present', id_column] + ([sign_column] if sign_column else [])
    seen = set()
    for line in _read(opener, required=tuple(required)):
        line['TF:TG'] = line.pop('#TF:TG')
        pair, signed = _pair(line, weights)
        if (pair['tf'], pair['target']) in seen:
            raise ValueError(f'CollecTRI2 lists {(pair["tf"], pair["target"])} twice')
        seen.add((pair['tf'], pair['target']))
        up, down = pair['up'], pair['down']
        if up == down:
            decision = PMID_TIE if up else NO_SIGNED_EVIDENCE
        else:
            majority = '1' if up > down else '-1'
            decision = PMID_MAJORITY if pair['weight'] == majority else CONTRADICTED
        yield {**pair, 'sign_decision': decision}
    if not seen:
        raise ValueError('CollecTRI2 has no pairs')
