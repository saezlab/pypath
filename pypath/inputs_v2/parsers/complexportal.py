"""Read the ComplexTAB files of the Complex Portal.

A ComplexTAB file is a tab-separated table with one row per complex of one
species (``9606.tsv``: manually curated human complexes; ``9606_predicted.tsv``:
machine-learning predicted human complexes). Quote characters occur inside the
text columns and are data. A missing value is ``-``.

Each row is yielded as a dictionary with these fields:

- ``accession``: the complex accession (``CPX-1``)
- ``name``: the recommended name; None where the file has ``-``
- ``synonyms``: the aliases, in file order
- ``taxon``: the NCBI taxonomy ID of the complex
- ``members``: the molecules of the complex, in file order. A molecule is
  ``{'accession': ..., 'set': None, 'stoichiometry': n}``; a defined set of
  molecules (written ``[P0DP23,P0DP24,P0DP25](1)``, for paralogs that cannot be
  told apart) is ``{'accession': None, 'set': [...], 'stoichiometry': n}``.
  The accession is a UniProt accession (with an isoform or ``-PRO_`` chain
  suffix when the file has one), a ChEBI ID, an RNAcentral ID or the accession
  of another complex. The stoichiometry is the number as the file writes it;
  0 means that it is not known.
- ``evidence_code`` and ``evidence_label``: the ECO code of the complex and the
  name written after it
- ``experimental_evidence``: references to the experiments (``intact:EBI-...``,
  ``wwpdb:...``, ``emdb:...``) as ``(database, accession)`` pairs
- ``go_terms``: the GO term accessions
- ``xrefs``: the cross references as ``{'database', 'accession', 'qualifier'}``
- ``description``: the free-text description; None where the file has ``-``
- ``source``: the database that curated the complex, as ``{'name', 'mi'}``
  (``psi-mi:"MI:0469"(IntAct)``); ``mi`` is None where the file writes ``-``

The columns that are free text (properties, comment, disease, ligand, agonist,
antagonist) and the assembly type are not read. The "Expanded participant
list" is derived from the members and is not read either.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

COLUMNS = (
    '#Complex ac',
    'Recommended name',
    'Aliases for complex',
    'Taxonomy identifier',
    'Identifiers (and stoichiometry) of molecules in complex',
    'Evidence Code',
    'Experimental evidence',
    'Go Annotations',
    'Cross references',
    'Description',
    'Source',
)
_MISSING = '-'
_MEMBER = re.compile(r'(?P<accession>.+)\((?P<stoichiometry>\d+)\)')
_SET = re.compile(r'\[(?P<members>[^\[\]]+)\]')
_EVIDENCE = re.compile(r'(?P<code>ECO:\d{7})\((?P<label>.*)\)')
_GO_TERM = re.compile(r'(?P<accession>GO:\d{7})\(.*\)')
_XREF = re.compile(r'(?P<database>[^:]+):(?P<accession>.+)\((?P<qualifier>[^()]*)\)')
_SOURCE = re.compile(r'psi-mi:"(?P<mi>[^"]*)"\((?P<name>[^()]*)\)')


def _value(text: str) -> str | None:
    text = text.strip()
    return None if text in {'', _MISSING} else text


def _items(text: str) -> list[str]:
    value = _value(text)
    return value.split('|') if value else []


def _member(text: str, accession: str) -> dict[str, Any]:
    match = _MEMBER.fullmatch(text)
    if not match:
        raise ValueError(f'Complex Portal {accession}: malformed molecule {text!r}')
    molecule = match['accession']
    members = _SET.fullmatch(molecule)
    return {
        'accession': None if members else molecule,
        'set': [item.strip() for item in members['members'].split(',')] if members else None,
        'stoichiometry': int(match['stoichiometry']),
    }


def _pairs(text: str, accession: str) -> list[tuple[str, str]]:
    pairs = []
    for item in _items(text):
        database, colon, reference = item.partition(':')
        if not colon:
            raise ValueError(f'Complex Portal {accession}: malformed reference {item!r}')
        pairs.append((database, reference))
    return pairs


def _parsed(row: dict[str, Any]) -> dict[str, Any]:
    accession = row['#Complex ac'].strip()
    evidence = _EVIDENCE.fullmatch(row['Evidence Code'].strip())
    if not evidence:
        raise ValueError(f'Complex Portal {accession}: malformed evidence code {row["Evidence Code"]!r}')
    source = _SOURCE.fullmatch(row['Source'].strip())
    if not source:
        raise ValueError(f'Complex Portal {accession}: malformed source {row["Source"]!r}')
    go_terms = []
    for item in _items(row['Go Annotations']):
        match = _GO_TERM.fullmatch(item)
        if not match:
            raise ValueError(f'Complex Portal {accession}: malformed GO annotation {item!r}')
        go_terms.append(match['accession'])
    xrefs = []
    for item in _items(row['Cross references']):
        match = _XREF.fullmatch(item)
        if not match:
            raise ValueError(f'Complex Portal {accession}: malformed cross reference {item!r}')
        xrefs.append(match.groupdict())
    return {
        'accession': accession,
        'name': _value(row['Recommended name']),
        'synonyms': _items(row['Aliases for complex']),
        'taxon': row['Taxonomy identifier'].strip(),
        'members': [
            _member(item, accession)
            for item in _items(row['Identifiers (and stoichiometry) of molecules in complex'])
        ],
        'evidence_code': evidence['code'],
        'evidence_label': evidence['label'],
        'experimental_evidence': _pairs(row['Experimental evidence'], accession),
        'go_terms': go_terms,
        'xrefs': xrefs,
        'description': _value(row['Description']),
        'source': {'name': source['name'], 'mi': _value(source['mi'])},
    }


def iter_complextab(opener, **kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The complexes of a ComplexTAB file; a changed header or an empty file is an error."""
    count = 0
    for row in iter_tsv(opener, quotechar=None, strict=True, **kwargs):
        missing = [column for column in COLUMNS if column not in row]
        if missing:
            raise ValueError(f'Complex Portal ComplexTAB lacks columns {missing}')
        count += 1
        yield _parsed(row)
    if not count:
        raise ValueError('Complex Portal ComplexTAB has no complexes')
