"""Read the InterPro entries, their hierarchy, their GO terms and the entries of human proteins.

Four files, all official:

- ``entry.list`` of an InterPro release: tab-separated ``ENTRY_AC``,
  ``ENTRY_TYPE`` (``Active_site``, ``Binding_site``, ``Conserved_site``,
  ``Domain``, ``Family``, ``Homologous_superfamily``, ``PTM``, ``Repeat``) and
  ``ENTRY_NAME``, one line per entry.
- ``ParentChildTreeFile.txt``: the hierarchy of the entries. A line is an
  entry (``IPR000008::C2 domain::``) and the number of ``--`` before it is its
  depth; the entries below it at the next depth are its children. An entry
  with several parents appears under each of them.
- ``interpro2go``: the GO terms of the entries, as lines such as
  ``InterPro:IPR000003 Retinoid X receptor/HNF4 > GO:DNA binding ; GO:0003677``
  after comment lines that start with ``!``.
- UniProt's cross-references to InterPro for the reviewed human proteins
  (``accession`` and ``xref_interpro``): the InterPro entries that match each
  protein. The positions of the matches are not in this file.
"""

from __future__ import annotations

import re
from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.base import Download
from pypath.inputs_v2.parsers.base import iter_tsv, read_text

_ENTRY = re.compile(r'^((?:--)*)(IPR\d{6})::(.*)::$')
_GO = re.compile(r'^InterPro:(IPR\d{6}) .* > GO:.* ; (GO:\d{7})$')
_ACCESSION = re.compile(r'^IPR\d{6}$')
ENTRY_COLUMNS = ('ENTRY_AC', 'ENTRY_TYPE', 'ENTRY_NAME')


def _hierarchy(tree_download: Download, force_refresh: bool = False, **kwargs: Any) -> dict[str, tuple[str, list[str]]]:
    """Name and parents of each entry in the tree file, from its indented lines."""
    names: dict[str, str] = {}
    parents: dict[str, list[str]] = {}
    stack: list[str] = []
    for line in read_text(tree_download.open(force_refresh=force_refresh, **kwargs)).splitlines():
        match = _ENTRY.match(line.rstrip('\r'))
        if not match:
            raise ValueError(f'Malformed line in the InterPro tree file: {line!r}')
        depth, accession = len(match[1]) // 2, match[2]
        if depth > len(stack):
            raise ValueError(f'Tree depth jumps at {line!r}')
        del stack[depth:]
        if names.setdefault(accession, match[3]) != match[3]:
            raise ValueError(f'InterPro entry {accession} has two names in the tree file')
        parents.setdefault(accession, [])
        if stack and stack[-1] not in parents[accession]:
            parents[accession].append(stack[-1])
        stack.append(accession)
    return {accession: (names[accession], parents[accession]) for accession in names}


def iter_entries(
    opener, *, tree_download: Download, force_refresh: bool = False, **kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per InterPro entry of ``entry.list``.

    ``accession``, ``type`` and ``name`` are as ``entry.list`` writes them and
    ``parents`` the accessions of the entry's parents in the tree file (none
    for an entry outside the hierarchy or at its top). An entry of the tree
    file that ``entry.list`` lacks, one whose name differs between the files
    and a malformed accession are errors.
    """
    tree = _hierarchy(tree_download, force_refresh=force_refresh, **kwargs)
    listed = set()
    for row in iter_tsv(opener, required=ENTRY_COLUMNS):
        accession = row['ENTRY_AC']
        if not _ACCESSION.match(accession or '') or not row['ENTRY_TYPE']:
            raise ValueError(f'Malformed InterPro entry row {row!r}')
        listed.add(accession)
        name, parents = tree.get(accession, (row['ENTRY_NAME'], []))
        if name != row['ENTRY_NAME']:
            raise ValueError(
                f'InterPro entry {accession} is named {row["ENTRY_NAME"]!r} in '
                f'entry.list and {name!r} in the tree file'
            )
        yield {
            'accession': accession,
            'type': row['ENTRY_TYPE'],
            'name': row['ENTRY_NAME'],
            'parents': parents,
        }
    unknown = sorted(set(tree) - listed)
    if unknown:
        raise ValueError(f'Entries of the InterPro tree file not in entry.list: {unknown[:5]}')


def iter_go_terms(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per InterPro entry of ``interpro2go``.

    ``accession`` and ``go_terms``, the GO accessions (``GO:0003677``) of the
    entry in the order of the file. The file's GO names are not read.
    """
    terms: dict[str, list[str]] = {}
    for line in read_text(opener).splitlines():
        line = line.strip()
        if not line or line.startswith('!'):
            continue
        match = _GO.match(line)
        if not match:
            raise ValueError(f'Malformed line in interpro2go: {line!r}')
        if match[2] not in terms.setdefault(match[1], []):
            terms[match[1]].append(match[2])
    if not terms:
        raise ValueError('interpro2go has no mappings')
    for accession, go_terms in terms.items():
        yield {'accession': accession, 'go_terms': go_terms}


def iter_protein_entries(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per protein of UniProt's table of InterPro cross-references.

    ``accession`` is the UniProt accession and ``entries`` the InterPro
    accessions of the protein in the order UniProt lists them (empty for a
    protein without any).
    """
    for row in iter_tsv(opener, required=('Entry', 'InterPro')):
        entries = [item for item in (row['InterPro'] or '').split(';') if item]
        bad = [item for item in entries if not _ACCESSION.match(item)]
        if not row['Entry'] or bad:
            raise ValueError(f'Malformed UniProt InterPro row {row!r}')
        yield {'accession': row['Entry'], 'entries': entries}

