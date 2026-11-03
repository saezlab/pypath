"""
Parsers of the MSigDB gene set files.

MSigDB distributes each collection of gene sets as a JSON file
(``h.all.v2026.1.Hs.json``) next to the GMT files, which hold only the gene
symbols. The JSON file is one object: the key is the gene set name
(``HALLMARK_ADIPOGENESIS``) and the value has the fields

    ``collection``       the code of the collection, as the module names it
    ``systematicName``   MSigDB's stable accession of the set (``M5890``)
    ``geneSymbols``      the HGNC symbols of the genes (the GMT columns, in the
                         same order)
    ``msigdbURL``        the page of the set at the MSigDB site (the GMT
                         description column)
    ``pmid``             the PubMed identifier of the publication the set comes
                         from, or an empty string
    ``exactSource``      what the set was made from, or an empty string: a
                         GO, HPO, Reactome, WikiPathways or KEGG MEDICUS
                         accession, a UniProt accession (``C3:TFT:GTRD``), or
                         text (a table or figure of the publication, a data
                         set or a description of the comparison)
    ``externalDetailsURL``, ``filteredBySimilarity``,
    ``externalNamesForSimilarTerms``   not read

The JSON files hold the sets and symbols of the GMT files (the same names and
the same gene lists in the same order, checked on release 2026.1.Hs), so they
are read instead of them. A JSON object has no order; the sets of a file are
read in the order of their names, which is the order of the GMT files. A gene
set is in one collection (and sub-collection) file; the "all" files of a
collection repeat the sets of its sub-collections and are not read.

The module passes the collections as a mapping of the collection code (MSigDB's
own, ``C2:CP:REACTOME``) to its parts:

    ``title``, ``description``   its name and what it is, for the vocabulary
    ``parent``                   the code of the collection it is a part of, or None
    ``download``                 the ``Download`` of its JSON file, or None for a
                                 collection that only groups sub-collections

Gene set names and systematic names are unique across the files read; one that
occurs twice raises, as does a set that is missing a field, is in a file of
another collection than its ``collection`` or has a malformed systematic name or
PubMed identifier.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict
from collections.abc import Generator, Mapping
from typing import Any

from pypath.inputs_v2._publications import is_pmid
from pypath.inputs_v2.parsers.base import read_text

# The fields every set of a file has; a file without one is not the file the
# module reads.
FIELDS = ('collection', 'systematicName', 'geneSymbols', 'msigdbURL', 'pmid', 'exactSource')
# ``M5890``; the first set of release 2026.1.Hs is ``M0``.
SYSTEMATIC_NAME = re.compile(r'M\d+')


def _check(collections: Mapping[str, Mapping[str, Any]]) -> None:
    for code, collection in collections.items():
        parent = collection.get('parent')
        if parent is not None and parent not in collections:
            raise ValueError(f'MSigDB collection {code} has the unknown parent {parent}')


def _file(code: str, download: Any, force_refresh: bool) -> dict[str, dict[str, Any]]:
    """The object of the JSON file of a collection: gene set name -> its fields."""
    try:
        sets = json.loads(read_text(download.open(force_refresh=force_refresh)))
    except json.JSONDecodeError as error:
        raise ValueError(f'MSigDB collection {code}: expected a JSON file ({error})') from error
    if not isinstance(sets, dict):
        raise ValueError(f'MSigDB collection {code}: expected one JSON object of gene sets')
    return sets


def _gene_sets(
    collections: Mapping[str, Mapping[str, Any]],
    force_refresh: bool,
) -> Generator[tuple[str, str, dict[str, Any]], None, None]:
    """The sets of each file: the code of the collection, the name and the fields.

    The sets of a file are in the order of their names.
    """
    seen: dict[str, str] = {}
    systematic: dict[str, str] = {}
    for code, collection in collections.items():
        download = collection.get('download')
        if download is None:
            continue
        sets = _file(code, download, force_refresh)
        for name in sorted(sets):
            entry = sets[name]
            if not isinstance(entry, dict) or any(field not in entry for field in FIELDS):
                raise ValueError(
                    f'MSigDB gene set {name} of the collection {code} lacks one of {FIELDS}'
                )
            if entry['collection'] != code:
                raise ValueError(
                    f'MSigDB gene set {name} is in the file of the collection {code} '
                    f'but names the collection {entry["collection"]!r}'
                )
            if name in seen:
                raise ValueError(
                    f'MSigDB gene set {name} is in the collections {seen[name]} and {code}'
                )
            seen[name] = code
            accession = entry['systematicName']
            if not isinstance(accession, str) or not SYSTEMATIC_NAME.fullmatch(accession):
                raise ValueError(
                    f'MSigDB gene set {name} has the malformed systematic name {accession!r}'
                )
            if accession in systematic:
                raise ValueError(
                    f'MSigDB systematic name {accession} is of the gene sets '
                    f'{systematic[accession]} and {name}'
                )
            systematic[accession] = name
            pmid = entry['pmid']
            if pmid != '' and not is_pmid(pmid):
                raise ValueError(f'MSigDB gene set {name} has the malformed PubMed id {pmid!r}')
            if not isinstance(entry['geneSymbols'], list):
                raise ValueError(f'MSigDB gene set {name} has no list of gene symbols')
            yield code, name, entry
        if not sets:
            raise ValueError(f'MSigDB collection {code} has no gene sets')


def iter_terms(
    opener=None,
    *,
    collections: Mapping[str, Mapping[str, Any]],
    force_refresh: bool = False,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """The collections, then the gene sets of every file.

    A collection row has ``kind`` ``collection``, ``id`` (the code), ``name``
    (the title), ``description`` and ``parent`` (code or None). A gene set row
    has ``kind`` ``gene_set``, ``id`` (the name), ``systematic_name`` (``M5890``),
    ``url`` (the page of the set), ``pmid`` (the PubMed identifier as written, or
    None), ``exact_source`` (the ``exactSource`` field as written, or None) and
    ``collection`` (the code of its file's collection).
    """
    _check(collections)
    for code, collection in collections.items():
        yield {
            'kind': 'collection',
            'id': code,
            'name': collection['title'],
            'description': collection.get('description'),
            'parent': collection.get('parent'),
        }
    for code, name, entry in _gene_sets(collections, force_refresh):
        yield {
            'kind': 'gene_set',
            'id': name,
            'systematic_name': entry['systematicName'],
            'url': entry['msigdbURL'] or None,
            'pmid': entry['pmid'] or None,
            'exact_source': entry['exactSource'] or None,
            'collection': code,
        }


def iter_genes(
    opener=None,
    *,
    collections: Mapping[str, Mapping[str, Any]],
    force_refresh: bool = False,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """One row for each gene, with the names of all the gene sets it is a member of.

    Fields: ``gene`` (the symbol as written) and ``gene_sets`` (the names, in
    the order of the files and, in a file, of the names). A gene listed twice in
    a set raises.
    """
    _check(collections)
    members: dict[str, list[str]] = defaultdict(list)
    for _, name, entry in _gene_sets(collections, force_refresh):
        genes = [gene for gene in entry['geneSymbols'] if gene.strip()]
        if len(set(genes)) != len(genes):
            raise ValueError(f'MSigDB gene set {name} lists a gene twice')
        for gene in genes:
            members[gene].append(name)
    for gene, gene_sets in members.items():
        yield {'gene': gene, 'gene_sets': gene_sets}
