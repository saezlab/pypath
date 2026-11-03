"""
Readers for the CancerSEA download page and its gene sets.

The page (``goDownload``) has a table of the 14 functional states, with the
description of each, the number of genes and the sources of the gene set, and a
link to the state's gene set, a tab separated text file with the Ensembl gene
identifier and the gene symbol of each member (``EnsembleID``, ``GeneName``). The
shared HTML table and TSV readers read them; the functions here find the table
of states, open the gene set of each and check it against the page.
"""

from __future__ import annotations

from collections.abc import Generator, Mapping
import re
from typing import Any

from pypath.inputs_v2.base import Download
from pypath.inputs_v2.parsers.base import iter_html_tables, iter_tsv, read_text

__all__ = ['iter_gene_sets', 'iter_states']

# The header of the table of states.
STATE_HEADER = ['Signature', 'Description', 'No. genes', 'Source', 'Download']
GENE_COLUMNS = ('EnsembleID', 'GeneName')
_FILE = re.compile(r'download/signature/([A-Za-z_]+)\.txt$')


def _states(opener) -> list[dict[str, Any]]:
    """The rows of the table of states; a page without it raises."""
    for table in iter_html_tables(read_text(opener)):
        rows = table['rows']
        if rows and [cell['text'] for cell in rows[0]] == STATE_HEADER:
            break
    else:
        raise ValueError('CancerSEA download page has no table of states')
    states = []
    for row in rows[1:]:
        if len(row) != len(STATE_HEADER) or len(row[4]['links']) != 1:
            raise ValueError(f'CancerSEA state row of unexpected shape: {[cell["text"] for cell in row]!r}')
        link = _FILE.search(row[4]['links'][0])
        if link is None or not row[2]['text'].isdigit():
            raise ValueError(f'CancerSEA state row of unexpected shape: {[cell["text"] for cell in row]!r}')
        states.append({
            'state': row[0]['text'],
            'description': row[1]['text'],
            'n_genes': int(row[2]['text']),
            'sources': row[3]['text'],
            'file': link[1],
        })
    if not states:
        raise ValueError('CancerSEA download page has no states')
    return states


def iter_states(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Yield one record per functional state of the page.

    Output fields: ``state`` (the name as written, ``Cell Cycle``),
    ``description``, ``n_genes`` (the number of genes the page gives) and
    ``sources`` (what the page lists as the sources of the gene set: PubMed
    identifiers and database names).
    """
    for state in _states(opener):
        yield {key: state[key] for key in ('state', 'description', 'n_genes', 'sources')}


def iter_gene_sets(
    opener,
    *,
    downloads: Mapping[str, Download],
    force_refresh: bool = False,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield one record per gene of each functional state's gene set.

    ``downloads`` maps the file name of a state (``Cell_Cycle``) to its
    download. A state of the page without a download, a download for a state
    that is not on the page, a gene set without genes and a gene set whose
    size is not the one on the page raise ``ValueError``.

    Output fields: ``state`` and ``description`` as in :func:`iter_states`,
    ``ensembl_id`` and ``symbol`` (``EnsembleID`` and ``GeneName``).
    """
    states = _states(opener)
    files = {state['file'] for state in states}
    if files != set(downloads):
        raise ValueError(
            f'CancerSEA states and gene set downloads differ: {sorted(files ^ set(downloads))}'
        )
    for state in states:
        genes = list(
            iter_tsv(
                downloads[state['file']].open(force_refresh=force_refresh, **kwargs),
                required=GENE_COLUMNS,
            )
        )
        if len(genes) != state['n_genes']:
            raise ValueError(
                f'CancerSEA {state["state"]}: the page says {state["n_genes"]} genes, the file has {len(genes)}'
            )
        for gene in genes:
            yield {
                'state': state['state'],
                'description': state['description'],
                'ensembl_id': gene['EnsembleID'].strip(),
                'symbol': gene['GeneName'].strip(),
            }
