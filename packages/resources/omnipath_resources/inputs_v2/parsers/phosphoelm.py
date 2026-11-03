"""Read the Phospho.ELM dump and the kinase list of the Phospho.ELM site.

The dump is the member ``phosphoELM_vertebrate_<date>.dump`` of the archive
``phosphoELM_vertebrate_latest.dump.tgz``: a tab separated table with a header
row and one line for each phosphorylation site, kinase and publication.

    acc     sequence  position  code  pmids     kinases  source  species       entry_date
    P00519  MGQQP...  412       Y     10964922  Abl      LTP     Homo sapiens  2005-03-14 ...

``sequence`` is the whole sequence of the protein ``acc`` names, ``position``
the place of the residue in it, ``code`` the residue, ``pmids`` one PubMed ID
(or ``N.N.`` for none), ``kinases`` one kinase name or nothing, ``source``
whether the site was found in a low-throughput (LTP) or high-throughput (HTP)
study, ``species`` the Latin name of the organism of the protein.

The kinases are Phospho.ELM's own names ("Abl", "Aurora A", "CK2_group"), not
gene symbols. The kinase list (http://phospho.elm.eu.org/kinases.html) is a
table of these names, each with the UniProt accession Phospho.ELM assigns to
it, for the names that have one.
"""

from __future__ import annotations

from collections.abc import Collection, Generator
import re
from types import SimpleNamespace
from typing import Any

from pypath.inputs_v2._publications import is_pmid
from pypath.inputs_v2.base import Download
from pypath.inputs_v2.parsers.base import iter_html_tables, iter_tsv, read_text

# Columns of the dump, in source order.
COLUMNS = (
    'acc',
    'sequence',
    'position',
    'code',
    'pmids',
    'kinases',
    'source',
    'species',
    'entry_date',
)
# The member of the archive: the file carries the date of the release.
_MEMBER_REGEX = re.compile(r'^phosphoELM_vertebrate_[\w.-]+\.dump$')
# The "source" column.
THROUGHPUT_LABELS = frozenset({'LTP', 'HTP'})
_UNIPROT = (
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|'
    r'[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
)
# The accession of a protein: UniProt with the number of the isoform where the
# protein is one, or an Ensembl protein.
_UNIPROT_REGEX = re.compile(rf'^{_UNIPROT}(?:-\d+)?$')
_ENSEMBL_PROTEIN_REGEX = re.compile(r'^ENS[A-Z]*P\d+$')
_RESIDUE_REGEX = re.compile(r'^[A-Z]$')
# The reference of a line: a PubMed ID, "N.N." (nomen nescio) or nothing. Of
# the PubMed IDs, 17 lines have six digits (11 distinct, all of 1975-1979). At
# NCBI nine are titled on the substrate's phosphorylation, kinase or protein;
# 176960 and 934118 are not (a pituitary study, a paper on childhood
# constipation). All are kept, as before.
_NO_PMID = frozenset({'N.N.', ''})
_DATE_REGEX = re.compile(r'^\d{4}-\d\d-\d\d \d\d:\d\d:\d\d(?:\.\d+)?\+\d\d$')
# A name of the kinase list that stands for a family of kinases.
GROUP_SUFFIX = '_group'
_KINASE_LIST_COLUMNS = (
    'Name',
    '#Instances',
    'Uniprot',
    'HUGO Id',
    'Classification',
    'Description',
)


def kinase_accessions(opener) -> dict[str, str]:
    """The UniProt accession the kinase list of Phospho.ELM gives to each kinase name.

    ``opener`` is the download of the page. The kinase list is its table with
    the id ``kinaselist``; names without an accession are not in the result.
    A page without the table, with other columns, a repeated name or an
    accession of another form is an error.
    """
    tables = [
        table
        for table in iter_html_tables(read_text(opener))
        if table['attrs'].get('id') == 'kinaselist'
    ]
    if len(tables) != 1:
        raise ValueError(
            f'Phospho.ELM kinase page has {len(tables)} tables "kinaselist", expected 1'
        )
    header, *rows = tables[0]['rows']
    if tuple(cell['text'] for cell in header) != _KINASE_LIST_COLUMNS:
        raise ValueError(
            f'Phospho.ELM kinase list has other columns: {[c["text"] for c in header]}'
        )
    if not rows:
        raise ValueError('Phospho.ELM kinase list has no rows')
    names: set[str] = set()
    accessions: dict[str, str] = {}
    for cells in rows:
        if len(cells) != len(_KINASE_LIST_COLUMNS):
            raise ValueError(f'Phospho.ELM kinase list row has {len(cells)} cells')
        name, accession = cells[0]['text'], cells[2]['text']
        if not name or name in names:
            raise ValueError(f'Phospho.ELM kinase list: name {name!r} is empty or repeated')
        names.add(name)
        if accession:
            if not re.fullmatch(_UNIPROT, accession):
                raise ValueError(
                    f'Phospho.ELM kinase list: accession {accession!r} of {name}'
                )
            accessions[name] = accession
    return accessions


def _dump_handle(opener):
    """The one member of the archive that is the dump, as an opener."""
    members = opener.result if opener else None
    if not isinstance(members, dict):
        raise ValueError('Phospho.ELM download is not an archive')
    dumps = [name for name in members if _MEMBER_REGEX.match(name)]
    if len(dumps) != 1:
        raise ValueError(
            f'Phospho.ELM archive has {len(dumps)} dump files, expected 1: {list(members)}'
        )
    return SimpleNamespace(result=members[dumps[0]])


def iter_phosphoelm(
    opener,
    *,
    organisms: Collection[str],
    kinases_download: Download | None = None,
    kinase_required: bool = True,
    force_refresh: bool = False,
    **_kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per line of the dump.

    ``opener`` is the archive. A line with a species that is not in
    ``organisms``, an accession, residue, reference, throughput or date of
    another form, or a residue that is not the one at ``position`` in
    ``sequence`` is an error, as is an archive without the dump, with other
    columns or without lines.

    By default a line without a kinase is not yielded: it is a phosphorylation
    site, not a kinase-substrate relation (``kinase_required=False`` yields it
    as well, for counting). The lines of one site and kinase for several
    publications, and repeated lines, are all yielded as they are.

    ``kinases_download`` is the ``Download`` of the kinase list page; with it,
    a kinase name that has an accession there gets it in ``enzyme_uniprot``.

    Each dictionary has the columns as written, except ``sequence``, which is
    left out, and:

    - ``position``: as an integer.
    - ``accession``: ``acc``, with the isoform suffix if it has one.
    - ``accession_database``: ``uniprot`` or ``ensembl_protein``.
    - ``pmid``: the PubMed ID of ``pmids``; None for ``N.N.`` and for no value.
    - ``enzyme``: the kinase name of ``kinases``; None if empty.
    - ``enzyme_is_group``: whether the name ends in ``_group``.
    - ``enzyme_uniprot``: the UniProt accession of the kinase list; None if
      the kinase has none there or ``kinases_download`` is not given.
    """
    accessions = (
        kinase_accessions(kinases_download.open(force_refresh=force_refresh))
        if kinases_download
        else {}
    )
    count = 0
    for row in iter_tsv(_dump_handle(opener), strict=True, quotechar=None):
        if tuple(row) != COLUMNS:
            raise ValueError(f'Phospho.ELM dump has columns {list(row)}')
        count += 1
        sequence = row.pop('sequence')
        position, residue = row['position'], row['code']
        if not position.isdigit() or int(position) < 1:
            raise ValueError(f'Phospho.ELM position not understood: {position!r}')
        position = int(position)
        if not _RESIDUE_REGEX.match(residue) or sequence[position - 1 : position] != residue:
            raise ValueError(
                f'Phospho.ELM {row["acc"]}: residue {residue!r} is not at position {position}'
            )
        accession = row['acc']
        if _UNIPROT_REGEX.match(accession):
            database = 'uniprot'
        elif _ENSEMBL_PROTEIN_REGEX.match(accession):
            database = 'ensembl_protein'
        else:
            raise ValueError(f'Phospho.ELM accession not understood: {accession!r}')
        if not (row['pmids'].isdigit() or row['pmids'] in _NO_PMID):
            raise ValueError(f'Phospho.ELM reference not understood: {row["pmids"]!r}')
        if row['source'] not in THROUGHPUT_LABELS:
            raise ValueError(f'Unknown Phospho.ELM throughput: {row["source"]!r}')
        if row['species'] not in organisms:
            raise ValueError(f'Unknown Phospho.ELM species: {row["species"]!r}')
        if not _DATE_REGEX.match(row['entry_date']):
            raise ValueError(f'Phospho.ELM date not understood: {row["entry_date"]!r}')
        enzyme = row['kinases'] or None
        if enzyme is None and kinase_required:
            continue
        yield {
            **row,
            'position': position,
            'accession': accession,
            'accession_database': database,
            'pmid': row['pmids'] if is_pmid(row['pmids']) else None,
            'enzyme': enzyme,
            'enzyme_is_group': bool(enzyme and enzyme.endswith(GROUP_SUFFIX)),
            'enzyme_uniprot': accessions.get(enzyme),
        }
    if not count:
        raise ValueError('Phospho.ELM dump has no rows')
