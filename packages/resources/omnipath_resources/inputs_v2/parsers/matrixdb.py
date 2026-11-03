"""
Parsers of the MatrixDB downloads.

MatrixDB (https://matrixdb.expasy.org/download) offers, among others:

- ``matrixdb_CORE_4_0.tab.zip``: the manually curated interactions of MatrixDB as a
  PSI-MI TAB 2.7 file, ``matrixdb_CORE.tab``. The header has 42 named columns and
  370 empty ones after them, and every row ends with empty cells. Interactors are
  proteins (UniProt), chemicals (ChEBI), complexes (Complex Portal) and the
  biomolecules that MatrixDB identifies itself, such as glycosaminoglycans
  (``matrixdb:GAG_1``). ``matrixdb:`` entries in the alternative identifier column
  give the MatrixDB identifier of the interactor.
- ``ecm_proteins.csv``: the extracellular matrix (ECM) and secreted proteins of
  UniProtKB, selected by the Gene Ontology cellular component terms and UniProt
  keywords of its own criteria table, with columns ``protein_id``,
  ``species (NCBI)``, ``go_terms``, ``uniprot_keywords``, ``dataset`` and
  ``interaction_count``.
"""

from __future__ import annotations

import re
from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_csv, iter_tsv

# The 42 columns of the file, in order, and the field each is read into.
INTERACTION_COLUMNS = (
    ('#ID(s) interactor A', 'id_a'),
    ('ID(s) interactor B', 'id_b'),
    ('Alt. ID(s) interactor A', 'alt_ids_a'),
    ('Alt. ID(s) interactor B', 'alt_ids_b'),
    ('Alias(es) interactor A', 'aliases_a'),
    ('Alias(es) interactor B', 'aliases_b'),
    ('Interaction detection method(s)', 'detection_methods'),
    ('Publication 1st author(s)', 'first_authors'),
    ('Publication Identifier(s)', 'publications'),
    ('Taxid interactor A', 'taxid_a'),
    ('Taxid interactor B', 'taxid_b'),
    ('Interaction type(s)', 'interaction_types'),
    ('Source database(s)', 'source_databases'),
    ('Interaction identifier(s)', 'interaction_ids'),
    ('Confidence value(s)', 'confidence'),
    ('Expansion method(s)', 'expansion_methods'),
    ('Biological role(s) interactor A', 'biological_role_a'),
    ('Biological role(s) interactor B', 'biological_role_b'),
    ('Experimental role(s) interactor A', 'experimental_role_a'),
    ('Experimental role(s) interactor B', 'experimental_role_b'),
    ('Type(s) interactor A', 'types_a'),
    ('Type(s) interactor B', 'types_b'),
    ('Xref(s) interactor A', 'xrefs_a'),
    ('Xref(s) interactor B', 'xrefs_b'),
    ('Interaction Xref(s)', 'interaction_xrefs'),
    ('Annotation(s) interactor A', 'annotations_a'),
    ('Annotation(s) interactor B', 'annotations_b'),
    ('Interaction annotation(s)', 'interaction_annotations'),
    ('Host organism(s)', 'host_organisms'),
    ('Interaction parameter(s)', 'parameters'),
    ('Creation date', 'creation_date'),
    ('Update date', 'update_date'),
    ('Checksum(s) interactor A', 'checksums_a'),
    ('Checksum(s) interactor B', 'checksums_b'),
    ('Interaction Checksum(s)', 'interaction_checksums'),
    ('Negative', 'negative'),
    ('Feature(s) interactor A', 'features_a'),
    ('Feature(s) interactor B', 'features_b'),
    ('Stoichiometry(s) interactor A', 'stoichiometry_a'),
    ('Stoichiometry(s) interactor B', 'stoichiometry_b'),
    ('Identification method participant A', 'identification_method_a'),
    ('Identification method participant B', 'identification_method_b'),
)
# The columns the module reads; the others are in the row but unused.
INTERACTION_FIELDS = tuple(field for _, field in INTERACTION_COLUMNS)


def _unquote(cell: str) -> str:
    """
    Undo the quoting of a cell of the interaction file.

    The file was written with CSV quoting but is read as tab separated text: a cell
    that has a double quote is either wrapped in quotes, with its own quotes
    doubled (``"psi-mi:""MI:0407""(direct interaction)"``), or has its quotes
    doubled without the wrapping (``psi-mi:""MI:0407""(direct interaction)``).
    Both mean ``psi-mi:"MI:0407"(direct interaction)``.
    """
    cell = cell.strip()
    if len(cell) > 1 and cell[0] == '"' and cell[-1] == '"':
        cell = cell[1:-1]
    return cell.replace('""', '"')


def iter_interactions(opener, **_kwargs: Any) -> Generator[dict[str, str], None, None]:
    """
    One dictionary per row of the interaction file.

    The fields are named in ``INTERACTION_COLUMNS``; every cell is text with its
    quoting undone (see ``_unquote``), ``-`` standing for no value. The cells after
    the 42nd must be empty. One row of the 4.0 release (the sixth, the first of the
    Ptprc interactions of Ori et al. 2011) stops after the interaction annotation
    column: the cells from the host organism on, among them ``negative``, are empty
    instead of ``-`` or ``false``. An empty cell is read as ``-``, so that row has
    no negative flag; the 1659 other rows say ``false``. A different header, a row
    with other content after the 42nd cell, a negative column other than ``false``,
    ``true`` or none, or an empty file is an error.
    """
    expected = [name for name, _ in INTERACTION_COLUMNS]
    count = 0
    for number, row in enumerate(
        iter_tsv(opener, quotechar=None, allow_repeated=True), start=2
    ):
        names = [name for name in row if name is not None]
        # The empty names of the padding are one key, the last of them.
        if names[: len(expected)] != expected or names[len(expected) :] not in ([], ['']):
            raise ValueError(f'MatrixDB interaction file has columns {names[:44]}')
        extra = [row.get('')] + (row.get(None) or [])
        if any(cell and cell.strip() for cell in extra):
            raise ValueError(f'MatrixDB interaction row {number} has content after the 42nd cell')
        record = {field: _unquote(row[name] or '') or '-' for name, field in INTERACTION_COLUMNS}
        if record['negative'] not in ('false', 'true', '-'):
            raise ValueError(
                f'MatrixDB interaction row {number}: negative is {record["negative"]!r}'
            )
        count += 1
        yield record
    if not count:
        raise ValueError('MatrixDB interaction file has no rows')


ECM_COLUMNS = (
    'protein_id',
    'species (NCBI)',
    'go_terms',
    'uniprot_keywords',
    'dataset',
    'interaction_count',
)
_ACCESSION = re.compile(
    r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
)


def iter_ecm_proteins(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """
    One dictionary per protein of ``ecm_proteins.csv``.

    - ``accession``: the UniProt accession.
    - ``taxon``: the NCBI taxon id, digits.
    - ``go_terms``, ``keywords``: the names of the GO cellular component terms and
      UniProt keywords of the criteria that select the protein, as the file writes
      them (the cells separate several by ``|``; empty gives an empty list).
    - ``dataset``: ``Swiss-Prot`` or ``UniProtKB unreviewed (TrEMBL)``.
    - ``interaction_count``: the number of MatrixDB interactions of the protein, an
      integer.

    A protein listed twice, an accession or taxon that is not one, a count that is
    not a whole number or an empty file is an error.
    """
    seen = set()
    for row in iter_csv(opener, required=ECM_COLUMNS):
        accession = (row['protein_id'] or '').strip()
        taxon = (row['species (NCBI)'] or '').strip()
        count = (row['interaction_count'] or '').strip()
        if (
            not _ACCESSION.fullmatch(accession)
            or not taxon.isdigit()
            or not count.isdigit()
            or accession in seen
        ):
            raise ValueError(f'MatrixDB ECM protein row not understood: {row!r}')
        seen.add(accession)
        yield {
            'accession': accession,
            'taxon': taxon,
            'go_terms': [t for t in (row['go_terms'] or '').split('|') if t],
            'keywords': [t for t in (row['uniprot_keywords'] or '').split('|') if t],
            'dataset': (row['dataset'] or '').strip(),
            'interaction_count': int(count),
        }
    if not seen:
        raise ValueError('MatrixDB ECM protein file has no rows')


def iter_list_terms(opener, **kwargs: Any) -> Generator[dict[str, str], None, None]:
    """The one list that ``ecm_proteins.csv`` is, after checking that it has rows."""
    for _ in iter_ecm_proteins(opener, **kwargs):
        break
    yield {'list': 'ecm_proteins'}
