"""
Parser of the protein types table of the Cell-Cell Interaction Database.

The Bader lab's table has a header line and two columns:

    hgnc_symbol   curated_protein_type
    SBSPON        ECM/Receptor/Ligand
    A2M           Ligand

``curated_protein_type`` is the manually curated classification of the gene as
ligand, receptor and extracellular matrix (ECM) protein. A gene in more than one
of the three lists has all its roles in the cell, joined by ``/``.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_tsv

COLUMNS = ('hgnc_symbol', 'curated_protein_type')


def iter_genes(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """
    One dictionary per gene of the table.

    - ``symbol``: the HGNC symbol, as written.
    - ``roles``: the roles of the gene, in the order of the cell
      (``ECM``, ``Ligand``, ``Receptor``).
    - ``protein_type``: the cell as written.

    A symbol that occurs twice, a row without a symbol or a role, or a role that
    is repeated within a cell, is an error. The role names are checked by the
    mapper.
    """
    seen = set()
    for row in iter_tsv(opener, quotechar=None, strict=True, required=COLUMNS):
        symbol = (row['hgnc_symbol'] or '').strip()
        protein_type = (row['curated_protein_type'] or '').strip()
        roles = [role.strip() for role in protein_type.split('/')]
        if not symbol or not all(roles) or len(set(roles)) != len(roles):
            raise ValueError(
                f'Cell-Cell Interaction Database row not understood: {row!r}'
            )
        if symbol in seen:
            raise ValueError(
                f'Cell-Cell Interaction Database lists {symbol} twice'
            )
        seen.add(symbol)
        yield {'symbol': symbol, 'roles': roles, 'protein_type': protein_type}
    if not seen:
        raise ValueError('Cell-Cell Interaction Database table has no rows')


def iter_roles(opener, **kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """The distinct roles of the table, in order of first appearance."""
    seen = {}
    for gene in iter_genes(opener, **kwargs):
        for role in gene['roles']:
            seen.setdefault(role, None)
    for role in seen:
        yield {'role': role}
