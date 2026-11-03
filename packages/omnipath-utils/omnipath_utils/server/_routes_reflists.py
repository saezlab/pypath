"""Reference list endpoints."""

from __future__ import annotations

from litestar import Controller, get
from sqlalchemy import text
from sqlalchemy.orm import Session
from litestar.params import Parameter

from omnipath_utils.db._connection import SCHEMA


class ReflistController(Controller):
    """Reference list endpoints."""

    path = '/reflists'

    @get('/list-names')
    async def list_names(self) -> list[str]:
        """List available reference list names."""
        return ['swissprot', 'trembl', 'uniprot']

    @get('/{list_name:str}')
    async def get_reflist(
        self,
        session: Session,
        list_name: str,
        ncbi_tax_id: int = Parameter(default=9606),
    ) -> dict:
        """Get a reference list from the database."""

        rows = session.execute(
            text(
                f'SELECT identifier FROM {SCHEMA}.reflist '
                'WHERE list_name = :name AND ncbi_tax_id = :tax'
            ),
            {'name': list_name, 'tax': ncbi_tax_id},
        ).fetchall()

        ids = sorted(row[0] for row in rows)

        return {
            'list_name': list_name,
            'ncbi_tax_id': ncbi_tax_id,
            'count': len(ids),
            'identifiers': ids,
        }
