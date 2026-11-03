"""Taxonomy endpoints."""

from __future__ import annotations

from litestar import Controller, get
from sqlalchemy.orm import Session
from litestar.params import Parameter

from omnipath_utils.taxonomy import (
    ensure_oma_code,
    ensure_kegg_code,
    ensure_latin_name,
    ensure_common_name,
    ensure_ncbi_tax_id,
    ensure_ensembl_name,
    ensure_mirbase_name,
)


class TaxonomyController(Controller):
    """Taxonomy / organism resolution endpoints."""

    path = '/taxonomy'

    @get('/resolve')
    async def resolve(
        self,
        organism: str = Parameter(
            description='Organism name, code, or NCBI Taxonomy ID',
        ),
    ) -> dict:
        """Resolve an organism identifier to all name forms."""
        taxid = ensure_ncbi_tax_id(organism)

        if taxid is None:
            return {
                'error': f'Unknown organism: {organism}',
                'ncbi_tax_id': None,
            }

        return {
            'ncbi_tax_id': taxid,
            'common_name': ensure_common_name(taxid),
            'latin_name': ensure_latin_name(taxid),
            'ensembl_name': ensure_ensembl_name(taxid),
            'kegg_code': ensure_kegg_code(taxid),
            'mirbase_code': ensure_mirbase_name(taxid),
            'oma_code': ensure_oma_code(taxid),
        }

    @get('/organisms')
    async def organisms(
        self,
        session: Session,
        has_data: bool = Parameter(
            default=False,
            description='Only organisms with mapping data in the database',
        ),
        limit: int = Parameter(
            default=0,
            description='Max results (0 = no limit)',
        ),
    ) -> list[dict]:
        """List organisms from the database.

        By default returns all organisms. Use has_data=true to show only
        those with mapping data (protein/gene tables loaded).
        """
        from sqlalchemy import text

        from omnipath_utils.db._connection import SCHEMA

        if has_data:
            query = text(f"""
                SELECT DISTINCT o.ncbi_tax_id, o.latin_name, o.common_name,
                       o.ensembl_name, o.kegg_code, o.mirbase_code, o.oma_code
                FROM {SCHEMA}.organism o
                JOIN {SCHEMA}.id_mapping m ON o.ncbi_tax_id = m.ncbi_tax_id
                WHERE m.ncbi_tax_id > 0
                ORDER BY o.ncbi_tax_id
            """)
        else:
            lim = f' LIMIT {limit}' if limit > 0 else ''
            query = text(f"""
                SELECT ncbi_tax_id, latin_name, common_name,
                       ensembl_name, kegg_code, mirbase_code, oma_code
                FROM {SCHEMA}.organism
                ORDER BY ncbi_tax_id
                {lim}
            """)

        rows = session.execute(query).fetchall()

        return [
            {
                'ncbi_tax_id': row[0],
                'latin_name': row[1],
                'common_name': row[2],
                'ensembl_name': row[3],
                'kegg_code': row[4],
                'mirbase_code': row[5],
                'oma_code': row[6],
            }
            for row in rows
        ]
