"""Orthology endpoints."""

from __future__ import annotations

from litestar import Controller, get
from litestar.params import Parameter


class OrthologyController(Controller):
    path = '/orthology'

    @get('/translate')
    async def translate_orthologs(
        self,
        identifiers: str = Parameter(description='Comma-separated identifiers'),
        source: int = Parameter(default=9606, description='Source organism'),
        target: int = Parameter(default=10090, description='Target organism'),
        id_type: str = Parameter(default='genesymbol', description='ID type'),
        resource: str | None = Parameter(
            default=None, required=False, description='Force resource'
        ),
        min_sources: int = Parameter(
            default=1, description='Min supporting DBs (HCOP)'
        ),
        raw: bool = Parameter(
            default=False, description='Skip post-processing'
        ),
    ) -> dict:
        """Translate identifiers to orthologs in another organism."""
        from omnipath_utils.orthology import translate

        id_list = [i.strip() for i in identifiers.split(',') if i.strip()]

        result = translate(
            id_list,
            source=source,
            target=target,
            id_type=id_type,
            resource=resource,
            min_sources=min_sources,
            raw=raw,
        )

        mapped = {k: sorted(v) for k, v in result.items() if v}
        unmapped = [i for i in id_list if not result.get(i)]

        return {
            'results': mapped,
            'unmapped': unmapped,
            'meta': {
                'source': source,
                'target': target,
                'id_type': id_type,
                'resource': resource,
                'min_sources': min_sources,
                'total_input': len(id_list),
                'total_mapped': len(mapped),
            },
        }

    @get('/table')
    async def get_orthology_table(
        self,
        source: int = Parameter(default=9606),
        target: int = Parameter(default=10090),
        id_type: str = Parameter(default='genesymbol'),
        resource: str | None = Parameter(default=None, required=False),
    ) -> dict:
        """Get a full orthology table."""

        from omnipath_utils.orthology import get_table

        result = get_table(
            source=source,
            target=target,
            id_type=id_type,
            resource=resource,
        )

        mapped = {k: sorted(v) for k, v in result.items() if v}

        return {
            'table': mapped,
            'meta': {
                'source': source,
                'target': target,
                'id_type': id_type,
                'resource': resource,
                'count': len(mapped),
            },
        }

