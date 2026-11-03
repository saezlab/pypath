"""KEGG organism code access."""

from __future__ import annotations

import logging
import collections

from pypath.share.downloads import dm

_log = logging.getLogger(__name__)

KEGG_LIST_URL = 'https://rest.kegg.jp/list/organism'

KeggOrganism = collections.namedtuple(
    'KeggOrganism',
    ('kegg_id', 'code', 'name', 'common_name', 'lineage'),
)


def kegg_organisms():
    """Download KEGG organism list.

    Yields KeggOrganism named tuples.
    """
    path = dm.download(KEGG_LIST_URL)
    if not path:
        _log.error('Failed to download KEGG organism list')
        return

    with open(path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue

            kegg_id = parts[0]
            code = parts[1]
            name_field = parts[2]
            lineage = parts[3] if len(parts) > 3 else ''

            # Parse "Homo sapiens (human)"
            if '(' in name_field:
                name = name_field.split('(')[0].strip()
                common = name_field.split('(')[1].rstrip(')').strip()
            else:
                name = name_field.strip()
                common = ''

            yield KeggOrganism(
                kegg_id=kegg_id,
                code=code,
                name=name,
                common_name=common,
                lineage=lineage,
            )


def kegg_code_to_taxid():
    """Build a KEGG code -> name mapping for cross-referencing."""
    return {org.code: org for org in kegg_organisms()}
