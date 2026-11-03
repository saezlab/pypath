"""Cumulative-build idempotency guard (007 T005).

``_populate_unichem`` and ``_populate_ramp`` COPY into ``id_mapping``; on a
cumulative rebuild they must DELETE their slice first so re-running does not
append duplicate rows. This checks the resulting invariant directly against a
built utils DB (gated on ``OMNIPATH_UTILS_DB_URL``): no backend has duplicate
``(source_type_id, target_type_id, ncbi_tax_id, source_id, target_id)`` rows.

The full "load twice → count unchanged" check runs during the cumulative build
(007 T021); this fast invariant guards every commit.
"""

from __future__ import annotations

import os

import pytest
from sqlalchemy import create_engine, text

DB_URL = os.environ.get('OMNIPATH_UTILS_DB_URL')

pytestmark = pytest.mark.skipif(
    not DB_URL,
    reason='OMNIPATH_UTILS_DB_URL not set (needs a built utils DB)',
)


@pytest.mark.parametrize('backend', ['unichem', 'ramp'])
def test_no_duplicate_rows_per_backend(backend):
    """No exact-duplicate id_mapping rows for the append-prone backends."""
    engine = create_engine(DB_URL)
    with engine.connect() as conn:
        bk = conn.execute(
            text(
                'SELECT id FROM omnipath_utils.backend WHERE name = :n'
            ),
            {'n': backend},
        ).scalar()

        if bk is None:
            pytest.skip(f'backend {backend!r} not present in this DB')

        dup = conn.execute(
            text(
                'SELECT count(*) FROM ('
                '  SELECT 1 FROM omnipath_utils.id_mapping'
                '  WHERE backend_id = :bk'
                '  GROUP BY source_type_id, target_type_id, ncbi_tax_id,'
                '           source_id, target_id'
                '  HAVING count(*) > 1'
                '  LIMIT 1'
                ') d'
            ),
            {'bk': bk},
        ).scalar()

    assert dup == 0, (
        f'{backend}: found duplicate id_mapping rows — the DELETE-slice-'
        f'before-COPY fix is missing or a load ran twice without it'
    )
