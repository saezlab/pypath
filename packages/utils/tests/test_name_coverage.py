"""C1/C2/C3/C5 -- per-resource name coverage & no-regression (DB-backed).

Runs against a built instance via OMNIPATH_UTILS_TEST_DB. Encodes the
acceptance SQL in contracts/name-coverage.md.
"""

import os

import pytest
from sqlalchemy import text

pytestmark = pytest.mark.skipif(
    not os.environ.get('OMNIPATH_UTILS_TEST_DB'),
    reason='set OMNIPATH_UTILS_TEST_DB to a built instance to run',
)

# SC-008/C3 baseline (utils2, 2026-06-23; see specs baseline-counts.txt).
BASELINE_ID_MAPPING = 213_220_981
BASELINE_ID_MAPPING_FTP = 744_242_420

NAME_TYPES = ('name', 'synonym', 'iupac', 'traditional_iupac')


@pytest.fixture(scope='module')
def session():
    from sqlalchemy.orm import Session

    from omnipath_utils.db._connection import get_engine

    engine = get_engine(os.environ['OMNIPATH_UTILS_TEST_DB'])
    with Session(engine) as s:
        yield s


def _name_rows_by_backend(session):
    rows = session.execute(text("""
        SELECT b.name, count(*) FROM omnipath_utils.id_mapping_long m
        JOIN omnipath_utils.backend b ON b.id = m.backend_id
        JOIN omnipath_utils.id_type s ON s.id = m.source_type_id
        WHERE s.name IN ('name','synonym','iupac','traditional_iupac')
        GROUP BY b.name
    """)).fetchall()
    return {r[0]: r[1] for r in rows}


@pytest.mark.parametrize(
    'backend', ['chebi', 'chembl', 'hmdb', 'ramp', 'kegg_compound'],
)
def test_c1_resource_has_name_rows(session, backend):
    counts = _name_rows_by_backend(session)
    assert counts.get(backend, 0) > 0, f'{backend} has no name rows'


def test_c2_name_types_reach_two_targets_incl_chebi(session):
    rows = session.execute(text("""
        SELECT s.name, count(DISTINCT t.name) AS n,
               string_agg(DISTINCT t.name, ',') AS targets
        FROM omnipath_utils.id_mapping_long m
        JOIN omnipath_utils.id_type s ON s.id = m.source_type_id
        JOIN omnipath_utils.id_type t ON t.id = m.target_type_id
        WHERE s.name IN ('name','synonym','iupac','traditional_iupac')
        GROUP BY s.name
    """)).fetchall()
    assert rows, 'no name source types present'
    for name_type, n_targets, targets in rows:
        assert n_targets >= 2, f'{name_type} reaches <2 targets: {targets}'
        assert 'chebi' in targets.split(','), f'{name_type} misses chebi'


def test_c3_database_id_tables_unchanged(session):
    n = session.execute(
        text('SELECT count(*) FROM omnipath_utils.id_mapping')
    ).scalar()
    assert n == BASELINE_ID_MAPPING, 'id_mapping row count changed (regression)'
    nftp = session.execute(
        text('SELECT count(*) FROM omnipath_utils.id_mapping_ftp')
    ).scalar()
    assert nftp == BASELINE_ID_MAPPING_FTP, 'id_mapping_ftp changed (regression)'


def test_c5_idempotent_slice_reload(session):
    """Re-running one slice leaves an equal id_mapping_long count (R6)."""
    from omnipath_utils.db._build import DatabaseBuilder

    builder = DatabaseBuilder(db_url=os.environ['OMNIPATH_UTILS_TEST_DB'])
    before = session.execute(text("""
        SELECT count(*) FROM omnipath_utils.id_mapping_long m
        JOIN omnipath_utils.backend b ON b.id=m.backend_id
        JOIN omnipath_utils.id_type s ON s.id=m.source_type_id
        JOIN omnipath_utils.id_type t ON t.id=m.target_type_id
        WHERE b.name='chebi' AND s.name='name' AND t.name='chebi'
    """)).scalar()
    data = {'taurine': {'CHEBI:15891'}}
    # reload a 1-row probe slice on a throwaway backend, then re-run -> equal
    builder._populate_long_slice(data, 'name', 'chebi', 'kegg_compound')
    a = builder._populate_long_slice(data, 'name', 'chebi', 'kegg_compound')
    b = builder._populate_long_slice(data, 'name', 'chebi', 'kegg_compound')
    assert a == b == 1
    session.commit()
    after = session.execute(text("""
        SELECT count(*) FROM omnipath_utils.id_mapping_long m
        JOIN omnipath_utils.backend b ON b.id=m.backend_id
        JOIN omnipath_utils.id_type s ON s.id=m.source_type_id
        JOIN omnipath_utils.id_type t ON t.id=m.target_type_id
        WHERE b.name='chebi' AND s.name='name' AND t.name='chebi'
    """)).scalar()
    assert after == before  # the chebi slice is untouched by the probe
