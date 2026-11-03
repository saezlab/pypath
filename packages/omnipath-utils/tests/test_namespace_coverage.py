"""Namespace coverage & SwissLipids structure-key well-formedness (spec 011
T084/T085, User Story 5).

Runs against a built instance via ``OMNIPATH_UTILS_TEST_DB``. Every chemical
namespace the build declares (``id_types.yaml``, ``entity_type:
small_molecule`` with a non-empty ``backends`` mapping) must hold rows
somewhere (``id_mapping`` or ``id_mapping_long``, either direction) or have a
reviewed exemption (``namespace_exemption`` -- research.md R7: "the
exemption is data, not silence").
"""

import os
import re

import pytest
from sqlalchemy import text

pytestmark = pytest.mark.skipif(
    not os.environ.get('OMNIPATH_UTILS_TEST_DB'),
    reason='set OMNIPATH_UTILS_TEST_DB to a built instance to run',
)

_INCHIKEY_PATTERN = re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$')


@pytest.fixture(scope='module')
def session():
    from sqlalchemy.orm import Session

    from omnipath_utils.db._connection import get_engine

    engine = get_engine(os.environ['OMNIPATH_UTILS_TEST_DB'])
    with Session(engine) as s:
        yield s


def _declared_chemical_namespaces() -> list[str]:
    """Every namespace id_types.yaml declares reachable via some loader."""
    from omnipath_utils.mapping._id_types import IdTypeRegistry

    registry = IdTypeRegistry.get()
    names = []
    for name in registry.all_names():
        info = registry.info(name)
        if info.get('entity_type') == 'small_molecule' and info.get('backends'):
            names.append(name)
    return sorted(names)


def _has_rows(session, namespace: str) -> bool:
    row = session.execute(
        text("""
            SELECT
                EXISTS(
                    SELECT 1 FROM omnipath_utils.id_mapping m
                    JOIN omnipath_utils.id_type it
                      ON it.id IN (m.source_type_id, m.target_type_id)
                    WHERE it.name = :ns
                )
                OR EXISTS(
                    SELECT 1 FROM omnipath_utils.id_mapping_long m
                    JOIN omnipath_utils.id_type it
                      ON it.id IN (m.source_type_id, m.target_type_id)
                    WHERE it.name = :ns
                )
        """),
        {'ns': namespace},
    ).scalar()
    return bool(row)


def _is_exempt(session, namespace: str) -> bool:
    row = session.execute(
        text(
            'SELECT 1 FROM omnipath_utils.namespace_exemption WHERE namespace = :ns'
        ),
        {'ns': namespace},
    ).first()
    return row is not None


@pytest.mark.parametrize('namespace', _declared_chemical_namespaces())
def test_namespace_has_rows_or_a_recorded_exemption(session, namespace):
    has_rows = _has_rows(session, namespace)
    exempt = _is_exempt(session, namespace) if not has_rows else False
    assert has_rows or exempt, (
        f'{namespace!r} has no id_mapping/id_mapping_long rows and no '
        f'recorded namespace_exemption -- a namespace may never fail silently'
    )


def test_every_exemption_has_a_real_reason(session):
    rows = session.execute(
        text('SELECT namespace, reason FROM omnipath_utils.namespace_exemption')
    ).fetchall()
    assert rows, 'no exemptions recorded at all'
    for namespace, reason in rows:
        assert reason and len(reason) > 20, (
            f'{namespace!r} exemption has no real reason: {reason!r}'
        )


def test_an_exempt_namespace_is_not_also_secretly_loaded(session):
    """An exemption that turns out to have rows anyway is stale, not honest."""
    rows = session.execute(
        text('SELECT namespace FROM omnipath_utils.namespace_exemption')
    ).fetchall()
    for (namespace,) in rows:
        assert not _has_rows(session, namespace), (
            f'{namespace!r} is recorded exempt but actually has rows now -- '
            f'remove the exemption instead of leaving it stale'
        )


# --------------------------------------------------------- T085: SwissLipids


def test_swisslipids_structure_keys_are_well_formed(session):
    """Every SwissLipids -> InChIKey value in id_mapping is a bare, 27-char
    key -- no leftover 'InChIKey=' format prefix (T086's fix)."""
    rows = session.execute(
        text("""
            SELECT m.target_id FROM omnipath_utils.id_mapping m
            JOIN omnipath_utils.id_type s ON s.id = m.source_type_id
              AND s.name = 'swisslipids'
            JOIN omnipath_utils.id_type t ON t.id = m.target_type_id
              AND t.name = 'inchikey'
        """)
    ).fetchall()
    assert rows, 'no swisslipids -> inchikey rows to check'
    malformed = [v for (v,) in rows if not _INCHIKEY_PATTERN.match(v)]
    assert not malformed, (
        f'{len(malformed)} malformed SwissLipids structure keys, e.g. '
        f'{malformed[:5]!r}'
    )


def test_swisslipids_reaches_a_real_structure(session):
    """A known SwissLipids id translates to a real, well-formed InChIKey."""
    row = session.execute(
        text("""
            SELECT m.target_id FROM omnipath_utils.id_mapping m
            JOIN omnipath_utils.id_type s ON s.id = m.source_type_id
              AND s.name = 'swisslipids'
            JOIN omnipath_utils.id_type t ON t.id = m.target_type_id
              AND t.name = 'inchikey'
            LIMIT 1
        """)
    ).first()
    assert row is not None
    assert _INCHIKEY_PATTERN.match(row[0])
