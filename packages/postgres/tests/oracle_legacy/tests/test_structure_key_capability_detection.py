"""spec 011 T082 -- structure_key_computation reads a real signal instead
of a hardcoded 'unavailable'.

The signal is the utils build's own build_info row (spec 011, omnipath-utils
DatabaseBuilder.record_structure_key_capability). Unit-tested with a fake
psycopg2 connection -- the DB-backed integration check lives in
test_build_capability.py.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

from omnipath_build.db.resources import _detect_structure_key_computation


def _fake_connect(row):
    """A fake psycopg2.connect(...) whose cursor's fetchone() returns row."""
    conn = MagicMock()
    cur = conn.cursor.return_value.__enter__.return_value
    cur.fetchone.return_value = row
    return conn


class TestNoUtilsUrl:
    def test_no_url_reports_unavailable_with_reason(self):
        result = _detect_structure_key_computation(None)
        assert result['capability'] == 'structure_key_computation'
        assert result['available'] is False
        assert 'OMNIPATH_BUILD_UTILS_PG_URL' in result['reason']


class TestRealSignal:
    def test_available_when_build_recorded_it(self):
        with patch('psycopg2.connect', return_value=_fake_connect(('available',))):
            result = _detect_structure_key_computation('postgresql://fake/db')
        assert result['available'] is True

    def test_unavailable_when_build_recorded_it(self):
        with patch(
            'psycopg2.connect', return_value=_fake_connect(('unavailable',)),
        ):
            result = _detect_structure_key_computation('postgresql://fake/db')
        assert result['available'] is False
        assert 'chemistry' in result['reason'].lower()

    def test_unavailable_when_never_recorded(self):
        # No row at all -- an older utils build that predates T082.
        with patch('psycopg2.connect', return_value=_fake_connect(None)):
            result = _detect_structure_key_computation('postgresql://fake/db')
        assert result['available'] is False
        assert 'never run' in result['reason'].lower() or (
            'not recorded' in result['reason'].lower()
        )

    def test_connection_failure_reports_unavailable_not_raises(self):
        with patch('psycopg2.connect', side_effect=OSError('unreachable')):
            result = _detect_structure_key_computation('postgresql://fake/db')
        assert result['available'] is False
        assert result['reason']
