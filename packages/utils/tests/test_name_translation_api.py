"""spec 011 T063/T067/T068/T070 -- the /mapping/translate HTTP layer for the
generic 'name' axis: the ``name_scope`` parameter, ``matched_as`` in the
response, and no on-demand loading / retry note on a plain miss.

Mocks ``translate_ids`` (same pattern as ``TestMappingEndpoints`` in
test_server.py) so these run without a live DB. The DB-backed contract
tests (C1-C4, C9, C10) live in ``test_name_translation.py``.
"""

from unittest.mock import patch

import pytest

litestar = pytest.importorskip('litestar', reason='litestar not installed')
from sqlalchemy import text as sa_text, create_engine  # noqa: E402
from sqlalchemy.orm import sessionmaker  # noqa: E402
from sqlalchemy.pool import StaticPool  # noqa: E402
from litestar.testing import TestClient  # noqa: E402

from omnipath_utils.server._app import create_app  # noqa: E402


@pytest.fixture
def app():
    engine = create_engine(
        'sqlite://',
        echo=False,
        connect_args={'check_same_thread': False},
        poolclass=StaticPool,
    )
    with engine.connect() as conn:
        for ddl in [
            'CREATE TABLE IF NOT EXISTS id_type (id INTEGER PRIMARY KEY, name TEXT)',
            'CREATE TABLE IF NOT EXISTS backend (id INTEGER PRIMARY KEY, name TEXT)',
            'CREATE TABLE IF NOT EXISTS organism (ncbi_tax_id INTEGER PRIMARY KEY)',
            'CREATE TABLE IF NOT EXISTS id_mapping (id INTEGER PRIMARY KEY)',
            'CREATE TABLE IF NOT EXISTS reflist (id INTEGER PRIMARY KEY)',
            'CREATE TABLE IF NOT EXISTS build_info (id INTEGER PRIMARY KEY, source_type TEXT, target_type TEXT, ncbi_tax_id INTEGER, backend TEXT, row_count INTEGER, built_at TEXT, duration_secs REAL, status TEXT)',
        ]:
            conn.execute(sa_text(ddl))
        conn.commit()

    TestSessionLocal = sessionmaker(bind=engine)
    with patch('omnipath_utils.server._app.create_session_factory') as mock_sf:
        mock_sf.return_value = TestSessionLocal
        _app = create_app(db_url='sqlite://')
        yield _app


@pytest.fixture
def client(app):
    return TestClient(app)


class TestNameScopeParameter:
    @patch('omnipath_utils.server._routes_mapping.translate_ids')
    def test_name_scope_reaches_translate_ids(self, mock_translate, client):
        mock_translate.return_value = ({'GABA': {'CHEBI:16865'}}, {'chebi'})
        client.get(
            '/mapping/translate',
            params={
                'identifiers': 'GABA', 'id_type': 'name',
                'target_id_type': 'chebi', 'name_scope': 'synonym',
            },
        )
        assert mock_translate.call_args.kwargs['name_scope'] == 'synonym'

    @patch('omnipath_utils.server._routes_mapping.translate_ids')
    def test_default_name_scope_when_omitted(self, mock_translate, client):
        mock_translate.return_value = ({}, set())
        client.get(
            '/mapping/translate',
            params={
                'identifiers': 'GABA', 'id_type': 'name',
                'target_id_type': 'chebi',
            },
        )
        assert mock_translate.call_args.kwargs['name_scope'] == 'default'


class TestMatchedAsInResponse:
    @patch('omnipath_utils.server._routes_mapping.translate_ids')
    def test_matched_as_present_for_name_axis(self, mock_translate, client):
        def fake_translate(session, ids, src, tgt, tax, **kw):
            if kw.get('name_meta') is not None:
                kw['name_meta']['GABA'] = 'synonym'
            return {'GABA': {'CHEBI:16865'}}, {'chebi'}

        mock_translate.side_effect = fake_translate
        resp = client.get(
            '/mapping/translate',
            params={
                'identifiers': 'GABA', 'id_type': 'name',
                'target_id_type': 'chebi',
            },
        )
        meta = resp.json()['meta']
        assert meta['matched_as'] == {'GABA': 'synonym'}
        assert meta['name_scope'] == 'default'

    @patch('omnipath_utils.server._routes_mapping.translate_ids')
    def test_matched_as_absent_for_a_non_name_query(self, mock_translate, client):
        mock_translate.return_value = ({'TP53': {'P04637'}}, {'uniprot'})
        resp = client.get(
            '/mapping/translate',
            params={
                'identifiers': 'TP53', 'id_type': 'genesymbol',
                'target_id_type': 'uniprot',
            },
        )
        assert 'matched_as' not in resp.json()['meta']
        assert 'name_scope' not in resp.json()['meta']


class TestNameAxisHonestMiss:
    @patch('omnipath_utils.server._routes_mapping.translate_ids')
    def test_a_name_miss_never_triggers_loading(self, mock_translate, client):
        # C4 / FR: the public service performs no on-demand loading for
        # names, and a miss reports loading:false with no retry note.
        mock_translate.return_value = ({}, set())
        resp = client.get(
            '/mapping/translate',
            params={
                'identifiers': 'not-a-real-name', 'id_type': 'name',
                'target_id_type': 'chebi', 'raw': 'true',
            },
        )
        meta = resp.json()['meta']
        assert meta['loading'] is False
        assert 'loading_note' not in meta
