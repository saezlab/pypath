"""Tests for the web service."""

from unittest.mock import patch

import pytest

litestar = pytest.importorskip('litestar', reason='litestar not installed')
from sqlalchemy import (  # noqa: E402
    text as sa_text,
    create_engine,
)
from sqlalchemy.orm import sessionmaker  # noqa: E402
from sqlalchemy.pool import StaticPool  # noqa: E402
from litestar.testing import TestClient  # noqa: E402

from omnipath_utils.server._app import create_app  # noqa: E402


@pytest.fixture
def app():
    """Create test app with an in-memory SQLite session."""

    engine = create_engine(
        'sqlite://',
        echo=False,
        connect_args={'check_same_thread': False},
        poolclass=StaticPool,
    )

    # Create tables that the health endpoint queries
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
    """Create test client."""
    return TestClient(app)


class TestHealthEndpoint:
    def test_health(self, client):
        resp = client.get('/health')
        assert resp.status_code == 200
        data = resp.json()
        assert data['status'] == 'ok'
        assert 'database' in data
        assert 'stats' in data
        assert 'backends' in data  # empty list when schema unavailable
        assert 'builds' in data

    def test_health_id_mapping_ftp_graceful(self, client):
        # The id_mapping_ftp stat uses Postgres-only reltuples/to_regclass; on
        # the SQLite test DB the guard must skip it without breaking /health.
        resp = client.get('/health')
        assert resp.status_code == 200
        stats = resp.json()['stats']
        assert 'id_mapping' in stats          # curated table is counted
        assert 'id_mapping_ftp' not in stats  # absent -> gracefully skipped


class TestTaxonomyEndpoints:
    def test_resolve_human(self, client):
        resp = client.get('/taxonomy/resolve', params={'organism': 'human'})
        assert resp.status_code == 200
        data = resp.json()
        assert data['ncbi_tax_id'] == 9606
        assert data['common_name'] == 'human'
        assert data['kegg_code'] == 'hsa'

    def test_resolve_by_taxid(self, client):
        resp = client.get('/taxonomy/resolve', params={'organism': '10090'})
        assert resp.status_code == 200
        assert resp.json()['common_name'] == 'mouse'

    def test_resolve_unknown(self, client):
        resp = client.get('/taxonomy/resolve', params={'organism': 'alien'})
        assert resp.status_code == 200
        assert resp.json()['ncbi_tax_id'] is None

    def test_organisms_list(self, client):
        resp = client.get('/taxonomy/organisms')
        assert resp.status_code == 200
        data = resp.json()
        assert len(data) >= 20
        taxids = {o['ncbi_tax_id'] for o in data}
        assert 9606 in taxids


class TestMappingEndpoints:
    def test_id_types(self, client):
        resp = client.get('/mapping/id-types')
        assert resp.status_code == 200
        data = resp.json()
        assert len(data) > 50
        names = {t['name'] for t in data}
        assert 'uniprot' in names

    @patch('omnipath_utils.server._routes_mapping.translate_ids')
    def test_translate_get(self, mock_translate, client):
        mock_translate.return_value = ({'TP53': {'P04637'}}, {'uniprot'})
        resp = client.get(
            '/mapping/translate',
            params={
                'identifiers': 'TP53',
                'id_type': 'genesymbol',
                'target_id_type': 'uniprot',
            },
        )
        assert resp.status_code == 200
        assert 'results' in resp.json()

    @patch('omnipath_utils.server._routes_mapping.translate_ids')
    def test_translate_post(self, mock_translate, client):
        mock_translate.return_value = ({'TP53': {'P04637'}}, {'uniprot'})
        resp = client.post(
            '/mapping/translate',
            json={
                'identifiers': ['TP53'],
                'id_type': 'genesymbol',
                'target_id_type': 'uniprot',
            },
        )
        assert resp.status_code in (200, 201)


class TestReflistEndpoints:
    def test_list_names(self, client):
        resp = client.get('/reflists/list-names')
        assert resp.status_code == 200
        assert 'swissprot' in resp.json()




class TestIdentifyEndpoint:
    @patch('omnipath_utils.server._routes_mapping.identify_ids')
    def test_identify(self, mock_identify, client):
        mock_identify.return_value = {
            'P04637': [
                {'id_type': 'uniprot', 'role': 'source', 'mappings_count': 5},
            ],
        }
        resp = client.get(
            '/mapping/identify',
            params={'identifiers': 'P04637'},
        )
        assert resp.status_code == 200
        data = resp.json()
        assert 'results' in data
        assert 'P04637' in data['results']
        assert data['meta']['total_input'] == 1

    @patch('omnipath_utils.server._routes_mapping.identify_ids')
    def test_identify_multiple(self, mock_identify, client):
        mock_identify.return_value = {
            'P04637': [
                {'id_type': 'uniprot', 'role': 'source', 'mappings_count': 5},
            ],
            'HMDB0000001': [
                {'id_type': 'hmdb', 'role': 'source', 'mappings_count': 3},
            ],
        }
        resp = client.get(
            '/mapping/identify',
            params={'identifiers': 'P04637,HMDB0000001'},
        )
        assert resp.status_code == 200
        data = resp.json()
        assert data['meta']['total_input'] == 2

    @patch('omnipath_utils.server._routes_mapping.identify_ids')
    def test_identify_custom_taxid(self, mock_identify, client):
        mock_identify.return_value = {'Trp53': []}
        resp = client.get(
            '/mapping/identify',
            params={'identifiers': 'Trp53', 'ncbi_tax_id': 10090},
        )
        assert resp.status_code == 200
        assert resp.json()['meta']['ncbi_tax_id'] == 10090


class TestAllMappingsEndpoint:
    @patch('omnipath_utils.server._routes_mapping.get_all_mappings')
    def test_all_mappings(self, mock_all, client):
        mock_all.return_value = {
            'P04637': {'genesymbol': ['TP53'], 'entrez': ['7157']},
        }
        resp = client.get(
            '/mapping/all',
            params={
                'identifiers': 'P04637',
                'id_type': 'uniprot',
            },
        )
        assert resp.status_code == 200
        data = resp.json()
        assert 'results' in data
        assert 'P04637' in data['results']
        assert 'genesymbol' in data['results']['P04637']
        assert data['meta']['id_type'] == 'uniprot'
        assert data['meta']['total_input'] == 1

    @patch('omnipath_utils.server._routes_mapping.get_all_mappings')
    def test_all_mappings_alias_resolution(self, mock_all, client):
        mock_all.return_value = {}
        resp = client.get(
            '/mapping/all',
            params={
                'identifiers': 'P04637',
                'id_type': 'gene_symbol',
            },
        )
        assert resp.status_code == 200


class TestRESTParams:
    @patch('omnipath_utils.server._routes_mapping.translate_ids')
    def test_raw_parameter(self, mock_translate, client):
        """Test that raw=true is accepted and reflected in meta."""
        mock_translate.return_value = ({'TP53': {'P04637'}}, {'uniprot'})
        resp = client.get(
            '/mapping/translate',
            params={
                'identifiers': 'TP53',
                'id_type': 'genesymbol',
                'target_id_type': 'uniprot',
                'raw': 'true',
            },
        )
        assert resp.status_code == 200
        data = resp.json()
        assert data['meta']['raw'] is True

    @patch('omnipath_utils.server._routes_mapping.translate_ids')
    def test_backend_parameter(self, mock_translate, client):
        """Test that backend parameter is accepted and reflected in meta."""
        mock_translate.return_value = ({'TP53': {'P04637'}}, {'uniprot'})
        resp = client.get(
            '/mapping/translate',
            params={
                'identifiers': 'TP53',
                'id_type': 'genesymbol',
                'target_id_type': 'uniprot',
                'backend': 'biomart',
            },
        )
        assert resp.status_code == 200
        data = resp.json()
        assert data['meta']['backend'] == ['uniprot']  # reports actual DB backends used

    @patch('omnipath_utils.server._routes_mapping.translate_ids')
    def test_raw_skips_fallbacks(self, mock_translate, client):
        """When raw=true, fallbacks should not be applied."""

        # Only return result for exact match; lowercase should miss
        def side_effect(session, ids, src, tgt, tax):
            return ({i: {'P04637'} for i in ids if i == 'TP53'}, {'uniprot'})

        mock_translate.side_effect = side_effect

        resp = client.get(
            '/mapping/translate',
            params={
                'identifiers': 'tp53',
                'id_type': 'genesymbol',
                'target_id_type': 'uniprot',
                'raw': 'true',
            },
        )
        assert resp.status_code == 200
        data = resp.json()
        # In raw mode, lowercase "tp53" should remain unmapped
        assert 'tp53' in data['unmapped']

    @patch('omnipath_utils.server._routes_mapping.translate_ids')
    def test_non_raw_applies_fallbacks(self, mock_translate, client):
        """When raw=false (default), fallbacks should be applied."""
        call_count = {'n': 0}

        def side_effect(session, ids, src, tgt, tax):
            call_count['n'] += 1
            result = {}
            for i in ids:
                if i == 'TP53':
                    result[i] = {'P04637'}
            return (result, {'uniprot'})

        mock_translate.side_effect = side_effect

        resp = client.get(
            '/mapping/translate',
            params={
                'identifiers': 'tp53',
                'id_type': 'genesymbol',
                'target_id_type': 'uniprot',
            },
        )
        assert resp.status_code == 200
        # Should have called translate_ids multiple times (fallbacks)
        assert call_count['n'] > 1
