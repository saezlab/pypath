"""spec 011 T076/T077/T082 -- the /mapping/translate HTTP layer for
structure-to-anything translation: a route reports its
``structure_key_computation`` capability, unavailable with a reason when
the chemistry toolkit is absent.

Mocks ``translate_ids`` (same pattern as ``test_name_translation_api.py``)
so these run without a live DB or rdkit.
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


class TestStructureCapabilityReported:
    @patch('omnipath_utils.server._routes_mapping.translate_ids')
    def test_smiles_query_reports_capability_available(self, mock_translate, client):
        mock_translate.return_value = ({'CCO': {'CHEBI:1'}}, {'chebi'})
        with patch(
            'omnipath_utils.server._routes_mapping.chemistry_available',
            return_value=True,
        ):
            resp = client.get(
                '/mapping/translate',
                params={
                    'identifiers': 'CCO', 'id_type': 'smiles',
                    'target_id_type': 'chebi',
                },
            )
        cap = resp.json()['meta']['structure_key_computation']
        assert cap == {'available': True, 'reason': None}

    @patch('omnipath_utils.server._routes_mapping.translate_ids')
    def test_smiles_query_reports_unavailable_with_reason(
        self, mock_translate, client,
    ):
        # T077: the route reports itself unavailable, with a reason, when
        # the chemistry extra is absent -- not a silent empty result.
        mock_translate.return_value = ({}, set())
        with patch(
            'omnipath_utils.server._routes_mapping.chemistry_available',
            return_value=False,
        ):
            resp = client.get(
                '/mapping/translate',
                params={
                    'identifiers': 'CCO', 'id_type': 'smiles',
                    'target_id_type': 'chebi', 'raw': 'true',
                },
            )
        cap = resp.json()['meta']['structure_key_computation']
        assert cap['available'] is False
        assert cap['reason']

    @patch('omnipath_utils.server._routes_mapping.translate_ids')
    def test_inchi_to_inchikey_reports_capability(self, mock_translate, client):
        # source in (inchi, smiles) with target inchikey also depends on
        # the toolkit (the live-compute fallback, C6).
        mock_translate.return_value = ({}, set())
        with patch(
            'omnipath_utils.server._routes_mapping.chemistry_available',
            return_value=False,
        ):
            resp = client.get(
                '/mapping/translate',
                params={
                    'identifiers': 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',
                    'id_type': 'inchi', 'target_id_type': 'inchikey',
                    'raw': 'true',
                },
            )
        assert resp.json()['meta']['structure_key_computation']['available'] is False

    @patch('omnipath_utils.server._routes_mapping.translate_ids')
    def test_capability_absent_for_a_non_structure_query(
        self, mock_translate, client,
    ):
        mock_translate.return_value = ({'TP53': {'P04637'}}, {'uniprot'})
        resp = client.get(
            '/mapping/translate',
            params={
                'identifiers': 'TP53', 'id_type': 'genesymbol',
                'target_id_type': 'uniprot',
            },
        )
        assert 'structure_key_computation' not in resp.json()['meta']

    @patch('omnipath_utils.server._routes_mapping.translate_ids')
    def test_inchikey_lookup_unaffected_by_missing_toolkit(
        self, mock_translate, client,
    ):
        # Plain inchikey<->id translation is stored-data lookup, no live
        # rdkit needed -- capability is not reported for it (research R14:
        # "identifier-to-identifier translation is unaffected").
        mock_translate.return_value = ({'CHEBI:1': {'XXXXXXXXXXXXXX-UHFFFAOYSA-N'}}, {'chebi'})
        resp = client.get(
            '/mapping/translate',
            params={
                'identifiers': 'CHEBI:1', 'id_type': 'chebi',
                'target_id_type': 'inchikey',
            },
        )
        assert 'structure_key_computation' not in resp.json()['meta']
