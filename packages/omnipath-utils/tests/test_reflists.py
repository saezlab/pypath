"""Tests for reference lists module."""

from unittest.mock import patch

from omnipath_utils.reflists._manager import ReferenceListManager


class TestReferenceListManager:
    def setup_method(self):
        ReferenceListManager._instance = None

    def test_singleton(self):
        r1 = ReferenceListManager.get()
        r2 = ReferenceListManager.get()
        assert r1 is r2

    def test_unknown_list(self):
        mgr = ReferenceListManager()
        result = mgr.get_reflist('nonexistent')
        assert result == set()

    @patch.object(ReferenceListManager, '_load_swissprot')
    def test_swissprot_cached(self, mock_load):
        mock_load.return_value = {'P04637', 'P00533'}
        mgr = ReferenceListManager()

        result1 = mgr.get_reflist('swissprot', 9606)
        result2 = mgr.get_reflist('swissprot', 9606)

        assert result1 == {'P04637', 'P00533'}
        assert result1 is result2  # same object, cached
        mock_load.assert_called_once()

    @patch.object(ReferenceListManager, '_load_swissprot')
    def test_is_swissprot(self, mock_load):
        mock_load.return_value = {'P04637', 'P00533'}
        mgr = ReferenceListManager()

        assert mgr.is_swissprot('P04637')
        assert not mgr.is_swissprot('FAKEID')

    @patch.object(ReferenceListManager, '_load_swissprot')
    @patch.object(ReferenceListManager, '_load_trembl')
    def test_all_uniprot_combines(self, mock_trembl, mock_swiss):
        mock_swiss.return_value = {'P04637'}
        mock_trembl.return_value = {'A0A024R1R8'}
        mgr = ReferenceListManager()

        result = mgr.get_reflist('uniprot', 9606)
        assert result == {'P04637', 'A0A024R1R8'}

    @patch.object(ReferenceListManager, '_load_swissprot')
    def test_different_organisms_separate(self, mock_load):
        mock_load.side_effect = [{'P04637'}, {'Q9Z0V2'}]
        mgr = ReferenceListManager()

        human = mgr.get_reflist('swissprot', 9606)
        mouse = mgr.get_reflist('swissprot', 10090)

        assert human != mouse
        assert mock_load.call_count == 2
