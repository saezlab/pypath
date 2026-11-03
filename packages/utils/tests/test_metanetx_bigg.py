"""Tests for MetaNetX, BiGG mapping backends and HMDB normalisation."""

from unittest.mock import patch, MagicMock

from omnipath_utils.mapping.backends._metanetx import MetaNetXBackend
from omnipath_utils.mapping.backends._bigg import BiggBackend
from omnipath_utils.mapping._special import normalise_hmdb


# ---- MetaNetX backend tests ----


class TestMetaNetXBackend:
    """Tests for the MetaNetX mapping backend."""

    def test_backend_registered(self):
        """MetaNetX backend is discoverable."""
        from omnipath_utils.mapping.backends import get_backend

        backend = get_backend("metanetx")
        assert backend is not None
        assert backend.name == "metanetx"

    def test_read_bigg_to_chebi(self):
        """Mock pypath and verify read returns dict[str, set[str]]."""
        with patch("pypath.inputs.metanetx.metanetx_mapping") as mock:
            mock.return_value = {
                "atp": {"CHEBI:30616"},
                "gtp": {"CHEBI:15996"},
            }
            backend = MetaNetXBackend()
            result = backend.read("bigg", "chebi", 0)
            assert result == {
                "atp": {"CHEBI:30616"},
                "gtp": {"CHEBI:15996"},
            }
            mock.assert_called_once_with("bigg", "chebi")

    def test_read_empty_when_no_data(self):
        """Returns empty dict when metanetx_mapping returns nothing."""
        with patch("pypath.inputs.metanetx.metanetx_mapping") as mock:
            mock.return_value = {}
            backend = MetaNetXBackend()
            result = backend.read("bigg", "chebi", 0)
            assert result == {}

    def test_fallback_on_import_error(self):
        """Returns empty dict when pypath not available."""
        with patch.object(
            MetaNetXBackend, "_read_via_pypath", side_effect=ImportError
        ):
            backend = MetaNetXBackend()
            result = backend.read("bigg", "chebi", 0)
            assert result == {}

    def test_read_handles_exception(self):
        """Returns empty dict when metanetx_mapping raises."""
        with patch("pypath.inputs.metanetx.metanetx_mapping") as mock:
            mock.side_effect = ValueError("download failed")
            backend = MetaNetXBackend()
            result = backend.read("bigg", "chebi", 0)
            assert result == {}


# ---- BiGG backend tests ----


class TestBiggBackend:
    """Tests for the BiGG mapping backend."""

    def test_backend_registered(self):
        """BiGG backend is discoverable."""
        from omnipath_utils.mapping.backends import get_backend

        backend = get_backend("bigg")
        assert backend is not None
        assert backend.name == "bigg"

    def test_read_bigg_to_chebi(self):
        """bigg->chebi reads forward mapping."""
        with patch("pypath.inputs.bigg.bigg_metabolite_mapping") as mock:
            mock.return_value = {
                "atp": {"CHEBI:30616", "CHEBI:15422"},
                "gtp": {"CHEBI:15996"},
            }
            backend = BiggBackend()
            result = backend.read("bigg", "chebi", 0)
            assert "atp" in result
            assert "CHEBI:30616" in result["atp"]
            mock.assert_called_once_with("chebi")

    def test_read_chebi_to_bigg_reverses(self):
        """chebi->bigg should reverse the mapping."""
        with patch("pypath.inputs.bigg.bigg_metabolite_mapping") as mock:
            mock.return_value = {"atp": {"CHEBI:30616"}}
            backend = BiggBackend()
            result = backend.read("chebi", "bigg", 0)
            assert "CHEBI:30616" in result
            assert "atp" in result["CHEBI:30616"]

    def test_unsupported_pair_returns_empty(self):
        """Neither source nor target is bigg -> empty dict."""
        backend = BiggBackend()
        result = backend.read("hmdb", "kegg", 0)
        assert result == {}

    def test_read_empty_when_no_data(self):
        """Returns empty dict when pypath returns nothing."""
        with patch("pypath.inputs.bigg.bigg_metabolite_mapping") as mock:
            mock.return_value = {}
            backend = BiggBackend()
            result = backend.read("bigg", "chebi", 0)
            assert result == {}

    def test_read_handles_exception(self):
        """Returns empty dict when bigg_metabolite_mapping raises."""
        with patch("pypath.inputs.bigg.bigg_metabolite_mapping") as mock:
            mock.side_effect = ValueError("download failed")
            backend = BiggBackend()
            result = backend.read("bigg", "chebi", 0)
            assert result == {}

    def test_reverse_multi_target(self):
        """Reverse mapping with multiple targets per bigg ID."""
        with patch("pypath.inputs.bigg.bigg_metabolite_mapping") as mock:
            mock.return_value = {
                "atp": {"CHEBI:30616", "CHEBI:15422"},
            }
            backend = BiggBackend()
            result = backend.read("chebi", "bigg", 0)
            assert "CHEBI:30616" in result
            assert "CHEBI:15422" in result
            assert result["CHEBI:30616"] == {"atp"}
            assert result["CHEBI:15422"] == {"atp"}


# ---- HMDB normalisation tests ----


class TestNormaliseHmdb:
    """Tests for HMDB ID normalisation."""

    def test_old_5digit_format(self):
        """Old 5-digit HMDB IDs are zero-padded to 7 digits."""
        assert normalise_hmdb("HMDB00001") == "HMDB0000001"

    def test_new_7digit_unchanged(self):
        """Already correct 7-digit HMDB IDs pass through."""
        assert normalise_hmdb("HMDB0000001") == "HMDB0000001"

    def test_non_hmdb_unchanged(self):
        """Non-HMDB strings are returned as-is."""
        assert normalise_hmdb("CHEBI:12345") == "CHEBI:12345"

    def test_longer_digits_unchanged(self):
        """7-digit IDs with larger numbers pass through."""
        assert normalise_hmdb("HMDB0062694") == "HMDB0062694"

    def test_case_insensitive(self):
        """Lowercase prefix is handled; output uses hardcoded HMDB prefix."""
        result = normalise_hmdb("hmdb00001")
        assert result == "HMDB0000001"

    def test_empty_string(self):
        """Empty string returns empty string."""
        assert normalise_hmdb("") == ""

    def test_old_4digit_format(self):
        """Very short old IDs are also zero-padded."""
        assert normalise_hmdb("HMDB0001") == "HMDB0000001"

    def test_six_digit_format(self):
        """6-digit IDs are zero-padded to 7."""
        assert normalise_hmdb("HMDB000001") == "HMDB0000001"


# ---- HMDB normalisation in translate_core ----


class TestHmdbNormalisationInTranslate:
    """Verify HMDB normalisation is applied in translate_core."""

    @patch("omnipath_utils.mapping._mapper.Mapper._load_table")
    def test_old_hmdb_format_normalised_in_translate(self, mock_load):
        """translate_core normalises old HMDB IDs before lookup."""
        from omnipath_utils.mapping._table import MappingTable
        from omnipath_utils.mapping._mapper import Mapper
        from omnipath_utils.mapping._translate import translate_core

        mock_load.return_value = None
        Mapper._instance = None
        mapper = Mapper()

        table = MappingTable(
            data={"HMDB0000001": {"CHEBI:18367"}},
            id_type="hmdb",
            target_id_type="chebi",
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        # Pass old-format ID: should be normalised before lookup
        result = translate_core(
            ["HMDB00001"], "hmdb", "chebi", 9606, raw=True,
        )
        assert "HMDB0000001" in result
        assert result["HMDB0000001"] == {"CHEBI:18367"}

    @patch("omnipath_utils.mapping._mapper.Mapper._load_table")
    def test_new_hmdb_format_works_in_translate(self, mock_load):
        """translate_core passes through already-normalised HMDB IDs."""
        from omnipath_utils.mapping._table import MappingTable
        from omnipath_utils.mapping._mapper import Mapper
        from omnipath_utils.mapping._translate import translate_core

        mock_load.return_value = None
        Mapper._instance = None
        mapper = Mapper()

        table = MappingTable(
            data={"HMDB0000001": {"CHEBI:18367"}},
            id_type="hmdb",
            target_id_type="chebi",
            ncbi_tax_id=9606,
        )
        mapper.tables[table.key] = table
        Mapper._instance = mapper

        result = translate_core(
            ["HMDB0000001"], "hmdb", "chebi", 9606, raw=True,
        )
        assert "HMDB0000001" in result
        assert result["HMDB0000001"] == {"CHEBI:18367"}
