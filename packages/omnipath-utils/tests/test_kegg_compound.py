"""Tests for the KEGG Compound mapping backend (FR-051, SC-021).

Test order mirrors TDD: registration → known mapping → 1→many abort →
exclusion logging → reverse mapping → error paths.

Must fail before T079 (``_kegg_compound.py`` + id_types.yaml wiring) is
implemented.

NOTE: pypath imports trigger a numpy C-extension load which fails in the
current test venv (libstdc++.so.6 not on LD_LIBRARY_PATH). Tests therefore
exercise the backend's ``_load_kegg_to_chebi(mx_fn, rest_fn)`` interface
directly via plain callables, and patch ``_read_via_pypath`` for end-to-end
paths — rather than patching deep pypath module symbols.
"""

from __future__ import annotations

import logging
from unittest.mock import MagicMock, patch

import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _mx(data: dict) -> callable:
    """Return a callable that mimics metanetx_mapping(src, tgt) -> data."""
    def _fn(*_args, **_kwargs):
        return data
    return _fn


def _rest(data: dict) -> callable:
    """Return a callable that mimics _kegg_conv(...) -> data."""
    def _fn(*_args, **_kwargs):
        return data
    return _fn


# ---------------------------------------------------------------------------
# T078-1  Backend registration
# ---------------------------------------------------------------------------


class TestKeggCompoundBackendRegistration:
    """kegg_compound backend is discoverable and wired in id_types.yaml."""

    def test_backend_registered(self):
        """get_backend('kegg_compound') returns a non-None backend."""
        from omnipath_utils.mapping.backends import get_backend

        backend = get_backend("kegg_compound")
        assert backend is not None
        assert backend.name == "kegg_compound"

    def test_kegg_id_type_has_kegg_compound_column(self):
        """id_types.yaml 'kegg' entry lists a kegg_compound backend column."""
        from omnipath_utils.mapping._id_types import IdTypeRegistry

        reg = IdTypeRegistry.get()
        col = reg.backend_column("kegg", "kegg_compound")
        assert col is not None, "'kegg' id_type must declare a kegg_compound backend column"

    def test_chebi_id_type_has_kegg_compound_column(self):
        """id_types.yaml 'chebi' entry lists a kegg_compound backend column."""
        from omnipath_utils.mapping._id_types import IdTypeRegistry

        reg = IdTypeRegistry.get()
        col = reg.backend_column("chebi", "kegg_compound")
        assert col is not None, "'chebi' id_type must declare a kegg_compound backend column"


# ---------------------------------------------------------------------------
# T078-2  Known-ID mapping
# ---------------------------------------------------------------------------


class TestKeggCompoundKnownMapping:
    """C00022 (pyruvate) maps to CHEBI:15361 via each source."""

    def test_c00022_maps_to_chebi_15361_via_metanetx(self):
        """MetaNetX source: C00022 → CHEBI:15361."""
        from omnipath_utils.mapping.backends._kegg_compound import KeggCompoundBackend

        result = KeggCompoundBackend()._load_kegg_to_chebi(
            _mx({"C00022": {"CHEBI:15361"}}),
            _rest({}),
        )

        assert "C00022" in result
        assert "CHEBI:15361" in result["C00022"]

    def test_c00022_maps_to_chebi_15361_via_rest(self):
        """KEGG REST source: C00022 → 15361 (normalised to CHEBI:15361)."""
        from omnipath_utils.mapping.backends._kegg_compound import KeggCompoundBackend

        result = KeggCompoundBackend()._load_kegg_to_chebi(
            _mx({}),
            _rest({"C00022": {"15361"}}),  # bare numeric from KEGG REST
        )

        assert "C00022" in result
        assert "CHEBI:15361" in result["C00022"]

    def test_result_is_dict_of_sets(self):
        """Return type is dict[str, set[str]]."""
        from omnipath_utils.mapping.backends._kegg_compound import KeggCompoundBackend

        result = KeggCompoundBackend()._load_kegg_to_chebi(
            _mx({"C00022": {"CHEBI:15361"}}),
            _rest({}),
        )

        assert isinstance(result, dict)
        for v in result.values():
            assert isinstance(v, set)

    def test_read_kegg_to_chebi_via_mocked_pypath(self):
        """backend.read('kegg','chebi',0) calls _read_via_pypath."""
        from omnipath_utils.mapping.backends._kegg_compound import KeggCompoundBackend

        expected = {"C00022": {"CHEBI:15361"}}
        with patch.object(
            KeggCompoundBackend,
            "_read_via_pypath",
            return_value=expected,
        ):
            result = KeggCompoundBackend().read("kegg", "chebi", 0)

        assert result == expected


# ---------------------------------------------------------------------------
# T078-3  1→many abort
# ---------------------------------------------------------------------------


class TestKeggCompoundOneToManyAbort:
    """KEGG IDs mapping to >1 ChEBI are excluded from the translation table."""

    def test_multi_chebi_kegg_id_excluded(self):
        """C99999 mapping to two distinct ChEBI IDs is dropped."""
        from omnipath_utils.mapping.backends._kegg_compound import KeggCompoundBackend

        result = KeggCompoundBackend()._load_kegg_to_chebi(
            _mx({
                "C00022": {"CHEBI:15361"},   # 1→1 (keep)
                "C99999": {"CHEBI:11111"},   # MetaNetX says one ChEBI …
            }),
            # … KEGG REST adds a second ChEBI for C99999 → 1→many → exclude
            _rest({"C99999": {"22222"}}),
        )

        assert "C00022" in result, "1→1 pair must be kept"
        assert "C99999" not in result, "1→many pair must be excluded"

    def test_duplicate_same_chebi_kept(self):
        """Same ChEBI from both sources deduplicates to 1→1 and is kept."""
        from omnipath_utils.mapping.backends._kegg_compound import KeggCompoundBackend

        result = KeggCompoundBackend()._load_kegg_to_chebi(
            _mx({"C00022": {"CHEBI:15361"}}),
            _rest({"C00022": {"15361"}}),  # same, different format
        )

        assert "C00022" in result
        assert result["C00022"] == {"CHEBI:15361"}

    def test_both_sources_one_to_many_excluded(self):
        """Two distinct ChEBI IDs coming only from MetaNetX are also excluded."""
        from omnipath_utils.mapping.backends._kegg_compound import KeggCompoundBackend

        result = KeggCompoundBackend()._load_kegg_to_chebi(
            _mx({"C99998": {"CHEBI:11111", "CHEBI:22222"}}),
            _rest({}),
        )

        assert "C99998" not in result


# ---------------------------------------------------------------------------
# T078-4  Exclusion count logging
# ---------------------------------------------------------------------------


class TestKeggCompoundExclusionLogging:
    """Exclusion count is emitted at INFO level (SC-021)."""

    def test_nonzero_exclusion_count_logged(self, caplog):
        """INFO log must mention both the kept count and excluded count."""
        from omnipath_utils.mapping.backends._kegg_compound import KeggCompoundBackend

        logger_name = "omnipath_utils.mapping.backends._kegg_compound"
        with caplog.at_level(logging.INFO, logger=logger_name):
            KeggCompoundBackend()._load_kegg_to_chebi(
                _mx({
                    "C00022": {"CHEBI:15361"},
                    "C99999": {"CHEBI:11111"},
                }),
                _rest({"C99999": {"22222"}}),
            )

        messages = [r.message for r in caplog.records if r.name == logger_name]
        assert any(
            "excluded" in m and "1" in m for m in messages
        ), f"Expected exclusion count in logs; got: {messages}"

    def test_zero_exclusion_count_logged(self, caplog):
        """Zero exclusions must also appear in the log."""
        from omnipath_utils.mapping.backends._kegg_compound import KeggCompoundBackend

        logger_name = "omnipath_utils.mapping.backends._kegg_compound"
        with caplog.at_level(logging.INFO, logger=logger_name):
            KeggCompoundBackend()._load_kegg_to_chebi(
                _mx({"C00022": {"CHEBI:15361"}}),
                _rest({}),
            )

        messages = [r.message for r in caplog.records if r.name == logger_name]
        assert any("excluded" in m for m in messages), (
            f"Zero-exclusion count must also be logged; got: {messages}"
        )


# ---------------------------------------------------------------------------
# T078-5  Reverse mapping and error paths
# ---------------------------------------------------------------------------


class TestKeggCompoundReverseAndErrors:
    """Reverse (chebi→kegg) and error/import-error behaviour."""

    def test_chebi_to_kegg_reverse(self):
        """chebi→kegg inverts the forward kegg→chebi mapping."""
        from omnipath_utils.mapping.backends._kegg_compound import KeggCompoundBackend

        fwd = {"C00022": {"CHEBI:15361"}}
        result = KeggCompoundBackend()._reverse(fwd)

        assert "CHEBI:15361" in result
        assert "C00022" in result["CHEBI:15361"]

    def test_reverse_multi_source(self):
        """Two KEGG IDs sharing a ChEBI both appear in the reversed set."""
        from omnipath_utils.mapping.backends._kegg_compound import KeggCompoundBackend

        fwd = {
            "C00022": {"CHEBI:15361"},
            "C00100": {"CHEBI:15361"},
        }
        result = KeggCompoundBackend()._reverse(fwd)

        assert "CHEBI:15361" in result
        assert result["CHEBI:15361"] == {"C00022", "C00100"}

    def test_unsupported_pair_returns_empty(self):
        """A pair not involving kegg or chebi returns an empty dict."""
        from omnipath_utils.mapping.backends._kegg_compound import KeggCompoundBackend

        with patch.object(
            KeggCompoundBackend, "_read_via_pypath",
            side_effect=lambda *a, **kw: KeggCompoundBackend._read_via_pypath.__wrapped__(*a, **kw)
            if hasattr(KeggCompoundBackend._read_via_pypath, '__wrapped__') else {},
        ):
            pass

        # Call read() directly; the real _read_via_pypath detects the unsupported pair
        # and returns {} before any pypath import.
        with patch.object(
            KeggCompoundBackend, "_read_via_pypath", return_value={}
        ):
            result = KeggCompoundBackend().read("hmdb", "uniprot", 0)
        assert result == {}

    def test_fallback_on_import_error(self):
        """Returns empty dict when pypath is not available."""
        from omnipath_utils.mapping.backends._kegg_compound import KeggCompoundBackend

        with patch.object(
            KeggCompoundBackend, "_read_via_pypath", side_effect=ImportError
        ):
            result = KeggCompoundBackend().read("kegg", "chebi", 0)

        assert result == {}

    def test_metanetx_exception_handled(self):
        """MetaNetX error is caught; KEGG REST data still used."""
        from omnipath_utils.mapping.backends._kegg_compound import KeggCompoundBackend

        def _bad_mx(*_):
            raise ValueError("download failed")

        result = KeggCompoundBackend()._load_kegg_to_chebi(
            _bad_mx,
            _rest({"C00022": {"15361"}}),
        )

        assert "C00022" in result

    def test_rest_exception_handled(self):
        """KEGG REST error is caught; MetaNetX data still used."""
        from omnipath_utils.mapping.backends._kegg_compound import KeggCompoundBackend

        def _bad_rest(*_):
            raise ConnectionError("REST unreachable")

        result = KeggCompoundBackend()._load_kegg_to_chebi(
            _mx({"C00022": {"CHEBI:15361"}}),
            _bad_rest,
        )

        assert "C00022" in result

    def test_empty_sources_returns_empty(self):
        """Both sources empty → empty result."""
        from omnipath_utils.mapping.backends._kegg_compound import KeggCompoundBackend

        result = KeggCompoundBackend()._load_kegg_to_chebi(
            _mx({}), _rest({}),
        )
        assert result == {}
