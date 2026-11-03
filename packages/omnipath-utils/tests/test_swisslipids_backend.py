"""Unit test for T086's SwissLipids format-prefix fix -- no DB, no download.

The live end-to-end path (``populate_mapping('swisslipids', 'inchikey', ...)``)
is blocked in this sandbox by a pre-existing, unrelated third-party bug
(``cachedir/_open.py``'s ``UnboundLocalError: cannot access local variable
'ext'``, hit re-fetching SwissLipids' raw file -- see deferred-items.md, not
caused by this change). This test isolates the actual fix -- the
``_read_via_pypath`` override that strips the raw ``InChIKey=`` format prefix
-- by mocking the base class's read instead of hitting the network.
"""

from __future__ import annotations

from unittest.mock import patch

from omnipath_utils.mapping.backends._swisslipids import SwissLipidsBackend


def test_swisslipids_strips_the_inchikey_format_prefix():
    backend = SwissLipidsBackend()
    raw = {
        'SLM:000000003': {'InChIKey=LZKPPSAEINBHRP-KORIGIIASA-O'},
        'SLM:000000006': {'InChIKey=VFYHIOGWELBGIK-DLBZAZTESA-O'},
    }
    with patch(
        'omnipath_utils.mapping.backends._inputs_v2_adapter'
        '.InputsV2Backend._read_via_pypath',
        return_value=raw,
    ):
        result = backend._read_via_pypath(
            'swisslipids', 'inchikey', 0,
            src_col='Lipid ID', tgt_col='InChI key (pH7.3)',
        )
    assert result == {
        'SLM:000000003': {'LZKPPSAEINBHRP-KORIGIIASA-O'},
        'SLM:000000006': {'VFYHIOGWELBGIK-DLBZAZTESA-O'},
    }


def test_swisslipids_drops_placeholder_and_unparseable_keys():
    backend = SwissLipidsBackend()
    raw = {
        'SLM:000000001': {'InChIKey=none'},
        'SLM:000000002': {'InChIKey=LZKPPSAEINBHRP-KORIGIIASA-O', 'garbage'},
    }
    with patch(
        'omnipath_utils.mapping.backends._inputs_v2_adapter'
        '.InputsV2Backend._read_via_pypath',
        return_value=raw,
    ):
        result = backend._read_via_pypath(
            'swisslipids', 'inchikey', 0,
            src_col='Lipid ID', tgt_col='InChI key (pH7.3)',
        )
    assert 'SLM:000000001' not in result
    assert result == {'SLM:000000002': {'LZKPPSAEINBHRP-KORIGIIASA-O'}}


def test_swisslipids_leaves_non_inchikey_targets_untouched():
    """The prefix-strip only applies when target_id_type == 'inchikey'."""
    backend = SwissLipidsBackend()
    raw = {'SLM:000000003': {'LMSL01000001'}}
    with patch(
        'omnipath_utils.mapping.backends._inputs_v2_adapter'
        '.InputsV2Backend._read_via_pypath',
        return_value=raw,
    ):
        result = backend._read_via_pypath(
            'swisslipids', 'lipidmaps', 0,
            src_col='Lipid ID', tgt_col='LIPID MAPS ID',
        )
    assert result == raw
