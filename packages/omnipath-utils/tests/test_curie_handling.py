"""CURIE normalisation: all mapping lookups strip/normalise CURIE prefixes.

``chebi:17612``, ``CHEBI:17612`` and ``17612`` must all resolve identically; the
response echoes the caller's original input.
"""

from collections import defaultdict
from unittest.mock import MagicMock

import pytest

import omnipath_utils.db._query as q


class TestStripCurie:
    def test_chebi_banana_forms_collapse(self):
        for v in ('chebi:17612', 'CHEBI:17612', 'Chebi:17612', '17612'):
            assert q.strip_curie('chebi', v) == 'CHEBI:17612'

    def test_hmdb_prefix_stripped(self):
        assert q.strip_curie('hmdb', 'hmdb:HMDB0000094') == 'HMDB0000094'
        assert q.strip_curie('hmdb', 'HMDB0000094') == 'HMDB0000094'

    def test_kegg_prefix_stripped(self):
        assert q.strip_curie('kegg', 'kegg:C00031') == 'C00031'
        assert q.strip_curie('kegg', 'kegg.compound:C00031') == 'C00031'
        assert q.strip_curie('kegg', 'C00031') == 'C00031'

    def test_pubchem_prefix_stripped(self):
        assert q.strip_curie('pubchem', 'pubchem.compound:5793') == '5793'
        assert q.strip_curie('pubchem', 'pubchem:5793') == '5793'

    def test_unknown_prefix_left_intact(self):
        # a colon that is not a known prefix for this type is preserved
        assert q.strip_curie('genesymbol', 'foo:bar') == 'foo:bar'


class TestLookupKey:
    def test_name_lowercased(self):
        assert q._lookup_key('name', '  Taurine ') == 'taurine'

    def test_structure_verbatim(self):
        inchi = 'InChI=1S/H2O/h1H2'
        assert q._lookup_key('inchi', inchi) == inchi

    def test_chebi_curie_normalised(self):
        assert q._lookup_key('chebi', 'chebi:17612') == 'CHEBI:17612'

    def test_inchikey_prefixed_lowercase_normalised(self):
        # C8: a prefixed, lower-case structure key translates. InChIKey's own
        # convention uses '=' (InChIKey=...), which strip_curie's ':'-only
        # split never catches -- this must go through the real normalizer.
        key = 'lfqscwfljhtthz-uhfffaoysa-n'
        assert q._lookup_key('inchikey', f'InChIKey={key}') == key.upper()
        assert q._lookup_key('inchikey', f'inchikey={key}') == key.upper()

    def test_inchikey_bare_lowercase_uppercased(self):
        assert q._lookup_key(
            'inchikey', 'lfqscwfljhtthz-uhfffaoysa-n',
        ) == 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'

    def test_inchikey_bare_uppercase_unchanged(self):
        assert q._lookup_key(
            'inchikey', 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N',
        ) == 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'

    def test_smiles_canonicalized_when_chemistry_available(self, monkeypatch):
        # C7: two equivalent SMILES spellings must produce the same key.
        from omnipath_utils.mapping import _chemistry

        monkeypatch.setattr(_chemistry, 'chemistry_available', lambda: True)
        monkeypatch.setattr(
            _chemistry, 'canonicalize_smiles', lambda s: 'CANONICAL',
        )
        assert q._lookup_key('smiles', 'OCC') == 'CANONICAL'
        assert q._lookup_key('smiles', 'CCO') == 'CANONICAL'

    def test_smiles_falls_back_to_verbatim_without_chemistry(self, monkeypatch):
        from omnipath_utils.mapping import _chemistry

        monkeypatch.setattr(_chemistry, 'chemistry_available', lambda: False)
        assert q._lookup_key('smiles', '  CCO  ') == 'CCO'


class TestTranslateCurie:
    def test_chebi_curie_variants_resolve_and_rekey(self, monkeypatch):
        # the DB stores CHEBI:17612; a fake table returns a hit for that key.
        def fake_query_table(session, table, ids, src, tgt, tax):
            res = defaultdict(set)
            if ids and 'CHEBI:17612' in ids:
                res['CHEBI:17612'].add('HMDB0000094')
            return res, {'chebi'}

        monkeypatch.setattr(q, '_query_table', fake_query_table)
        monkeypatch.setattr(q, '_ftp_types', lambda s: frozenset())
        monkeypatch.setattr(q, '_ftp_table_exists', lambda s: False)

        for variant in ('chebi:17612', 'CHEBI:17612', '17612'):
            out, _ = q.translate_ids(MagicMock(), [variant], 'chebi', 'hmdb', 0)
            # response echoes the caller's original input form
            assert out == {variant: {'HMDB0000094'}}
