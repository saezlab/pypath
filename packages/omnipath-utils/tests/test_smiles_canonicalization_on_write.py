"""spec 011 T079 -- canonical SMILES stored beside the source form.

``_canonicalize_smiles_in_data`` runs on the plain ``{source: {targets}}``
dict *before* ``_long_rows`` projects it into COPY rows -- ``_long_rows``
itself stays byte-verbatim for structures (test_long_copy_helper.py's own
documented, deliberate contract: "'C' != 'c': folding would corrupt it").
Canonicalization is chemistry-aware (rdkit), not naive folding, so it does
not corrupt that same aromaticity/stereochemistry information -- it is a
different transform applied at a different layer.
"""

from omnipath_utils.db._build import DatabaseBuilder

_canon = DatabaseBuilder._canonicalize_smiles_in_data


class TestSourceSideSmiles:
    def test_two_spellings_of_one_molecule_collapse_to_one_key(self):
        # ethanol: two equivalent SMILES spellings.
        data = {'CCO': {'CHEBI:1'}, 'OCC': {'CHEBI:1'}}
        out = _canon(data, 'smiles', 'chebi')
        assert len(out) == 1

    def test_garbage_keeps_its_raw_form(self):
        data = {'not-a-real-smiles': {'CHEBI:1'}}
        out = _canon(data, 'smiles', 'chebi')
        assert out == {'not-a-real-smiles': {'CHEBI:1'}}


class TestTargetSideSmiles:
    def test_two_spellings_of_one_molecule_collapse_to_one_target(self):
        data = {'CHEBI:1': {'CCO'}, 'CHEBI:2': {'OCC'}}
        out = _canon(data, 'chebi', 'smiles')
        targets = out['CHEBI:1'] | out['CHEBI:2']
        assert len(targets) == 1


class TestNonSmilesUnaffected:
    def test_neither_side_smiles_passes_through_unchanged(self):
        data = {'CHEBI:1': {'HMDB0000094'}}
        out = _canon(data, 'chebi', 'hmdb')
        assert out == data

    def test_inchi_side_unaffected_even_when_smiles_is_the_other_side(self):
        inchi = 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
        data = {inchi: {'CCO'}}
        out = _canon(data, 'inchi', 'smiles')
        assert set(out.keys()) == {inchi}  # inchi key untouched


class TestDegradedWithoutChemistry:
    def test_passes_through_unchanged_when_unavailable(self, monkeypatch):
        from omnipath_utils.mapping import _chemistry

        monkeypatch.setattr(_chemistry, 'chemistry_available', lambda: False)
        data = {'CCO': {'CHEBI:1'}, 'OCC': {'CHEBI:1'}}
        out = _canon(data, 'smiles', 'chebi')
        assert out == data
