"""spec 011 T078-T080 -- the chemistry toolkit binding (an optional extra).

``rdkit`` is installed via ``uv run --with rdkit`` when unavailable locally,
same as any other test dependency -- these are not gated on a special env
var, they just skip individual assertions when rdkit genuinely cannot be
imported (matching the capability's own "degrades in a declared way",
research R14).
"""

from __future__ import annotations

import pytest

from omnipath_utils.mapping import _chemistry as chem

rdkit = pytest.importorskip('rdkit', reason='rdkit not installed')


def test_chemistry_available_reports_true_when_rdkit_importable():
    assert chem.chemistry_available() is True


def test_compute_inchikey_from_inchi():
    # Ethanol -- contracts/translation-api.md's own worked example.
    inchi = 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
    assert chem.compute_inchikey(inchi=inchi) == 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'


def test_compute_inchikey_from_smiles():
    assert chem.compute_inchikey(smiles='CCO') == 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'


def test_compute_inchikey_prefers_inchi_when_both_given():
    # A deliberately mismatched smiles must not override a valid inchi.
    inchi = 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
    assert chem.compute_inchikey(inchi=inchi, smiles='c1ccccc1') == (
        'LFQSCWFLJHTTHZ-UHFFFAOYSA-N'
    )


def test_compute_inchikey_returns_none_for_garbage():
    assert chem.compute_inchikey(inchi='not a real inchi') is None
    assert chem.compute_inchikey(smiles='not a real smiles') is None


def test_compute_inchikey_returns_none_for_no_input():
    assert chem.compute_inchikey() is None


def test_canonicalize_smiles_normalizes_equivalent_spellings():
    # C7: two equivalent SMILES spellings for ethanol must canonicalize
    # identically.
    a = chem.canonicalize_smiles('CCO')
    b = chem.canonicalize_smiles('OCC')
    assert a == b
    assert a is not None


def test_canonicalize_smiles_returns_none_for_garbage():
    assert chem.canonicalize_smiles('not a real smiles') is None


class TestDegradedWithoutRdkit:
    """chemistry_available() reports False and every function degrades to
    None instead of raising, when rdkit genuinely cannot be imported
    (research R14: "structure-to-anything routes are unavailable and say
    so; identifier-to-identifier translation is unaffected").
    """

    def test_functions_return_none_when_unavailable(self, monkeypatch):
        monkeypatch.setattr(chem, 'chemistry_available', lambda: False)
        assert chem.compute_inchikey(inchi='InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3') is None
        assert chem.canonicalize_smiles('CCO') is None
