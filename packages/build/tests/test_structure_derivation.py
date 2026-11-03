from omnipath_resolver.canonical import structures
from omnipath_resolver.canonical.match import votes_for
from omnipath_resolver.canonical.policy import CHEMICAL_POLICY
from omnipath_build.extract.observations import RawEntityObservation


def test_standard_structure_identity_preserves_stereo_salts_and_isotopes():
    ethanol = structures.derive("CCO")
    assert ethanol["inchikey"] == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    assert ethanol["inchi"].startswith("InChI=1S/")
    assert (
        structures.derive("C[C@H](O)C(=O)O")["inchikey"]
        != structures.derive("C[C@@H](O)C(=O)O")["inchikey"]
    )
    assert structures.derive("CCO.[Na+]")["inchikey"] != ethanol["inchikey"]
    assert structures.derive("[13CH3]CO")["inchikey"] != ethanol["inchikey"]
    assert structures.derive("*CC")["status"] == "generic_structure"
    assert structures.derive("not a smiles")["status"] == "invalid_smiles"


def test_explicit_key_wins_and_cached_derivation_retains_provenance(tmp_path, monkeypatch):
    monkeypatch.setenv("OMNIPATH_STRUCTURE_CACHE", str(tmp_path))
    obs = RawEntityObservation(
        "e", "chemical_entity", "name", "ethanol", identifiers=[{"ns": "smiles", "id": "CCO"}]
    )
    votes, observed = votes_for(obs, CHEMICAL_POLICY)
    assert any(v.ns == "inchikey" and v.id == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N" for v in votes)
    assert observed["smiles"] == ["CCO"]
    assert obs.structure_derivations[0]["status"] == "derived"
    structures._cached.cache_clear()
    monkeypatch.setattr(
        structures, "derive", lambda _: (_ for _ in ()).throw(AssertionError("cache missed"))
    )
    assert structures.cached_derivation("CCO")["inchikey"] == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    obs.identifiers.append({"ns": "inchikey", "id": "AAAAAAAAAAAAAA-BBBBBBBBSA-N"})
    votes, _ = votes_for(obs, CHEMICAL_POLICY)
    assert [v.id for v in votes if v.ns == "inchikey"] == ["AAAAAAAAAAAAAA-BBBBBBBBSA-N"]
