"""Global reference defaults and authoritative PubChem identity regression tests."""

from library_fixture import build_fixture_library
from omnipath_resolver import LibraryMatcher
from omnipath_build.silver import RawEntityObservation


def test_conflicting_authoritative_identities_remain_unresolved(tmp_path):
    # The shared "external-cid-claim" reference (see library_fixture) has an
    # external resource claim CID 962 for ASPIRIN although HMDB assigns it to
    # WATER. Deliberately conflicting PubChem evidence verifies ownership, not
    # real-world chemistry. It is a complete production build, copied per test.
    library = build_fixture_library(tmp_path, "external-cid-claim")
    matcher = LibraryMatcher(library)
    observation = RawEntityObservation(
        entity_key="bindingdb-example",
        entity_type="chemical_entity",
        namespace="pubchem",
        identifier="962",
        taxon=None,
        identifiers=[{"ns": "chebi", "id": "CHEBI:15377"}],
    )
    try:
        mouse = RawEntityObservation(
            entity_key="mouse",
            entity_type="protein",
            namespace="uniprot",
            identifier="P02340",
            taxon="10090",
            identifiers=[],
        )
        matches = matcher.match({observation.entity_key: observation, "mouse": mouse})
        result = matches[observation.entity_key]
        assert matches["mouse"].matched
        assert matches["mouse"].taxon == "10090"
    finally:
        matcher.close()
    assert result.canonical_identifier == "962"
    assert not result.matched
    assert "CHEBI:15377" in result.aliases["chebi"]
