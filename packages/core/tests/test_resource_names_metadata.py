"""The names of a resource in the metadata the API reads: slug, short, full and synonyms."""

import re

from omnipath_core import resource_metadata as metadata


def catalog():
    return metadata._catalog()


def test_every_resource_carries_a_short_and_a_full_name():
    missing = [slug for slug, entry in catalog().items() if not entry.get("short") or not entry.get("full_name")]
    assert not missing, missing


def test_every_slug_is_lowercase_letters_and_digits_and_every_short_name_has_no_space_or_underscore():
    assert [s for s in catalog() if not re.fullmatch(r"[a-z0-9]+", s)] == []
    assert [s for s, e in catalog().items() if re.search(r"[\s_]", e["short"])] == []
    assert [s for s, e in catalog().items() if "_" in e["full_name"]] == []


def test_the_renamed_resources_have_their_new_slugs_only():
    slugs = set(catalog())
    assert {"ramp", "omnipathontology", "phenolexplorer", "psimi", "switcheselm"} <= slugs
    assert not {"rampdb", "omnipath_ontology", "phenol_explorer", "psi_mi", "switches_elm"} & slugs


def test_a_slug_a_short_name_and_a_synonym_each_find_the_resource_whatever_their_case():
    slug = next(s for s, e in catalog().items() if e.get("synonyms"))
    entry = catalog()[slug]
    for spelling in (slug, slug.upper(), entry["short"], entry["short"].lower(), entry["synonyms"][0]):
        assert metadata.resolve_resource_name(spelling) == slug, spelling


def test_a_name_no_resource_has_resolves_to_nothing():
    assert metadata.resolve_resource_name("nonexistentresource") is None
    assert metadata.resolve_resource_name("") is None


def test_the_metadata_of_a_resource_carries_its_names_even_for_an_older_build_snapshot():
    found = metadata.resource_metadata("reactome", {"resource_metadata": {"website": "https://x"}})
    assert found["short"] and found["full_name"] and isinstance(found["synonyms"], list)
