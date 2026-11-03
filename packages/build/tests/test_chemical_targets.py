import json

import pytest

from library_fixture import ASPIRIN, VARIANT, WATER, build_fixture_library
from omnipath_resolver.resolver import EntityResolver
from omnipath_build.silver import RawEntityObservation, SilverExtractor
from writer_fixture import write_observations

MISSING = "ZZZZZZZZZZZZZZ-UHFFFAOYSA-N"


@pytest.fixture(scope="module")
def chemical_library(tmp_path_factory):
    # Authoritative identities are built by the same identity builders as real
    # resources (shared "exact-structures" reference: PubChem
    # CID 1 is WATER, ChEBI:1 a WATER-connectivity VARIANT, CHEMBL1 ASPIRIN);
    # modifying retired Parquet exports cannot affect lookup.
    return build_fixture_library(
        tmp_path_factory.mktemp("chemical-target-reference"), "exact-structures"
    )


@pytest.fixture
def resolver(chemical_library):
    instance = EntityResolver(library_dir=chemical_library)
    try:
        yield instance
    finally:
        instance.close()


def obs(extra):
    return RawEntityObservation(
        entity_key="x",
        entity_type="chemical_entity",
        namespace="pubchem",
        identifier="1",
        taxon=None,
        identifiers=[dict(ns=ns, id=ident) for ns, ident in extra],
    )


def test_supplied_key_overrides_conflicting_ids_and_does_not_expand(resolver):
    r = resolver
    try:
        targets = r.resolve_entity_targets(
            {"x": obs([("inchikey", ASPIRIN), ("chebi", "CHEBI:1")])}
        )["x"]
        assert [t.canonical_identifier for t in targets] == [ASPIRIN]
        assert targets[0].matched
    finally:
        r.close()


def test_absent_supplied_key_keeps_structure_identity(resolver):
    r = resolver
    try:
        targets = r.resolve_entity_targets({"x": obs([("inchikey", MISSING)])})["x"]
        assert [t.canonical_identifier for t in targets] == [MISSING]
        assert not targets[0].matched
    finally:
        r.close()


def test_conflicting_structures_do_not_project_by_connectivity(resolver):
    r = resolver
    try:
        targets = r.resolve_entity_targets({"x": obs([("chebi", "CHEBI:1")])})["x"]
        assert len(targets) == 1 and not targets[0].matched
        assert r.resolution_stats()["input_entities"] == 1
        different_connectivity = r.resolve_entity_targets({"y": obs([("chembl", "CHEMBL1")])})["y"]
        assert len(different_connectivity) == 1 and not different_connectivity[0].matched
        targets = r.resolve_entity_targets(
            {"z": obs([("inchikey", WATER), ("inchikey", VARIANT)])}
        )["z"]
        assert len(targets) == 1 and not targets[0].matched
        targets = r.resolve_entity_targets(
            {"bad": obs([("inchikey", WATER), ("inchikey", MISSING)])}
        )["bad"]
        assert len(targets) == 1 and not targets[0].matched
    finally:
        r.close()


def test_projected_relations_keep_payloads(resolver):
    r = resolver
    x = SilverExtractor("fixture", "test")
    x.process_record(
        {
            "subject": {
                "type": "chemical_entity",
                "identifiers": [
                    {"type": "pubchem", "value": "1"},
                    {"type": "chebi", "value": "CHEBI:1"},
                ],
            },
            "predicate": "interacts_with",
            "object": {"type": "protein", "identifiers": [{"type": "uniprot", "value": "P04637"}]},
        },
        {"measurement": "same evidence"},
        "row1",
        0,
    )
    try:
        tables = write_observations(r, x.entities, x.relations, x.payloads)
        entities, relations, payloads = (
            tables["entity"],
            tables["relation"],
            tables["evidence_payloads"],
        )
        by_key = {e["entity_key"]: e for e in entities}
        assert {(e["entity_type"], e["namespace"], e["identifier"]) for e in entities} == {
            ("chemical_entity", "pubchem", "1"),
            ("protein", "entrez", "7157"),
            ("protein", "uniprot", "P04637"),
        }
        assert len(relations) == len(payloads) == 1
        relation = relations[0]
        assert relation["evidence_count"] == 1
        assert relation["subject_entity_key"] in by_key and relation["object_entity_key"] in by_key
        (evidence,) = tables.children("relation_evidence", relation)
        forms = [evidence[side + "_molecular_form"] for side in ("subject", "object")]
        assert sum(form is None for form in forms) == 1
        product = by_key[next(form for form in forms if form)["protein_entity_key"]]
        assert (product["namespace"], product["identifier"]) == ("uniprot", "P04637")
        assert payloads[0]["relation_key"] == relation["relation_key"]
        assert json.loads(payloads[0]["payload_json"]) == {"measurement": "same evidence"}
    finally:
        r.close()


def test_compiled_exact_structure_identities_remain_distinct(resolver):
    observations = {
        key: RawEntityObservation(
            entity_key=key,
            entity_type="chemical_entity",
            namespace=namespace,
            identifier=identifier,
            taxon=None,
            identifiers=[],
        )
        for key, namespace, identifier in [
            ("water", "pubchem", "1"),
            ("variant", "chebi", "CHEBI:1"),
            ("aspirin", "chembl", "CHEMBL1"),
        ]
    }
    expected = {"water": WATER, "variant": VARIANT, "aspirin": ASPIRIN}
    targets = resolver.resolve_entity_targets(observations)
    assert WATER.split("-")[0] == VARIANT.split("-")[0]
    assert WATER != VARIANT
    for key, values in targets.items():
        assert len(values) == 1
        target = values[0]
        assert target.matched
        assert target.entity_type == "chemical_entity"
        assert (target.canonical_namespace, target.canonical_identifier) == (
            "inchikey",
            expected[key],
        )
        assert target.protein_identifier is None
