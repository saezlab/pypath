"""Gene references preserve source types, explicit products and occurrence counts."""

import pytest

from library_fixture import build_fixture_library
from omnipath_resolver import EntityResolver
from test_resolver import obs
from test_partition_contract import run_rows


@pytest.fixture
def library(tmp_path):
    return build_fixture_library(tmp_path / "reference")


def test_gene_targets_and_accession_identity(library):
    resolver = EntityResolver(library_dir=library, defer_aliases=True)
    try:
        observations = {
            "gene": obs("gene", "gene", "entrez", "7157"),
            "shared": obs("shared", "gene", "entrez", "801"),
            "shared2": obs("shared2", "gene", "entrez", "805"),
            "orphan": obs("orphan", "gene", "entrez", "55"),
            "protein": obs("protein", "protein", "uniprot", "A0A0U1RQF1"),
            "secondary": obs("secondary", "protein", "uniprot", "Q15086"),
        }
        targets = resolver.resolve_entity_targets(observations)

        expected = {
            "gene": "7157",
            "shared": "801",
            "shared2": "805",
            "orphan": "55",
            "protein": "7157",
            "secondary": "7157",
        }
        for key, identifier in expected.items():
            assert len(targets[key]) == 1
            target = targets[key][0]
            assert (target.canonical_namespace, target.canonical_identifier) == (
                "entrez",
                identifier,
            )
            assert target.entity_type == observations[key].entity_type
            assert target.gene_candidates == (f"entrez:{identifier}",)
        assert targets["protein"][0].protein_identifier == "A0A0U1RQF1"
        assert targets["secondary"][0].protein_identifier == "P04637"
        assert all(
            targets[key][0].protein_identifier is None
            for key in ("gene", "shared", "shared2", "orphan")
        )
        assert resolver.resolution_stats()["resolved_entities"] == len(observations)
        assert resolver.resolution_stats()["input_entities"] == len(observations)
    finally:
        resolver.close()


def entity(kind, ns, value):
    return {"type": kind, "identifiers": [{"type": ns, "value": value}]}


def test_writer_gene_references_are_partition_invariant(library, tmp_path):
    gene = entity("gene", "entrez", "7157")
    other_gene = entity("gene", "entrez", "999")
    phenotype = entity("ontology_class", "hpo", "HP:0000001")
    rows = [
        (
            "hpo",
            {
                "subject": gene,
                "predicate": "associated_with",
                "object": phenotype,
                "annotations": [{"term": "description", "value": "gene evidence"}],
            },
        ),
        ("network", {"subject": gene, "predicate": "affects", "object": other_gene}),
        (
            "protein",
            {
                "subject": entity("protein", "uniprot", "P04637"),
                "predicate": "associated_with",
                "object": phenotype,
            },
        ),
    ]
    tables = run_rows(tmp_path / "all", rows, 3, library)
    assert run_rows(tmp_path / "split", rows[::-1], 1, library).comparable() == tables.comparable()
    entities, relations, payloads = (
        tables["entity"],
        tables["relation"],
        tables["evidence_payloads"],
    )
    assert {(e["entity_type"], e["namespace"], e["identifier"]) for e in entities} == {
        ("gene", "entrez", "7157"),
        ("gene", "entrez", "999"),
        ("protein", "entrez", "7157"),
        ("protein", "uniprot", "P04637"),
        ("ontology_class", "hpo", "HP:0000001"),
    }
    assert len(relations) == len(payloads) == 3
    assert sum(r["evidence_count"] for r in relations) == 3
    products = {e["entity_key"]: e for e in entities if e["namespace"] == "uniprot"}
    for relation in relations:
        evidence = tables.children("relation_evidence", relation)[0]
        if relation["subject_type"] == "gene":
            assert evidence["subject_molecular_form"] is None
        else:
            assert (
                products[evidence["subject_molecular_form"]["protein_entity_key"]]["identifier"]
                == "P04637"
            )
        assert relation["subject_reference_entity_key"] == "entrez:7157"
    assert len({r["subject_entity_key"] for r in relations}) == 2
    for ev in tables["relation_evidence"]:
        if ev["row_id"] == "hpo":
            assert len(ev["annotations"]) == 1
    for e in entities:
        if e["identifier"] == "P04637":
            assert not any(
                i["ns"] == "uniprot" and i["id"] == "A0A0U1RQF1"
                for i in tables.children("entity_identifier", e)
            )


def test_symmetric_gene_reference_keeps_one_occurrence(library, tmp_path):
    gene = entity("gene", "entrez", "999")
    rows = [
        (
            "symmetric",
            {
                "subject": gene,
                "predicate": "interacts_with",
                "object": gene,
                "annotations": [{"term": "description", "value": "once"}],
            },
        )
    ]
    tables = run_rows(tmp_path / "symmetric", rows, 1, library)
    entities, relations, payloads = (
        tables["entity"],
        tables["relation"],
        tables["evidence_payloads"],
    )
    assert len(entities) == len(relations) == len(payloads) == 1
    assert entities[0]["entity_type"] == "gene"
    assert relations[0]["subject_entity_key"] == relations[0]["object_entity_key"]
    (evidence,) = tables["relation_evidence"]
    assert evidence["subject_molecular_form"] is None
    assert relations[0]["evidence_count"] == 1
    assert len(evidence["annotations"]) == 1


def test_identity_library_without_kv_store_requires_rebuild(library):
    import shutil

    from omnipath_resolver.identity_kv import KvMissing

    shutil.rmtree(library / "kv")
    with pytest.raises(KvMissing, match="build-identity-kv"):
        EntityResolver(library_dir=library)


def test_directory_that_is_not_an_identity_library_is_rejected(library):
    (library / "manifest.json").unlink()  # A read-only hard link to the shared template.
    (library / "manifest.json").write_text('{"format": "omnipath-full-two-index-msgpack-zstd-v2"}')
    with pytest.raises(ValueError, match="identity directory is required"):
        EntityResolver(library_dir=library)


def test_shared_protein_and_unmapped_gene_keep_all_evidence(library, tmp_path):
    phenotype = entity("ontology_class", "hpo", "HP:0000001")
    rows = [
        (
            gene,
            {
                "subject": entity("gene", "entrez", gene),
                "predicate": "associated_with",
                "object": phenotype,
            },
        )
        for gene in ["801", "805", "55"]
    ]
    tables = run_rows(tmp_path / "shared", rows, 3, library)
    entities, relations, payloads = (
        tables["entity"],
        tables["relation"],
        tables["evidence_payloads"],
    )
    assert {(e["entity_type"], e["identifier"]) for e in entities} == {
        ("gene", "801"),
        ("gene", "805"),
        ("gene", "55"),
        ("ontology_class", "HP:0000001"),
    }
    assert sorted(r["evidence_count"] for r in relations) == [1, 1, 1]
    assert len(payloads) == 3


def test_protein_source_with_gene_identifier_does_not_assert_product(library):
    resolver = EntityResolver(library_dir=library)
    try:
        targets = resolver.resolve_entity_targets({"p": obs("p", "protein", "entrez", "7157")})["p"]
        assert len(targets) == 1
        assert targets[0].matched
        assert targets[0].canonical_identifier == "7157"
        assert targets[0].entity_type == "protein"
        assert targets[0].protein_identifier is None
    finally:
        resolver.close()


def test_product_review_status_does_not_expand_gene_observation(library):
    # Gene 4242 only has unreviewed products (entry name P3333x_HUMAN carries the
    # accession itself), both in the shared fixture hubs.
    resolver = EntityResolver(library)
    try:
        targets = resolver.resolve_entity_targets({"g": obs("g", "gene", "entrez", "4242")})["g"]
        assert len(targets) == 1
        assert targets[0].canonical_identifier == "4242"
        assert targets[0].entity_type == "gene"
        assert targets[0].protein_identifier is None
    finally:
        resolver.close()


def test_ambiguous_and_conflicting_gene_links_preserve_asserted_product(library):
    resolver = EntityResolver(library)
    try:
        observations = {
            "multi": obs("multi", "protein", "uniprot", "P0DP23"),
            "conflict": obs("conflict", "protein", "uniprot", "P04637", [("entrez", "805")]),
        }
        targets = resolver.resolve_entity_targets(observations)
        multi, conflict = targets["multi"][0], targets["conflict"][0]
        assert all(len(values) == 1 for values in targets.values())
        assert (multi.canonical_namespace, multi.canonical_identifier) == ("uniprot", "P0DP23")
        assert multi.gene_mapping_status == "ambiguous"
        assert set(multi.gene_candidates) == {"entrez:801", "entrez:805", "entrez:808"}
        assert (conflict.canonical_namespace, conflict.canonical_identifier) == (
            "uniprot",
            "P04637",
        )
        assert conflict.gene_mapping_status == "conflict"
        assert set(conflict.gene_candidates) == {"entrez:7157", "entrez:805"}
        assert conflict.protein_gene_candidates == ("entrez:7157",)
        assert conflict.matched and multi.matched
    finally:
        resolver.close()
