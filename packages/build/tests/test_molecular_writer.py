"""Molecular occurrence identity survives gene grouping and Parquet reduction."""

import json

import pytest

from omnipath_resolver.resolver import ResolvedEntityTarget
from omnipath_build.silver import SilverExtractor
from omnipath_build.writer import ParquetWriter
from omnipath_core.keys import entity_key
from tables_fixture import read_tables


class MappedResolver:
    """Fixed gene/product assignments; resolver policy is tested separately."""

    library_dir = None

    def resolve_entity_targets(self, entities, **kwargs):
        result = {}
        for key, raw in entities.items():
            product = raw.identifier.split("-", 1)[0] if raw.entity_type == "protein" else None
            gene = "2" if product == "P67890" else "1"
            result[key] = [
                ResolvedEntityTarget(
                    canonical_namespace="entrez",
                    canonical_identifier=gene,
                    label=f"GENE{gene}",
                    entity_type=raw.entity_type,
                    taxon="9606",
                    aliases={},
                    protein_namespace="uniprot" if product else None,
                    protein_identifier=product,
                    protein_gene_candidates=(f"entrez:{gene}",) if product else (),
                    gene_candidates=(f"entrez:{gene}",),
                    gene_mapping_status="resolved",
                )
            ]
        return result


def protein(identifier, *, variants=(), modifications=()):
    return {
        "type": "protein",
        "identifiers": [{"type": "uniprot", "value": identifier}],
        "molecular_form": {
            "isoform_identifier": {"ns": "uniprot", "id": identifier}
            if "-" in identifier
            else None,
            "variants": [{"description": value} for value in variants],
            "modifications": list(modifications),
        },
    }


def build(path, rows, *, shards=False):
    writer = ParquetWriter(path)
    resolver = MappedResolver()
    for i, record in enumerate(rows):
        extractor = SilverExtractor("fixture", "molecular")
        extractor.process_record(record, record, str(i), i)
        if shards:
            shard = ParquetWriter(path / f"shard-{i}")
            shard.append_observations(extractor, resolver)
            writer.import_observation_shard(shard.seal_observation_shard())
        else:
            writer.append_observations(extractor, resolver)
    return read_tables(writer.close()["files"])


def test_different_products_and_forms_share_gene_without_losing_occurrences(tmp_path):
    rows = [
        {
            "subject": protein("P12345-2", variants=("X", "Y")),
            "predicate": "affects",
            "object": protein("P67890"),
        },
        {
            "subject": protein("Q12345", variants=("X", "Z")),
            "predicate": "affects",
            "object": protein("P67890"),
        },
        {
            "subject": {"type": "gene", "identifiers": [{"type": "entrez", "value": "1"}]},
            "predicate": "affects",
            "object": protein("P67890"),
        },
    ]
    tables = build(tmp_path / "direct", rows)
    entities, relations = tables["entity"], tables["relation"]
    assert len(relations) == 2
    by_key = {e["entity_key"]: e for e in entities}
    protein_relation = next(r for r in relations if r["subject_type"] == "protein")
    gene_relation = next(r for r in relations if r["subject_type"] == "gene")
    assert gene_relation["subject_entity_key"] != protein_relation["subject_entity_key"]
    assert (
        gene_relation["subject_reference_entity_key"]
        == protein_relation["subject_reference_entity_key"]
        == "entrez:1"
    )
    assert tables.children("relation_evidence", gene_relation)[0]["subject_molecular_form"] is None
    assert protein_relation["evidence_count"] == 2
    observed = {}
    for evidence in tables.children("relation_evidence", protein_relation):
        form = evidence["subject_molecular_form"]
        product = by_key[form["protein_entity_key"]]
        assert product["reference_entity_key"] == "entrez:1"
        assert product["entity_type"] == "protein"
        observed[product["identifier"]] = [v["description"] for v in form["variants"]]
    assert observed == {"P12345": ["X", "Y"], "Q12345": ["X", "Z"]}
    anchor = by_key[protein_relation["subject_entity_key"]]
    assert anchor["namespace"] == "entrez" and anchor["entity_type"] == "protein"
    assert not any(i["id"] == "P12345-2" for i in tables.children("entity_identifier", anchor))
    assert (
        next(
            e for e in tables.children("relation_evidence", protein_relation) if e["row_id"] == "0"
        )["subject_molecular_form"]["isoform_identifier"]["id"]
        == "P12345-2"
    )
    sharded = build(tmp_path / "sharded", rows, shards=True)
    assert sharded.comparable() == tables.comparable()


def test_explicit_forms_do_not_collide_before_resolution(tmp_path):
    extractor = SilverExtractor("fixture", "molecular")
    for i, combination in enumerate((("X", "Y"), ("X", "Z"))):
        item = {
            "subject": protein("P12345", variants=combination),
            "predicate": "affects",
            "object": protein("P67890"),
        }
        extractor.process_record(item, item, str(i), i)
    assert len(extractor.entities) == 3
    writer = ParquetWriter(tmp_path)
    writer.append_observations(extractor, MappedResolver())
    tables = read_tables(writer.close()["files"])
    (relation,) = tables["relation"]
    combinations = {
        tuple(v["description"] for v in e["subject_molecular_form"]["variants"])
        for e in tables.children("relation_evidence", relation)
    }
    assert combinations == {("X", "Y"), ("X", "Z")}


def test_standalone_forms_and_unknown_coordinates_remain_evidence(tmp_path):
    rows = [
        protein(
            "P12345-2", modifications=({"term": "phosphorylation", "residue": "S", "position": 15},)
        ),
        protein("P12345-2", variants=("X",)),
    ]
    tables = build(tmp_path, rows)
    entities, payloads = tables["entity"], tables["evidence_payloads"]
    assert tables["relation"] == []
    anchor = next(e for e in entities if e["namespace"] == "entrez")
    assert len(tables.children("entity_evidence", anchor)) == 2
    first = tables.children("entity_evidence", anchor)[0]["molecular_form"]
    assert first["modifications"][0]["coordinate_reference"] == {
        "identifier": None,
        "coordinate_system": "unknown",
        "position_base": None,
    }
    assert first["variants"] is None
    assert tables.children("entity_evidence", anchor)[1]["molecular_form"]["modifications"] is None
    assert all(p["entity_key"] == anchor["entity_key"] for p in payloads)


def test_symmetric_relation_keeps_each_form_with_its_endpoint(tmp_path):
    a, b = protein("P12345-2", variants=("A",)), protein("P67890", variants=("B",))
    rows = [
        {"subject": a, "predicate": "interacts_with", "object": b},
        {"subject": b, "predicate": "interacts_with", "object": a},
    ]
    tables = build(tmp_path, rows)
    (relation,) = tables["relation"]
    assert relation["evidence_count"] == 2
    sides = {
        "entrez:1": ("A", entity_key("protein", "uniprot", "P12345")),
        "entrez:2": ("B", entity_key("protein", "uniprot", "P67890")),
    }
    for side in ("subject", "object"):
        label, product = sides[relation[f"{side}_reference_entity_key"]]
        for occurrence in tables.children("relation_evidence", relation):
            form = occurrence[f"{side}_molecular_form"]
            assert form["variants"][0]["description"] == label
            assert form["protein_entity_key"] == product


def test_source_isoform_is_captured_even_without_explicit_form():
    record = {"type": "protein", "identifiers": [{"type": "uniprot", "value": "P12345-2"}]}
    extractor = SilverExtractor("fixture", "molecular")
    extractor.process_record(record, record, "1", 1)
    (raw,) = extractor.entities.values()
    assert raw.molecular_form["isoform_identifier"] == {"ns": "uniprot", "id": "P12345-2"}
    json.dumps(raw.molecular_form)


@pytest.mark.parametrize(
    "source_type,namespace,identifier",
    [
        ("protein", "uniprot", "P12345-2"),
        ("transcript", "enst", "ENST00000000001.3"),
    ],
)
def test_primary_sequence_survives_an_explicit_feature_form(
    tmp_path, source_type, namespace, identifier
):
    item = {
        "type": source_type,
        "identifiers": [{"type": namespace, "value": identifier}],
        "molecular_form": {"variants": [{"description": "reported variant"}]},
    }
    tables = build(tmp_path, [item])
    entities = tables["entity"]
    anchor = next(e for e in entities if e["namespace"] == "entrez")
    form = tables.children("entity_evidence", anchor)[0]["molecular_form"]
    assert {"ns": namespace, "id": identifier} in form["sequence_identifiers"]
    assert form["variants"][0]["description"] == "reported variant"
    if source_type == "protein":
        assert form["isoform_identifier"] == {"ns": namespace, "id": identifier}


def test_explicit_gene_variant_is_retained_without_inferred_product(tmp_path):
    item = {
        "type": "gene",
        "identifiers": [{"type": "entrez", "value": "1"}],
        "molecular_form": {
            "variants": [
                {
                    "reference": "C",
                    "alternate": "T",
                    "position": 100,
                    "coordinate_reference": {
                        "identifier": {"ns": "refseq", "id": "NC_000001.11"},
                        "coordinate_system": "genomic",
                        "position_base": 1,
                    },
                }
            ]
        },
    }
    tables = build(tmp_path, [item])
    assert tables["relation"] == []
    (entity,) = tables["entity"]
    assert entity["entity_type"] == "gene"
    assert entity["reference_entity_key"] == "entrez:1"
    form = tables.children("entity_evidence", entity)[0]["molecular_form"]
    assert form["protein_entity_key"] is None
    assert form["transcript_entity_key"] is None
    assert form["variants"][0]["coordinate_reference"]["identifier"]["id"] == "NC_000001.11"


def test_complexes_with_distinct_member_forms_do_not_merge(tmp_path):
    from omnipath_resolver.resolver import EntityResolver

    class ComplexResolver(MappedResolver):
        def resolve_entity_targets(self, entities, **kwargs):
            with_native = EntityResolver()
            try:
                result = super().resolve_entity_targets(
                    {k: v for k, v in entities.items() if v.entity_type == "protein"}
                )
                result.update(
                    with_native.resolve_entity_targets(
                        {k: v for k, v in entities.items() if v.entity_type != "protein"}
                    )
                )
                return result
            finally:
                with_native.close()

    writer = ParquetWriter(tmp_path)
    extractor = SilverExtractor("fixture", "complexes")
    for i, accession in enumerate(("P12345-1", "P12345-2")):
        item = {
            "type": "macromolecular_complex",
            "identifiers": [{"type": "signor", "value": f"C{i}"}],
            "membership": [{"member": protein(accession), "predicate": "has_member"}],
        }
        extractor.process_record(item, item, str(i), i)
    writer.append_observations(extractor, ComplexResolver())
    tables = read_tables(writer.close()["files"])
    complexes = [row for row in tables["entity"] if row["entity_type"] == "macromolecular_complex"]
    assert len(complexes) == 2
    assert len({row["entity_key"] for row in complexes}) == 2
    assert all(row["namespace"] == "complex" for row in complexes)


def test_conflicting_source_keeps_native_reference_and_occurrence_diagnostics(tmp_path):
    class ConflictResolver:
        library_dir = None

        def resolve_entity_targets(self, entities, **kwargs):
            return {
                key: [
                    ResolvedEntityTarget(
                        canonical_namespace="uniprot",
                        canonical_identifier="P04637",
                        label="TP53",
                        entity_type=raw.entity_type,
                        taxon="9606",
                        aliases={},
                        protein_namespace="uniprot",
                        protein_identifier="P04637",
                        protein_gene_candidates=("entrez:7157",),
                        gene_candidates=("entrez:55", "entrez:7157"),
                        gene_mapping_status="conflict",
                    )
                ]
                for key, raw in entities.items()
            }

    item = protein("P04637-2")
    extractor = SilverExtractor("fixture", "conflicting")
    extractor.process_record(item, item, "entity", 0)
    relation = {"subject": item, "predicate": "affects", "object": item}
    extractor.process_record(relation, relation, "relation", 1)
    writer = ParquetWriter(tmp_path)
    writer.append_observations(extractor, ConflictResolver())
    tables = read_tables(writer.close()["files"])
    (entity,) = tables["entity"]
    (relation,) = tables["relation"]
    assert entity["reference_entity_key"] == "uniprot:P04637"
    assert entity["gene_reference_keys"] == ["entrez:7157"]
    assert relation["subject_reference_entity_key"] == entity["reference_entity_key"]
    assert relation["object_reference_entity_key"] == entity["reference_entity_key"]
    annotations = tables.children("entity_evidence", entity)[0]["annotations"]
    assert {a["value"] for a in annotations if a["term"] == "omnipath:gene_mapping_candidate"} == {
        "entrez:55",
        "entrez:7157",
    }
    assert any(
        a["term"] == "omnipath:gene_mapping_status" and a["value"] == "conflict"
        for a in annotations
    )
    annotations = tables.children("relation_evidence", relation)[0]["annotations"]
    assert {a["scope"] for a in annotations if a["term"] == "omnipath:gene_mapping_status"} == {
        "subject",
        "object",
    }


@pytest.mark.parametrize("with_gene", [False, True])
@pytest.mark.parametrize(
    "namespace,identifier,product_namespace",
    [
        ("uniprot", "Q9Y6K9", "uniprot"),
        ("uniprot", "Q9Y6K9-2", "uniprot"),
        ("uniprot", "Q9Y6K9-2-PRO_000001", "uniprot"),
        ("ensembl", "ENSP999999999999.73", "ensp"),
        ("refseq", "NP_999999999999.73", "refseq"),
        ("refseq_protein", "NP_999999999999.73", "refseq_protein"),
    ],
)
def test_reported_native_product_is_retained_without_catalogue_links(
    tmp_path, with_gene, namespace, identifier, product_namespace
):
    from library_fixture import build_fixture_library
    from omnipath_resolver.resolver import EntityResolver

    library = build_fixture_library(tmp_path / "reference")
    resolver = EntityResolver(library_dir=library)
    ids = [{"type": namespace, "value": identifier}]
    if with_gene:
        ids.append({"type": "entrez", "value": "7157"})
    item = {"type": "protein", "identifiers": ids}
    extractor = SilverExtractor("fixture", "native-product")
    extractor.process_record(item, item, "standalone", 0)
    relation = {
        "subject": item,
        "predicate": "affects",
        "object": {"type": "gene", "identifiers": [{"type": "entrez", "value": "55"}]},
    }
    extractor.process_record(relation, relation, "relation", 1)
    writer = ParquetWriter(tmp_path / "output", library_dir=library)
    try:
        writer.append_observations(extractor, resolver)
        tables = read_tables(writer.close()["files"])
    finally:
        resolver.close()
    entities = {row["entity_key"]: row for row in tables["entity"]}
    (relation,) = tables["relation"]
    product_key = entity_key("protein", product_namespace, identifier)
    product = entities[product_key]
    assert product["reference_entity_key"] == f"{product_namespace}:{identifier}"
    assert product["gene_reference_keys"] == []
    assert product["label"] == identifier
    form = tables.children("relation_evidence", relation)[0]["subject_molecular_form"]
    assert form["protein_entity_key"] == product_key
    assert any(
        annotation["term"] == "omnipath:protein_mapping_status"
        and annotation["value"] == "reported"
        and annotation["scope"] == "subject"
        for annotation in tables.children("relation_evidence", relation)[0]["annotations"]
    )
    source_endpoint = entities[relation["subject_entity_key"]]
    assert any(
        annotation["term"] == "omnipath:protein_mapping_status"
        and annotation["value"] == "reported"
        for annotation in tables.children("entity_evidence", source_endpoint)[0]["annotations"]
    )
    assert any(
        alias["ns"] == product_namespace and alias["id"] == identifier and alias["source"] == "raw"
        for alias in tables.children("entity_identifier", product)
    )
    if with_gene:
        assert not any(
            alias["source"] == "resolver" for alias in tables.children("entity_identifier", product)
        )
    if namespace == "uniprot" and identifier.endswith("-2"):
        assert form["isoform_identifier"]["id"] == identifier
    if namespace != "uniprot" or "-" in identifier:
        assert identifier in {value["id"] for value in form["sequence_identifiers"]}
    if with_gene:
        anchor = entities[relation["subject_entity_key"]]
        assert anchor["namespace"] == "entrez"
        assert anchor["identifier"] == "7157"
        assert relation["subject_entity_key"] != product_key
        assert (
            tables.children("entity_evidence", anchor)[0]["molecular_form"]["protein_entity_key"]
            == product_key
        )
        assert product["label"] != anchor["label"]


@pytest.mark.parametrize(
    "entity_type,namespace,identifier",
    [
        ("gene", "uniprot", "Q9Y6K9"),
        ("physical_entity", "uniprot", "Q9Y6K9"),
        ("protein", "genesymbol", "Q9Y6K9"),
        ("protein", "entrez", "7157"),
    ],
)
def test_gene_symbol_or_generic_source_does_not_create_native_product(
    tmp_path, entity_type, namespace, identifier
):
    from omnipath_resolver.resolver import EntityResolver

    item = {"type": entity_type, "identifiers": [{"type": namespace, "value": identifier}]}
    extractor = SilverExtractor("fixture", "negative-native")
    extractor.process_record(item, item, "entity", 0)
    resolver = EntityResolver(library_dir=tmp_path / "absent")
    writer = ParquetWriter(tmp_path / "output")
    try:
        writer.append_observations(extractor, resolver)
        tables = read_tables(writer.close()["files"])
    finally:
        resolver.close()
    (entity,) = tables["entity"]
    assert tables.children("entity_evidence", entity)[0]["molecular_form"] is None
    assert not any(
        annotation["term"] == "omnipath:protein_mapping_status"
        for annotation in tables.children("entity_evidence", entity)[0]["annotations"]
    )


@pytest.mark.parametrize("known_accession", ["P04637", "Q99999"])
def test_reported_protein_diagnostics_are_endpoint_scoped_and_catalogue_aliases_unchanged(
    tmp_path, known_accession
):
    from library_fixture import build_fixture_library
    from omnipath_resolver.resolver import EntityResolver

    library = build_fixture_library(tmp_path / "reference")
    resolver = EntityResolver(library_dir=library)
    known = {"type": "protein", "identifiers": [{"type": "uniprot", "value": known_accession}]}
    reported = {"type": "protein", "identifiers": [{"type": "uniprot", "value": "Q9Y6K9"}]}
    extractor = SilverExtractor("fixture", "protein-provenance")
    extractor.process_record(known, known, "known", 0)
    extractor.process_record(reported, reported, "reported", 1)
    relation = {"subject": known, "predicate": "affects", "object": reported}
    extractor.process_record(relation, relation, "relation", 2)
    writer = ParquetWriter(tmp_path / "output", library_dir=library)
    try:
        writer.append_observations(extractor, resolver)
        tables = read_tables(writer.close()["files"])
    finally:
        resolver.close()
    entities = {row["entity_key"]: row for row in tables["entity"]}
    (relation,) = tables["relation"]
    reported_status = [
        annotation
        for annotation in tables.children("relation_evidence", relation)[0]["annotations"]
        if annotation["term"] == "omnipath:protein_mapping_status"
    ]
    assert [(annotation["value"], annotation["scope"]) for annotation in reported_status] == [
        ("reported", "object")
    ]
    known_endpoint = entities[relation["subject_entity_key"]]
    assert not any(
        annotation["term"] == "omnipath:protein_mapping_status"
        for evidence in tables.children("entity_evidence", known_endpoint)
        for annotation in evidence["annotations"]
    )
    known_product = entities[entity_key("protein", "uniprot", known_accession)]
    assert any(
        alias["source"] == "resolver"
        for alias in tables.children("entity_identifier", known_product)
    )
    reported_endpoint = entities[relation["object_entity_key"]]
    assert any(
        annotation["term"] == "omnipath:protein_mapping_status"
        and annotation["value"] == "reported"
        for evidence in tables.children("entity_evidence", reported_endpoint)
        for annotation in evidence["annotations"]
    )
