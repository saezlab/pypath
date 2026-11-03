"""Molecular specificity, unknown coordinates and Parquet table contracts."""

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from omnipath_core import (
    ENTITY_EVIDENCE_TABLE,
    ENTITY_TABLE,
    RELATION_EVIDENCE_TABLE,
    SILVER_ENTITY_SCHEMA,
    SILVER_RELATION_SCHEMA,
    Entity,
    EntityRef,
    Identifier,
    MolecularForm,
    MolecularModification,
    MolecularRegion,
    MolecularVariant,
    SequenceIdentifier,
    CoordinateReference,
    entity_key,
    molecular_form_from_identifiers,
    normalize_molecular_form,
)


def test_identifier_retention_is_specificity_driven():
    assert molecular_form_from_identifiers([Identifier("uniprot", "P04637")]) is None
    form = molecular_form_from_identifiers(
        [
            Identifier("uniprot", "P04637-2"),
            Identifier("uniprot", "P04637-PRO_0000000001"),
            Identifier("refseq_protein", "NP_000537.3"),
            Identifier("enst", "ENST00000269305.9"),
            Identifier("entrez", "7157"),
        ]
    )
    assert form["isoform_identifier"] == {"ns": "uniprot", "id": "P04637-2"}
    assert len(form["sequence_identifiers"]) == 4
    assert form["sequence_identifiers"][2]["id"] == "NP_000537.3"
    assert form["protein_entity_key"] is None
    assert form["modifications"] is None
    multiple = molecular_form_from_identifiers(
        [
            Identifier("uniprot", "P04637-2"),
            Identifier("uniprot", "P04637-3"),
        ]
    )
    assert multiple["isoform_identifier"] is None
    assert len(multiple["sequence_identifiers"]) == 2


def test_empty_is_unspecified_and_product_keys_are_output_only():
    assert normalize_molecular_form(MolecularForm()) is None
    assert normalize_molecular_form({"modifications": [], "variants": []}) is None
    with pytest.raises(ValueError, match="resolution"):
        normalize_molecular_form({"protein_entity_key": "resolved"}, allow_resolved=False)
    with pytest.raises(ValueError, match="Unknown molecular form fields"):
        normalize_molecular_form({"entity_type": "protein"})


def test_structured_features_keep_coordinate_reference_and_unknown():
    identifier = SequenceIdentifier("uniprot", "P04637-2")
    form = normalize_molecular_form(
        MolecularForm(
            isoform_identifier=identifier,
            modifications=[
                MolecularModification(
                    term="MOD:00696",
                    residue="S",
                    position=15,
                    coordinate_reference=CoordinateReference(identifier, "protein", 1),
                )
            ],
            variants=[MolecularVariant(reference="R", alternate="H", position=175)],
        )
    )
    assert form["modifications"][0]["coordinate_reference"] == {
        "identifier": {"ns": "uniprot", "id": "P04637-2"},
        "coordinate_system": "protein",
        "position_base": 1,
    }
    assert form["variants"][0]["coordinate_reference"] == {
        "identifier": None,
        "coordinate_system": "unknown",
        "position_base": None,
    }
    # A resolved protein does not retroactively become a coordinate reference.
    form["protein_entity_key"] = entity_key("protein", "uniprot", "P04637")
    assert (
        normalize_molecular_form(form)["variants"][0]["coordinate_reference"]["identifier"] is None
    )


@pytest.mark.parametrize(
    "feature",
    [
        {"position": -1},
        {"position": True},
        {"position": 20, "end_position": 19},
        {"position": 1, "coordinate_reference": {"position_base": 2}},
        {"position": 0, "coordinate_reference": {"position_base": 1}},
        {"coordinate_reference": {"coordinate_system": "canonical"}},
    ],
)
def test_invalid_coordinate_values_fail_instead_of_becoming_unknown(feature):
    with pytest.raises(ValueError):
        normalize_molecular_form({"variants": [feature]})


def test_source_type_is_independent_of_gene_reference():
    gene_reference = "entrez:7157"
    protein = Entity("protein", [Identifier("entrez", "7157")])
    assert protein.type == "protein"
    assert protein.molecular_form is None
    assert EntityRef("protein", "entrez", "7157").molecular_form is None
    table = pa.Table.from_pylist(
        [
            {
                "entity_key": entity_key("protein", "entrez", "7157"),
                "entity_type": "protein",
                "namespace": "entrez",
                "identifier": "7157",
                "reference_entity_key": gene_reference,
                "gene_reference_keys": [gene_reference],
            },
            {
                "entity_key": entity_key("gene", "entrez", "7157"),
                "entity_type": "gene",
                "namespace": "entrez",
                "identifier": "7157",
                "reference_entity_key": gene_reference,
                "gene_reference_keys": [gene_reference],
            },
        ],
        schema=ENTITY_TABLE,
    )
    assert len(set(table["entity_key"].to_pylist())) == 2
    assert table["entity_type"].to_pylist() == ["protein", "gene"]


def test_occurrence_pairing_and_standalone_forms_round_trip(tmp_path):
    x = molecular_form_from_identifiers([Identifier("uniprot", "P04637-2")])
    y = molecular_form_from_identifiers([Identifier("uniprot", "P38398-3")])
    relations = pa.Table.from_pylist(
        [
            {"relation_id": 0, "ordinal": 0, "row_id": "xy",
             "subject_molecular_form": x, "object_molecular_form": y},
            {"relation_id": 0, "ordinal": 1, "row_id": "yx",
             "subject_molecular_form": y, "object_molecular_form": x},
        ],
        schema=RELATION_EVIDENCE_TABLE,
    )  # fmt: skip
    path = tmp_path / "relation_evidence.parquet"
    pq.write_table(relations, path)
    evidence = pq.read_table(path).to_pylist()
    assert evidence[0]["subject_molecular_form"] == evidence[1]["object_molecular_form"]
    assert evidence[0]["object_molecular_form"] == evidence[1]["subject_molecular_form"]
    entities = pa.Table.from_pylist(
        [{"entity_id": 0, "ordinal": 0, "row_id": "standalone", "molecular_form": x}],
        schema=ENTITY_EVIDENCE_TABLE,
    )
    path = tmp_path / "entity_evidence.parquet"
    pq.write_table(entities, path)
    assert pq.read_table(path).to_pylist()[0]["molecular_form"] == x
    assert (
        SILVER_ENTITY_SCHEMA.field("molecular_form").type
        == ENTITY_EVIDENCE_TABLE.field("molecular_form").type
    )
    assert SILVER_RELATION_SCHEMA.field("subject").type.field("molecular_form").type == (
        RELATION_EVIDENCE_TABLE.field("subject_molecular_form").type
    )


def test_regions_keep_partial_extent_and_round_trip(tmp_path):
    reference = CoordinateReference(SequenceIdentifier("uniprot", "P01308"), "protein", 1)
    form = normalize_molecular_form(
        MolecularForm(
            sequence_identifiers=[SequenceIdentifier("uniprot", "P01308-PRO_0000015819")],
            regions=[
                MolecularRegion(
                    type="Chain",
                    identifier=SequenceIdentifier("uniprot_feature", "PRO_0000015819"),
                    position=25,
                    end_position=54,
                    coordinate_reference=reference,
                    description="Insulin B chain",
                ),
                # '?..30': an unknown start stays unknown.
                MolecularRegion(type="Chain", end_position=30, coordinate_reference=reference),
            ],
        ),
        allow_resolved=False,
    )
    assert form["regions"][0]["identifier"] == {"ns": "uniprot_feature", "id": "PRO_0000015819"}
    assert form["regions"][1]["position"] is None
    assert form["regions"][1]["end_position"] == 30
    assert form["modifications"] is None and form["variants"] is None
    unreferenced = normalize_molecular_form({"regions": [{"type": "Fragment", "position": 3}]})
    assert unreferenced["regions"][0]["coordinate_reference"]["coordinate_system"] == "unknown"
    with pytest.raises(ValueError, match="end_position"):
        normalize_molecular_form({"regions": [{"position": 9, "end_position": 2}]})
    with pytest.raises(ValueError, match="Unknown molecular feature fields"):
        normalize_molecular_form({"regions": [{"start": 1}]})
    table = pa.Table.from_pylist(
        [{"entity_id": 0, "ordinal": 0, "row_id": "chain", "molecular_form": form}],
        schema=ENTITY_EVIDENCE_TABLE,
    )
    path = tmp_path / "entity_evidence.parquet"
    pq.write_table(table, path)
    assert pq.read_table(path).to_pylist()[0]["molecular_form"] == form
