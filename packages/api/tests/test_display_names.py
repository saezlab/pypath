from omnipath_api.queries.shape.display_name import preferred_name
from omnipath_api.queries.shape.entity import shape_entity_summary
from omnipath_api.queries.entity_details import page_details


def test_resolver_label_is_the_display_name_on_every_page():
    entity = shape_entity_summary(
        dict(
            entity_key="a",
            entity_type="chemical_entity",
            namespace="inchikey",
            identifier="QNAYBMKLOCPYGJ-UQEXSWPGSA-N",
            label="alanine-d7",
            identifiers=[dict(ns="name", id=f"Long chemical name {i}") for i in range(25)]
            + [dict(ns="name", id="Alanine")],
        )
    )
    first = page_details(entity)
    assert first["displayName"] == "alanine-d7"
    assert first["identifiersNextCursor"] == "20"
    assert page_details(entity, offset=20)["displayName"] == "alanine-d7"


def test_numeric_chemical_labels_fall_back_to_a_named_identifier():
    entity = shape_entity_summary(
        dict(
            entity_key="a",
            entity_type="small_molecule",
            namespace="pubchem",
            identifier="1234",
            label="1234",
            identifiers=[dict(ns="chebi", id="CHEBI:15377")],
        )
    )
    assert entity["displayName"] == "CHEBI:15377"


def test_group_names_prefer_plain_names():
    assert preferred_name(["Alanine", "alanine-d7", "CID_1234"]) == "Alanine"
