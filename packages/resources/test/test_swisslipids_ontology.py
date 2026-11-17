"""The ontology relations of a SwissLipids lipid: one per parent and class it names."""

from omnipath_resources.inputs_v2 import swisslipids


def objects(row):
    return [
        (relation.predicate.name if hasattr(relation.predicate, "name") else str(relation.predicate),
         relation.object.identifier)
        for relation in swisslipids._ontology_relations(row)
        if "subclass_of" in str(relation.predicate)
    ]


def test_a_parent_and_a_class_each_give_a_relation():
    found = objects(
        {"Lipid ID": "SLM:000000010", "Parent": "SLM:000000001", "Lipid class*": "SLM:000000002"}
    )
    assert [identifier for _, identifier in found] == ["SLM:000000001", "SLM:000000002"]


def test_several_parents_in_one_cell_each_give_a_relation():
    found = objects(
        {
            "Lipid ID": "SLM:000000010",
            "Parent": "SLM:000000001 | SLM:000000002",
            "Lipid class*": "SLM:000000003|SLM:000000004",
        }
    )
    assert [identifier for _, identifier in found] == [
        "SLM:000000001",
        "SLM:000000002",
        "SLM:000000003",
        "SLM:000000004",
    ]


def test_an_id_named_twice_or_naming_the_lipid_itself_gives_no_extra_relation():
    found = objects(
        {
            "Lipid ID": "SLM:000000010",
            "Parent": "SLM:000000001 | SLM:000000010",
            "Lipid class*": "SLM:000000001 | SLM:000000002",
        }
    )
    assert [identifier for _, identifier in found] == ["SLM:000000001", "SLM:000000002"]


def test_empty_and_placeholder_entries_are_dropped():
    found = objects(
        {"Lipid ID": "SLM:000000010", "Parent": " | None | SLM:000000001 | ", "Lipid class*": "nan"}
    )
    assert [identifier for _, identifier in found] == ["SLM:000000001"]


def test_no_relation_carries_a_pipe_in_its_identifier():
    relations = swisslipids._ontology_relations(
        {
            "Lipid ID": "SLM:000000010",
            "Parent": "SLM:000000001 | SLM:000000002",
            "Lipid class*": "SLM:000000003 | SLM:000000004",
            "Components*": "SLM:000000418 (sn2) / SLM:000000419",
        }
    )
    assert relations and all("|" not in relation.object.identifier for relation in relations)
