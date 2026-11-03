"""Ontology projection uses declared statement context, not resource/type lists."""

from omnipath_api.queries.parquet import ParquetQueries
from omnipath_core.fixtures import write_resource


def test_composition_is_not_an_ontology_but_axioms_are(tmp_path):
    entities = [
        dict(entity_key=key, identifier=identifier, label=label, namespace="arbitrary")
        for key, identifier, label in [
            ("food", "FOOD00005", "Allium"),
            ("compound", "FDB1", "Compound"),
            ("term", "TERM", "Term"),
            ("class", "CLASS", "Class"),
            ("unrelated", "OTHER", "Other"),
        ]
    ]
    labels = {e["entity_key"]: e["label"] for e in entities}
    relations = [
        dict(
            relation_key=f"{subject}-{obj}",
            subject_entity_key=subject,
            subject_label=labels[subject],
            object_entity_key=obj,
            object_label=labels[obj],
            predicate=predicate,
            statement_kind=kind,
        )
        for subject, obj, predicate, kind in [
            ("food", "compound", "has_part", "relation"),
            ("class", "term", "has_part", "ontology"),
            ("term", "unrelated", "related_to", "ontology"),
        ]
    ]
    write_resource(tmp_path / "resources/arbitrary_source/1", entities, relations)
    engine = ParquetQueries(data_root=tmp_path)
    food = engine.get_ontology_tree(["FOOD00005"])["root"]
    assert food is None
    assert engine.get_ontology_children("FOOD00005")["children"] == []
    tree = engine.get_ontology_tree(["TERM"])["root"]
    assert tree["id"] == "CLASS"
    assert [item["id"] for item in tree["children"]] == ["TERM"]
    assert [item["id"] for item in engine.get_ontology_children("CLASS")["children"]] == ["TERM"]
