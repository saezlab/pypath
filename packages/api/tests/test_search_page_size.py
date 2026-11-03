"""An entity found by search is the same entity whatever the page size."""

from omnipath_core.fixtures import fixture_key as key, write_resource
from omnipath_api.queries.parquet import ParquetQueries


def _entity(identifier, label, **extra):
    return dict(
        entity_key=key(identifier),
        identifier=identifier,
        label=label,
        namespace="test",
        entity_type="protein",
        **extra,
    )


def _relations(resource, identifier, count):
    return [
        dict(
            relation_key=key(f"{resource}-{identifier}-{i}"),
            subject_entity_key=key(identifier),
            object_entity_key=key(identifier),
            predicate="associated_with",
        )
        for i in range(count)
    ]


def test_relation_count_and_label_do_not_depend_on_the_page_size(tmp_path):
    # X matches "abc" exactly in one resource and only by a long alias in the other.
    # Enough better matches lie between the two that a small page stops reading hits
    # before the second one.
    fillers = [_entity(f"F{i}", f"abcde{i:02d}") for i in range(70)]
    write_resource(
        tmp_path / "resources/first/1",
        [_entity("X1", "abc"), _entity("Y1", "abcd"), *fillers],
        _relations("first", "X1", 1),
    )
    write_resource(
        tmp_path / "resources/second/1",
        [
            _entity(
                "X1",
                "Something else",
                identifiers=[dict(ns="alias", id="abcdefghijklmnopqrstuvwxyz", source="second")],
            )
        ],
        _relations("second", "X1", 5),
    )
    engine = ParquetQueries(tmp_path)
    try:
        found = {
            limit: engine.search_entities_api("abc", limit=limit)["entities"][0]
            for limit in (1, 2, 20, 100)
        }
        assert {entity["entityPk"] for entity in found.values()} == {key("X1")}
        assert {entity["relationCount"] for entity in found.values()} == {6}
        assert len({entity["label"] for entity in found.values()}) == 1
        # The count is the one the relation search finds for the entity.
        assert engine.search_relations_api({"entity_pks": [key("X1")]})["total"] == 6
        # Filters still apply to every copy.
        filtered = engine.search_entities_api("abc", filters={"entity_types": ["gene"]}, limit=1)
        assert filtered["entities"] == []
    finally:
        engine.close()
