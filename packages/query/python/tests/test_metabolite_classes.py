"""The relation filter metabolite_classes: a relation is kept when one of its ends has a
row of the release's metabolite_class subset in a bucket named.

The release is small: two proteins and three chemicals, one in the bucket ``hmdb``, one
in ``chebi_role`` and one in none. The table ``metabolite_class`` is written straight
into the database, as the server copies a subset's Parquet in."""

import json
import shutil

import duckdb
import pytest

from omnipath_client import steps
from omnipath_client.database import build_release_database
from omnipath_core.fixtures import fixture_key as key, write_resource

ENTITIES = [
    ("P00001", "PROT1", "uniprot", "protein"),
    ("P00002", "PROT2", "uniprot", "protein"),
    ("HMDB0000063", "cortisol", "hmdb", "chemical_entity"),
    ("CHEMBL0001", "unclaimed", "chembl", "chemical_entity"),
    ("CHEBI:00003", "role metabolite", "chebi", "chemical_entity"),
]
# name: subject, object, predicate
RELATIONS = {
    "r_hmdb": ("HMDB0000063", "P00001", "affects"),
    "r_unclaimed": ("P00001", "CHEMBL0001", "interacts_with"),
    "r_proteins": ("P00001", "P00002", "interacts_with"),
    "r_role": ("CHEBI:00003", "P00002", "affects"),
    "r_hmdb_object": ("P00002", "HMDB0000063", "affects"),
}
CLASSES = {"HMDB0000063": "hmdb", "CHEBI:00003": "chebi_role"}
QUERIES = ("relation_page", "relation_counts", "relation_filter_options", "relation_facets")


def _label(ident):
    return next(label for i, label, _, _ in ENTITIES if i == ident)


def _kind(ident):
    return next(kind for i, _, _, kind in ENTITIES if i == ident)


def _write_release(root):
    entities = [
        dict(entity_key=key(ident), identifier=ident, label=label, namespace=ns,
             entity_type=kind, taxon="9606" if kind == "protein" else None)
        for ident, label, ns, kind in ENTITIES
    ]
    relations = [
        dict(relation_key=key(name), subject_entity_key=key(s), subject_label=_label(s),
             object_entity_key=key(o), object_label=_label(o), subject_type=_kind(s),
             object_type=_kind(o), predicate=predicate, category="interaction", taxon="9606",
             sources=["testres"])
        for name, (s, o, predicate) in RELATIONS.items()
    ]
    write_resource(root / "resources" / "testres" / "1", entities, relations)
    (root / "releases").mkdir()
    (root / "releases" / "1.json").write_text(
        json.dumps(dict(schema_version=1, version="1", resources={"testres": "1"}))
    )


@pytest.fixture(scope="module")
def databases(tmp_path_factory):
    """The release's database without the subset table, and a copy with it."""
    tmp = tmp_path_factory.mktemp("metabolite_classes")
    _write_release(tmp / "data")
    bare = tmp / "bare.duckdb"
    build_release_database(tmp / "data", "1", bare, memory_limit="1GB", threads=1)
    classed = tmp / "classed.duckdb"
    shutil.copy(bare, classed)
    with duckdb.connect(str(classed)) as db:
        db.execute("""CREATE TABLE metabolite_class (
            entity_key VARCHAR NOT NULL, label VARCHAR, inchikey VARCHAR, hmdb VARCHAR,
            metabolite_class VARCHAR NOT NULL, match_kind VARCHAR NOT NULL,
            via_entity_key VARCHAR, supported_by VARCHAR[] NOT NULL, build_id VARCHAR NOT NULL)""")
        for ident, bucket in CLASSES.items():
            db.execute(
                "INSERT INTO metabolite_class VALUES (?, ?, NULL, NULL, ?, 'entity', NULL, ['testres'], 'b1')",
                [key(ident), _label(ident), bucket],
            )
    return {"bare": bare, "classed": classed}


def _run(path, name, **params):
    with duckdb.connect(str(path), read_only=True) as con:
        steps.apply_settings(con)
        return json.loads(steps.run(con, name, params, output="json"))


def _page(path, **params):
    response = _run(path, "relation_page", limit=100, records=True, **params)
    found = {(r["subject_entity_key"], r["object_entity_key"]) for r in response["records"] or []}
    names = {name for name, (s, o, _) in RELATIONS.items() if (key(s), key(o)) in found}
    assert response["total"] == len(names)
    return names


def _exported(path, **params):
    with duckdb.connect(str(path), read_only=True) as con:
        steps.apply_settings(con)
        rows = steps.run(con, "relations", params).to_pylist()
    return {name for name, (s, o, _) in RELATIONS.items()
            for row in rows if (row["subject_entity_key"], row["object_entity_key"]) == (key(s), key(o))}


# Selections that take each path of the relation queries: no endpoint filter, one
# entity (a range of the endpoint table), scope entities (candidates), a source (scope).
SELECTIONS = [
    {},
    {"entities": [key("P00001")]},
    {"scope_entity_ids": [key("P00001")]},
    {"sources": ["testres"]},
]


@pytest.mark.parametrize("selection", SELECTIONS)
def test_the_filter_unset_changes_no_answer(databases, selection):
    for name in QUERIES:
        bare = _run(databases["bare"], name, **selection)
        assert _run(databases["classed"], name, **selection) == bare, name
        assert _run(databases["classed"], name, metabolite_classes=None, **selection) == bare, name
    assert _page(databases["classed"], **selection) == _page(databases["bare"], **selection)


@pytest.mark.parametrize("selection", SELECTIONS)
def test_a_bucket_keeps_the_relations_of_its_chemicals_at_either_end(databases, selection):
    every = _page(databases["classed"], **selection)
    kept = _page(databases["classed"], metabolite_classes=["hmdb"], **selection)
    assert kept == every & {"r_hmdb", "r_hmdb_object"}
    assert kept and "r_unclaimed" not in kept and "r_proteins" not in kept
    assert _exported(databases["classed"], metabolite_classes=["hmdb"], **selection) == kept
    both = _page(databases["classed"], metabolite_classes=["hmdb", "chebi_role"], **selection)
    assert both == every & {"r_hmdb", "r_hmdb_object", "r_role"}
    counts = _run(databases["classed"], "relation_counts", metabolite_classes=["hmdb"], **selection)
    assert sum(counts["categories"].values()) == len(kept)


def test_the_filter_narrows_the_filter_options(databases):
    options = _run(databases["classed"], "relation_filter_options", metabolite_classes=["chebi_role"])
    # only r_role, a chemical affecting a protein, is left
    assert options["predicatesByCategory"] == {"interaction": ["affects"]}
    entities = _run(databases["classed"], "relation_filter_options",
                    entities=[key("P00001")], metabolite_classes=["chebi_role"])
    assert entities["predicatesByCategory"] == {}


@pytest.mark.parametrize("selection", SELECTIONS)
def test_an_unknown_bucket_keeps_no_relation(databases, selection):
    assert _page(databases["classed"], metabolite_classes=["hmbd"], **selection) == set()
    assert _exported(databases["classed"], metabolite_classes=["hmbd"], **selection) == set()
    counts = _run(databases["classed"], "relation_counts", metabolite_classes=["hmbd"], **selection)
    assert sum(counts["categories"].values()) == 0
    options = _run(databases["classed"], "relation_filter_options", metabolite_classes=["hmbd"], **selection)
    assert options["predicatesByCategory"] == {} and options["sources"] == []


@pytest.mark.parametrize("name", [*QUERIES, "relations"])
def test_the_filter_on_a_release_without_the_subset_is_an_error(databases, name):
    with duckdb.connect(str(databases["bare"]), read_only=True) as con:
        steps.apply_settings(con)
        with pytest.raises(steps.QueryError, match="no metabolite_class subset"):
            steps.run(con, name, {"metabolite_classes": ["hmdb"]},
                      output="table" if name == "relations" else "json")
