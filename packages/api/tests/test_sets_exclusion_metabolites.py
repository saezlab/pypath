"""``/sets/exclusion-metabolites``: the exclusion table of a served release, as a table."""

import pytest
from fastapi.testclient import TestClient

from metsigdb_served import key, serve
from omnipath_api.queries.release import ReleaseQueries
from omnipath_api.release_database import build_served_database, publish_newest
from omnipath_api.server import create_app
from omnipath_core import subset_layout
from omnipath_core.fixtures import (
    metsigdb_exclusion_rows,
    metsigdb_report_rows,
    write_dataset,
    write_resource,
)

PATH = "/sets/exclusion-metabolites"
DEFAULT = [
    "entity",
    "metabolite_label",
    "inchikey",
    "hmdb",
    "pubchem",
    "chebi",
    "kegg",
    "exclusion_class",
    "assignment_source",
    "origin",
    "match_kind",
]
ROWS = metsigdb_exclusion_rows()


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    root = serve(tmp_path_factory.mktemp("data"))
    engine = ReleaseQueries(root)
    with TestClient(create_app(queries=engine)) as client:
        yield client
    engine.close()


def get(client, path=PATH, **params):
    return client.get(path, params=params)


def rows(client, **params):
    response = get(client, **params)
    assert response.status_code == 200, response.text
    return response.json()["rows"]


def test_the_table_is_returned_with_the_default_fields(client):
    body = get(client).json()
    assert body["count"] == len(body["rows"]) == len(ROWS)
    assert body["has_more"] is False and body["total"] is None
    assert list(body["rows"][0]) == DEFAULT


def test_a_class_filter_returns_that_class_only(client):
    found = rows(client, exclusion_class="Cofactors")
    assert {r["exclusion_class"] for r in found} == {"cofactors"}
    assert {r["entity"] for r in found} == {key("atp")}


def test_a_list_of_classes_is_the_union(client):
    found = rows(client, exclusion_class="ions,small_molecules")
    assert {r["entity"] for r in found} == {key("sodium"), key("water")}


def test_the_entity_may_be_named_by_its_key(client):
    found = rows(client, entity=key("water"))
    assert {r["entity"] for r in found} == {key("water")}


def test_a_metabolite_outside_every_class_has_no_row(client):
    assert rows(client, entity=key("glucose")) == []


def test_paging_neither_repeats_nor_skips_a_row(client):
    seen = []
    for offset in range(len(ROWS) + 1):
        seen += [(r["entity"], r["exclusion_class"]) for r in rows(client, limit=1, offset=offset)]
    assert seen == sorted(seen, key=lambda pair: (pair[1], pair[0]))
    assert len(seen) == len(set(seen)) == len(ROWS)


def test_total_is_the_size_of_the_whole_result(client):
    body = get(client, limit=1, total="true").json()
    assert (body["count"], body["has_more"], body["total"]) == (1, True, len(ROWS))


def test_all_fields_adds_the_assignment(client):
    row = rows(client, fields="all", limit=1)[0]
    assert {"seed_namespace", "seed_identifier", "evidence", "via_entity", "build_id"} <= set(row)


def test_a_parameter_the_route_does_not_publish_is_refused_by_name(client):
    response = get(client, hmdb="HMDB0000538")
    assert response.status_code == 400 and "hmdb" in response.text


@pytest.mark.parametrize(
    "params",
    [{"exclusion_class": "xenon"}, {"origin": "guessed"}, {"fields": "nope"}, {"entity": "??"}],
)
def test_a_value_the_route_does_not_publish_is_refused(client, params):
    assert get(client, **params).status_code == 400


def test_the_report_lists_what_the_build_measured(client):
    body = get(client, f"{PATH}/report", total="true").json()
    assert body["total"] == len(metsigdb_report_rows()) == len(body["rows"])
    assert {r["kind"] for r in body["rows"]} == {"class_size"}
    assert isinstance(body["rows"][0]["detail"], dict)


def test_the_report_filters_by_class(client):
    found = get(client, f"{PATH}/report", exclusion_class="ions").json()["rows"]
    assert [r["exclusion_class"] for r in found] == ["ions"]


def test_a_release_without_the_subset_answers_503(tmp_path):
    root = serve(tmp_path, subsets=False)
    engine = ReleaseQueries(root)
    with TestClient(create_app(queries=engine)) as bare:
        for path in (PATH, f"{PATH}/report"):
            response = bare.get(path)
            assert response.status_code == 503
            assert "exclusion metabolite" in response.json()["detail"]
    engine.close()


def test_the_metabolite_is_named_and_found_by_the_identifiers_the_release_holds(tmp_path):
    atp = key("atp")
    chemical = dict(
        entity_key=atp,
        identifier="HMDB0000538",
        label="ATP",
        namespace="hmdb",
        entity_type="chemical_entity",
        identifiers=[
            dict(ns="hmdb", id="HMDB00538", source="chem"),
            dict(ns="chebi", id="CHEBI:30616", source="chem"),
            dict(ns="kegg", id="C00002", source="chem"),
            dict(ns="inchikey", id="ZKHQWZAMYRWXGA-KQYNXXCUSA-J", source="chem"),
        ],
    )
    root = write_dataset(tmp_path / "data")
    write_resource(root / "resources" / "chem" / "1", [chemical])
    publish_newest(root, "1")
    subset_layout.build_product(
        root,
        "1",
        "exclusion_metabolites",
        _writer(
            {
                "exclusion_metabolite": [r for r in ROWS if r["entity_key"] == atp],
                "exclusion_metabolite_report": metsigdb_report_rows(),
            }
        ),
    )
    build_served_database(root, "1", keep=None, memory_limit="1GB", threads=1)
    engine = ReleaseQueries(root)
    with TestClient(create_app(queries=engine)) as http:
        for name in ("HMDB0000538", "hmdb00538", "CHEBI:30616", "C00002", "ZKHQWZAMYRWXGA-KQYNXXCUSA-J"):
            found = rows(http, entity=name)
            assert {r["entity"] for r in found} == {atp}, name
        row = rows(http, entity="C00002", exclusion_class="cofactors")[0]
        assert (row["metabolite_label"], row["hmdb"], row["chebi"], row["kegg"]) == (
            "ATP",
            "HMDB0000538",
            "30616",
            "C00002",
        )
        assert rows(http, entity="C00001") == []
    engine.close()


def _writer(tables):
    import pyarrow as pa
    import pyarrow.parquet as pq

    def write(directory, context):
        for name, content in tables.items():
            pq.write_table(pa.Table.from_pylist(content), directory / f"{name}.parquet")
        return {"rows": sum(len(content) for content in tables.values())}

    return write
