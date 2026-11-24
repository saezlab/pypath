from fastapi.testclient import TestClient
from served import served_queries
from omnipath_api.server import create_app
from omnipath_core.fixtures import write_resource


def _client(tmp_path):
    folder = tmp_path / "resources" / "test" / "1"
    folder.mkdir(parents=True)
    rows = [
        dict(
            entity_key="tp53-protein",
            entity_type="protein",
            namespace="uniprot",
            identifier="P04637",
            label="TP53",
            identifiers=[
                dict(ns="uniprot", id="P04637", source="test", is_canonical=True),
                dict(ns="genesymbol", id="TP53", source="test", is_canonical=False),
            ],
            annotations=[],
        ),
        dict(
            # the accession is not its canonical identifier, label or key
            entity_key="tp53-gene",
            entity_type="gene",
            namespace="ncbigene",
            identifier="7157",
            label="tumor protein p53 gene",
            identifiers=[
                dict(ns="ncbigene", id="7157", source="test", is_canonical=True),
                dict(ns="uniprot", id="P04637", source="test", is_canonical=False),
            ],
            annotations=[],
        ),
        dict(
            entity_key="egfr-protein",
            entity_type="protein",
            namespace="uniprot",
            identifier="P00533",
            label="EGFR",
            identifiers=[dict(ns="uniprot", id="P00533", source="test", is_canonical=True)],
            annotations=[],
        ),
    ]
    write_resource(folder, rows, [])
    return TestClient(create_app(queries=served_queries(data_root=tmp_path)))


def test_resolve_matches_each_identifier_exactly_and_in_order(tmp_path):
    client = _client(tmp_path)
    identifiers = ["P04637", " p04637 ", "uniprot:P04637", "P0463", "EGFR", "tp53-gene", "", "P04637"]
    response = client.post("/entities/resolve", json={"identifiers": identifiers})
    assert response.status_code == 200
    body = response.json()
    matches = {(m["identifier"], tuple(sorted(m["entityIds"]))) for m in body["matches"]}
    # one match per non-empty identifier, in order, repeats kept
    assert [m["identifier"] for m in body["matches"]] == [
        "P04637", "p04637", "uniprot:P04637", "P0463", "EGFR", "tp53-gene", "P04637"]
    both = ("ncbigene|7157", "uniprot|P04637")
    assert ("P04637", both) in matches  # also the gene, whose canonical id is not the accession
    assert ("p04637", both) in matches  # case-insensitive, trimmed
    assert ("uniprot:P04637", ("uniprot|P04637",)) in matches  # a CURIE keeps its namespace
    assert ("P0463", ()) in matches  # a prefix is not an identifier
    assert ("EGFR", ("uniprot|P00533",)) in matches
    assert ("tp53-gene", ("ncbigene|7157",)) in matches  # an entity key
    assert sorted(e["entityPk"] for e in body["entities"]) == ["egfr-protein", "tp53-gene", "tp53-protein"]


def test_resolve_falls_back_to_search_only_when_nothing_resolves(tmp_path):
    client = _client(tmp_path)
    body = client.post("/entities/resolve", json={"identifiers": ["tumor protein"]}).json()
    assert [m["identifier"] for m in body["matches"]] == ["tumor protein"]
    found = {e["entityPk"] for e in body["entities"]}
    assert found and set(body["matches"][0]["entityIds"]) == {
        f"{e['canonicalIdentifierType']}|{e['canonicalIdentifier']}" for e in body["entities"]}
    assert client.post("/entities/resolve", json={"identifiers": [" "]}).json() == {
        "entities": [], "matches": []}
