"""Network-independent integration tests: a mocked catalog and file host, real Parquet
written by the core fixtures and real DuckDB; the prebuilt queries are compared with
the API's queries (omnipath_api.queries) on the same files."""

import hashlib
import json

import duckdb
import httpx
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from omnipath_client import ClientError, available, connect, init
from omnipath_client.installation import DEFAULT_API_URL
from omnipath_core.fixtures import fixture_key as key, write_resource
from omnipath_core.schema import PUBLISHED_TABLES, SERVING_TABLES
from omnipath_api.queries import ReleaseQueries
from omnipath_api.release_database import serve_newest

TAGS = {
    "alpha": [("topic", "signaling"), ("molecule", "proteins")],
    "beta": [("topic", "metabolism"), ("molecule", "proteins")],
    "gamma": [("topic", "metabolism"), ("molecule", "metabolites")],
}


def entities():
    return [
        dict(
            entity_key=key(ident),
            identifier=ident,
            label=label,
            namespace=ns,
            entity_type=kind,
            taxon=taxon,
            identifiers=[dict(ns="hgnc", id="11998", source="test")] if label == "TP53" else [],
            reference_entity_key="entrez:7157" if label == "TP53" else None,
            gene_reference_keys=["entrez:7157"] if label == "TP53" else [],
        )
        for ident, label, ns, kind, taxon in [
            ("P04637", "TP53", "uniprot", "protein", "9606"),
            ("CHEBI:27732", "caffeine", "chebi", "small_molecule", None),
            ("quote", "O'Brien", "uniprot", "protein", "9606"),
        ]
    ]


def relations(resource):
    return [
        dict(
            relation_key=key(name),
            subject_entity_key=key(subject),
            subject_label=subject_label,
            object_entity_key=key(object_),
            object_label=object_label,
            predicate=predicate,
            category="interaction",
            taxon="9606",
            sources=[resource],
        )
        for name, subject, subject_label, object_, object_label, predicate in [
            ("r1", "P04637", "TP53", "CHEBI:27732", "caffeine", "affects"),
            ("r2", "quote", "O'Brien", "P04637", "TP53", "interacts_with"),
        ]
    ]


class Service:
    """A catalog API and a static file host over fixture resources."""

    def __init__(self, root):
        self.root, self.requests, self.versions = root, [], {}
        self.truncated, self.without_url = set(), set()
        self.taxonomy = None  # the release's taxon names (a Parquet file), when it has them
        for name in TAGS:
            self.publish(name, "1.0")

    def publish(self, name, version):
        folder = write_resource(
            self.root / name / version,
            entities(),
            relations(name) if name != "gamma" else [],
            [dict(relation_key=key("r1"), source=name, row_id="0", payload_json="{}")],
        )
        # Union-by-name must tolerate additive schema changes across resources.
        if name == "beta":
            table = pq.read_table(folder / "entity.parquet")
            pq.write_table(
                table.append_column("new_field", pa.array(["x"] * len(table))),
                folder / "entity.parquet",
            )
        self.versions[name] = version

    def record(self, name, evidence=False):
        version = self.versions[name]
        files = []
        for table in PUBLISHED_TABLES:
            if table == "evidence_payloads" and not evidence:
                continue
            path = self.root / name / version / f"{table}.parquet"
            url = f"https://data.example/resources/{name}/{version}/{table}.parquet"
            files.append(
                dict(
                    name=path.name,
                    url=None if url in self.without_url else url,
                    size_bytes=path.stat().st_size,
                    rows=pq.read_metadata(path).num_rows,
                )
            )
        tags = [dict(id=tag, dimension=dimension) for dimension, tag in TAGS[name]]
        return dict(resource_id=name, version=version, files=files, tags=tags, license="CC0")

    def respond(self, request):
        self.requests.append(request)
        url, path = str(request.url), request.url.path
        if request.url.host == "data.example":
            name, version, file = path.split("/")[2:]
            data = (self.root / name / version / file).read_bytes()
            return httpx.Response(200, content=data[:10] if url in self.truncated else data)
        if path == "/api/resources":
            references = {}
            if self.taxonomy is not None:
                data = self.taxonomy.read_bytes()
                version = hashlib.sha256(data).hexdigest()
                references["taxonomy"] = dict(version=version, size_bytes=len(data))
            resources = [self.record(n) for n in self.versions]
            return httpx.Response(200, json={"resources": resources, "references": references})
        if path == "/api/references/taxonomy":
            return httpx.Response(200, content=self.taxonomy.read_bytes())
        if path.startswith("/api/resources/"):
            return httpx.Response(200, json=self.record(path.split("/")[3], evidence=True))
        raise AssertionError(f"Unexpected request: {request.url}")

    def downloads(self):
        return [r.url.path.split("/", 2)[2] for r in self.requests if r.url.host == "data.example"]


@pytest.fixture
def service(tmp_path):
    service = Service(tmp_path / "server")
    with httpx.Client(transport=httpx.MockTransport(service.respond)) as http:
        service.http = http
        service.data = tmp_path / "data"

        def install(*args, **kwargs):
            return init(
                *args,
                data_dir=service.data,
                api_url="https://api.example/api",
                http_client=http,
                **kwargs,
            )

        service.init = install
        yield service


def names(client):
    return [r["resource_id"] for r in client.resources()]


def test_init_selects_by_name_and_tag(service):
    listed = available(api_url="https://api.example/api", http_client=service.http)
    assert [r["resource_id"] for r in listed] == ["alpha", "beta", "gamma"]
    service.init(topic="signaling", release="2026.09")
    assert service.requests[1].url.params["release"] == "2026.09"
    with connect(service.data) as db:
        assert names(db) == ["alpha"]
        assert db.release == "2026.09"
        assert db.resources()[0]["license"] == "CC0"
    # One tag of every given dimension; several tags of one dimension are alternatives.
    service.init(topic=["signaling", "metabolism"], molecule="metabolites")
    with connect(service.data) as db:
        assert names(db) == ["alpha", "gamma"]
    with pytest.raises(
        ValueError, match="Unknown topic.*available: \\['metabolism', 'signaling'\\]"
    ):
        service.init(topic="astrology")
    with pytest.raises(ValueError, match="Unknown resources"):
        service.init("missing")
    with pytest.raises(ValueError, match="at least one"):
        service.init([])
    with pytest.raises(ValueError, match="No resource matches"):
        service.init("alpha", topic="metabolism")
    service.init()
    with connect(service.data) as db:
        assert len(db.resources()) == 3


def test_init_is_incremental_and_follows_the_release(service):
    service.init("alpha")
    assert len(service.downloads()) == 7
    assert not any("evidence_payloads" in path for path in service.downloads())
    service.requests.clear()
    service.init("alpha")
    assert service.downloads() == []
    service.init("beta", include_evidence=True)
    assert {path.split("/")[0] for path in service.downloads()} == {"beta"}
    assert len(service.downloads()) == 8
    # Raw evidence alone is added to a complete version.
    service.requests.clear()
    service.init("alpha", include_evidence=True)
    assert service.downloads() == ["alpha/1.0/evidence_payloads.parquet"]
    # A new version replaces the old one, also for resources not named this time,
    # and keeps the evidence that was installed.
    service.publish("alpha", "2.0")
    service.requests.clear()
    service.init("beta")
    assert len(service.downloads()) == 8
    assert [p.name for p in (service.data / "resources/alpha").iterdir()] == ["2.0"]
    with connect(service.data) as db:
        assert {r["resource_id"]: r["version"] for r in db.resources()} == {
            "alpha": "2.0",
            "beta": "1.0",
        }
        assert db.sql("SELECT * FROM resources ORDER BY ALL").fetchall() == [
            ("alpha", "2.0", True),
            ("beta", "1.0", True),
        ]
        assert db.sql("SELECT count(*) FROM evidence_payloads WHERE _resource = 'alpha'").fetchone() == (1,)


def test_derived_serving_tables_match_the_build(service):
    service.init("alpha")
    for table in SERVING_TABLES:
        derived = pq.read_table(service.data / f"resources/alpha/1.0/{table}.parquet")
        built = pq.read_table(service.root / f"alpha/1.0/{table}.parquet")
        assert derived.schema.names == built.schema.names
        assert derived.num_rows

        def rows(table):
            return sorted(tuple(row.values()) for row in table.to_pylist())

        assert rows(derived) == rows(built)


def queries(db):
    """What a user sees through SQL: the same in both database modes."""
    return dict(
        tables=sorted(
            db.sql(
                "SELECT table_schema, table_name FROM information_schema.tables "
                "WHERE table_schema <> 'main' OR table_name NOT IN ('database_info') "
                "AND NOT starts_with(table_name, 'op_')"
            ).fetchall()
        ),
        entities=db.sql("SELECT label, resources FROM entity ORDER BY ALL").fetchall(),
        quoted=db.sql("SELECT identifier FROM entity WHERE label = ?", ["O'Brien"]).fetchall(),
        one=db.sql("SELECT count(*) FROM relation WHERE list_contains(resources, 'alpha')").fetchone(),
        joined=sorted(
            db.sql(
                """SELECT e.label, unnest(r.resources) FROM entity e
                JOIN relation r ON r.subject_n = e.entity_n WHERE e.identifier = ?""",
                ["P04637"],
            ).fetchall()
        ),
        endpoints=db.sql(
            """SELECT count(*) FROM relation_endpoint p JOIN entity e USING (entity_n)
            WHERE e.entity_key = ? AND list_contains(p.resources, 'alpha')""",
            [key("P04637")],
        ).fetchone(),
        terms=db.sql("SELECT count(*) FROM entity_term WHERE term = 'tp53'").fetchone(),
    )


def test_sql_on_native_tables(service, monkeypatch):
    service.init(["alpha", "beta"])
    database = service.data / "omnipath.duckdb"

    def forbidden(*args, **kwargs):
        raise AssertionError("A query attempted HTTP")

    with monkeypatch.context() as patched:
        patched.setattr(httpx.Client, "get", forbidden)
        patched.setattr(httpx.Client, "stream", forbidden)
        with connect(service.data) as db:
            assert db.path == database
            tables = queries(db)
            # One merged table each, no tables per resource.
            assert ("main", "relation") in tables["tables"] and ("main", "entity_term") in tables["tables"]
            assert {schema for schema, _ in tables["tables"]} == {"main"}
            assert ("main", "evidence_payloads") not in tables["tables"]
            assert tables["entities"] == [
                ("O'Brien", ["alpha", "beta"]), ("TP53", ["alpha", "beta"]), ("caffeine", ["alpha", "beta"])
            ]
            assert tables["quoted"] == [("quote",)]
            assert tables["one"] == (2,)
            assert tables["joined"] == [("TP53", "alpha"), ("TP53", "beta")]
            assert tables["endpoints"] == (2,) and tables["terms"] == (1,)
            # The database is read-only; temporary tables are the session's own.
            with pytest.raises(duckdb.Error):
                db.sql("CREATE TABLE mine AS SELECT 1")
            db.connection.execute("CREATE TEMP TABLE mine AS SELECT label FROM entity")
            assert db.sql("SELECT count(*) FROM mine").fetchone() == (3,)
            types = db.sql(
                "SELECT table_type FROM information_schema.tables "
                "WHERE table_name = 'entity' AND table_schema = 'main'"
            ).fetchall()
            assert types == [("BASE TABLE",)]
            # The database format's indexes.
            assert sorted(r[0] for r in db.sql("SELECT index_name FROM duckdb_indexes()").fetchall()) == [
                "entity_annotation_n", "entity_evidence_n", "entity_key", "entity_n",
                "relation_evidence_n", "relation_n",
            ]
    # A new version and a new resource: the database is written anew.
    service.publish("alpha", "2.0")
    before = database.stat().st_ino
    service.init("gamma")
    assert database.stat().st_ino != before
    with connect(service.data) as db:
        assert db.sql(
            "SELECT r, count(*) FROM (SELECT unnest(resources) AS r FROM entity) GROUP BY ALL ORDER BY ALL"
        ).fetchall() == [("alpha", 3), ("beta", 3), ("gamma", 3)]
        assert db.sql("SELECT * FROM resources ORDER BY ALL").fetchall() == [
            ("alpha", "2.0", False),
            ("beta", "1.0", False),
            ("gamma", "1.0", False),
        ]
    # Raw evidence for one resource.
    service.init("gamma", include_evidence=True)
    with connect(service.data) as db:
        assert db.sql("SELECT _resource FROM evidence_payloads").fetchall() == [("gamma",)]
        assert db.sql("SELECT count(*) FROM entity").fetchone() == (3,)
        assert json.loads((service.data / "installed.json").read_text())["materialize"] is True
    # Nothing changed: the database stays.
    before = database.stat().st_ino
    service.init()
    assert database.stat().st_ino == before


def test_prebuilt_queries_answer_like_the_api(service):
    service.init(["alpha", "beta"], include_evidence=True)
    with connect(service.data) as db:
        found = db.search("TP53")
        assert found.columns == [
            "entity_key",
            "label",
            "entity_type",
            "namespace",
            "identifier",
            "taxon",
            "relation_count",
            "reference_entity_key",
            "gene_reference_keys",
        ]
        # One merged entity for the two resource copies, as the explorer shows it.
        row = dict(zip(found.columns, found.fetchone()))
        assert found.count("*").fetchone() == (1,)
        assert row["entity_key"] == key("P04637") and row["label"] == "TP53"
        assert row["reference_entity_key"] == "entrez:7157"
        assert db.search("11998").project("label").fetchall() == [("TP53",)]
        assert db.search("TP53", entity_types="small_molecule").fetchall() == []
        assert db.search("no such thing").fetchall() == []

        relations_ = db.relations("TP53")
        assert {"relation_key", "annotations", "evidence", "sources"} <= set(relations_.columns)
        # The two resources state the same two relations: one row each, with both sources.
        assert relations_.count("*").fetchone() == (2,)
        assert db.relations("TP53", sources="alpha", predicates=["affects"]).project(
            "subject_label, object_label, sources"
        ).fetchall() == [("TP53", "caffeine", ["alpha", "beta"])]
        with pytest.raises(ValueError, match="unavailable"):
            db.relations("TP53", sources="gamma")
        with pytest.raises(ValueError, match="Unknown Biolink predicate: nonsense"):
            db.relations("TP53", predicates="nonsense")
        assert db.relations(["caffeine", "O'Brien"]).count("*").fetchone() == (2,)
        assert db.relations(key("quote")).count("*").fetchone() == (1,)
        assert db.relations("gene:entrez:7157").count("*").fetchone() == (2,)
        assert db.relations(limit=1).count("*").fetchone() == (1,)
        # Results are DuckDB relations like any other.
        assert db.connection.sql(
            "SELECT count(*) FROM found JOIN relations_ ON entity_key = subject_entity_key"
        ).fetchone() == (1,)
        with pytest.raises(TypeError, match="Unknown filters: \\['source'\\]"):
            db.relations("TP53", source="alpha")
        with pytest.raises(TypeError, match="once"):
            db.relations("TP53", entity_ids=["TP53"])
        with pytest.raises(ValueError, match="limit"):
            db.relations(limit=100_001)

        details = db.entity(key("P04637"))
        assert details["entity"]["label"] == "TP53"
        assert db.entity("TP53") is None
        # The queries by name.
        assert db.query("relations", entities=["TP53"]).count("*").fetchone() == (2,)
        assert db.sql(
            "SELECT count(*) FROM evidence_payloads WHERE relation_key = ?", [key("r1")]
        ).fetchone() == (2,)
    # Without raw evidence the queries work as well.
    service.data = service.data.parent / "plain"
    service.init("alpha")
    with connect(service.data) as db:
        assert db.search("TP53").count("*").fetchone() == (1,)
        assert db.relations("TP53").count("*").fetchone() == (2,)


def test_failed_download_installs_nothing(service):
    service.truncated.add("https://data.example/resources/beta/1.0/relation.parquet")
    with pytest.raises(ClientError, match="Truncated"):
        service.init(["alpha", "beta"])
    # The resource finished before the failure stays, and the database describes it;
    # the failed one leaves no files.
    assert [p.name for p in (service.data / "resources/beta").iterdir()] == []
    with connect(service.data) as db:
        assert names(db) == ["alpha"]
        assert db.sql("SELECT resource_id FROM resources").fetchall() == [("alpha",)]
    service.truncated.clear()
    service.init("beta")
    with connect(service.data) as db:
        assert db.sql("SELECT count(DISTINCT r) FROM (SELECT unnest(resources) AS r FROM relation)").fetchone() == (2,)


def test_missing_data_and_damaged_files(service, tmp_path):
    with pytest.raises(ClientError, match="run omnipath_client.init"):
        connect(tmp_path / "empty")
    service.init("alpha")
    path = service.data / "resources/alpha/1.0/entity.parquet"
    original = path.read_bytes()
    path.write_bytes(original[:-1])
    with pytest.raises(ClientError, match="Missing or truncated"):
        connect(service.data)
    # init repairs it.
    service.init("alpha")
    assert path.read_bytes() == original
    (service.data / "resources/alpha/1.0/entity_term.parquet").unlink()
    with pytest.raises(ClientError, match="Missing or truncated"):
        connect(service.data)
    service.init("alpha")
    # A database that does not describe the installation is written again on connect.
    (service.data / "omnipath.duckdb").unlink()
    with connect(service.data) as db:
        assert db.sql("SELECT count(*) FROM entity").fetchone() == (3,)
    manifest = service.data / "installed.json"
    data = json.loads(manifest.read_text())
    data["resources"][0]["resource_id"] = "../escape"
    manifest.write_text(json.dumps(data))
    with pytest.raises(ClientError, match="Invalid resource"):
        connect(service.data)


def test_lifecycle(service):
    service.init("alpha")
    db = connect(service.data)
    relation = db.sql("SELECT * FROM entity")
    assert db.search("TP53").fetchall()
    db.close()
    db.close()
    for use in (lambda: db.sql("SELECT 1"), lambda: db.entity(key("P04637")), lambda: db.search("TP53")):
        with pytest.raises(ClientError, match="closed"):
            use()
    with pytest.raises(duckdb.Error):
        relation.fetchall()
    # Two sessions read the same database at once.
    with connect(service.data) as one, connect(service.data) as two:
        assert one.sql("SELECT count(*) FROM entity").fetchone() == (3,)
        assert two.sql("SELECT count(*) FROM entity").fetchone() == (3,)


def test_missing_url_and_http_failure(service, tmp_path):
    service.without_url.add("https://data.example/resources/alpha/1.0/entity.parquet")
    with pytest.raises(ClientError, match="static download URL"):
        service.init("alpha")
    assert not (service.data / "installed.json").exists()
    assert not (service.data / "omnipath.duckdb").exists()
    with httpx.Client(transport=httpx.MockTransport(lambda _: httpx.Response(503))) as http:
        with pytest.raises(httpx.HTTPStatusError):
            init(data_dir=tmp_path / "failed", http_client=http)
    assert not (tmp_path / "failed").exists()


def test_data_directory_survives_move_and_rejects_symlinks(service, tmp_path):
    service.init("alpha")
    moved = tmp_path / "moved"
    service.data.rename(moved)
    # The views named the old location; connecting writes them for the new one.
    with connect(moved) as db:
        assert db.sql("SELECT count(*) FROM entity").fetchone() == (3,)
    artifact = moved / "resources/alpha/1.0/entity.parquet"
    outside = tmp_path / "outside.parquet"
    artifact.rename(outside)
    artifact.symlink_to(outside)
    with pytest.raises(ClientError, match="escapes"):
        connect(moved)


def test_locations_fall_back_to_environment(service, monkeypatch, tmp_path):
    seen = []

    def respond(request):
        seen.append(str(request.url).split("?")[0])
        return httpx.Response(200, json={"resources": []})

    with httpx.Client(transport=httpx.MockTransport(respond)) as http:
        monkeypatch.setenv("OMNIPATH_API_URL", "https://example.invalid/api/")
        available(http_client=http)
        available(api_url="https://explicit.invalid/api", http_client=http)
        monkeypatch.delenv("OMNIPATH_API_URL")
        available(http_client=http)
    assert seen == [
        "https://example.invalid/api/resources",
        "https://explicit.invalid/api/resources",
        f"{DEFAULT_API_URL}/resources",
    ]
    service.init("alpha")
    monkeypatch.setenv("OMNIPATH_DATA_DIR", str(service.data))
    with connect() as db:
        assert db.data_dir == service.data


def api_search(api, query, **filters):
    """The API's search in the client's columns."""
    found = api.search_entities_api(query, filters=filters, limit=20)["entities"]
    return [
        (e["entityPk"], e.get("label"), e.get("entityType"),
         getattr(e.get("canonicalIdentifierType"), "value", e.get("canonicalIdentifierType")),
         e.get("canonicalIdentifier"), e.get("taxonomyId"), e.get("relationCount"),
         e.get("referenceEntityKey"), e.get("geneReferenceKeys") or [])
        for e in found
    ]


def canonical(table):
    """Rows comparable as a multiset. Product records carry the entity columns of every
    installed resource, also under a source filter, so empty fields are left out."""
    rows = []
    for row in table.to_pylist():
        row["referenced_product_records"] = [
            {k: v for k, v in record.items() if v is not None}
            for record in row["referenced_product_records"]
        ]
        rows.append(json.dumps(row, sort_keys=True, default=str))
    return sorted(rows)


def test_prebuilt_queries_match_the_engine(service):
    """The macros give the API's answers."""
    service.init(["alpha", "beta", "gamma"], include_evidence=True)
    # The API's engine over a release of the same files.
    serve_newest(service.data, memory_limit="1GB", threads=1)
    api = ReleaseQueries(service.data)
    try:
        with connect(service.data) as db:
            for filters in [
                {"entity_ids": ["TP53"]},
                {"entity_ids": ["TP53"], "sources": ["alpha"]},
                {"entity_ids": ["TP53"], "gene_mode": "source"},
                {"entity_ids": ["caffeine", "O'Brien"], "gene_mode": "both"},
                {"entity_ids": ["gene:entrez:7157"]},
                {"entity_ids": ["uniprot:P04637"]},
                {"reference_entity_keys": ["entrez:7157"]},
                {"predicates": ["biolink:interacts_with"]},
                {"interaction_types": ["protein"], "partner_types": ["small_molecule"]},
                {"taxonomy_ids": ["9606"], "relation_categories": ["interaction"]},
                {"sources": ["beta"]},
                {},
            ]:
                data, _ = api.export_slice(filters, None)
                expected = pq.read_table(pa.BufferReader(data))
                rest = dict(filters)
                actual = db.relations(rest.pop("entity_ids", None), **rest).to_arrow_table()
                assert actual.schema.names == expected.schema.names, filters
                assert canonical(actual) == canonical(expected), filters
            for query, filters in [
                ("TP53", {}), ("tp", {}), ("caffeine", {}), ("CHEBI:27732", {}),
                ("uniprot:P04637", {}), ("11998", {}), ("brien", {}), ("", {}),
                (key("P04637"), {}), ("TP53", {"sources": ["beta"]}),
                ("TP53", {"entity_types": ["protein"]}), ("TP53", {"taxonomy_ids": ["9606"]}),
            ]:
                got = [tuple(row) for row in db.search(query, **filters).fetchall()]
                assert got == api_search(api, query, **filters), (query, filters)
            for ident in ("P04637", "CHEBI:27732", "quote"):
                want = json.loads(json.dumps(api.get_entity_details(key(ident)), default=str))
                assert db.entity(key(ident)) == want, ident
            assert db.entity("0" * 64) is None
    finally:
        api.close()


def test_api_mode_answers_like_local_mode(service):
    """The hosted API (here in-process over the same files) gives the local answers."""
    from fastapi.testclient import TestClient
    from omnipath_api.server import create_app

    service.init(["alpha", "beta", "gamma"], include_evidence=True)
    serve_newest(service.data, memory_limit="1GB", threads=1)
    engine = ReleaseQueries(service.data)
    try:
        with TestClient(create_app(queries=engine)) as http, connect(service.data) as local, connect(
            mode="api", api_url="http://testserver", http_client=http
        ) as remote:
            assert (local.mode, remote.mode) == ("local", "api")
            for filters in [
                {"entity_ids": ["TP53"]},
                {"entity_ids": ["TP53"], "sources": ["alpha"]},
                {"entity_ids": ["TP53"], "gene_mode": "source"},
                {"entity_ids": ["gene:entrez:7157"]},
                {"reference_entity_keys": ["entrez:7157"]},
                {"predicates": ["biolink:interacts_with"]},
                {"sources": ["beta"]},
                {},
            ]:
                rest = dict(filters)
                entities = rest.pop("entity_ids", None)
                want = local.relations(entities, **rest).to_arrow_table()
                got = remote.relations(entities, **rest).to_arrow_table()
                assert got.schema.names == want.schema.names, filters
                assert canonical(got) == canonical(want), filters
            for query, filters in [("TP53", {}), ("tp", {}), ("caffeine", {}), ("", {}),
                                   ("TP53", {"sources": ["beta"]}), ("TP53", {"entity_types": ["protein"]})]:
                assert remote.search(query, **filters).fetchall() == local.search(query, **filters).fetchall()
            for ident in ("P04637", "CHEBI:27732", "quote"):
                assert remote.entity(key(ident)) == local.entity(key(ident)), ident
            assert remote.entity("0" * 64) is None
            with pytest.raises(TypeError, match="Unknown filters"):
                remote.relations("TP53", no_such_filter=1)
            with pytest.raises(ClientError, match="local installation"):
                remote.sql("SELECT 1")
            assert {r["resource_id"] for r in remote.resources()} == {"alpha", "beta", "gamma"}
    finally:
        engine.close()


def write_taxonomy(path, name):
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(
        pa.table(dict(
            taxon_id=["9606", "10090"],
            current_taxon_id=["9606", "10090"],
            scientific_name=["Homo sapiens", "Mus musculus"],
            common_name=[name, "house mouse"],
        )),
        path,
    )
    return path


def test_init_downloads_the_release_taxon_names(service):
    service.init("alpha")
    with connect(service.data) as db:
        assert db.sql("SELECT count(*) FROM op_taxonomy").fetchone()[0] == 0
    service.taxonomy = write_taxonomy(service.root / "taxonomy.parquet", "human")
    service.init("alpha")
    tp53 = key("P04637")
    with connect(service.data) as db:
        assert db.entity(tp53)["entity"]["taxonomyName"] == "human"
    downloads = [r for r in service.requests if r.url.path == "/api/references/taxonomy"]
    assert len(downloads) == 1
    # The same names are not downloaded or written again; new ones write the database
    # anew, since the merged entities' summaries name taxa.
    written = (service.data / "omnipath.duckdb").stat().st_ino
    service.init("alpha")
    assert len([r for r in service.requests if r.url.path == "/api/references/taxonomy"]) == 1
    assert (service.data / "omnipath.duckdb").stat().st_ino == written
    write_taxonomy(service.taxonomy, "people")
    service.init("alpha")
    with connect(service.data) as db:
        assert db.entity(tp53)["entity"]["taxonomyName"] == "people"
        assert db.search("TP53").fetchone() is not None
    assert len(list((service.data / "references/taxonomy").iterdir())) == 1


def test_query_layer_of_another_client_version_is_rewritten(service):
    service.init("alpha")
    database = service.data / "omnipath.duckdb"
    with duckdb.connect(str(database)) as db:
        db.execute("UPDATE database_info SET query_layer = 'older'")
        db.execute("DROP MACRO op_entity_ok")
    before = database.stat().st_ino
    with connect(service.data) as db:
        assert db.search("TP53").count("*").fetchone() == (1,)
        assert db.sql("SELECT query_layer <> 'older' FROM database_info").fetchone() == (True,)
    # The tables stayed: only the query layer was written.
    assert database.stat().st_ino == before


def test_shared_contract_is_up_to_date():
    import importlib.util
    from pathlib import Path

    path = Path(__file__).resolve().parents[2] / "scripts" / "generate_contract.py"
    spec = importlib.util.spec_from_file_location("generate_contract", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.main(check=True) == 0, "run packages/query/scripts/generate_contract.py"


