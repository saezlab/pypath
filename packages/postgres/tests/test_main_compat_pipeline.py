"""Durability contract around faithful main product helpers."""

from types import SimpleNamespace

import pytest

from omnipath_subsets import scientific as pipeline


class Cursor:
    def __init__(self, conn):
        self.connection = conn

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return None

    def execute(self, query, params=None):
        self.connection.queries.append((query, params))


class Connection:
    def __init__(self):
        self.pending = []
        self.durable = ["base", "metsigdb"]
        self.commit_count = 0
        self.rollback_count = 0
        self.queries = []

    def cursor(self):
        return Cursor(self)

    def commit(self):
        self.durable.extend(self.pending)
        self.pending.clear()
        self.commit_count += 1

    def rollback(self):
        self.pending.clear()
        self.rollback_count += 1


@pytest.mark.parametrize("product", ["network_views", "cosmos"])
def test_main_internal_commits_wait_for_complete_product_checkpoint(monkeypatch, product):
    conn = Connection()

    def build(facade, *args, **kwargs):
        # The main helpers historically commit at intermediate publication
        # boundaries. Both writes must survive or roll back as one product.
        facade.pending.append(f"{product}.tables")
        facade.commit()
        facade.pending.append(f"{product}.indexes")
        facade.commit()
        if product == "cosmos":
            assert kwargs["use_published_identifiers"] is True
            assert kwargs["utils_db_url"] is None
        return {"rows": 4}

    if product == "network_views":
        monkeypatch.setattr(pipeline.network_views, "apply_all", build)
    else:
        monkeypatch.setattr(pipeline.cosmos, "build_cosmos_projection", build)
    result = pipeline.run_product(conn, product, schema="scratch")
    assert result["result"] == {"rows": 4}
    assert conn.commit_count == 0
    assert conn.durable == ["base", "metsigdb"]
    conn.pending.append(f"{product}.metadata")
    conn.commit()
    assert conn.durable == [
        "base",
        "metsigdb",
        f"{product}.tables",
        f"{product}.indexes",
        f"{product}.metadata",
    ]


def test_cancellation_keeps_prior_products_and_rolls_back_current_product(monkeypatch):
    conn = Connection()

    def build(facade, *args, **kwargs):
        facade.pending.append("cosmos.edges")
        facade.commit()
        facade.pending.append("cosmos.labels")
        raise KeyboardInterrupt

    monkeypatch.setattr(pipeline.cosmos, "build_cosmos_projection", build)
    with pytest.raises(KeyboardInterrupt):
        pipeline.run_product(conn, "cosmos", schema="scratch")
    conn.rollback()
    assert conn.durable == ["base", "metsigdb"]
    assert conn.pending == []
    assert conn.commit_count == 0
    assert conn.rollback_count == 1


def test_commit_facade_keeps_real_rollback_and_cursor_connection():
    conn = Connection()
    facade = pipeline.CommitDeferredConnection(conn)
    facade.pending.append("unfinished")
    facade.commit()
    assert conn.pending == ["unfinished"]
    with facade.cursor() as cur:
        assert cur.connection is conn
    facade.rollback()
    assert conn.pending == []
    assert conn.rollback_count == 1


def test_shared_pipeline_keeps_main_order_without_online_resolution(monkeypatch):
    from omnipath_postgres.relational import pipeline

    calls = []
    functions = {
        "populate_identifier_authority": "identifier_authority",
        "rebuild_derived_tables": "derived_tables",
        "rebuild_chemical_resolution_levels": "chemical_resolution_levels",
        "rebuild_chemical_ambiguous_name_candidates": "chemical_ambiguous_name_candidates",
        "classify_chemical_class": "chemical_class",
        "classify_interaction_class": "interaction_class",
        "rebuild_interaction_tables": "interactions",
        "populate_entity_labels": "entity_labels",
        "populate_entity_name": "entity_name",
        "populate_chemical_labels": "chemical_labels",
        "classify_metabolic_domain": "metabolic_domain",
        "rebuild_bitmap_tables": "bitmaps",
        "rebuild_resource_overlap_summary": "resource_overlap",
        "sync_resources_table": "resources",
        "sync_data_source_licenses": "licenses",
        "emit_build_manifest": "build_manifest",
    }

    def step(name):
        def call(*args, **kwargs):
            calls.append(name)
            if name == "derived_tables":
                assert kwargs["interactions"] is False
                assert kwargs["use_external_mappings"] is False
            if name == "build_manifest":
                assert kwargs["utils_db_url"] is None
            if name == "interactions":
                return pipeline.InteractionDeriveStats(step_seconds={"interaction_header": 1.5})
            if name == "build_manifest":
                assert kwargs["derive_cost"]["interaction_header"]["seconds"] == 1.5
                assert kwargs["derive_cost"]["interaction_header"]["rows"] == 0
                assert kwargs["deferral_cost"] is None
            return 0 if name == "identifier_authority" else {}

        return call

    for function, name in functions.items():
        monkeypatch.setattr(pipeline, function, step(name))
    monkeypatch.setattr(
        pipeline,
        "_metadata_inputs",
        lambda records: {
            "reactome": [SimpleNamespace(call=SimpleNamespace(config=SimpleNamespace(mints=[])))],
        },
    )
    monkeypatch.setattr(pipeline, "_availability", lambda conn, schema: [])
    result = pipeline.run_main_derivations(
        Connection(), schema="scratch", published_metadata={"reactome": {}}
    )
    assert calls == list(functions.values())
    assert list(result["phase_seconds"]) == calls
    assert result["results"]["identifier_authority"]["declarations_available"] is False
