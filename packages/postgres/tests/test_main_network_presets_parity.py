"""Exact main preset registry/build oracle, without a new serving layer.

Main registers metadata; its API consumer executes scopes and collapse. These
checks compare complete descriptors and stored rows, exercise literal loaded
resource names and legacy migration, and ensure apply/refresh preserve bounded
human/mouse, signed/unsigned and known/unknown-license fact data. They do not
claim scientific execution coverage for a preset query that main does not ship.
"""

from __future__ import annotations

from collections import Counter
from contextlib import closing
from dataclasses import asdict, replace
from datetime import datetime
import ast
import importlib
import json
from pathlib import Path
import uuid

import psycopg2
import pytest
from oracles import read_legacy

import test_main_parity_contract as reference_contract

MAIN = reference_contract.MAIN
ORACLE_PREFIX = reference_contract.ORACLE_PREFIX
_catalogue = reference_contract._catalogue
oracle = reference_contract.oracle
port = reference_contract.port
paired_schemas = reference_contract.paired_schemas


@pytest.fixture(scope="module")
def parity_connection(postgres_dsn):
    # Conftest owns a disposable local cluster. This file does not inspect an
    # externally supplied parity DSN or connect to a deployment.
    with closing(psycopg2.connect(postgres_dsn)) as connection:
        try:
            yield connection
        finally:
            connection.rollback()


FROZEN_DIGESTS = {
    "__init__.py": "d3f0b6f4791e31eb8eeabeccce26905e163e25b024e3ac2ceecc5769cdb7f034",
    "_definitions.py": "f714d8a8938e52d5cdcb1415ade551f99282954b4d8194b65cd16e0da1c854d4",
    "_framework.py": "11468b717324ebac81cb01b113dd8f540517ce4c93fec5fe35dd56ddbd13688e",
    "sql/metalinksdb.sql": "fdd6628666246e55b074a9d25da7c83d8514c259bf1df9b99b6d185eb8d5b135",
    "sql/metalinksdb_annotations.sql": "435b95b724c89ce2f75d6289a0510ed4b51bda70d638fce81c6aadbca3bf4437",
}


def _modules():
    return (
        importlib.import_module(ORACLE_PREFIX + ".network_views"),
        importlib.import_module("omnipath_subsets.network_views"),
    )


# The partner and predicate gates are this port's own: one chemical end let a
# resource's annotation links of its chemicals into the dataset (2026-10-19).
_PORT_ONLY_FILTERS = {"partner_types", "predicates"}
# Resources this port adds to a preset; main has no module for them (2026-10-25).
_PORT_ONLY_SOURCES = {"metachatdb"}


def _normalize(value):
    """Only the declared concrete chemical type gate may differ."""
    if isinstance(value, dict):
        result = {}
        for key, item in value.items():
            if key in _PORT_ONLY_FILTERS:
                continue
            if key in ("included_sources", "resources") and isinstance(item, (tuple, list)):
                result[key] = [s for s in item if s not in _PORT_ONLY_SOURCES]
                continue
            if key == "entity_types" and item in (
                ["Chemical:OM:0037"],
                ["chemical_entity", "small_molecule"],
            ):
                result[key] = ["published_chemical_gate"]
            elif key == "entity_type_gate" and item in ("Chemical:OM:0037", "chemical_entity"):
                result[key] = "published_chemical_gate"
            else:
                result[key] = _normalize(item)
        return result
    if isinstance(value, (tuple, list)):
        return [_normalize(item) for item in value]
    return value


def test_preset_reference_framework_and_all_descriptors_match_main(oracle, port):
    for relative, digest in FROZEN_DIGESTS.items():
        read_legacy("network_views/" + relative, digest)
    # This implementation needs no vocabulary adaptation at all.
    from omnipath_subsets import network_views

    current_framework = Path(network_views.__file__).parent / "_framework.py"
    assert ast.dump(ast.parse(current_framework.read_text())) == ast.dump(
        ast.parse(read_legacy("network_views/_framework.py"))
    )
    expected, actual = _modules()
    assert [_normalize(asdict(d)) for d in expected.NETWORKS] == [
        _normalize(asdict(d)) for d in actual.NETWORKS
    ]
    assert [d.name for d in actual.NETWORKS] == ["metalinksdb", "liana", "reactions"]
    assert actual.METALINKSDB.composition["components"][0]["parameters"]["filters"][
        "entity_types"
    ] == ["chemical_entity", "small_molecule"]
    assert actual.METALINKSDB.curation["entity_type_gate"] == "chemical_entity"


def test_main_exports_registry_builds_and_no_preset_query_or_fold(oracle):
    expected, actual = _modules()
    assert expected.__all__ == actual.__all__
    assert not {"query", "iter_records", "collapse", "fold"} & set(expected.__all__)
    for module in (expected, actual):
        assert all(
            d.is_preset and not d.sql_files and not d.matviews and d.sql_text() == ""
            for d in module.NETWORKS
        )
        assert all(d.schema is None and d.combined_relation is None for d in module.NETWORKS)
        assert all(d.license_scope is None and d.evidence_scope is None for d in module.NETWORKS)
        assert [d.collapse_mode for d in module.NETWORKS] == ["endpoints", "endpoints", "none"]
        assert [d.grain for d in module.NETWORKS] == ["interaction", "interaction", "participant"]
        assert module.LIANA.included_sources == ("connectomedb2025",)
        assert module.REACTIONS.included_sources == ("kegg", "rhea", "metatlas", "recon3d")
        assert "reactome" not in module.REACTIONS.included_sources
        assert "humangem" not in module.METALINKSDB.included_sources
        assert module.METALINKSDB.labels["resources"] == {"metatlas": "humangem"}
        assert [s["operation"] for s in module.METALINKSDB.composition["steps"]] == [
            "exclude",
            "collapse",
            "annotate",
        ]
        # Organisms are a serving dimension. The current preset supplies no
        # human-only gate, rather than storing a hidden 9606 restriction.
        assert "organism_scope" not in module.NetworkDefinition.__dataclass_fields__
        assert "9606" not in json.dumps([asdict(d) for d in module.NETWORKS])
        assert "api-service executes the algebra" in module.NetworkDefinition.__doc__


def _registry_rows(conn, schema):
    from psycopg2 import sql

    with conn.cursor() as cur:
        cur.execute(
            sql.SQL("SELECT * FROM {}.network_registry ORDER BY name").format(
                sql.Identifier(schema)
            )
        )
        names = [d.name for d in cur.description]
        rows, stamps = {}, {}
        for values in cur.fetchall():
            row = dict(zip(names, values, strict=True))
            stamp = row.pop("built_at")
            assert isinstance(stamp, datetime) and stamp.tzinfo is not None
            stamps[row["name"]] = stamp
            rows[row["name"]] = _normalize(row)
        return rows, stamps


def _registry_ddl(conn, schema):
    catalogue = _catalogue(conn, schema)
    return {
        grain: Counter(
            {row: count for row, count in catalogue[grain].items() if row[0] == "network_registry"}
        )
        for grain in ("columns", "constraints", "indexes")
    }


def _base_rows(conn, schema):
    from psycopg2 import sql

    with conn.cursor() as cur:
        result = {}
        for table in ("data_source", "data_source_license", "interaction_fact_resource"):
            cur.execute(
                sql.SQL("SELECT * FROM {}.{}").format(sql.Identifier(schema), sql.Identifier(table))
            )
            result[table] = Counter(
                json.dumps(row, sort_keys=True, default=str) for row in cur.fetchall()
            )
        cur.execute(
            "SELECT c.relname,c.relkind FROM pg_class c JOIN pg_namespace n ON n.oid=c.relnamespace WHERE n.nspname=%s ORDER BY 1",
            [schema],
        )
        result["objects"] = Counter(
            row for row in cur.fetchall() if not row[0].startswith("network_registry")
        )
        return result


def test_legacy_registry_migration_defaults_and_data_match_main(paired_schemas, oracle):
    from psycopg2 import sql

    conn, schemas = paired_schemas
    expected, actual = _modules()
    try:
        for schema in schemas:
            with conn.cursor() as cur:
                cur.execute(
                    sql.SQL("""
                    CREATE TABLE {}.network_registry (
                      name text PRIMARY KEY,kind text NOT NULL,
                      schema_name text NOT NULL,combined_relation text NOT NULL,
                      included_sources text[] NOT NULL,
                      built_at timestamptz NOT NULL DEFAULT now())
                """).format(sql.Identifier(schema))
                )
                cur.execute(
                    sql.SQL(
                        "INSERT INTO {}.network_registry(name,kind,schema_name,combined_relation,included_sources) VALUES ('retired','legacy','retired_network','relations',ARRAY['connectomedb','humangem'])"
                    ).format(sql.Identifier(schema))
                )
            conn.commit()
        for module, schema in zip((expected, actual), schemas, strict=True):
            module.ensure_network_registry(conn, registry_schema=schema)
        rows = [_registry_rows(conn, schema)[0] for schema in schemas]
        assert rows[0] == rows[1]
        row = rows[0]["retired"]
        assert row["schema_name"] == "retired_network" and row["combined_relation"] == "relations"
        assert row["included_sources"] == ["connectomedb", "humangem"]
        assert row["collapse_mode"] == "endpoints" and row["grain"] == "interaction"
        assert row["composition"] is None and row["license_scope"] is None
        assert row["evidence_scope"] is None and row["default_attributes"] is None
        assert _registry_ddl(conn, schemas[0]) == _registry_ddl(conn, schemas[1])
    finally:
        conn.rollback()
        with conn.cursor() as cur:
            for schema in schemas:
                cur.execute(
                    sql.SQL("DROP TABLE IF EXISTS {}.network_registry").format(
                        sql.Identifier(schema)
                    )
                )
        conn.commit()


def _substrate(conn, schema, module, *, biolink):
    from psycopg2 import sql

    prefix = "omnipath_postgres.relational" if biolink else module.__name__.rsplit(".", 1)[0]
    # Seed the real append-only class vocabulary before inserting finished
    # facts; no relation classification or resource parsing is synthesized.
    classifier = importlib.import_module(prefix + ".classify.interaction_class")
    classifier.classify_interaction_class(conn, schema=schema)

    def q(text):
        return sql.SQL(text).format(s=sql.Identifier(schema))

    # Actual main slugs are present except the two exact-name checks below.
    loaded = [n for n in module.METALINKSDB.included_sources if n != "metatlas"] + [
        "kegg",
        "connectomedb",
        "humangem",
    ]
    sources = {name: 610001 + i for i, name in enumerate(loaded)}
    with conn.cursor() as cur:
        for name, sid in sources.items():
            cur.execute(
                q("INSERT INTO {s}.data_source(source_id,name) VALUES (%s,%s)"), [sid, name]
            )
            # Unknown even though the levels are maximally permissive; these
            # metadata rows must survive registration without reinterpretation.
            cur.execute(
                q(
                    "INSERT INTO {s}.data_source_license(source_id,license_name,purpose_level,sharing_level,attrib_level,is_known) VALUES (%s,%s,%s,25,10,%s)"
                ),
                [
                    sid,
                    "fixture_unknown" if name in {"bindingdb", "connectomedb"} else "fixture_known",
                    25 if name in {"bindingdb", "connectomedb"} else 20,
                    name not in {"bindingdb", "connectomedb"},
                ],
            )
        for name in (
            "protein" if biolink else "Protein:MI:0326",
            "chemical_entity" if biolink else "Chemical:OM:0037",
        ):
            cur.execute(
                q(
                    "INSERT INTO {s}.vocab_entity_type(name) VALUES (%s) ON CONFLICT(name) DO NOTHING"
                ),
                [name],
            )
        cur.execute(q("SELECT name,entity_type_id FROM {s}.vocab_entity_type"))
        types = dict(cur.fetchall())
        entities = {
            name: str(uuid.UUID(int=70000 + i))
            for i, name in enumerate(
                (
                    "ligand_h",
                    "receptor_h",
                    "ligand_m",
                    "receptor_m",
                    "chemical",
                    "protein_h",
                    "protein_m",
                )
            )
        }
        for name, eid in entities.items():
            if name == "chemical":
                typ = "chemical_entity" if biolink else "Chemical:OM:0037"
            else:
                typ = "protein" if biolink else "Protein:MI:0326"
            cur.execute(
                q(
                    "INSERT INTO {s}.entity(entity_id,entity_type_id,canonical_identifier,resolution_status_id) VALUES (%s,%s,%s,2)"
                ),
                [eid, types[typ], "preset_" + name],
            )
        cur.execute(q("SELECT name,interaction_class_id FROM {s}.vocab_interaction_class"))
        classes = dict(cur.fetchall())
        records = (
            ("connectomedb", "ligand_h", "receptor_h", "ligand_receptor", 9606, True, None),
            ("connectomedb", "ligand_m", "receptor_m", "ligand_receptor", 10090, None, False),
            ("humangem", "chemical", "protein_h", "transport", 9606, None, None),
            ("bindingdb", "chemical", "protein_h", "signaling", 9606, True, None),
            ("stitch", "chemical", "protein_m", "signaling", 10090, None, None),
        )
        for i, (source, subject, obj, cls, organism, stimulation, inhibition) in enumerate(records):
            cur.execute(
                q(
                    "INSERT INTO {s}.interaction_fact_resource(interaction_fact_resource_id,subject_entity_id,object_entity_id,interaction_class_id,source_id,is_directed,is_stimulation,is_inhibition,subject_organism,object_organism,reference_pubmed_ids,curation_flags) VALUES (%s,%s,%s,%s,%s,TRUE,%s,%s,%s,%s,%s,%s)"
                ),
                [
                    str(uuid.UUID(int=71000 + i)),
                    entities[subject],
                    entities[obj],
                    classes[cls],
                    sources[source],
                    stimulation,
                    inhibition,
                    organism,
                    organism,
                    [str(100 + i)],
                    ["fixture_claim"],
                ],
            )
    conn.commit()
    return sources


def _full_scope(module):
    # Use the exact public fields demonstrated by main's own round-trip tests.
    # This tests storage, not execution of a hypothetical custom preset query.
    return module.NetworkDefinition(
        name="fixture_scoped",
        kind="ligand_receptor",
        included_sources=("connectomedb2025", "cellphonedb"),
        interaction_class_scope=("ligand_receptor", "transport"),
        evidence_scope={
            "evidence_type": ["experimental", "curated"],
            "predicate": ["binds"],
            "min_confidence": 2,
        },
        default_attributes=("endpoints", "references"),
        mandatory_attributes=("label", "evidence"),
        labels={"preset": "Scoped fixture", "columns": {"label": "Interaction label"}},
        curation={"moa_only": True, "affinity_cutoff": 6.0, "metabolite_class_gate": True},
        attribute_sources={"protein_localization": {"stage": "interim", "source": "uniprot"}},
        collapse_mode="assertion",
        license_scope={"purpose": 15, "sharing": 0, "attrib": 0},
    )


def test_complete_registry_literal_sources_apply_refresh_and_upsert_match_main(
    paired_schemas, oracle
):
    from psycopg2 import sql

    conn, (main, adapted) = paired_schemas
    expected, actual = _modules()
    sources = _substrate(conn, main, expected, biolink=False)
    assert sources == _substrate(conn, adapted, actual, biolink=True)
    baseline = {schema: _base_rows(conn, schema) for schema in (main, adapted)}
    assert baseline[main] == baseline[adapted]
    for prefix, module, schema in (
        (ORACLE_PREFIX, expected, main),
        ("omnipath_subsets", actual, adapted),
    ):
        framework = importlib.import_module(prefix + ".network_views._framework")
        assert framework._warn_unloaded(conn, module.LIANA, registry_schema=schema) == (
            "connectomedb2025",
        )
        assert framework._warn_unloaded(conn, module.METALINKSDB, registry_schema=schema) == (
            "metatlas",
        )
        assert framework._warn_unloaded(conn, module.REACTIONS, registry_schema=schema) == (
            "metatlas",
        )
        stats = module.apply_all(conn, module.NETWORKS, registry_schema=schema)
        assert stats.applied == ("metalinksdb", "liana", "reactions")
        module.register_network(conn, _full_scope(module), registry_schema=schema)
    original, original_stamps = _registry_rows(conn, main)
    current, current_stamps = _registry_rows(conn, adapted)
    assert original == current
    assert set(current) == {"metalinksdb", "liana", "reactions", "fixture_scoped"}
    assert current["liana"]["included_sources"] == ["connectomedb2025"]
    assert current["metalinksdb"]["labels"]["resources"] == {"metatlas": "humangem"}
    assert current["fixture_scoped"]["license_scope"] == {"purpose": 15, "sharing": 0, "attrib": 0}
    assert current["fixture_scoped"]["evidence_scope"]["min_confidence"] == 2
    assert _registry_ddl(conn, main) == _registry_ddl(conn, adapted)
    for schema in (main, adapted):
        assert _base_rows(conn, schema) == baseline[schema]
    # Build/refresh are registration operations for these current presets.
    for module, schema in ((expected, main), (actual, adapted)):
        assert module.refresh_all(conn, module.NETWORKS, registry_schema=schema).applied == (
            "metalinksdb",
            "liana",
            "reactions",
        )
        assert module.apply_all(conn, module.NETWORKS, registry_schema=schema).applied == (
            "metalinksdb",
            "liana",
            "reactions",
        )
        assert _registry_rows(conn, schema)[0] == original
        assert _base_rows(conn, schema) == baseline[schema]
    # A display/module spelling is not the loaded identity. Adding those exact
    # names clears the warnings without rewriting the old source's fact rows.
    for prefix, module, schema in (
        (ORACLE_PREFIX, expected, main),
        ("omnipath_subsets", actual, adapted),
    ):
        with conn.cursor() as cur:
            cur.executemany(
                sql.SQL("INSERT INTO {}.data_source(source_id,name) VALUES (%s,%s)")
                .format(sql.Identifier(schema))
                .as_string(conn),
                [(620001, "connectomedb2025"), (620002, "metatlas")],
            )
        conn.commit()
        framework = importlib.import_module(prefix + ".network_views._framework")
        assert all(
            framework._warn_unloaded(conn, d, registry_schema=schema) == () for d in module.NETWORKS
        )
        updated = replace(
            _full_scope(module),
            included_sources=("stitch",),
            default_attributes=("evidence",),
            mandatory_attributes=("role",),
            evidence_scope=None,
            labels={"preset": "Updated scoped fixture"},
            curation=None,
            attribute_sources=None,
            license_scope=None,
            grain="participant",
            collapse_mode="none",
            composition={
                "operation": "union",
                "components": [{"preset": "liana"}, {"preset": "reactions"}],
                "steps": [{"operation": "collapse"}],
            },
        )
        module.register_network(conn, updated, registry_schema=schema)
    original_after, original_after_stamps = _registry_rows(conn, main)
    current_after, current_after_stamps = _registry_rows(conn, adapted)
    assert original_after == current_after
    assert len(current_after) == 4 and current_after["fixture_scoped"]["included_sources"] == [
        "stitch"
    ]
    assert current_after["fixture_scoped"]["license_scope"] is None
    assert (
        current_after["fixture_scoped"]["grain"] == "participant"
        and current_after["fixture_scoped"]["collapse_mode"] == "none"
    )
    assert [c["preset"] for c in current_after["fixture_scoped"]["composition"]["components"]] == [
        "liana",
        "reactions",
    ]
    assert original_after_stamps["fixture_scoped"] >= original_stamps["fixture_scoped"]
    assert current_after_stamps["fixture_scoped"] >= current_stamps["fixture_scoped"]
    for schema in (main, adapted):
        remaining = _base_rows(conn, schema)
        assert (
            remaining["interaction_fact_resource"] == baseline[schema]["interaction_fact_resource"]
        )
        assert remaining["data_source_license"] == baseline[schema]["data_source_license"]
        assert remaining["objects"] == baseline[schema]["objects"]
    # Defaults do not execute a source, organism or license filter at build.
    assert all(
        original_after[name]["license_scope"] is None
        and original_after[name]["evidence_scope"] is None
        for name in ("liana", "metalinksdb", "reactions")
    )
