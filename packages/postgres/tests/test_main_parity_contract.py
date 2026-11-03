"""Independent main oracle checks, with an opt-in bounded PostgreSQL fixture.

The reference is imported under a private namespace from the frozen snapshot
(or the identical legacy copy). Only import paths change: its algorithms and
DDL are not patched. The SQL fixture declares fixed identities and assertions,
so these checks do not call any resolver, parser, cache or resource download.

Set OMNIPATH_MAIN_PARITY_DSN to an isolated migration database to exercise the
catalogue and duplicate/null/orientation fixture. Never set it to production.
"""

from __future__ import annotations

import ast
from collections import Counter
import importlib
import importlib.abc
import importlib.util
import json
import os
from pathlib import Path
import re
import sys
import types
import uuid

from oracles import LEGACY_ROOT, read_legacy, read_frozen

import pytest

REPO = Path(__file__).resolve().parents[3]
MAIN = LEGACY_ROOT.parent
ORACLE_PREFIX = "_omnipath_frozen_main_oracle"
CV_PREFIX = "_omnipath_frozen_main_cv_terms"
CV_ROOT = Path(__file__).parent / "oracle_cv_terms"
REFERENCE_SHA256 = {
    "metsigdb/sql/extract_classyfire.sql": "b1cf3eebb7e7d221a7c2ca8b68860e665be226847365600e8fe1f243ea0f4fd9",
    "labels/entity_labels.py": "6e0329767a0497272fbb8d42d4022f12c8d0741bbc65eb94e693279d309ed27c",
    "db/schema.py": "59c5593a52b4f42622ab6348e0b9316b6d5355fbe4f5d70e571269ab37a4589b",
    "db/derived_tables.py": "7a31e540727f3956302614ceaab9b1e964ffc2e5bc90fbf303b33093436a802e",
    "db/indexes.py": "6c9bd32d95db5d4a30d221fef75bd2fef40ea7382ee0f331babd08787bcc1457",
    "shared_interaction_schema.py": "7ea50c14b51c429732aea22c5263f9e6a3f21a9569fbb8a6b74be769cc230f5b",
    "classify/interaction_class.yaml": "3a6f691093602af11da0608b829a025d65d2df28f5a22386c63ae45cd3ab81a2",
    "resolver/identifier_types.py": "8c5fbbf679ff761e00553c201ee2391adf85bc5a309dc98d1ce7484f428e96a0",
    "classify/interaction_class.py": "b7898f8f06f41cf01ac25f6efdce51b5543ea67348d1f88d7071b8879a5e7118",
    "cosmos/build.py": "38b10025e5b54e78bf86fead14570e6e766783a91dd578fe57cae49a38080db0",
    "cosmos/translate.py": "a68a3e53470bd40d7d02926afc6b0bef042e157e709a6c498cd1148abfb596d1",
    "cosmos/sql/project_edges.sql": "f67915b1f693b21e930a40cce3d5cf956ce09b26124d752b977cb24e2e77d1c1",
    "cosmos/sql/project_connectors.sql": "1c92b28f6881bdc810ef85aaafeb6fdf327ce1408cfb3745695316b555f73588",
    "cosmos/sql/cosmos_edge_table.sql": "852dc2b0f84b30194a81842654ee0f014174ee43b33f6a0e013776d70c24581b",
}


class _ReferenceImports(ast.NodeTransformer):
    """Namespace imports while leaving reference expressions and SQL untouched."""

    def visit_ImportFrom(self, node):
        if node.module and (
            node.module == "omnipath_build" or node.module.startswith("omnipath_build.")
        ):
            node.module = node.module.replace("omnipath_build", ORACLE_PREFIX, 1)
        elif node.module and (
            node.module == "pypath.internals.cv_terms"
            or node.module.startswith("pypath.internals.cv_terms.")
        ):
            node.module = node.module.replace("pypath.internals.cv_terms", CV_PREFIX, 1)
        return node

    def visit_Import(self, node):
        for name in node.names:
            if name.name == "omnipath_build" or name.name.startswith("omnipath_build."):
                original = name.name
                name.name = name.name.replace("omnipath_build", ORACLE_PREFIX, 1)
                if name.asname is None and original == "omnipath_build":
                    name.asname = "omnipath_build"
        return node


class _ReferenceLoader(importlib.abc.Loader):
    def __init__(self, path):
        self.path = path

    def create_module(self, spec):
        return None

    def exec_module(self, module):
        tree = _ReferenceImports().visit(ast.parse(read_frozen(self.path), str(self.path)))
        module.__file__ = str(self.path)
        exec(compile(ast.fix_missing_locations(tree), str(self.path), "exec"), module.__dict__)


class _ReferenceFinder(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == CV_PREFIX or fullname.startswith(CV_PREFIX + "."):
            relative = fullname[len(CV_PREFIX) :].lstrip(".").replace(".", "/")
            base = CV_ROOT / relative
            file = base / "__init__.py" if base.is_dir() else base.with_suffix(".py")
            return (
                importlib.util.spec_from_loader(
                    fullname, _ReferenceLoader(file), is_package=base.is_dir()
                )
                if file.is_file()
                else None
            )
        if fullname != ORACLE_PREFIX and not fullname.startswith(ORACLE_PREFIX + "."):
            return None
        relative = fullname[len(ORACLE_PREFIX) :].lstrip(".").replace(".", "/")
        base = MAIN / "omnipath_build" / relative
        if base.is_dir():
            # Safe algorithm packages re-export entry points used by main.
            # Root/db/resolver/ingest initialization can configure old sessions
            # or ingestion, and is intentionally bypassed by direct imports.
            initializer = base / "__init__.py"
            if (
                relative in {"classify", "labels", "cosmos", "metsigdb", "network_views"}
                and initializer.is_file()
            ):
                return importlib.util.spec_from_loader(
                    fullname, _ReferenceLoader(initializer), is_package=True
                )
            return importlib.util.spec_from_loader(fullname, loader=None, is_package=True)
        file = base.with_suffix(".py")
        return (
            importlib.util.spec_from_loader(fullname, _ReferenceLoader(file))
            if file.is_file()
            else None
        )


@pytest.fixture(scope="module")
def oracle():
    if not (MAIN / "omnipath_build/db/schema.py").is_file():
        pytest.fail("Frozen main source is required for independent parity checks")
    for relative, digest in REFERENCE_SHA256.items():
        read_legacy(relative, digest)
    finder = _ReferenceFinder()
    sys.meta_path.insert(0, finder)
    try:
        yield types.SimpleNamespace(
            schema=importlib.import_module(ORACLE_PREFIX + ".db.schema"),
            derive=importlib.import_module(ORACLE_PREFIX + ".db.derived_tables"),
            indexes=importlib.import_module(ORACLE_PREFIX + ".db.indexes"),
        )
    finally:
        sys.meta_path.remove(finder)


@pytest.fixture(scope="module")
def port():
    from omnipath_postgres.relational.db import schema, derived_tables, indexes

    return types.SimpleNamespace(schema=schema, derive=derived_tables, indexes=indexes)


@pytest.mark.parametrize(
    "name",
    [
        "SOURCE_PARTITIONED_TABLES",
        "SOURCE_PARTITION_DROP_ORDER",
        "CONTENT_PRIMARY_KEYS",
    ],
)
def test_main_partition_and_primary_key_contracts_are_retained(oracle, port, name):
    assert getattr(port.schema, name) == getattr(oracle.schema, name)


@pytest.mark.parametrize(
    "function,arguments",
    [
        (
            "interaction_content_uuid_sql",
            {"participants": "ARRAY[a, a, b]", "interaction_class": "class_name"},
        ),
        (
            "interaction_record_uuid_sql",
            {
                "subject_entity_id": "subject",
                "object_entity_id": "object",
                "interaction_class": "class_name",
                "source": "source_name",
                "is_directed": "directed",
                "is_stimulation": "stim",
                "is_inhibition": "inhib",
            },
        ),
    ],
)
def test_identity_sql_matches_frozen_main_exactly(oracle, port, function, arguments):
    assert getattr(port.derive, function)(**arguments) == getattr(oracle.derive, function)(
        **arguments
    )


def test_stable_identifier_registry_matches_frozen_main(oracle, port):
    reference = importlib.import_module(ORACLE_PREFIX + ".resolver.identifier_types")
    current = importlib.import_module("omnipath_postgres.relational.identifier_types")
    for name, number in reference.IDENTIFIER_TYPE_IDS.items():
        assert current.IDENTIFIER_TYPE_IDS[name] == number


def _ddl_literals(path):
    result = Counter()
    for node in ast.walk(
        ast.parse(read_frozen(path) if path.is_relative_to(LEGACY_ROOT) else path.read_text())
    ):
        if not isinstance(node, ast.Constant) or not isinstance(node.value, str):
            continue
        text = node.value
        if not re.search(
            r"\bCREATE\s+(?:UNIQUE\s+)?(?:TABLE|INDEX|STATISTICS|EXTENSION)\b", text, re.I
        ):
            continue
        # The only allowed type-literal variation in a structural DDL expression.
        text = text.replace("Gene:MI:0250", "gene").replace("Protein:MI:0326", "protein")
        text = " ".join(text.split())
        # Status 5 truthfully marks published resolution; original statuses,
        # reason/FK fields and all other checks retain their main definitions.
        text = text.replace("status_id IN (1, 2, 3, 5)", "status_id IN (1, 2, 3)")
        result[text] += 1
    return result


@pytest.mark.parametrize("relative", ["db/indexes.py", "db/schema.py"])
def test_main_ddl_definitions_are_retained(oracle, port, relative):
    original = _ddl_literals(MAIN / "omnipath_build" / relative)
    current = _ddl_literals(Path(port.schema.__file__).parent.parent / relative)
    # Additive provenance/quantity tables are permitted; losing a main object is not.
    missing = original - current
    assert not missing, (
        f"Main DDL missing or changed beyond the Biolink type adaptation: {list(missing)}"
    )


@pytest.fixture(scope="module")
def parity_connection():
    dsn = os.environ.get("OMNIPATH_MAIN_PARITY_DSN")
    if not dsn:
        pytest.skip("Set OMNIPATH_MAIN_PARITY_DSN to an isolated migration PostgreSQL")
    import psycopg2

    connection = psycopg2.connect(dsn)
    try:
        yield connection
    finally:
        connection.rollback()
        connection.close()


def _catalogue(conn, schema):
    """Definitions, not object counts; include index operator classes and statistics."""
    with conn.cursor() as cur:
        cur.execute(
            """
            SELECT c.relname, c.relkind, a.attname, format_type(a.atttypid,a.atttypmod),
                   a.attnotnull, pg_get_expr(d.adbin,d.adrelid)
            FROM pg_class c JOIN pg_namespace n ON n.oid=c.relnamespace
            JOIN pg_attribute a ON a.attrelid=c.oid AND a.attnum>0 AND NOT a.attisdropped
            LEFT JOIN pg_attrdef d ON d.adrelid=c.oid AND d.adnum=a.attnum
            WHERE n.nspname=%s AND c.relkind IN ('r','p','v','m')
            ORDER BY c.relname,a.attnum
        """,
            [schema],
        )
        columns = cur.fetchall()
        cur.execute(
            """
            SELECT c.relname, co.contype, pg_get_constraintdef(co.oid, true), co.convalidated
            FROM pg_constraint co JOIN pg_class c ON c.oid=co.conrelid
            JOIN pg_namespace n ON n.oid=c.relnamespace
            WHERE n.nspname=%s ORDER BY c.relname,co.contype,pg_get_constraintdef(co.oid,true)
        """,
            [schema],
        )
        constraints = cur.fetchall()
        cur.execute(
            """
            SELECT t.relname, pg_get_indexdef(i.indexrelid),i.indisvalid,i.indisready
            FROM pg_index i JOIN pg_class t ON t.oid=i.indrelid
            JOIN pg_namespace n ON n.oid=t.relnamespace
            WHERE n.nspname=%s ORDER BY t.relname,pg_get_indexdef(i.indexrelid)
        """,
            [schema],
        )
        indexes = cur.fetchall()
        cur.execute(
            """
            SELECT c.relname, pg_get_partkeydef(c.oid), pg_get_expr(c.relpartbound,c.oid)
            FROM pg_class c JOIN pg_namespace n ON n.oid=c.relnamespace
            WHERE n.nspname=%s AND (c.relkind='p' OR c.relispartition) ORDER BY c.relname
        """,
            [schema],
        )
        partitions = cur.fetchall()
        cur.execute(
            """
            SELECT pg_get_statisticsobjdef(s.oid) FROM pg_statistic_ext s
            JOIN pg_namespace n ON n.oid=s.stxnamespace WHERE n.nspname=%s
            ORDER BY pg_get_statisticsobjdef(s.oid)
        """,
            [schema],
        )
        statistics = cur.fetchall()
        cur.execute(
            """
            SELECT c.relname, pg_get_viewdef(c.oid, true)
            FROM pg_class c JOIN pg_namespace n ON n.oid=c.relnamespace
            WHERE n.nspname=%s AND c.relkind IN ('v','m') ORDER BY c.relname
        """,
            [schema],
        )
        views = cur.fetchall()
        cur.execute(
            """
            SELECT p.proname, pg_get_function_identity_arguments(p.oid),
                   pg_get_functiondef(p.oid)
            FROM pg_proc p JOIN pg_namespace n ON n.oid=p.pronamespace
            WHERE n.nspname=%s AND p.prokind IN ('f','p')
            ORDER BY p.proname,pg_get_function_identity_arguments(p.oid)
        """,
            [schema],
        )
        functions = cur.fetchall()

    def normalize(value):
        if isinstance(value, str):
            value = value.replace(schema + ".", "<schema>.").replace(
                '"' + schema + '".', "<schema>."
            )
            value = value.replace("Gene:MI:0250", "gene").replace("Protein:MI:0326", "protein")
            return re.sub(r"(status_id = ANY \(ARRAY\[)1, 2, 3, 5(\]\))", r"\g<1>1, 2, 3\2", value)
        return value

    return {
        key: Counter(tuple(normalize(v) for v in row) for row in rows)
        for key, rows in {
            "columns": columns,
            "constraints": constraints,
            "indexes": indexes,
            "partitions": partitions,
            "statistics": statistics,
            "views": views,
            "functions": functions,
        }.items()
    }


@pytest.fixture(scope="module")
def paired_schemas(parity_connection, oracle, port):
    from psycopg2 import sql

    schemas = ["parity_main_" + uuid.uuid4().hex, "parity_port_" + uuid.uuid4().hex]
    try:
        with parity_connection.cursor() as cur:
            cur.execute("CREATE EXTENSION IF NOT EXISTS pg_trgm")
        for module, schema in zip((oracle.schema, port.schema), schemas, strict=True):
            module.ensure_schema(parity_connection, schema=schema)
            with parity_connection.cursor() as cur:
                # Definitions built during derive are part of the public contract too.
                derive = oracle.derive if module is oracle.schema else port.derive
                derive._create_derived_tables(cur, schema)
                derive._create_derived_indexes(cur, schema)
                # The companion preserves structured published quantities; it
                # adds no invented old annotations or resolver diagnostics.
                cur.execute(
                    sql.SQL("""
                    CREATE TABLE IF NOT EXISTS {}.annotation_quantity (
                      annotation_key uuid PRIMARY KEY REFERENCES {}.annotation(annotation_key),
                      has_numeric_value double precision, has_unit text,
                      has_unit_prefix text, has_binary_relation text,
                      source_field text, comparator text, published_value text)
                """).format(sql.Identifier(schema), sql.Identifier(schema))
                )
            parity_connection.commit()
        yield parity_connection, tuple(schemas)
    finally:
        parity_connection.rollback()
        with parity_connection.cursor() as cur:
            for schema in schemas:
                cur.execute(
                    sql.SQL("DROP SCHEMA IF EXISTS {} CASCADE").format(sql.Identifier(schema))
                )
        parity_connection.commit()


def test_executed_main_catalogue_matches(paired_schemas):
    conn, (main, adapted) = paired_schemas
    original, current = _catalogue(conn, main), _catalogue(conn, adapted)
    for kind in original:
        assert not original[kind] - current[kind], (
            f"Missing or changed main {kind}: {original[kind] - current[kind]}"
        )


ENTITY_IDS = {name: str(uuid.UUID(int=index)) for index, name in enumerate("abcdef", 1)}


def _binary_substrate(conn, schema, *, biolink):
    """Six entities, three sources, seven assertions; all identities are fixed."""
    from psycopg2 import sql

    def q(text):
        return sql.SQL(text).format(s=sql.Identifier(schema))

    sources = [(9001, "parity_a"), (9002, "parity_b"), (9003, "parity_silent")]
    with conn.cursor() as cur:
        cur.executemany(
            q("INSERT INTO {s}.data_source(source_id,name) VALUES (%s,%s)").as_string(conn), sources
        )
        cur.executemany(
            q("INSERT INTO {s}.dataset(dataset_id,source_id,name) VALUES (%s,%s,'test')").as_string(
                conn
            ),
            [(s, s) for s, _ in sources],
        )
        cur.execute(
            q("INSERT INTO {s}.vocab_entity_type(name) VALUES (%s) RETURNING entity_type_id"),
            ["protein" if biolink else "Protein:MI:0326"],
        )
        type_id = cur.fetchone()[0]
        cur.executemany(
            q(
                "INSERT INTO {s}.entity(entity_id,entity_type_id,taxonomy_id,canonical_identifier,resolution_status_id) VALUES (%s,%s,9606,%s,2)"
            ).as_string(conn),
            [(eid, type_id, "fixture_" + name) for name, eid in ENTITY_IDS.items()],
        )
        cur.execute(
            q(
                "INSERT INTO {s}.vocab_relation_category(name) VALUES ('interaction') RETURNING relation_category_id"
            )
        )
        category = cur.fetchone()[0]
        predicate_names = (
            ["interacts_with", "regulates"]
            if biolink
            else ["interacts_with", "controls", "positively_regulates", "negatively_regulates"]
        )
        cur.executemany(
            q(
                "INSERT INTO {s}.vocab_relation_predicate(name) VALUES (%s) ON CONFLICT(name) DO NOTHING"
            ).as_string(conn),
            [(p,) for p in predicate_names],
        )
        cur.execute(q("SELECT name,relation_predicate_id FROM {s}.vocab_relation_predicate"))
        predicates = dict(cur.fetchall())
        # Main keeps signed verbs as separate graph triples; Biolink qualifiers
        # keep their distinct evidence claims under one regulates graph triple.
        claims = [
            ("a", "b", 9001, "positive"),
            ("a", "b", 9001, "negative"),
            ("a", "b", 9002, "negative"),
            ("a", "b", 9003, None),
            ("b", "a", 9001, None),
            ("c", "d", 9001, "unsigned"),
            ("e", "f", 9001, "ligand_receptor"),
        ]
        graphs = {}
        annotations = {}

        def annotation(term, value):
            key = str(uuid.uuid5(uuid.NAMESPACE_URL, json.dumps([term, value])))
            if key not in annotations:
                cur.execute(
                    q("INSERT INTO {s}.annotation(annotation_key,term,value) VALUES (%s,%s,%s)"),
                    [key, term, value],
                )
                annotations[key] = (term, value)
            return key

        for i, (subject, obj, source, kind) in enumerate(claims, 1):
            pred = (
                "interacts_with"
                if kind in ("unsigned", "ligand_receptor")
                else (
                    "regulates"
                    if biolink
                    else {
                        "positive": "positively_regulates",
                        "negative": "negatively_regulates",
                    }.get(kind, "controls")
                )
            )
            triple = (subject, pred, obj)
            if triple not in graphs:
                rid = str(uuid.uuid5(uuid.NAMESPACE_URL, json.dumps(triple)))
                graphs[triple] = rid
                cur.execute(
                    q(
                        "INSERT INTO {s}.relation(relation_id,subject_entity_id,predicate_id,object_entity_id,relation_category_id) VALUES (%s,%s,%s,%s,%s)"
                    ),
                    [rid, ENTITY_IDS[subject], predicates[pred], ENTITY_IDS[obj], category],
                )
            evidence = str(uuid.UUID(int=100 + i))
            se = oe = None
            if kind == "ligand_receptor" and not biolink:
                se, oe = str(uuid.UUID(int=200 + i)), str(uuid.UUID(int=300 + i))
                for endpoint, eid, role in [
                    (se, ENTITY_IDS[subject], "Ligand:OM:7777"),
                    (oe, ENTITY_IDS[obj], "Receptor:OM:7778"),
                ]:
                    cur.execute(
                        q(
                            "INSERT INTO {s}.entity_evidence(source_id,entity_evidence_id,dataset_id,row_id,entity_role_id,entity_type_id,taxonomy_id) VALUES (%s,%s,%s,%s,1,%s,9606)"
                        ),
                        [source, endpoint, source, i, type_id],
                    )
                    cur.execute(
                        q(
                            "INSERT INTO {s}.entity_evidence_resolution(source_id,entity_evidence_id,status_id,entity_id) VALUES (%s,%s,2,%s)"
                        ),
                        [source, endpoint, eid],
                    )
                    cur.execute(
                        q("INSERT INTO {s}.entity_evidence_annotation VALUES (%s,%s,%s)"),
                        [source, endpoint, annotation(role, None)],
                    )
            cur.execute(
                q(
                    "INSERT INTO {s}.relation_evidence(source_id,relation_evidence_id,dataset_id,row_id,subject_entity_evidence_id,subject_entity_id,predicate_id,object_entity_evidence_id,object_entity_id,relation_category_id) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
                ),
                [
                    source,
                    evidence,
                    source,
                    i,
                    se,
                    None if se else ENTITY_IDS[subject],
                    predicates[pred],
                    oe,
                    None if oe else ENTITY_IDS[obj],
                    category,
                ],
            )
            cur.execute(
                q(
                    "INSERT INTO {s}.relation_evidence_relation(source_id,relation_id,relation_evidence_id) VALUES (%s,%s,%s)"
                ),
                [source, graphs[triple], evidence],
            )
            attrs = []
            if kind in ("positive", "negative"):
                attrs.append(
                    (
                        "object_direction_qualifier",
                        "increased" if kind == "positive" else "decreased",
                        1,
                    )
                    if biolink
                    else ("MI:2235" if kind == "positive" else "MI:2240", None, 1)
                )
            if kind == "ligand_receptor" and biolink:
                attrs += [
                    ("connectomedb:participant_role", "ligand", 2),
                    ("connectomedb:participant_role", "receptor", 3),
                ]
            # Repeating a shared reference must not create duplicate dictionary
            # or evidence-link rows, and each source retains its own reference.
            attrs += (
                [("publications", "PMID:12345", 1)] * 2
                if biolink
                else [("Pubmed:MI:0446", "12345", 1)] * 2
            )
            for term, value, scope in attrs:
                cur.execute(
                    q(
                        "INSERT INTO {s}.relation_evidence_annotation VALUES (%s,%s,%s,%s) ON CONFLICT DO NOTHING"
                    ),
                    [source, evidence, annotation(term, value), scope],
                )
    conn.commit()


def _scientific_rows(conn, schema):
    from psycopg2 import sql

    queries = {
        "facts": "SELECT f.subject_entity_id::text,f.object_entity_id::text,c.name,d.name,f.is_directed,f.is_stimulation,f.is_inhibition,f.reference_pubmed_ids,f.reference_dois,f.interaction_fact_resource_id::text,f.interaction_id::text,f.affinity,f.pchembl,f.score,f.curation_flags,f.subject_organism,f.object_organism,f.attributes->'transport' FROM {s}.interaction_fact_resource f JOIN {s}.vocab_interaction_class c USING(interaction_class_id) JOIN {s}.data_source d USING(source_id)",
        "headers": "SELECT i.interaction_id::text,c.name,i.arity,i.sources FROM {s}.interaction i JOIN {s}.vocab_interaction_class c USING(interaction_class_id)",
        "parties": "SELECT p.interaction_id::text,p.entity_id::text,r.name,p.side,p.ordinal,p.stoichiometry,p.compartment,p.organism,p.role_flag FROM {s}.interaction_party p JOIN {s}.vocab_relation_role r ON r.relation_role_id=p.role_id",
    }
    result = {}
    with conn.cursor() as cur:
        for kind, query in queries.items():
            cur.execute(sql.SQL(query).format(s=sql.Identifier(schema)))
            result[kind] = Counter(
                json.dumps(row, sort_keys=True, default=str) for row in cur.fetchall()
            )
    return result


def _derive_pair(conn, main, adapted, oracle, port):
    # Each fixture can add new predicates. The real downstream pipeline
    # classifies the complete loaded vocabulary before interaction derivation;
    # main's defensive ensure shortcut is not an incremental classifier.
    for prefix, module, schema in [
        (ORACLE_PREFIX, oracle.derive, main),
        ("omnipath_postgres.relational", port.derive, adapted),
    ]:
        classify = importlib.import_module(prefix + ".classify")
        classify.classify_interaction_class(conn, schema=schema)
        conn.commit()
        module.rebuild_interaction_tables(conn, schema=schema)
        conn.commit()


def _assert_scientific_parity(expected, actual, aliases=None):
    """Bounded deltas only; every displayed value belongs to a fixed fixture."""
    names = {value: key for key, value in {**ENTITY_IDS, **(aliases or {})}.items()}

    def display(counter):
        rows = []
        for encoded, count in list(counter.items())[:3]:
            row = json.loads(encoded)
            rows.append(
                {"row": [names.get(v, v) if isinstance(v, str) else v for v in row], "count": count}
            )
        return rows

    deltas = {}
    for kind in expected:
        missing = expected[kind] - actual[kind]
        extra = actual[kind] - expected[kind]
        if missing or extra:
            deltas[kind] = {
                "missing_count": sum(missing.values()),
                "extra_count": sum(extra.values()),
                "missing": display(missing),
                "extra": display(extra),
            }
    assert not deltas, "Fixture-only main parity delta: " + json.dumps(
        deltas, sort_keys=True, default=str
    )


@pytest.fixture(scope="module")
def binary_pair(paired_schemas):
    conn, (main, adapted) = paired_schemas
    _binary_substrate(conn, main, biolink=False)
    _binary_substrate(conn, adapted, biolink=True)
    return paired_schemas


def test_duplicate_null_sign_and_orientation_match_main(binary_pair, oracle, port):
    conn, (main, adapted) = binary_pair
    _derive_pair(conn, main, adapted, oracle, port)
    expected, actual = _scientific_rows(conn, main), _scientific_rows(conn, adapted)
    _assert_scientific_parity(expected, actual)
    # These assertions prevent an empty or under-specified oracle from passing.
    facts = [json.loads(row) for row in actual["facts"]]
    signed = [r for r in facts if r[:2] == [ENTITY_IDS["a"], ENTITY_IDS["b"]]]
    assert {(r[3], r[4], r[5], r[6]) for r in signed} == {
        ("parity_a", True, True, None),
        ("parity_a", True, None, True),
        ("parity_b", True, None, True),
        ("parity_silent", True, None, None),
    }
    unsigned = [r for r in facts if r[:2] == [ENTITY_IDS["c"], ENTITY_IDS["d"]]]
    assert len(unsigned) == 1 and unsigned[0][4:7] == [None, None, None]
    assert any(r[:2] == [ENTITY_IDS["b"], ENTITY_IDS["a"]] for r in facts)
    lr = [r for r in facts if r[:2] == [ENTITY_IDS["e"], ENTITY_IDS["f"]]]
    assert len(lr) == 1 and lr[0][2] == "ligand_receptor" and lr[0][4:7] == [True, None, None]


def _reaction_substrate(conn, schema, *, biolink):
    """Fixed stars: merged chemistry, cargo on two sides, pathway and wide key."""
    from psycopg2 import sql

    def q(text):
        return sql.SQL(text).format(s=sql.Identifier(schema))

    entities = {
        name: str(uuid.UUID(int=1000 + i))
        for i, name in enumerate(
            [
                "x",
                "y",
                "z",
                "enzyme",
                "reaction_a",
                "reaction_b",
                "movement",
                "pathway",
                "wide",
                "orphan",
                *[f"wide_{i}" for i in range(71)],
            ]
        )
    }
    events = {"reaction_a", "reaction_b", "movement", "wide", "orphan"}
    with conn.cursor() as cur:
        for name in sorted(entities):
            if name in events:
                typ = (
                    "molecular_activity"
                    if biolink
                    else ("Transport:OM:0035" if name == "movement" else "Reaction:OM:0015")
                )
            elif name == "pathway":
                typ = "pathway" if biolink else "Pathway:OM:0109"
            elif name == "enzyme":
                typ = "protein" if biolink else "Protein:MI:0326"
            else:
                typ = "chemical_entity" if biolink else "Small Molecule:MI:0328"
            cur.execute(
                q(
                    "INSERT INTO {s}.vocab_entity_type(name) VALUES (%s) ON CONFLICT(name) DO NOTHING"
                ),
                [typ],
            )
            cur.execute(q("SELECT entity_type_id FROM {s}.vocab_entity_type WHERE name=%s"), [typ])
            type_id = cur.fetchone()[0]
            cur.execute(
                q(
                    "INSERT INTO {s}.entity(entity_id,entity_type_id,taxonomy_id,canonical_identifier,resolution_status_id) VALUES (%s,%s,9606,%s,2)"
                ),
                [entities[name], type_id, "reaction_fixture_" + name],
            )
        cur.execute(
            q(
                "SELECT relation_category_id FROM {s}.vocab_relation_category WHERE name='interaction'"
            )
        )
        category = cur.fetchone()[0]
        predicates = (
            ("has_input", "has_output", "enabled_by", "catalyzes")
            if biolink
            else ("has_participant", "controls")
        )
        cur.executemany(
            q(
                "INSERT INTO {s}.vocab_relation_predicate(name) VALUES (%s) ON CONFLICT(name) DO NOTHING"
            ).as_string(conn),
            [(p,) for p in predicates],
        )
        cur.execute(q("SELECT name,relation_predicate_id FROM {s}.vocab_relation_predicate"))
        predicate_ids = dict(cur.fetchall())
        graph = {}
        ann_keys = {}
        next_evidence = 5000

        def add(
            event,
            member,
            role,
            source,
            number=None,
            compartment=None,
            *,
            catalyst_form="enabled_by",
        ):
            nonlocal next_evidence
            catalyst = role == "enzyme"
            if catalyst:
                predicate = catalyst_form if biolink else "controls"
                subject, obj = (
                    (event, member)
                    if biolink and catalyst_form == "enabled_by"
                    else (member, event)
                )
            else:
                predicate = (
                    ("has_input" if role == "reactant" else "has_output")
                    if biolink
                    else "has_participant"
                )
                subject, obj = event, member
            triple = (subject, predicate, obj)
            if triple not in graph:
                graph[triple] = str(
                    uuid.uuid5(uuid.NAMESPACE_URL, json.dumps(["reaction", *triple]))
                )
                cur.execute(
                    q(
                        "INSERT INTO {s}.relation(relation_id,subject_entity_id,predicate_id,object_entity_id,relation_category_id) VALUES (%s,%s,%s,%s,%s)"
                    ),
                    [
                        graph[triple],
                        entities[subject],
                        predicate_ids[predicate],
                        entities[obj],
                        category,
                    ],
                )
            next_evidence += 1
            evidence_id = str(uuid.UUID(int=next_evidence))
            cur.execute(
                q(
                    "INSERT INTO {s}.relation_evidence(source_id,relation_evidence_id,dataset_id,row_id,subject_entity_id,predicate_id,object_entity_id,relation_category_id) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)"
                ),
                [
                    source,
                    evidence_id,
                    source,
                    next_evidence,
                    entities[subject],
                    predicate_ids[predicate],
                    entities[obj],
                    category,
                ],
            )
            cur.execute(
                q(
                    "INSERT INTO {s}.relation_evidence_relation(source_id,relation_id,relation_evidence_id) VALUES (%s,%s,%s)"
                ),
                [source, graph[triple], evidence_id],
            )
            attrs = []
            if not biolink and not catalyst:
                attrs.append(
                    ("Reactant:OM:0310" if role == "reactant" else "Product:OM:0311", None)
                )
            if number is not None:
                attrs.append(("stoichiometry" if biolink else "Stoichiometry:OM:1226", str(number)))
            if compartment is not None:
                attrs.append(
                    (
                        "biopax:cellularLocation" if biolink else "Subcellular Location:OM:0604",
                        compartment,
                    )
                )
            for term, value in attrs:
                if (term, value) not in ann_keys:
                    key = str(uuid.uuid5(uuid.NAMESPACE_URL, json.dumps([term, value])))
                    ann_keys[term, value] = key
                    cur.execute(
                        q(
                            "INSERT INTO {s}.annotation(annotation_key,term,value) VALUES (%s,%s,%s) ON CONFLICT DO NOTHING"
                        ),
                        [key, term, value],
                    )
                cur.execute(
                    q(
                        "INSERT INTO {s}.relation_evidence_annotation VALUES (%s,%s,%s,%s) ON CONFLICT DO NOTHING"
                    ),
                    [source, evidence_id, ann_keys[term, value], 1 if biolink else 3],
                )

        # Different event identifiers, same chemistry, unioned source provenance.
        for event, source, form in [
            ("reaction_a", 9001, "enabled_by"),
            ("reaction_b", 9002, "catalyzes"),
        ]:
            for member, role, number, where in [
                ("x", "reactant", 2, "c"),
                ("y", "reactant", 1, "c"),
                ("z", "product", 1, "m"),
            ]:
                add(event, member, role, source, number, where)
            add(event, "enzyme", "enzyme", source, catalyst_form=form)
        # The same cargo is consumed twice in c and produced once in e.
        # Main has one canonical has_participant triple plus two evidence rows;
        # Biolink correctly has distinct has_input and has_output triples.
        add("movement", "x", "reactant", 9001, 2, "c")
        add("movement", "x", "product", 9001, 1, "e")
        add("movement", "y", "reactant", 9001, 1, "c")
        add("movement", "enzyme", "enzyme", 9001)
        # Predicate alone may not turn a pathway into a reaction hyperedge.
        add("pathway", "y", "reactant", 9001)
        add("pathway", "z", "product", 9001)
        # The signature exceeds a BTree index entry; main uses hash indexes.
        for i in range(71):
            add("wide", f"wide_{i}", "product" if i == 70 else "reactant", 9001)
        add("wide", "enzyme", "enzyme", 9001)
        add("orphan", "x", "reactant", 9001, 1, "c")
        add("orphan", "z", "product", 9001, 1, "m")
    conn.commit()
    return entities


@pytest.fixture(scope="module")
def reaction_pair(binary_pair):
    conn, (main, adapted) = binary_pair
    main_entities = _reaction_substrate(conn, main, biolink=False)
    entities = _reaction_substrate(conn, adapted, biolink=True)
    assert entities == main_entities
    return conn, (main, adapted), entities


def test_reaction_multiset_roles_and_catalyst_orientation_match_main(reaction_pair, oracle, port):
    conn, (main, adapted), entities = reaction_pair
    _derive_pair(conn, main, adapted, oracle, port)
    expected, actual = _scientific_rows(conn, main), _scientific_rows(conn, adapted)
    _assert_scientific_parity(expected, actual, entities)
    headers = [json.loads(row) for row in actual["headers"]]
    parties = [json.loads(row) for row in actual["parties"]]
    assert any(row[2] == 72 for row in headers), (
        "Wide reaction must survive the main hash-index path"
    )
    merged = [row for row in headers if row[2] == 4 and set(row[3]) == {"parity_a", "parity_b"}]
    assert len(merged) == 1, "Equal chemistry must merge despite different resource event IDs"
    cargo = [p for p in parties if p[1] == entities["x"]]
    assert any(p[2] == "reactant" and p[5] == "2" and p[6] == "c" for p in cargo)
    assert any(p[2] == "product" and p[5] == "1" and p[6] == "e" for p in cargo)
    assert all(p[1] != entities["pathway"] or p[2] not in ("reactant", "product") for p in parties)


def _measurement_substrate(conn, schema, *, biolink):
    from psycopg2 import sql

    def q(text):
        return sql.SQL(text).format(s=sql.Identifier(schema))

    with conn.cursor() as cur:

        def attach(source, evidence, term, value, *, quantity=None, unit=None):
            key = str(
                uuid.uuid5(
                    uuid.NAMESPACE_URL, json.dumps([term, value, quantity, unit], sort_keys=True)
                )
            )
            cur.execute(
                q(
                    "INSERT INTO {s}.annotation(annotation_key,term,value,unit) VALUES (%s,%s,%s,%s) ON CONFLICT DO NOTHING"
                ),
                [key, term, value, unit],
            )
            if quantity is not None:
                cur.execute(
                    q(
                        "INSERT INTO {s}.annotation_quantity(annotation_key,has_numeric_value,has_unit,has_unit_prefix,has_binary_relation,source_field,comparator,published_value) VALUES (%s,%s,%s,%s,%s,%s,%s,%s) ON CONFLICT DO NOTHING"
                    ),
                    [
                        key,
                        quantity.get("number"),
                        quantity.get("unit"),
                        quantity.get("prefix"),
                        quantity.get("relation"),
                        quantity.get("field"),
                        quantity.get("comparator"),
                        json.dumps(quantity, sort_keys=True),
                    ],
                )
            cur.execute(
                q(
                    "INSERT INTO {s}.relation_evidence_annotation VALUES (%s,%s,%s,1) ON CONFLICT DO NOTHING"
                ),
                [source, str(uuid.UUID(int=evidence)), key],
            )

        # Source A's positive assertion: true pChEMBL and a different generic
        # quantity must not be mixed; the negative assertion stays silent.
        attach(
            9001,
            101,
            "BAO:0000192" if biolink else "Ki:MI:0643",
            "3",
            quantity={"number": 3, "unit": "nM", "comparator": "="} if biolink else None,
            unit="nM",
        )
        attach(
            9001,
            101,
            "has_quantitative_value" if biolink else "Pchembl Value:OM:0708",
            "7.2",
            quantity={"number": 7.2, "field": "pchembl_value"} if biolink else None,
        )
        if biolink:
            attach(
                9001,
                101,
                "has_quantitative_value",
                "37",
                quantity={"number": 37, "unit": "Cel", "field": "temperature"},
            )
            attach(
                9001, 101, "has_quantitative_value", "9.5", quantity={"number": 9.5, "field": "pH"}
            )
            attach(9001, 101, "BAO:0000479", "0.1", quantity={"number": 0.1, "field": "KOFF"})
            attach(
                9001,
                101,
                "BAO:0000192",
                "1",
                quantity={"number": 1, "unit": "nM", "comparator": "<", "relation": "less_than"},
            )
        else:
            # Original numeric guards reject a bound and non-affinity terms.
            attach(9001, 101, "Ki:MI:0643", "<1", unit="nM")
            attach(9001, 101, "Koff:MI:2221", "0.1")
        attach(
            9001,
            101,
            "publications" if biolink else "Doi:MI:0574",
            "doi:10.1/parity" if biolink else "10.1/parity",
        )
        # Another resource: an equivalent unit spelling must produce the main
        # normalized concentration, while references remain source-local.
        attach(
            9002,
            103,
            "BAO:0000034" if biolink else "Kd:MI:0646",
            "0.004" if biolink else "4",
            quantity={"number": 0.004, "unit": "uM", "comparator": "="} if biolink else None,
            unit="uM" if biolink else "nM",
        )
        attach(
            9002,
            103,
            "publications" if biolink else "Pubmed:MI:0446",
            "PMID:99" if biolink else "99",
        )
        # Known source provenance distinguishes STITCH's combined score from
        # action scores and ChEMBL's assay-confidence metadata.
        cur.executemany(
            q("INSERT INTO {s}.data_source(source_id,name) VALUES (%s,%s)").as_string(conn),
            [(9004, "chembl"), (9005, "stitch")],
        )
        cur.executemany(
            q("INSERT INTO {s}.dataset(dataset_id,source_id,name) VALUES (%s,%s,%s)").as_string(
                conn
            ),
            [(9004, 9004, "mechanisms"), (9005, 9005, "test")],
        )
        cur.execute(
            q(
                "SELECT relation_id,predicate_id,relation_category_id FROM {s}.relation WHERE subject_entity_id=%s AND object_entity_id=%s"
            ),
            [ENTITY_IDS["c"], ENTITY_IDS["d"]],
        )
        relation, predicate, category = cur.fetchone()
        for source, evidence in [(9004, 20000), (9005, 20001)]:
            cur.execute(
                q(
                    "INSERT INTO {s}.relation_evidence(source_id,relation_evidence_id,dataset_id,row_id,subject_entity_id,predicate_id,object_entity_id,relation_category_id) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)"
                ),
                [
                    source,
                    str(uuid.UUID(int=evidence)),
                    source,
                    evidence,
                    ENTITY_IDS["c"],
                    predicate,
                    ENTITY_IDS["d"],
                    category,
                ],
            )
            cur.execute(
                q(
                    "INSERT INTO {s}.relation_evidence_relation(source_id,relation_id,relation_evidence_id) VALUES (%s,%s,%s)"
                ),
                [source, relation, str(uuid.UUID(int=evidence))],
            )
        attach(
            9005,
            20001,
            "has_confidence_score" if biolink else "Confidence Value:OM:1201",
            "8",
            quantity={"number": 8, "field": "combined_score"} if biolink else None,
        )
        if biolink:
            attach(
                9005,
                20001,
                "has_confidence_score",
                "999",
                quantity={"number": 999, "field": "action_score"},
            )
            attach(9001, 101, "chembl_confidence_score", "12")
        else:
            attach(9001, 101, "Chembl Assay Confidence:OM:0762", "12")
        # A mechanisms dataset proves origin even without generic annotations.
        if not biolink:
            attach(9004, 20000, "Chembl Mechanism:OM:0227", "123")
    conn.commit()


def test_measurement_semantics_publications_and_curation_match_main(binary_pair, oracle, port):
    conn, (main, adapted) = binary_pair
    _measurement_substrate(conn, main, biolink=False)
    _measurement_substrate(conn, adapted, biolink=True)
    _derive_pair(conn, main, adapted, oracle, port)
    _assert_scientific_parity(_scientific_rows(conn, main), _scientific_rows(conn, adapted))
    facts = [json.loads(row) for row in _scientific_rows(conn, adapted)["facts"]]
    positive = [
        r
        for r in facts
        if r[:4] == [ENTITY_IDS["a"], ENTITY_IDS["b"], "signaling", "parity_a"] and r[5] is True
    ]
    assert len(positive) == 1 and positive[0][11:14] == [3.0, 7.2, None]
    other = [
        r for r in facts if r[:4] == [ENTITY_IDS["a"], ENTITY_IDS["b"], "signaling", "parity_b"]
    ]
    assert len(other) == 1 and other[0][11] == 4.0
    assert set(other[0][7]) == {"12345", "99"}
    assert all("99" not in (r[7] or []) for r in facts if r[3] != "parity_b")
    stitch = [r for r in facts if r[3] == "stitch"]
    assert len(stitch) == 1 and stitch[0][13] == 8.0
    mechanism = [r for r in facts if r[3] == "chembl"]
    assert len(mechanism) == 1 and mechanism[0][14] == ["mechanism_of_action"]


def _direction_substrate(conn, schema, entities, *, biolink):
    from psycopg2 import sql

    def q(text):
        return sql.SQL(text).format(s=sql.Identifier(schema))

    with conn.cursor() as cur:
        # Canonical labels and shared published chemical aliases, with no
        # external translation. A fallback namespace must stay truthful.
        cur.execute(q("SELECT name,identifier_type_id FROM {s}.vocab_identifier_type"))
        id_types = dict(cur.fetchall())
        for name, namespace, identifier in [
            ("x", "Chebi:MI:0474", "123"),
            ("y", "Pubchem Compound:OM:0002", "678"),
            ("z", "Hmdb:OM:0004", "HMDB0000001"),
            ("enzyme", "Uniprot:MI:1097", "P12345"),
        ]:
            cur.execute(
                q(
                    "UPDATE {s}.entity SET canonical_identifier_type_id=%s,canonical_identifier=%s WHERE entity_id=%s"
                ),
                [id_types[namespace], identifier, entities[name]],
            )
        for name, value in [("x", "123"), ("y", "456"), ("z", "789")]:
            identifier = str(uuid.uuid5(uuid.NAMESPACE_URL, json.dumps(["Chebi:MI:0474", value])))
            cur.execute(
                q(
                    "INSERT INTO {s}.identifier_evidence(identifier_id,identifier_type_id,value) VALUES (%s,%s,%s) ON CONFLICT DO NOTHING"
                ),
                [identifier, id_types["Chebi:MI:0474"], value],
            )
            for source in (9001, 9002, 9001):
                cur.execute(
                    q(
                        "INSERT INTO {s}.entity_identifier(source_id,entity_id,identifier_id) VALUES (%s,%s,%s) ON CONFLICT DO NOTHING"
                    ),
                    [source, entities[name], identifier],
                )
        cur.execute(q("CREATE TABLE IF NOT EXISTS {s}.build_manifest (build_id text NOT NULL)"))
        cur.execute(q("INSERT INTO {s}.build_manifest(build_id) VALUES ('parity_fixture_v1')"))
        for index, (event, direction) in enumerate(
            [("reaction_a", "REVERSIBLE"), ("movement", "LEFT-TO-RIGHT")], 30000
        ):
            term = "biopax:conversionDirection" if biolink else "Conversion Direction:OM:1211"
            key = str(uuid.uuid5(uuid.NAMESPACE_URL, json.dumps([term, direction])))
            cur.execute(
                q(
                    "INSERT INTO {s}.annotation(annotation_key,term,value) VALUES (%s,%s,%s) ON CONFLICT DO NOTHING"
                ),
                [key, term, direction],
            )
            if biolink:
                cur.execute(
                    q(
                        "SELECT source_id,relation_evidence_id FROM {s}.relation_evidence WHERE source_id=9001 AND subject_entity_id=%s ORDER BY row_id LIMIT 1"
                    ),
                    [entities[event]],
                )
                source, evidence = cur.fetchone()
                cur.execute(
                    q("INSERT INTO {s}.relation_evidence_annotation VALUES (%s,%s,%s,1)"),
                    [source, evidence, key],
                )
            else:
                evidence = str(uuid.UUID(int=index))
                cur.execute(
                    q(
                        "INSERT INTO {s}.entity_evidence(source_id,entity_evidence_id,dataset_id,row_id,entity_role_id,entity_type_id,taxonomy_id) SELECT 9001,%s,9001,%s,1,entity_type_id,9606 FROM {s}.entity WHERE entity_id=%s"
                    ),
                    [evidence, index, entities[event]],
                )
                cur.execute(
                    q(
                        "INSERT INTO {s}.entity_evidence_resolution(source_id,entity_evidence_id,status_id,entity_id) VALUES (9001,%s,2,%s)"
                    ),
                    [evidence, entities[event]],
                )
                cur.execute(
                    q("INSERT INTO {s}.entity_evidence_annotation VALUES (9001,%s,%s)"),
                    [evidence, key],
                )
    conn.commit()


def _cosmos_rows(conn, schema):
    from psycopg2 import sql

    with conn.cursor() as cur:
        # Exclude only the insertion-order surrogate. All labels, identities,
        # compartments, sources and published direction remain in the oracle.
        cur.execute(
            sql.SQL("""
            SELECT build_id,source_label,target_label,
                source_entity_id::text,target_entity_id::text,
                source_type,target_type,source_compartment,target_compartment,
                mor,interaction_type,reaction_entity_id::text,interaction_id::text,
                reaction_index,orphan,reverse,direction,
                source_id_type,target_id_type,sources
            FROM {}.cosmos_edge
        """).format(sql.Identifier(schema))
        )
        return Counter(json.dumps(row, sort_keys=True, default=str) for row in cur.fetchall())


def test_cosmos_full_edges_labels_and_statistics_match_main(reaction_pair, oracle, port):
    conn, (main, adapted), entities = reaction_pair
    _direction_substrate(conn, main, entities, biolink=False)
    _direction_substrate(conn, adapted, entities, biolink=True)
    _derive_pair(conn, main, adapted, oracle, port)
    reference = importlib.import_module(ORACLE_PREFIX + ".cosmos")
    current = importlib.import_module("omnipath_subsets.cosmos")
    # The fixture sources are an explicit identical scope in both algorithms.
    # External mapping services are disabled; namespace fallbacks are tested.
    expected_stats = reference.build_cosmos_projection(
        conn, schema=main, reaction_sources=("parity_a", "parity_b"), utils_db_url=None
    )
    actual_stats = current.build_cosmos_projection(
        conn, schema=adapted, reaction_sources=("parity_a", "parity_b"), utils_db_url=None
    )
    expected, actual = _cosmos_rows(conn, main), _cosmos_rows(conn, adapted)
    _assert_scientific_parity({"cosmos": expected}, {"cosmos": actual}, entities)
    assert {k: v for k, v in vars(expected_stats).items() if k != "seconds"} == {
        k: v for k, v in vars(actual_stats).items() if k != "seconds"
    }
    rows = [json.loads(row) for row in actual]
    assert len(rows) > 0 and actual_stats.reactions == 4
    assert actual_stats.reversible_reactions == 1 and actual_stats.reverse_edges > 0
    assert actual_stats.orphan_reactions == 1
    assert actual_stats.connectors > 0
    assert any(r[10] == "transport" and r[11] == entities["movement"] for r in rows)
    assert all(r[16] == "reversible" for r in rows if r[15])
    # Re-run the structural comparison after actual scientific derivation so
    # late key/index restoration and resource-fact statistics are checked.
    original, current_catalogue = _catalogue(conn, main), _catalogue(conn, adapted)
    for kind in original:
        assert not original[kind] - current_catalogue[kind], (
            f"Main {kind} changed after scientific derivation: {original[kind] - current_catalogue[kind]}"
        )


def _ontology_substrate(conn, schema, *, biolink, definition_term):
    from psycopg2 import sql

    def q(text):
        return sql.SQL(text).format(s=sql.Identifier(schema))

    terms = {name: str(uuid.UUID(int=40000 + i)) for i, name in enumerate("abcd")}
    genes = {
        name: str(uuid.UUID(int=40100 + i))
        for i, name in enumerate(("attested", "shortest", "alphabetic"))
    }
    with conn.cursor() as cur:
        cur.executemany(
            q("INSERT INTO {s}.data_source(source_id,name) VALUES (%s,%s)").as_string(conn),
            [(9006, "chemont"), (9007, "ontology_alt")],
        )
        cur.executemany(
            q(
                "INSERT INTO {s}.dataset(dataset_id,source_id,name) VALUES (%s,%s,'ontology')"
            ).as_string(conn),
            [(9006, 9006), (9007, 9007)],
        )
        cur.execute(
            q(
                "INSERT INTO {s}.vocab_identifier_type(identifier_type_id,name) VALUES (9100,'Synonym:OM:0203') ON CONFLICT(name) DO NOTHING"
            )
        )
        cur.execute(q("SELECT name,identifier_type_id FROM {s}.vocab_identifier_type"))
        namespaces = dict(cur.fetchall())
        for typ in (
            "ontology_class" if biolink else "Cv Term:OM:0012",
            "gene" if biolink else "Gene:MI:0250",
        ):
            cur.execute(
                q(
                    "INSERT INTO {s}.vocab_entity_type(name) VALUES (%s) ON CONFLICT(name) DO NOTHING"
                ),
                [typ],
            )
        cur.execute(q("SELECT name,entity_type_id FROM {s}.vocab_entity_type"))
        types = dict(cur.fetchall())

        def identifier(namespace, value):
            key = str(uuid.uuid5(uuid.NAMESPACE_URL, json.dumps([namespace, value])))
            cur.execute(
                q(
                    "INSERT INTO {s}.identifier_evidence(identifier_id,identifier_type_id,value) VALUES (%s,%s,%s) ON CONFLICT DO NOTHING"
                ),
                [key, namespaces[namespace], value],
            )
            return key

        for name, eid in terms.items():
            value = "CHEMONT:" + name.upper()
            cur.execute(
                q(
                    "INSERT INTO {s}.entity(entity_id,entity_type_id,canonical_identifier_type_id,canonical_identifier,resolution_status_id) VALUES (%s,%s,%s,%s,2)"
                ),
                [
                    eid,
                    types["ontology_class" if biolink else "Cv Term:OM:0012"],
                    namespaces["Cv Term Accession:OM:0204"],
                    value,
                ],
            )
            primary = identifier("Cv Term Accession:OM:0204", value)
            for namespace, text in [
                ("Cv Term Accession:OM:0204", value),
                ("Cv Term Accession:OM:0204", value + "_ALT"),
                ("Name:OM:0202", name.upper() + " label"),
                ("Synonym:OM:0203", name.upper() + " synonym"),
            ]:
                key = identifier(namespace, text)
                cur.execute(
                    q(
                        "INSERT INTO {s}.entity_identifier VALUES (9006,%s,%s) ON CONFLICT DO NOTHING"
                    ),
                    [eid, key],
                )
            evidence = str(uuid.UUID(int=41000 + ord(name)))
            cur.execute(
                q(
                    "INSERT INTO {s}.entity_evidence(source_id,entity_evidence_id,dataset_id,row_id,entity_role_id,entity_type_id) VALUES (9006,%s,9006,%s,1,%s)"
                ),
                [evidence, ord(name), types["ontology_class" if biolink else "Cv Term:OM:0012"]],
            )
            cur.execute(
                q(
                    "INSERT INTO {s}.entity_evidence_resolution(source_id,entity_evidence_id,status_id,entity_id) VALUES (9006,%s,%s,%s)"
                ),
                [evidence, 5 if biolink else 2, eid],
            )
            cur.execute(
                q("INSERT INTO {s}.entity_evidence_identifier VALUES (9006,%s,%s)"),
                [evidence, primary],
            )
            ann = str(
                uuid.uuid5(uuid.NAMESPACE_URL, json.dumps([definition_term, name + " definition"]))
            )
            cur.execute(
                q("INSERT INTO {s}.annotation(annotation_key,term,value) VALUES (%s,%s,%s)"),
                [ann, definition_term, name + " definition"],
            )
            cur.execute(
                q("INSERT INTO {s}.entity_evidence_annotation VALUES (9006,%s,%s)"), [evidence, ann]
            )
        predicate = "subclass_of" if biolink else "is_a"
        for pred in (predicate, "part_of", "related_to"):
            cur.execute(
                q(
                    "INSERT INTO {s}.vocab_relation_predicate(name) VALUES (%s) ON CONFLICT(name) DO NOTHING"
                ),
                [pred],
            )
        cur.execute(q("SELECT name,relation_predicate_id FROM {s}.vocab_relation_predicate"))
        predicates = dict(cur.fetchall())
        for source, subject, pred, obj in [
            (9006, "a", predicate, "b"),
            (9006, "b", predicate, "c"),
            (9006, "a", predicate, "c"),
            (9007, "c", predicate, "d"),
            (9006, "a", "part_of", "d"),
            (9006, "b", "related_to", "d"),
        ]:
            cur.execute(
                q(
                    "INSERT INTO {s}.entity_ontology_relation VALUES (%s,%s,%s,%s,'chemont') ON CONFLICT DO NOTHING"
                ),
                [source, terms[subject], predicates[pred], terms[obj]],
            )
        for name, eid in genes.items():
            cur.execute(
                q(
                    "INSERT INTO {s}.entity(entity_id,entity_type_id,taxonomy_id,canonical_identifier_type_id,canonical_identifier,resolution_status_id) VALUES (%s,%s,9606,%s,%s,2)"
                ),
                [
                    eid,
                    types["gene" if biolink else "Gene:MI:0250"],
                    namespaces["Ensembl:MI:0476"],
                    "ENSG_" + name,
                ],
            )
        for name, source, symbol in [
            ("attested", 9001, "LONG_SYMBOL"),
            ("attested", 9002, "LONG_SYMBOL"),
            ("attested", 9001, "A"),
            ("shortest", 9001, "X"),
            ("shortest", 9002, "LONGER"),
            ("alphabetic", 9001, "ABA"),
            ("alphabetic", 9002, "AAB"),
        ]:
            key = identifier("Gene Name Primary:OM:0200", symbol)
            cur.execute(
                q("INSERT INTO {s}.entity_identifier VALUES (%s,%s,%s) ON CONFLICT DO NOTHING"),
                [source, genes[name], key],
            )
    conn.commit()
    return terms, genes


def _ontology_rows(conn, schema, module, sql_root):
    from psycopg2 import sql

    with conn.cursor() as cur:
        module._populate_entity_identifier_lookup(cur, schema)
        module._populate_entity_ontology_terms(cur, schema)
        # Execute the unchanged closure section of main's ClassyFire script.
        # It is bounded here by four terms, not by rebuilding the product.
        path = sql_root / "metsigdb/sql/extract_classyfire.sql"
        if path.is_relative_to(LEGACY_ROOT):
            script = read_frozen(path)
        else:
            from omnipath_subsets import metsigdb

            script = (Path(metsigdb.__file__).parent / "sql/extract_classyfire.sql").read_text()
        closure = script.split("-- The direct assignments,", 1)[0]
        cur.execute(sql.SQL("SET LOCAL search_path = {},pg_catalog").format(sql.Identifier(schema)))
        cur.execute(closure, {"hierarchy_source_id": 9006})
        cur.execute("SELECT node::text,ancestor::text,depth FROM metsigdb_chemont_ancestor")
        ancestors = cur.fetchall()
        cur.execute("DROP TABLE metsigdb_chemont_ancestor,metsigdb_chemont_edge")
        cur.execute(
            sql.SQL(
                "SELECT term_entity_id::text,term_id,ontology_prefix,label,definition,synonyms,synonyms_text,term_aliases,identifiers_text,ontology_id,sources,child_count FROM {}.entity_ontology_term"
            ).format(sql.Identifier(schema))
        )
        terms = cur.fetchall()
        cur.execute(
            sql.SQL(
                "SELECT d.name,r.subject_entity_id::text,CASE WHEN p.name='subclass_of' THEN 'is_a' ELSE p.name END,r.object_entity_id::text,r.ontology_id FROM {}.entity_ontology_relation r JOIN {}.data_source d USING(source_id) JOIN {}.vocab_relation_predicate p ON p.relation_predicate_id=r.predicate_id"
            ).format(*[sql.Identifier(schema)] * 3)
        )
        assertions = cur.fetchall()
    conn.commit()
    return {
        key: Counter(json.dumps(row, sort_keys=True, default=str) for row in rows)
        for key, rows in {
            "ontology_terms": terms,
            "ontology_assertions": assertions,
            "closure": ancestors,
        }.items()
    }


def test_ontology_term_labels_source_scope_and_closure_match_main(binary_pair, oracle, port):
    conn, (main, adapted) = binary_pair
    terms, genes = _ontology_substrate(
        conn, main, biolink=False, definition_term=oracle.derive.ONTOLOGY_DEFINITION_TERM
    )
    actual_terms, actual_genes = _ontology_substrate(
        conn, adapted, biolink=True, definition_term=port.derive.ONTOLOGY_DEFINITION_TERM
    )
    assert terms == actual_terms and genes == actual_genes
    expected = _ontology_rows(conn, main, oracle.derive, MAIN / "omnipath_build")
    current = _ontology_rows(conn, adapted, port.derive, Path(port.derive.__file__).parent.parent)
    _assert_scientific_parity(expected, current, {**terms, **genes})
    closure = [json.loads(row) for row in current["closure"]]
    assert [terms["a"], terms["c"], 1] in closure
    assert all(row[1] != terms["d"] for row in closure), (
        "Other source and non-is_a axioms must not enter this closure"
    )
    rows = [json.loads(row) for row in current["ontology_terms"]]
    a = [r for r in rows if r[0] == terms["a"]]
    assert len(a) == 1 and a[0][3:5] == ["A label", "a definition"]
    assert "A synonym" in a[0][5] and "CHEMONT:A_ALT" in a[0][7]
    c = [r for r in rows if r[0] == terms["c"]]
    assert set(c[0][10]) == {"chemont", "ontology_alt"}
    for prefix, schema in [(ORACLE_PREFIX, main), ("omnipath_postgres.relational", adapted)]:
        labels = importlib.import_module(prefix + ".labels.entity_labels")
        labels.populate_entity_labels(conn, schema=schema)
    from psycopg2 import sql

    def gene_rows(schema):
        with conn.cursor() as cur:
            cur.execute(
                sql.SQL(
                    "SELECT entity_id::text,label,label_rule FROM {}.entity WHERE entity_id=ANY(%s::uuid[])"
                ).format(sql.Identifier(schema)),
                [list(genes.values())],
            )
            return set(cur.fetchall())

    expected_labels = {
        (genes["attested"], "LONG_SYMBOL", "gene_symbol"),
        (genes["shortest"], "X", "gene_symbol"),
        (genes["alphabetic"], "AAB", "gene_symbol"),
    }
    assert gene_rows(main) == gene_rows(adapted) == expected_labels


def _maturation_substrate(conn, schema, *, biolink):
    from psycopg2 import sql

    def q(text):
        return sql.SQL(text).format(s=sql.Identifier(schema))

    entities = {
        name: str(uuid.UUID(int=45000 + i))
        for i, name in enumerate(
            (
                "pure_pre",
                "pure_mat",
                "wrong_pre",
                "wrong_mat",
                "protein_pre",
                "protein_mat",
                "shared_pre",
                "shared_mat",
            )
        )
    }
    with conn.cursor() as cur:
        cur.executemany(
            q("INSERT INTO {s}.data_source(source_id,name) VALUES (%s,%s)").as_string(conn),
            [(90008, "mirbase"), (90009, "maturation_other")],
        )
        cur.executemany(
            q(
                "INSERT INTO {s}.dataset(dataset_id,source_id,name) VALUES (%s,%s,'maturation')"
            ).as_string(conn),
            [(90008, 90008), (90009, 90009)],
        )
        cur.execute(q("SELECT name,identifier_type_id FROM {s}.vocab_identifier_type"))
        namespaces = dict(cur.fetchall())
        for typ in (
            "rna_product" if biolink else "Rna:MI:0320",
            "protein" if biolink else "Protein:MI:0326",
        ):
            cur.execute(
                q(
                    "INSERT INTO {s}.vocab_entity_type(name) VALUES (%s) ON CONFLICT(name) DO NOTHING"
                ),
                [typ],
            )
        cur.execute(q("SELECT name,entity_type_id FROM {s}.vocab_entity_type"))
        types = dict(cur.fetchall())
        for name, eid in entities.items():
            protein = name.startswith("protein")
            ns = (
                "Uniprot:MI:1097"
                if protein
                else (
                    "Mirbase Precursor:OM:0127"
                    if name.endswith("pre")
                    else "Mirbase Mature:OM:0128"
                )
            )
            typ = (
                ("protein" if biolink else "Protein:MI:0326")
                if protein
                else ("rna_product" if biolink else "Rna:MI:0320")
            )
            cur.execute(
                q(
                    "INSERT INTO {s}.entity(entity_id,entity_type_id,taxonomy_id,canonical_identifier_type_id,canonical_identifier,resolution_status_id) VALUES (%s,%s,9606,%s,%s,2)"
                ),
                [eid, types[typ], namespaces[ns], "maturation_" + name],
            )
        for predicate in ("derives_from", "OM:1257"):
            cur.execute(
                q(
                    "INSERT INTO {s}.vocab_relation_predicate(name) VALUES (%s) ON CONFLICT(name) DO NOTHING"
                ),
                [predicate],
            )
        cur.execute(q("SELECT name,relation_predicate_id FROM {s}.vocab_relation_predicate"))
        predicates = dict(cur.fetchall())
        cur.execute(
            q(
                "SELECT relation_category_id FROM {s}.vocab_relation_category WHERE name='interaction'"
            )
        )
        category = cur.fetchone()[0]
        graphs = {}
        for index, (case, source, mapped) in enumerate(
            [
                ("pure", 90008, True),
                ("wrong", 90009, False),
                ("protein", 90008, False),
                ("shared", 90008, True),
                ("shared", 90009, False),
            ],
            50000,
        ):
            predicate = "derives_from" if biolink or not mapped else "OM:1257"
            subject, obj = (
                (case + "_pre", case + "_mat")
                if mapped and not biolink
                else (case + "_mat", case + "_pre")
            )
            triple = (subject, predicate, obj)
            if triple not in graphs:
                graphs[triple] = str(
                    uuid.uuid5(uuid.NAMESPACE_URL, json.dumps(["maturation", *triple]))
                )
                cur.execute(
                    q(
                        "INSERT INTO {s}.relation(relation_id,subject_entity_id,predicate_id,object_entity_id,relation_category_id) VALUES (%s,%s,%s,%s,%s)"
                    ),
                    [
                        graphs[triple],
                        entities[subject],
                        predicates[predicate],
                        entities[obj],
                        category,
                    ],
                )
            evidence = str(uuid.UUID(int=index))
            cur.execute(
                q(
                    "INSERT INTO {s}.relation_evidence(source_id,relation_evidence_id,dataset_id,row_id,subject_entity_id,predicate_id,object_entity_id,relation_category_id) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)"
                ),
                [
                    source,
                    evidence,
                    source,
                    index,
                    entities[subject],
                    predicates[predicate],
                    entities[obj],
                    category,
                ],
            )
            cur.execute(
                q("INSERT INTO {s}.relation_evidence_relation VALUES (%s,%s,%s)"),
                [source, graphs[triple], evidence],
            )
    conn.commit()
    return entities


def test_mirbase_source_specific_maturation_and_generic_claims_match_main(
    binary_pair, oracle, port
):
    conn, (main, adapted) = binary_pair
    entities = _maturation_substrate(conn, main, biolink=False)
    assert entities == _maturation_substrate(conn, adapted, biolink=True)
    _derive_pair(conn, main, adapted, oracle, port)
    expected, actual = _scientific_rows(conn, main), _scientific_rows(conn, adapted)
    _assert_scientific_parity(expected, actual, entities)
    facts = [json.loads(row) for row in actual["facts"]]
    for case in ("pure", "shared"):
        mature = [
            r
            for r in facts
            if r[:4] == [entities[case + "_pre"], entities[case + "_mat"], "maturation", "mirbase"]
        ]
        assert len(mature) == 1 and mature[0][4:7] == [None, None, None]
    for case, source in [
        ("wrong", "maturation_other"),
        ("protein", "mirbase"),
        ("shared", "maturation_other"),
    ]:
        generic = [
            r
            for r in facts
            if r[:4] == [entities[case + "_mat"], entities[case + "_pre"], "other", source]
        ]
        assert len(generic) == 1 and generic[0][4:7] == [None, None, None]
