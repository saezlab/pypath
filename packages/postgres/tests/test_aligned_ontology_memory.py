"""Fresh-main ontology behavior from bounded synthetic published Parquets.

Frozen main's evidence writers leave ontology_terms_raw empty. Its serving
catalogue derives from ontology edges and relational identifiers/annotations.
Fixtures never call a real resolver, parser or resource build.
"""

import ast
from oracles import LEGACY_ROOT, read_legacy
from collections import Counter
from contextlib import closing
from copy import deepcopy
import hashlib
import importlib
import json
from pathlib import Path
from types import SimpleNamespace
import uuid

import duckdb
import pyarrow as pa
import psycopg2
from psycopg2 import sql
import pytest

from omnipath_core.fixtures import nested_rows
from omnipath_build.silver import SilverExtractor
from omnipath_build.writer import ParquetWriter
from omnipath_postgres import loader as aligned_loader, projection as aligned_projection
import test_main_parity_contract as reference_contract
from published_fixture import annotation, entity, fixture_rows, write_fixture


oracle = reference_contract.oracle
port = reference_contract.port
paired_schemas = reference_contract.paired_schemas


@pytest.fixture(scope="module")
def parity_connection(postgres_dsn):
    with closing(psycopg2.connect(postgres_dsn)) as connection:
        try:
            yield connection
        finally:
            connection.rollback()


TERM_COLUMNS = (
    "source_id",
    "term_entity_id",
    "term_id",
    "ontology_prefix",
    "label",
    "definition",
    "ontology_id",
    "synonyms",
    "synonyms_text",
    "sources",
)


def _identifier(value, *, namespace="synonym"):
    return dict(ns=namespace, id=value, is_canonical=False, source="declared-source")


def _description(value, *, term="description"):
    result = annotation(term, value=value)
    result.pop("scope")
    return result


def _term(key, identifier, *, namespace="go", label="Term label", entity_type="ontology_class"):
    item = entity(key, identifier, namespace=namespace, taxon=None, entity_type=entity_type)
    item["label"] = label
    return item


def _statement(subject, obj, *, kind="relation"):
    item = deepcopy(fixture_rows()[1][0])
    item.update(
        relation_key=subject["entity_key"],
        statement_kind=kind,
        subject_entity_key=subject["entity_key"],
        object_entity_key=obj["entity_key"],
        subject_type=subject["entity_type"],
        object_type=obj["entity_type"],
        predicate="associated_with",
        annotations=[annotation("description", value="relation noise")],
        evidence=[
            dict(
                source="noise",
                dataset="ontology",
                row_id="one",
                upstream_id="one",
                annotations=[annotation("description", value="evidence noise")],
            )
        ],
        evidence_count=1,
    )
    return item


def _fixtures(root, *, version="v1"):
    a = _term("term-a", "GO:0000001", label="Term α")
    a["identifiers"] = [
        _identifier("zeta"),
        _identifier("Alpha"),
        _identifier("zeta"),
        _identifier(None),
        _identifier("not a synonym", namespace="name"),
    ]
    # first(description ORDER BY ordinal) intentionally includes a NULL value.
    a["annotations"] = [
        _description("unrelated before", term="comment"),
        _description(None),
        _description("later description"),
        _description("unrelated after", term="publications"),
    ]
    b = _term("term-b", "GO:0000002", label=None)
    b.update(identifiers=None, annotations=None)
    c = _term("term-c", "GO:0000003", label="")
    c.update(identifiers=[], annotations=[])
    d = _term("term-d", "GO:0000004", label="Definition term")
    d["identifiers"] = [
        _identifier(None),
        _identifier(""),
        _identifier('a "quote"\\slash'),
        _identifier("last"),
        _identifier("last"),
    ]
    d["annotations"] = [
        _description("irrelevant", term="other"),
        _description("first description"),
        _description("later description"),
        _description("first description"),
    ]
    noise = _term("protein-noise", "P11111", namespace="uniprot", entity_type="protein")
    noise["identifiers"] = [_identifier("noise synonym")]
    noise["annotations"] = [
        _description("unscoped description"),
        _description("large quantity noise", term="has_quantitative_value"),
    ]
    noise["annotations"][1]["quantity"] = dict(
        has_numeric_value=20.0,
        has_unit="UO:0000001",
        has_unit_prefix="mega",
        has_binary_relation="greater_than",
        source_field="irrelevant field",
        comparator=">",
    )
    wrong_namespace = _term("upper-namespace", "GO:0000005", namespace="GO")
    wrong_namespace["annotations"] = [_description("namespace must match exactly")]
    upper = deepcopy(a)
    upper.update(
        label="Other resource label",
        identifiers=[_identifier("upper source")],
        annotations=[_description("other resource description")],
    )
    pathway = _term("reactome-pathway", "R-HSA-1", namespace="reactome", entity_type="pathway")
    pathway.update(
        identifiers=[_identifier("pathway alias")],
        annotations=[_description("pathway description")],
    )
    activity = _term(
        "reactome-activity", "R-HSA-2", namespace="reactome", entity_type="molecular_activity"
    )
    activity.update(
        identifiers=[_identifier("not pathway alias")], annotations=[_description("activity noise")]
    )
    gene = _term("reactome-gene", "R-HSA-3", namespace="reactome", entity_type="gene")
    gene["annotations"] = [_description("gene noise")]
    hpo = _term("hpo-term", "HP:0000001", namespace="hpo", label="HPO term")
    hpo.update(
        identifiers=[_identifier("HPO alias")], annotations=[_description("HPO description")]
    )
    unscoped = deepcopy(a)
    unscoped.update(
        label="Unscoped shared entity",
        identifiers=[_identifier("unscoped alias")],
        annotations=[_description("unscoped shared description")],
    )
    groups = [
        ("go", [a, b, c, d, noise, wrong_namespace], [_statement(a, b)]),
        ("GO", [upper], []),
        ("reactome", [pathway, activity, gene], [_statement(pathway, activity)]),
        ("hpo", [hpo], []),
        ("unscoped", [unscoped], []),
    ]
    selected = []
    for index, (resource, entities, statements) in enumerate(groups):
        # Keep case-distinct resources separate on macOS case-insensitive disks.
        directory = root / ("resource-" + str(index))
        write_fixture(directory, entities=entities, relations=statements, payloads=[])
        selected.append(SimpleNamespace(source=resource, version=version, directory=directory))
    assert sum(len(entities) + len(statements) for _resource, entities, statements in groups) == 14
    return tuple(selected)


def _entity_id(key):
    payload = json.dumps(["published-entity", key], ensure_ascii=False, separators=(",", ":"))
    return uuid.UUID(hex=hashlib.md5(payload.encode()).hexdigest())


@pytest.mark.parametrize("version", ["v1", "other-version"])
def test_base_ontology_terms_matches_fresh_main_and_keeps_relational_inputs(tmp_path, version):
    selected = _fixtures(tmp_path, version=version)
    observed = []
    with duckdb.connect() as connection:
        plan = aligned_projection.prepare_aligned_release(
            connection,
            selected,
            progress=lambda fields: aligned_loader._emit(
                lambda event, **fields: observed.append((event, fields)),
                "projection_progress",
                **fields,
            ),
        )
        copy = next(query for query in plan.queries if query.table == "ontology_terms")
        assert copy.columns == TERM_COLUMNS
        actual = Counter(connection.execute(copy.query).fetchall())
        actual_types = tuple(str(column[1]) for column in connection.description)
        assert actual == Counter()
        assert plan.counts["ontology_terms"] == 0
        assert actual_types == ("BIGINT", "UUID", *("VARCHAR",) * 8)
        # Empty base search rows do not discard the source-owned inputs used
        # by the actual main entity_ontology_term derivation.
        assert connection.execute(
            """SELECT value FROM ap_annotation_occurrence
            WHERE resource='go' AND version=? AND owner_kind='entity' AND owner_key='term-a'
            ORDER BY ordinal""",
            [version],
        ).fetchall() == [
            ("unrelated before",),
            (None,),
            ("later description",),
            ("unrelated after",),
        ]
        assert connection.execute(
            """SELECT item.id FROM ap_identifier_occurrence
            WHERE resource='go' AND version=? AND entity_key='term-a'
            AND item.ns='synonym' ORDER BY ordinal""",
            [version],
        ).fetchall() == [
            ("zeta",),
            ("Alpha",),
            ("zeta",),
            (None,),
        ]
        copied = [
            fields
            for event, fields in observed
            if event == "projection_progress"
            and fields["phase"] == "copy_prepared"
            and fields["table"] == "ontology_terms"
        ]
        assert copied == [dict(phase="copy_prepared", table="ontology_terms", rows=0)]


def _walk_plan(nodes):
    for node in nodes:
        yield node
        yield from _walk_plan(node.get("children", []))


def test_empty_base_ontology_query_scans_no_identifiers_or_annotations(tmp_path):
    selected = _fixtures(tmp_path)
    with duckdb.connect() as connection:
        aligned_projection.prepare_aligned_release(connection, selected)
        query = next(
            item.query for item in aligned_projection._queries() if item.table == "ontology_terms"
        )
        physical = json.loads(connection.execute("EXPLAIN (FORMAT JSON) " + query).fetchone()[1])
        nodes = list(_walk_plan(physical))
        assert any(node["name"] == "EMPTY_RESULT" for node in nodes)
        assert all("SCAN" not in node["name"] and "GROUP_BY" not in node["name"] for node in nodes)


def test_frozen_main_evidence_writers_concretely_leave_base_terms_empty():
    path = LEGACY_ROOT / "duckdb_load.py"
    source = read_legacy(
        "duckdb_load.py", "2a4f8e4a4ad6c8452b60213517940860f9d50c03d2e127014c4a3e4afd06959b"
    )
    # Execute exact isolated main writer definitions, preserving their function
    # bodies and schema declarations. No main runtime/import/resolver is used.
    wanted = {"_DuckDBEvidenceWriters", "_DuckDBRowWriter", "_create_duckdb_evidence_tables"}
    tree = ast.parse(source)
    definitions = [
        node
        for node in tree.body
        if (isinstance(node, (ast.ClassDef, ast.FunctionDef)) and node.name in wanted)
        or (
            isinstance(node, ast.Assign)
            and any(
                isinstance(target, ast.Name) and target.id.endswith("_SCHEMA")
                for target in node.targets
            )
        )
    ]
    namespace = {"duckdb": duckdb, "pa": pa}
    exec(compile(ast.Module(body=definitions, type_ignores=[]), str(path), "exec"), namespace)
    with duckdb.connect() as connection:
        namespace["_create_duckdb_evidence_tables"](connection)
        writers = namespace["_DuckDBEvidenceWriters"](connection, chunk_size=1)
        assert set(vars(writers)) == {
            "entity",
            "identifier",
            "entity_annotation",
            "relation_annotation",
            "annotation",
            "relation",
            "annotation_relation",
            "ontology_relation",
        }
        for writer in vars(writers).values():
            row = dict.fromkeys(writer.schema.names)
            if "source" in row:
                row["source"] = "synthetic-source"
            writer.write(row)
        writers.close()
        assert all(
            connection.execute("SELECT count(*) FROM " + writer.table).fetchone()[0] == 1
            for writer in vars(writers).values()
        )
        assert connection.execute("SELECT count(*) FROM ontology_terms_raw").fetchone()[0] == 0


def test_raw_taxonomy_capability_counts_match_frozen_distinct_contract(tmp_path):
    variants = [
        (None, "unknown"),
        ("", "unknown"),
        ("9606", "known"),
        ("", "known"),
        ("09606", "leading-zero"),
        ("9606", "leading-zero"),
        ("9606", "conflict"),
        ("10090", "conflict"),
    ]
    selected = []
    for index, (taxon, group) in enumerate(variants):
        item = entity("taxon:" + group, group, taxon=taxon)
        directory = tmp_path / str(index)
        write_fixture(directory, entities=[item], relations=[], payloads=[])
        selected.append(
            SimpleNamespace(
                source="taxon-resource-" + str(index), version="v1", directory=directory
            )
        )
    with duckdb.connect() as connection:
        plan = aligned_projection.prepare_aligned_release(connection, selected)
        conflicting = connection.execute("""SELECT count(*) FROM (
            SELECT entity_key FROM ap_entity_raw GROUP BY entity_key
            HAVING count(DISTINCT nullif(taxon,''))>1) q""").fetchone()[0]
        null_and_known = connection.execute("""SELECT count(*) FROM (
            SELECT entity_key FROM ap_entity_raw GROUP BY entity_key
            HAVING count(DISTINCT nullif(taxon,''))=1 AND bool_or(nullif(taxon,'') IS NULL)) q""").fetchone()[
            0
        ]
        assert conflicting == plan.compatibility["shared_identity_taxonomy_conflicts"] == 2
        assert null_and_known == plan.compatibility["shared_identity_taxonomy_null_and_known"] == 1


class _DeclaredTargets:
    """Supply fixture identities to the real writer without reference resolution."""

    def resolve_entity_targets(self, entities, *, progress):
        assert progress is False
        return {
            key: [
                SimpleNamespace(
                    entity_type=item.entity_type,
                    canonical_namespace=item.namespace,
                    canonical_identifier=item.identifier,
                    taxon=item.taxon,
                    label=item.label,
                    node_id=None,
                    reference_library=None,
                    matched=False,
                    aliases={},
                )
            ]
            for key, item in entities.items()
        }


def _writer_fixture(directory, source, namespace):
    extractor = SilverExtractor(source, "ontology")
    records = [
        dict(
            type="ontology_class",
            identifiers=[dict(type=namespace, value=namespace.upper() + ":" + name)],
            label=name + " label",
            ontology_relations=[],
        )
        for name in ("a", "b", "c", "d")
    ]
    by_name = dict(zip(("a", "b", "c", "d"), records, strict=True))
    facts = [
        ("a", "subclass_of", "b"),
        ("b", "subclass_of", "c"),
        ("a", "part_of", "d"),
        ("b", "related_to", "d"),
    ]
    for subject, predicate, obj in facts:
        by_name[subject]["ontology_relations"].append(
            dict(
                predicate=predicate,
                object=dict(
                    type="ontology_class",
                    identifier_type=namespace,
                    identifier=namespace.upper() + ":" + obj,
                ),
            )
        )
    for number, record in enumerate(records):
        extractor.process_record(record, record, str(number), number)
    writer = ParquetWriter(directory)
    writer.append_observations(extractor, _DeclaredTargets())
    writer.close()
    entities = nested_rows(directory / "entity.parquet")
    statements = nested_rows(directory / "relation.parquet")
    assert len(entities) + len(statements) == 8
    assert {row["statement_kind"] for row in statements} == {"ontology"}
    return SimpleNamespace(source=source, version="writer-v1", directory=directory), entities, facts


@pytest.mark.parametrize("legacy_kind", [False, True])
def test_real_publisher_ontology_kind_projects_edges_and_reports_unknown_scope(
    tmp_path, legacy_kind
):
    selected = []
    expected = []
    for source, namespace, scope in (
        ("chemont", "chemont", "chemont"),
        (
            "unknown-ontology",
            "unknown_ontology",
            "omnipath:published-scope:unknown-ontology:unknown_ontology",
        ),
    ):
        resource, entities, facts = _writer_fixture(tmp_path / source, source, namespace)
        selected.append(resource)
        keys = {item["identifier"].rsplit(":", 1)[1]: item["entity_key"] for item in entities}
        expected.extend(
            (source, _entity_id(keys[subject]), predicate, _entity_id(keys[obj]), scope)
            for subject, predicate, obj in facts
        )
        if legacy_kind:
            # The deprecated adapter spelling remains supported, while the
            # first parameter case checks unmodified real publisher outputs.
            rows = nested_rows(resource.directory / "relation.parquet")
            for row in rows:
                row["statement_kind"] = "ontology_axiom"
            write_fixture(resource.directory, entities=entities, relations=rows, payloads=[])
    with duckdb.connect() as connection:
        plan = aligned_projection.prepare_aligned_release(connection, selected)
        query = next(query for query in plan.queries if query.table == "entity_ontology_relation")
        predicates = dict(
            connection.execute("SELECT id,name FROM ap_vocab_relation_predicate").fetchall()
        )
        sources = dict(plan.sources)
        actual = Counter(
            (sources[source], subject, predicates[predicate], obj, scope)
            for source, subject, predicate, obj, scope in connection.execute(query.query).fetchall()
        )
        assert actual == Counter(expected)
        assert plan.counts["relation"] == 0
        assert plan.counts["entity_ontology_relation"] == 8
        assert plan.compatibility["unknown_ontology_scope_statements"] == 4
        assert plan.counts["ontology_terms"] == 0


def _closure_rows(connection, schema, script):
    closure = script.split("-- The direct assignments,", 1)[0]
    with connection.cursor() as cursor:
        cursor.execute(
            sql.SQL("SET LOCAL search_path={},pg_catalog").format(sql.Identifier(schema))
        )
        cursor.execute(closure, {"hierarchy_source_id": 9006})
        cursor.execute("SELECT node::text,ancestor::text,depth FROM metsigdb_chemont_ancestor")
        result = Counter(cursor.fetchall())
        cursor.execute("DROP TABLE metsigdb_chemont_ancestor,metsigdb_chemont_edge")
    return result


def test_real_publisher_copied_ontology_edges_feed_exact_frozen_main_closure(
    tmp_path, postgres_dsn
):
    resource, entities, facts = _writer_fixture(tmp_path / "chemont", "chemont", "chemont")
    unknown, _unknown_entities, _unknown_facts = _writer_fixture(
        tmp_path / "unknown", "unknown-ontology", "unknown_ontology"
    )
    keys = {item["identifier"].rsplit(":", 1)[1]: item["entity_key"] for item in entities}
    original = read_legacy(
        "metsigdb/sql/extract_classyfire.sql",
        "b1cf3eebb7e7d221a7c2ca8b68860e665be226847365600e8fe1f243ea0f4fd9",
    )
    from omnipath_subsets import metsigdb

    adapted_script = (Path(metsigdb.__file__).parent / "sql/extract_classyfire.sql").read_text()
    schemas = ("ontology_oracle_" + uuid.uuid4().hex, "ontology_copy_" + uuid.uuid4().hex)
    with closing(psycopg2.connect(postgres_dsn)) as postgres, duckdb.connect() as duck:
        # Minimal isolated main-contract tables suffice for the exact frozen
        # hierarchy section. This uses conftest's temporary local cluster only.
        try:
            with postgres.cursor() as cursor:
                for schema in schemas:
                    cursor.execute(sql.SQL("CREATE SCHEMA {}").format(sql.Identifier(schema)))
                    cursor.execute(
                        sql.SQL("""CREATE TABLE {}.vocab_relation_predicate
                        (relation_predicate_id BIGINT PRIMARY KEY,name TEXT NOT NULL)""").format(
                            sql.Identifier(schema)
                        )
                    )
                    cursor.execute(
                        sql.SQL("""CREATE TABLE {}.entity_ontology_relation
                        (source_id BIGINT,subject_entity_id UUID,predicate_id BIGINT,
                         object_entity_id UUID,ontology_id TEXT)""").format(sql.Identifier(schema))
                    )
                cursor.executemany(
                    sql.SQL("INSERT INTO {}.vocab_relation_predicate VALUES(%s,%s)")
                    .format(sql.Identifier(schemas[0]))
                    .as_string(postgres),
                    [(1, "is_a"), (2, "part_of"), (3, "related_to")],
                )
                old_predicates = {"subclass_of": 1, "part_of": 2, "related_to": 3}
                cursor.executemany(
                    sql.SQL("INSERT INTO {}.entity_ontology_relation VALUES(%s,%s,%s,%s,%s)")
                    .format(sql.Identifier(schemas[0]))
                    .as_string(postgres),
                    [
                        (
                            9006,
                            str(_entity_id(keys[subject])),
                            old_predicates[predicate],
                            str(_entity_id(keys[obj])),
                            "chemont",
                        )
                        for subject, predicate, obj in facts
                    ],
                )
            plan = aligned_projection.prepare_aligned_release(
                duck,
                (resource, unknown),
                dimension_rows={"data_source": ((9006, "chemont"), (9008, "unknown-ontology"))},
            )
            with postgres.cursor() as cursor:
                cursor.executemany(
                    sql.SQL("INSERT INTO {}.vocab_relation_predicate VALUES(%s,%s)")
                    .format(sql.Identifier(schemas[1]))
                    .as_string(postgres),
                    plan.dimensions["vocab_relation_predicate"],
                )
            query = next(item for item in plan.queries if item.table == "entity_ontology_relation")
            spool = tmp_path / "copy"
            spool.mkdir()
            assert aligned_loader._copy(postgres, schemas[1], query, duck, spool).rows == 8
            expected = _closure_rows(postgres, schemas[0], original)
            actual = _closure_rows(postgres, schemas[1], adapted_script)
            assert (
                actual
                == expected
                == Counter(
                    [
                        (str(_entity_id(keys["a"])), str(_entity_id(keys["b"])), 1),
                        (str(_entity_id(keys["b"])), str(_entity_id(keys["c"])), 1),
                        (str(_entity_id(keys["a"])), str(_entity_id(keys["c"])), 2),
                    ]
                )
            )
        finally:
            postgres.rollback()


def _serving_fixture(root):
    """Nineteen entities/statements: supported namespace forms and negatives."""
    groups = [
        (
            "chebi",
            "chebi",
            "chemical_entity",
            "inchikey",
            ["AAAAAAAAAAAAAA-BBBBBBBBBB-C", "CCCCCCCCCCCCCC-DDDDDDDDDD-E"],
        ),
        ("mondo", "mondo", "ontology_class", "mondo", ["0000001", "0000002"]),
        ("kegg", "kegg_pathways", "pathway", "kegg_pathway", ["map00010", "map00020"]),
        (
            "uniprot",
            "uniprot_keywords",
            "ontology_class",
            "uniprot_keyword",
            ["KW-0001", "KW-0002"],
        ),
        (
            "go",
            "gene_ontology",
            "ontology_class",
            "go",
            ["http://purl.obolibrary.org/obo/GO_0000001", "GO:0000002"],
        ),
    ]
    selected, declared, axioms = [], [], []
    positive_keys = set()
    for source, ontology, entity_type, namespace, values in groups:
        entities = []
        for index, value in enumerate(values):
            key = source + ":serving:" + str(index)
            item = _term(
                key, value, namespace=namespace, label="label " + key, entity_type=entity_type
            )
            item["identifiers"] = [
                _identifier(value, namespace=namespace),
                _identifier("label " + key, namespace="name"),
                _identifier("synonym " + key),
            ]
            item["annotations"] = [_description("description " + key)]
            if source == "go" and index == 0:
                # Main chooses MIN of supported descriptions, not the first
                # published annotation or an unrelated comment value.
                item["annotations"] = [
                    _description("zulu first description " + key),
                    _description("description " + key),
                    _description("aaa unrelated comment", term="comment"),
                ]
            if source == "chebi":
                chebi_ids = ["CHEBI:20", "CHEBI:30"] if index == 0 else ["CHEBI:40"]
                item["identifiers"].extend(
                    _identifier(value, namespace="chebi") for value in chebi_ids
                )
                if index == 0:
                    # This is an original published CV identifier, not an
                    # adapter alias. It must not promote chemical canonical
                    # InChIKey to priority0 above the genuine Chebi candidates.
                    item["identifiers"].append(_identifier(item["identifier"], namespace="cv_term"))
            if source == "mondo" and index == 0:
                # A canonical CV term must outrank secondary ChEBI identifiers.
                item["identifiers"].append(_identifier("CHEBI:999", namespace="chebi"))
            entities.append(item)
            declared.append((source, ontology, item))
            if entity_type in {"ontology_class", "pathway"}:
                positive_keys.add(key)
        statements = []
        main = _statement(entities[0], entities[1], kind="ontology")
        main.update(relation_key=source + ":ontology:0", predicate="subclass_of")
        statements.append(main)
        axioms.append(
            (source, entities[0]["entity_key"], "subclass_of", entities[1]["entity_key"], ontology)
        )
        if source == "go":
            for key, namespace, identifier, entity_type in (
                ("negative-protein", "uniprot", "P12345", "protein"),
                ("negative-name", "name", "random ontology name", "ontology_class"),
            ):
                item = _term(
                    key,
                    identifier,
                    namespace=namespace,
                    label="label " + key,
                    entity_type=entity_type,
                )
                item["identifiers"] = [
                    _identifier(identifier, namespace=namespace),
                    _identifier("label " + key, namespace="name"),
                ]
                if key == "negative-protein":
                    item["identifiers"].extend(
                        [
                            _identifier(identifier, namespace="cv_term"),
                            _identifier("CHEBI:9991", namespace="chebi"),
                        ]
                    )
                item["annotations"] = [_description("description " + key)]
                entities.append(item)
                declared.append((source, ontology, item))
                stmt = _statement(item, entities[0], kind="ontology")
                stmt.update(relation_key=key + ":ontology", predicate="related_to")
                statements.append(stmt)
                axioms.append((source, key, "related_to", entities[0]["entity_key"], ontology))
        directory = root / source
        write_fixture(directory, entities=entities, relations=statements, payloads=[])
        selected.append(SimpleNamespace(source=source, version="serving-v1", directory=directory))
    assert len(declared) + len(axioms) == 19
    return tuple(selected), declared, axioms, positive_keys


def _declared_main_serving_facts(
    connection, schema, oracle, declared, axioms, positive_keys, source_ids
):
    """Declare supported legacy main facts; no resolution or projection oracle."""
    cv = importlib.import_module(reference_contract.CV_PREFIX)
    ns_label = {
        namespace: cv.cv_term_label_accession(getattr(cv.IdentifierNamespaceCv, member))
        for namespace, member in (
            ("inchikey", "STANDARD_INCHI_KEY"),
            ("chebi", "CHEBI"),
            ("name", "NAME"),
            ("synonym", "SYNONYM"),
            ("uniprot", "UNIPROT"),
            ("cv_term", "CV_TERM_ACCESSION"),
        )
    }
    with connection.cursor() as cursor:

        def q(text):
            return sql.SQL(text).format(s=sql.Identifier(schema))

        cursor.executemany(
            q("INSERT INTO {s}.data_source(source_id,name) VALUES(%s,%s)").as_string(connection),
            [(identifier, source) for source, identifier in source_ids.items()],
        )
        cursor.executemany(
            q(
                "INSERT INTO {s}.dataset(dataset_id,source_id,name) VALUES(%s,%s,'ontology')"
            ).as_string(connection),
            [(identifier, identifier) for identifier in source_ids.values()],
        )
        cursor.execute(
            q(
                "INSERT INTO {s}.vocab_entity_type(name) VALUES('Cv Term:OM:0012') ON CONFLICT DO NOTHING"
            )
        )
        cursor.execute(
            q("SELECT entity_type_id FROM {s}.vocab_entity_type WHERE name='Cv Term:OM:0012'")
        )
        entity_type_id = cursor.fetchone()[0]
        namespaces = set(ns_label.values()) | {
            item["namespace"] for _source, _ontology, item in declared
        }
        cursor.execute(q("SELECT name,identifier_type_id FROM {s}.vocab_identifier_type"))
        namespace_ids = dict(cursor.fetchall())
        next_namespace = max(namespace_ids.values(), default=0)
        for namespace in sorted(namespaces - namespace_ids.keys()):
            next_namespace += 1
            cursor.execute(
                q("INSERT INTO {s}.vocab_identifier_type(identifier_type_id,name) VALUES(%s,%s)"),
                (next_namespace, namespace),
            )
        cursor.execute(q("SELECT name,identifier_type_id FROM {s}.vocab_identifier_type"))
        namespace_ids = dict(cursor.fetchall())

        def identifier(namespace, value):
            identity = str(uuid.uuid5(uuid.NAMESPACE_URL, json.dumps([namespace, value])))
            cursor.execute(
                q(
                    "INSERT INTO {s}.identifier_evidence(identifier_id,identifier_type_id,value) VALUES(%s,%s,%s) ON CONFLICT DO NOTHING"
                ),
                (identity, namespace_ids[namespace], value),
            )
            return identity

        for ordinal, (source, _ontology, item) in enumerate(declared):
            key = item["entity_key"]
            entity_id = str(_entity_id(key))
            evidence_id = str(uuid.uuid5(uuid.NAMESPACE_URL, "main-serving-evidence:" + key))
            canonical_type = (
                ns_label["cv_term"]
                if key in positive_keys
                else ns_label.get(item["namespace"], item["namespace"])
            )
            cursor.execute(
                q(
                    "INSERT INTO {s}.entity(entity_id,entity_type_id,canonical_identifier_type_id,canonical_identifier,resolution_status_id) VALUES(%s,%s,%s,%s,2)"
                ),
                (entity_id, entity_type_id, namespace_ids[canonical_type], item["identifier"]),
            )
            cursor.execute(
                q(
                    "INSERT INTO {s}.entity_evidence(source_id,entity_evidence_id,dataset_id,row_id,entity_role_id,entity_type_id) VALUES(%s,%s,%s,%s,1,%s)"
                ),
                (source_ids[source], evidence_id, source_ids[source], ordinal, entity_type_id),
            )
            cursor.execute(
                q(
                    "INSERT INTO {s}.entity_evidence_resolution(source_id,entity_evidence_id,status_id,entity_id) VALUES(%s,%s,2,%s)"
                ),
                (source_ids[source], evidence_id, entity_id),
            )
            for published in item["identifiers"]:
                namespace = ns_label.get(published["ns"], published["ns"])
                ident_id = identifier(namespace, published["id"])
                cursor.execute(
                    q("INSERT INTO {s}.entity_identifier VALUES(%s,%s,%s) ON CONFLICT DO NOTHING"),
                    (source_ids[source], entity_id, ident_id),
                )
                # Main's legacy ontology identifiers use CV_TERM_ACCESSION;
                # the new path keeps raw occurrence types and joins its exact
                # value compatibility lookup rather than fabricating evidence.
                if (
                    key in positive_keys
                    and published["ns"] == item["namespace"]
                    and published["id"] == item["identifier"]
                ):
                    ident_id = identifier(ns_label["cv_term"], published["id"])
                    cursor.execute(
                        q(
                            "INSERT INTO {s}.entity_identifier VALUES(%s,%s,%s) ON CONFLICT DO NOTHING"
                        ),
                        (source_ids[source], entity_id, ident_id),
                    )
                cursor.execute(
                    q(
                        "INSERT INTO {s}.entity_evidence_identifier VALUES(%s,%s,%s) ON CONFLICT DO NOTHING"
                    ),
                    (source_ids[source], evidence_id, ident_id),
                )
                if published["ns"] == "chebi" and published["id"].startswith("CHEBI:"):
                    bare = identifier(namespace, published["id"].split(":", 1)[1])
                    cursor.execute(
                        q(
                            "INSERT INTO {s}.entity_identifier VALUES(%s,%s,%s) ON CONFLICT DO NOTHING"
                        ),
                        (source_ids[source], entity_id, bare),
                    )
            for ann_ordinal, published in enumerate(item["annotations"]):
                ann_id = str(
                    uuid.uuid5(
                        uuid.NAMESPACE_URL,
                        "main-serving-description:" + key + ":" + str(ann_ordinal),
                    )
                )
                term = (
                    oracle.derive.ONTOLOGY_DEFINITION_TERM
                    if published["term"] == "description"
                    else published["term"]
                )
                cursor.execute(
                    q("INSERT INTO {s}.annotation(annotation_key,term,value) VALUES(%s,%s,%s)"),
                    (ann_id, term, published["value"]),
                )
                cursor.execute(
                    q("INSERT INTO {s}.entity_evidence_annotation VALUES(%s,%s,%s)"),
                    (source_ids[source], evidence_id, ann_id),
                )
        for predicate in ("is_a", "related_to"):
            cursor.execute(
                q(
                    "INSERT INTO {s}.vocab_relation_predicate(name) VALUES(%s) ON CONFLICT DO NOTHING"
                ),
                (predicate,),
            )
        cursor.execute(q("SELECT name,relation_predicate_id FROM {s}.vocab_relation_predicate"))
        predicates = dict(cursor.fetchall())
        cursor.executemany(
            q("INSERT INTO {s}.entity_ontology_relation VALUES(%s,%s,%s,%s,%s)").as_string(
                connection
            ),
            [
                (
                    source_ids[source],
                    str(_entity_id(subject)),
                    predicates["is_a" if predicate == "subclass_of" else predicate],
                    str(_entity_id(obj)),
                    ontology,
                )
                for source, subject, predicate, obj, ontology in axioms
            ],
        )
    connection.commit()


def _serving_rows(connection, schema, module):
    with connection.cursor() as cursor:
        module._populate_entity_identifier_lookup(cursor, schema)
        module._populate_entity_ontology_terms(cursor, schema)
        cursor.execute(
            sql.SQL("""SELECT term_entity_id::text,term_id,ontology_prefix,label,definition,
            synonyms,synonyms_text,term_aliases,identifiers_text,ontology_id,sources,child_count
            FROM {}.entity_ontology_term""").format(sql.Identifier(schema))
        )
        rows = cursor.fetchall()
    connection.commit()
    # Main's DISTINCT arrays/string aggregations do not define member order;
    # compare every value and duplicate token while normalizing that order.
    result = []
    for row in rows:
        row = list(row)
        row[5], row[7], row[10] = sorted(row[5]), sorted(row[7]), sorted(row[10])
        row[8] = " ".join(sorted(row[8].split()))
        result.append(json.dumps(row))
    return Counter(result)


def test_serving_catalogue_preserves_frozen_main_bare_uri_canonical_and_description_facts(
    tmp_path,
    paired_schemas,
    oracle,
    port,
):
    postgres, (main, adapted) = paired_schemas
    resources, declared, axioms, positive_keys = _serving_fixture(tmp_path)
    source_ids = {
        source: 9100 + index for index, source in enumerate(sorted({row[0] for row in declared}))
    }
    _declared_main_serving_facts(
        postgres, main, oracle, declared, axioms, positive_keys, source_ids
    )
    expected = _serving_rows(postgres, main, oracle.derive)
    assert len(expected) == 11
    with postgres.cursor() as cursor:
        cursor.executemany(
            sql.SQL("INSERT INTO {}.data_source(source_id,name) VALUES(%s,%s)")
            .format(sql.Identifier(adapted))
            .as_string(postgres),
            [(identifier, source) for source, identifier in source_ids.items()],
        )
    postgres.commit()
    with duckdb.connect() as duck:
        plan = aligned_projection.prepare_aligned_release(
            duck, resources, dimension_rows=aligned_loader._read_dimensions(postgres, adapted)
        )
        cv_label = "Cv Term Accession:OM:0204"
        cv_links = duck.execute(
            """SELECT p.entity_key,it.value FROM ap_copy_entity_identifier link
            JOIN ap_identifier it USING(identifier_id)
            JOIN ap_vocab_identifier_type ns ON ns.id=it.identifier_type_id
            JOIN ap_entity_canonical_input p ON ap_uuid(to_json(list_value('published-entity',p.entity_key)))=link.entity_id
            WHERE ns.name=?""",
            [cv_label],
        ).fetchall()
        original_cv = [
            (item["entity_key"], identifier["id"])
            for _source, _ontology, item in declared
            for identifier in item["identifiers"]
            if identifier["ns"] == "cv_term"
        ]
        generated_cv = [
            (item["entity_key"], item["identifier"])
            for _source, _ontology, item in declared
            if item["entity_key"] in positive_keys
        ]
        assert Counter(cv_links) == Counter([*original_cv, *generated_cv])
        assert Counter(
            duck.execute(
                """SELECT entity_key,alias_value FROM ap_ontology_lookup_alias
            WHERE namespace_name=?""",
                [cv_label],
            ).fetchall()
        ) == Counter(generated_cv)
        # Only the two genuine original CV occurrences may appear in EEI.
        # Compatibility aliases never become published identifier evidence.
        evidence_cv = duck.execute(
            """SELECT e.entity_key,it.value FROM ap_copy_entity_evidence_identifier link
            JOIN ap_identifier it USING(identifier_id)
            JOIN ap_vocab_identifier_type ns ON ns.id=it.identifier_type_id
            JOIN ap_entity_occurrence e USING(source_id,entity_evidence_id)
            WHERE ns.name=?""",
            [cv_label],
        ).fetchall()
        assert Counter(evidence_cv) == Counter(original_cv) and len(evidence_cv) == 2
        assert plan.compatibility["unknown_ontology_scope_statements"] == 0
        assert plan.counts["ontology_terms"] == plan.counts["relation"] == 0
        aligned_loader._write_dimensions(postgres, adapted, plan.dimensions)
        for _source_id, source in plan.sources:
            port.schema.ensure_source_partitions(postgres, schema=adapted, source=source)
        with postgres.cursor() as cursor:
            for statement in aligned_projection.companion_ddl(adapted):
                if '."annotation_quantity"' not in statement:
                    cursor.execute(statement)
        spool = tmp_path / "copy"
        spool.mkdir()
        for query in plan.queries:
            aligned_loader._copy(postgres, adapted, query, duck, spool)
        postgres.commit()
    actual = _serving_rows(postgres, adapted, port.derive)
    assert actual == expected
    rows = {json.loads(row)[0]: json.loads(row) for row in actual}
    assert len(rows) == 11
    for _source, _ontology, item in declared:
        key = item["entity_key"]
        if key == "negative-name":
            assert str(_entity_id(key)) not in rows
            continue
        row = rows[str(_entity_id(key))]
        expected_id = (
            "CHEBI:9991"
            if key == "negative-protein"
            else ("CHEBI:20" if key == "chebi:serving:0" else "CHEBI:40")
            if key.startswith("chebi:")
            else item["identifier"]
        )
        assert row[1] == expected_id
        assert row[3] == "label " + key
        assert row[4] == "description " + key
    # A second ontology root for one canonical term must not multiply the
    # eligible entity/value CV proof or change its main row ranking. The fact
    # is added equally to the two local schemas without any extra source build.
    with postgres.cursor() as cursor:
        for schema, predicate in ((main, "is_a"), (adapted, "subclass_of")):
            cursor.execute(
                sql.SQL(
                    "SELECT relation_predicate_id FROM {}.vocab_relation_predicate WHERE name=%s"
                ).format(sql.Identifier(schema)),
                (predicate,),
            )
            predicate_id = cursor.fetchone()[0]
            cursor.execute(
                sql.SQL("INSERT INTO {}.entity_ontology_relation VALUES(%s,%s,%s,%s,%s)").format(
                    sql.Identifier(schema)
                ),
                (
                    source_ids["go"],
                    str(_entity_id("mondo:serving:0")),
                    predicate_id,
                    str(_entity_id("go:serving:0")),
                    "gene_ontology",
                ),
            )
    postgres.commit()
    scoped_expected = _serving_rows(postgres, main, oracle.derive)
    scoped_actual = _serving_rows(postgres, adapted, port.derive)
    assert scoped_actual == scoped_expected and len(scoped_actual) == 12
    shared = [
        json.loads(row)
        for row in scoped_actual
        if json.loads(row)[0] == str(_entity_id("mondo:serving:0"))
    ]
    assert {row[9] for row in shared} == {"mondo", "gene_ontology"}
    assert all(row[1] == "0000001" for row in shared)
    # Reader-only adversarial lookup rows: an alias on another entity, and an
    # alias with a different value on this entity, must not make a raw protein
    # name identifier's description eligible. These do not alter published inputs
    # or their evidence links, and exercise the exact EXISTS key/value boundary.
    with postgres.cursor() as cursor:
        cursor.execute(
            sql.SQL("SELECT identifier_type_id FROM {}.vocab_identifier_type WHERE name=%s").format(
                sql.Identifier(adapted)
            ),
            ("Cv Term Accession:OM:0204",),
        )
        cv_type = cursor.fetchone()[0]
        cursor.execute(
            sql.SQL("SELECT count(*) FROM {}.entity_evidence_identifier").format(
                sql.Identifier(adapted)
            )
        )
        original_evidence_links = cursor.fetchone()[0]
        for key, value in (
            ("negative-name", "different-CV-value"),
            ("mondo:serving:0", "random ontology name"),
        ):
            alias_id = str(uuid.uuid5(uuid.NAMESPACE_URL, "adversarial lookup:" + value))
            cursor.execute(
                sql.SQL(
                    "INSERT INTO {}.identifier_evidence(identifier_id,identifier_type_id,value) VALUES(%s,%s,%s)"
                ).format(sql.Identifier(adapted)),
                (alias_id, cv_type, value),
            )
            cursor.execute(
                sql.SQL("INSERT INTO {}.entity_identifier VALUES(%s,%s,%s)").format(
                    sql.Identifier(adapted)
                ),
                (
                    source_ids["go" if key == "negative-name" else "mondo"],
                    str(_entity_id(key)),
                    alias_id,
                ),
            )
    postgres.commit()
    mutated = {
        json.loads(row)[0]: json.loads(row) for row in _serving_rows(postgres, adapted, port.derive)
    }
    assert mutated[str(_entity_id("negative-name"))][1] == "different-CV-value"
    assert mutated[str(_entity_id("negative-name"))][4] is None
    with postgres.cursor() as cursor:
        cursor.execute(
            sql.SQL("SELECT count(*) FROM {}.entity_evidence_identifier").format(
                sql.Identifier(adapted)
            )
        )
        assert cursor.fetchone()[0] == original_evidence_links


def test_scalar_label_lookup_fallback_is_owned_truthful_and_has_no_new_occurrences(tmp_path):
    swiss = []
    for index, label in enumerate(
        (
            "Genuine Swiss label",
            "SLM:000002",
            "scalar-key:2",
            "SWISS-INCHI-3",
            "Protein label",
            None,
            "Other fallback label",
        )
    ):
        namespace, value, entity_type = (
            ("uniprot", "P33333", "protein")
            if index == 4
            else ("inchikey", "SWISS-INCHI-" + str(index), "chemical_entity")
        )
        item = _term(
            "scalar-key:" + str(index),
            value,
            namespace=namespace,
            label=label,
            entity_type=entity_type,
        )
        item["identifiers"] = [_identifier(value, namespace=namespace)]
        if index != 4:
            item["identifiers"].append(
                _identifier("SLM:" + str(index + 1).zfill(6), namespace="swisslipids")
            )
        if index == 0:
            # Main filters these missing/empty identifier values before label
            # selection, so they cannot suppress a genuine scalar fallback.
            item["identifiers"].extend([_identifier(None), _identifier("")])
        swiss.append(item)
    swiss_statements = []
    for index, item in enumerate(swiss[1:], 1):
        statement = _statement(item, swiss[0], kind="ontology")
        statement.update(relation_key="swiss-scalar:" + str(index), predicate="related_to")
        swiss_statements.append(statement)
    shared = deepcopy(swiss[-1])
    shared["label"] = "Other scalar label"
    shared["identifiers"].append(_identifier("Existing source Name", namespace="name"))
    go = [
        _term("go-scalar:0", "GO:12345", label="Genuine GO label"),
        _term("go-scalar:1", "GO:23456", label="GO:23456"),
    ]
    for item in go:
        item["identifiers"] = [_identifier(item["identifier"], namespace="go")]
    go_statement = _statement(go[0], go[1], kind="ontology")
    go_statement.update(relation_key="go-scalar:ontology", predicate="subclass_of")
    unknown = _term("unknown-scalar", "GO:99999", label="Unowned scalar label")
    unknown["identifiers"] = [_identifier(unknown["identifier"], namespace="go")]
    groups = [
        ("swisslipids", swiss, swiss_statements),
        ("chebi", [shared], []),
        ("go", go, [go_statement]),
        ("unknown-source", [unknown], []),
    ]
    assert sum(len(entities) + len(statements) for _source, entities, statements in groups) == 18
    selected = []
    for source, entities, statements in groups:
        directory = tmp_path / source
        write_fixture(directory, entities=entities, relations=statements, payloads=[])
        selected.append(SimpleNamespace(source=source, version="label-v1", directory=directory))
    with duckdb.connect() as duck:
        plan = aligned_projection.prepare_aligned_release(
            duck, selected, retain_published_provenance=True
        )
        names = Counter(
            duck.execute("""SELECT link.entity_id,i.value
            FROM ap_copy_entity_identifier link JOIN ap_identifier i USING(identifier_id)
            JOIN ap_vocab_identifier_type ns ON ns.id=i.identifier_type_id
            WHERE ns.name='Name:OM:0202'""").fetchall()
        )
        assert names == Counter(
            [
                (_entity_id("scalar-key:0"), "Genuine Swiss label"),
                (_entity_id("scalar-key:6"), "Existing source Name"),
                (_entity_id("go-scalar:0"), "Genuine GO label"),
            ]
        )
        evidence_names = Counter(
            duck.execute("""SELECT e.entity_key,i.value
            FROM ap_copy_entity_evidence_identifier link JOIN ap_identifier i USING(identifier_id)
            JOIN ap_vocab_identifier_type ns ON ns.id=i.identifier_type_id
            JOIN ap_entity_occurrence e USING(source_id,entity_evidence_id)
            WHERE ns.name='Name:OM:0202'""").fetchall()
        )
        assert evidence_names == Counter([("scalar-key:6", "Existing source Name")])
        assert plan.counts["parquet_identifier_occurrence"] == sum(
            len(item["identifiers"])
            for _source, entities, _statements in groups
            for item in entities
        )
