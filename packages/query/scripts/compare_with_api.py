"""Compare the client's prebuilt queries with the API's queries on an installed directory.

    OP_THREADS=2 nice -n 19 uv run python packages/query/scripts/compare_with_api.py <data_dir> [relations|search|entity ...]

Runs each case through ``omnipath_client`` (the macros in the database) and through
``omnipath_api.queries`` on the same Parquet files: ``export_slice`` for ``relations``,
``search_entities_api`` as the client shaped it before for ``search``, and
``get_entity_details`` for ``entity``. Relations are compared as multisets (neither
orders them), search results in order, entity details as whole JSON values. Times are
the median of three warm calls (the API caches entity details, so for ``entity``
its first call is shown). ``OP_THREADS`` limits both to that many threads. Exits 1 when
an answer differs; a tie order that the API's page size decides is reported but not
counted (see the client README).
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import statistics
import sys
import time

import pyarrow as pa
import pyarrow.parquet as pq

from omnipath_api.queries import ParquetQueries
import omnipath_client as op

THREADS = int(os.environ.get("OP_THREADS", "0")) or None


def timed(run, repeats=3):
    start = time.perf_counter()
    result = run()
    first = time.perf_counter() - start
    rest = []
    for _ in range(repeats):
        start = time.perf_counter()
        run()
        rest.append(time.perf_counter() - start)
    return result, first, statistics.median(rest)


def canonical(table: pa.Table) -> list[str]:
    return sorted(json.dumps(row, sort_keys=True, default=str) for row in table.to_pylist())


def relation_cases(db) -> dict[str, dict]:
    def one(sql):
        return (db.sql(sql).fetchone() or [None])[0]

    tp53_reference = one(
        "SELECT reference_entity_key FROM entity WHERE label = 'TP53' "
        "AND reference_entity_key IS NOT NULL ORDER BY ALL LIMIT 1"
    )
    protein = one(
        "SELECT subject_molecular_form.protein_entity_key FROM relation_evidence "
        "WHERE subject_molecular_form.protein_entity_key IS NOT NULL "
        "GROUP BY 1 ORDER BY count(*) DESC, 1 LIMIT 1"
    )
    isoform = one(
        "SELECT f.isoform_identifier.ns || ':' || f.isoform_identifier.id FROM ("
        "SELECT object_molecular_form AS f FROM relation_evidence) "
        "WHERE f.isoform_identifier.id IS NOT NULL GROUP BY 1 ORDER BY count(*) DESC, 1 LIMIT 1"
    )
    installed = {r for (r,) in db.sql("SELECT resource_id FROM resources").fetchall()}
    cases = {
        "TP53": {"entity_ids": ["TP53"]},
        "TP53 as subject": {"entity_ids": ["TP53"], "gene_mode": "source"},
        "TP53 and MDM2, both ends": {"entity_ids": ["TP53", "MDM2"], "gene_mode": "both"},
        "uniprot:P04637": {"entity_ids": ["uniprot:P04637"]},
        "ATP (CHEBI:30616)": {"entity_ids": ["CHEBI:30616"]},
        "glucose": {"entity_ids": ["glucose"]},
        "EGFR": {"entity_ids": ["EGFR"]},
        "EGFR protein-chemical": {
            "entity_ids": ["EGFR"], "interaction_types": ["protein"], "partner_types": ["chemical_entity"],
        },
        "EGFR, biolink:interacts_with": {"entity_ids": ["EGFR"], "predicates": ["biolink:interacts_with"]},
        "imatinib": {"entity_ids": ["imatinib"]},
        "no match": {"entity_ids": ["no-such-entity-xyz"]},
    }
    if tp53_reference:
        cases["gene group of TP53"] = {"entity_ids": [f"gene:{tp53_reference}"]}
        cases["reference key of TP53"] = {"reference_entity_keys": [tp53_reference]}
    if protein:
        cases["molecular form (protein key)"] = {"protein_entity_keys": [protein]}
    if isoform:
        cases["molecular form (isoform), target"] = {
            "isoform_identifiers": [isoform], "molecular_endpoint_mode": "target",
        }
    for resource, extra in (
        ("signor", {}),
        ("signor", {"signs": [-1], "taxonomy_ids": ["9606"]}),
        ("signor", {"object_direction_qualifier": ["increased"]}),
        ("trrust", {"min_evidence_count": 2}),
        ("collectri", {"is_directed": True}),
        ("biogrid", {"entity_ids": ["TP53"]}),
        ("chembl", {"entity_ids": ["EGFR"], "taxonomy_ids": ["9606"]}),
    ):
        if resource in installed:
            label = f"{resource} " + ", ".join(f"{k}={v}" for k, v in extra.items())
            cases[label.strip()] = {"sources": [resource], **extra}
    return cases


def compare_relations(client, api) -> int:
    print(f"{'relations':42} {'rows':>7}  {'API s':>8}  {'client s':>8}  answer")
    failures = 0
    for name, filters in relation_cases(client).items():
        data, _, api_s = timed(lambda: api.export_slice(filters, None)[0])
        expected = pq.read_table(pa.BufferReader(data))
        if expected.num_rows >= op.client.RELATION_EXPORT_LIMIT:
            print(f"{name[:42]:42} skipped: the API's {expected.num_rows} rows are its limit, in no order")
            continue
        filters = dict(filters)
        entities = filters.pop("entity_ids", None)
        actual, _, client_s = timed(lambda: client.relations(entities, **filters).to_arrow_table())
        if actual.schema != expected.schema:
            verdict = "DIFFERENT COLUMNS"
        elif canonical(actual) != canonical(expected):
            verdict = f"DIFFERENT ROWS ({expected.num_rows} vs {actual.num_rows})"
        else:
            verdict = "same"
        failures += verdict != "same"
        print(f"{name[:42]:42} {expected.num_rows:>7}  {api_s:8.2f}  {client_s:8.2f}  {verdict}", flush=True)
    return failures


SEARCH_CASES = [
    ("TP53", {}), ("tp5", {}), ("EGFR", {}), ("glucose", {}), ("ATP", {}), ("CHEBI:30616", {}),
    ("uniprot:P04637", {}), ("entrez:7157", {}), ("kinase", {}), ("kinase", {"limit": 100}),
    ("a", {}), ("insulin receptor", {}), ("ase activity", {}), ("imatinib", {}),
    ("xyzzy-nothing", {}), ("", {}), ("gluc", {"entity_types": ["chemical_entity"]}),
    ("TP53", {"taxonomy_ids": ["9606"]}), ("TP53", {"entity_ids": ["MDM2", "TP53"]}),
]


def api_search(api, query, limit=20, **filters):
    """The API's search shaped as the client returned it before the macros."""
    found = api.search_entities_api(
        query, filters={k: v if isinstance(v, list) else [v] for k, v in filters.items()}, limit=limit
    )
    return [
        (e["entityPk"], e.get("label"), e.get("entityType"),
         # the API's namespace is a str enum; the client returns its value
         getattr(e.get("canonicalIdentifierType"), "value", e.get("canonicalIdentifierType")),
         e.get("canonicalIdentifier"), e.get("taxonomyId"), e.get("relationCount"),
         e.get("referenceEntityKey"), e.get("geneReferenceKeys") or [])
        for e in found["entities"]
    ]


def compare_search(client, api) -> int:
    print(f"{'search':42} {'rows':>7}  {'API s':>8}  {'client s':>8}  answer")
    failures = 0
    cases = list(SEARCH_CASES)
    key = client.sql("SELECT entity_key FROM entity WHERE label = 'TP53' ORDER BY relation_count DESC LIMIT 1").fetchone()
    if key:
        cases.append((key[0], {}))
    if "signor" in client._catalog:
        cases.append(("TP53", {"sources": ["signor"]}))
    for query, filters in cases:
        want, _, api_s = timed(lambda: api_search(api, query, **filters))
        got, _, client_s = timed(lambda: client.search(query, **filters).fetchall())
        got = [tuple(row) for row in got]
        if want == got:
            verdict = "same"
        elif sorted(map(repr, want)) == sorted(map(repr, got)):
            verdict = "same rows, other tie order"
        else:
            verdict = f"DIFFERENT ({len(want)} vs {len(got)} rows)"
            failures += 1
        name = repr(query if len(query) < 20 else query[:12] + "...") + (f" {filters}" if filters else "")
        print(f"{name[:42]:42} {len(want):>7}  {api_s:8.2f}  {client_s:8.2f}  {verdict}", flush=True)
    return failures


def compare_entity(client, api) -> int:
    print(f"{'entity':42} {'rel.':>7}  {'API s':>8}  {'client s':>8}  answer")
    failures = 0
    cases = []
    for query in ("TP53", "EGFR", "glucose", "ATP", "CHEBI:30616", "kinase", "insulin receptor", "imatinib"):
        row = client.search(query).fetchone()
        if row:
            cases.append((f"search {query!r}", row[0]))
    for name, sql in (
        ("most structural relations", "SELECT key FROM relation_endpoint WHERE key_kind = 'entity' "
         "AND predicate IN ('has_input', 'has_output', 'enabled_by', 'has_member', 'has_part') "
         "GROUP BY key ORDER BY count(*) DESC, key LIMIT 1"),
        ("most associations", "SELECT key FROM relation_endpoint WHERE key_kind = 'entity' "
         "AND category = 'association' GROUP BY key ORDER BY count(*) DESC, key LIMIT 1"),
        ("an ontology term", "SELECT entity_key FROM entity WHERE has_hierarchy ORDER BY child_count DESC, entity_key LIMIT 1"),
        ("most annotations", "SELECT entity_key FROM entity ORDER BY annotation_count DESC, entity_key LIMIT 1"),
        ("in most resources", "SELECT entity_key FROM entity GROUP BY entity_key ORDER BY count(*) DESC, entity_key LIMIT 1"),
    ):
        row = client.sql(sql).fetchone()
        if row:
            cases.append((name, row[0]))
    cases.append(("unknown key", "0" * 64))
    for name, key in cases:
        start = time.perf_counter()
        want = api.get_entity_details(key)
        api_s = time.perf_counter() - start
        want = json.loads(json.dumps(want, default=str)) if want is not None else None
        got, _, client_s = timed(lambda: client.entity(key))
        verdict = "same" if want == got else "DIFFERENT"
        failures += verdict != "same"
        relations = len(want["relationships"]) if want else 0
        print(f"{name[:42]:42} {relations:>7}  {api_s:8.2f}  {client_s:8.2f}  {verdict}", flush=True)
    return failures


def main(data_dir: str, parts: list[str]) -> int:
    root = Path(data_dir).resolve()
    client = op.connect(root, threads=THREADS)
    api = ParquetQueries(root, memory_limit="2GB", threads=THREADS or os.cpu_count() or 2)
    failures = 0
    try:
        for part in parts or ["relations", "search", "entity"]:
            failures += {"relations": compare_relations, "search": compare_search, "entity": compare_entity}[part](client, api)
            print()
    finally:
        api.close()
        client.close()
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2:]))
