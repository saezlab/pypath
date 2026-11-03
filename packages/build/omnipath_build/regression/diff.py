"""Compare two result sets that resolved the same stored observations.

Outputs (in the output directory):

``differences.parquet``  every observation that is not ``same``: ids, outcomes,
                         accepted entities (and labels), protein/gene fields,
                         its votes, occurrence count and example datasets
``crosstab.parquet``     (resource, library, outcome_a, outcome_b) counts
``change_types.parquet`` (resource, library, change_type) counts
``summary.md``           short human summary with examples
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import duckdb

from .store import list_resources, resource_dir, sample_select_sql

# change_type values that do not count as "accepted entities differ"
NOT_ENTITY_CHANGES = ("same", "same_entities_other_fields", "missing_result")

CHANGE_TYPES = (
    "same",
    "same_entities_other_fields",
    "resolved_to_none",
    "none_to_resolved",
    "b_adds_entities",
    "b_drops_entities",
    "overlap",
    "different_entities",
    "missing_result",
)

_JOINED = """
SELECT '{resource}' AS resource, q.input_id, q.library, q.entity_type, q.namespace,
  q.identifier, q.taxon,
  coalesce(occ.occurrences, 0) AS occurrences, occ.datasets,
  coalesce(a.outcome, 'MISSING') AS outcome_a, coalesce(b.outcome, 'MISSING') AS outcome_b,
  coalesce(list_sort(a.entities), []::VARCHAR[]) AS entities_a,
  coalesce(list_sort(b.entities), []::VARCHAR[]) AS entities_b,
  a.entity_labels AS labels_a, b.entity_labels AS labels_b,
  a.protein_entity_id AS protein_a, b.protein_entity_id AS protein_b,
  a.gene_mapping_status AS gene_status_a, b.gene_mapping_status AS gene_status_b,
  coalesce(list_sort(a.gene_candidates), []::VARCHAR[]) AS gene_candidates_a,
  coalesce(list_sort(b.gene_candidates), []::VARCHAR[]) AS gene_candidates_b
FROM (SELECT * FROM read_parquet('{queries}'){restrict}) q
LEFT JOIN read_parquet('{a}') a ON a.input_id = q.input_id
LEFT JOIN read_parquet('{b}') b ON b.input_id = q.input_id
LEFT JOIN (
  SELECT input_id, sum(n)::BIGINT AS occurrences,
         list_slice(list(dataset ORDER BY n DESC, dataset), 1, 5) AS datasets
  FROM {occurrences} GROUP BY input_id
) occ ON occ.input_id = q.input_id
"""

_CHANGE_TYPE = """
CASE
  WHEN outcome_a = 'MISSING' OR outcome_b = 'MISSING' THEN 'missing_result'
  WHEN entities_a = entities_b THEN
    CASE WHEN coalesce(protein_a, '') = coalesce(protein_b, '')
          AND coalesce(gene_status_a, '') = coalesce(gene_status_b, '')
          AND gene_candidates_a = gene_candidates_b
         THEN 'same' ELSE 'same_entities_other_fields' END
  WHEN len(entities_a) > 0 AND len(entities_b) = 0 THEN 'resolved_to_none'
  WHEN len(entities_a) = 0 THEN 'none_to_resolved'
  WHEN len(list_intersect(entities_a, entities_b)) = len(entities_a) THEN 'b_adds_entities'
  WHEN len(list_intersect(entities_a, entities_b)) = len(entities_b) THEN 'b_drops_entities'
  WHEN len(list_intersect(entities_a, entities_b)) > 0 THEN 'overlap'
  ELSE 'different_entities'
END
"""

_VOTE = """
ns || ':' || identifier
  || CASE WHEN scope <> '' THEN ' scope=' || scope ELSE '' END
  || CASE WHEN anchor <> '' THEN ' anchor' ELSE '' END
  || CASE WHEN route <> 1 THEN ' route=' || route::VARCHAR ELSE '' END
  || CASE WHEN gene_only THEN ' gene_only' ELSE '' END
"""


def _sql_path(path: Path) -> str:
    return str(path).replace("'", "''")


def compare(
    observations: str | Path,
    a_dir: str | Path,
    b_dir: str | Path,
    output: str | Path,
    *,
    resources: list[str] | None = None,
    memory_limit: str = "2GB",
    threads: int = 2,
    examples: int = 3,
    sample: int | None = None,
) -> dict[str, Any]:
    """Diff two result sets; writes the output files and returns the summary dict."""
    observations, a_dir, b_dir, output = map(Path, (observations, a_dir, b_dir, output))
    output.mkdir(parents=True, exist_ok=True)
    selected = resources or sorted(
        set(list_resources(observations))
        & {p.parent.name for p in a_dir.glob("*/results.parquet")}
        & {p.parent.name for p in b_dir.glob("*/results.parquet")}
    )
    db_path = output / ".diff.duckdb"
    db_path.unlink(missing_ok=True)
    con = duckdb.connect(str(db_path))
    con.execute(f"SET memory_limit='{memory_limit}'")
    con.execute(f"SET threads={int(threads)}")
    con.execute("SET preserve_insertion_order=false")
    created = False
    skipped: dict[str, str] = {}
    for resource in selected:
        obs = resource_dir(observations, resource)
        files = dict(
            queries=obs / "queries.parquet",
            votes=obs / "votes.parquet",
            a=resource_dir(a_dir, resource) / "results.parquet",
            b=resource_dir(b_dir, resource) / "results.parquet",
        )
        absent = [name for name, path in files.items() if not path.is_file()]
        if absent:
            skipped[resource] = "missing " + ", ".join(absent)
            continue
        occ_path = obs / "occurrences.parquet"
        occurrences = (
            f"read_parquet('{_sql_path(occ_path)}')"
            if occ_path.is_file()
            else "(SELECT NULL::VARCHAR AS input_id, NULL::VARCHAR AS dataset, NULL::BIGINT AS n "
            "WHERE false)"
        )
        joined = _JOINED.format(
            resource=resource.replace("'", "''"),
            queries=_sql_path(files["queries"]),
            a=_sql_path(files["a"]),
            b=_sql_path(files["b"]),
            occurrences=occurrences,
            restrict=(
                f" WHERE input_id IN ({sample_select_sql(files['queries'], sample)})"
                if sample
                else ""
            ),
        )
        con.execute(
            f"CREATE OR REPLACE TEMP TABLE joined AS SELECT *, {_CHANGE_TYPE} AS change_type FROM ({joined})"
        )
        votes = f"""
            SELECT input_id, list({_VOTE} ORDER BY ordinal) AS votes
            FROM read_parquet('{_sql_path(files["votes"])}')
            WHERE input_id IN (SELECT input_id FROM joined WHERE change_type <> 'same')
            GROUP BY input_id"""
        statements = {
            "differences": f"""SELECT j.*, coalesce(v.votes, []::VARCHAR[]) AS votes
                FROM joined j LEFT JOIN ({votes}) v ON v.input_id = j.input_id
                WHERE j.change_type <> 'same'""",
            "crosstab": """SELECT resource, library, outcome_a, outcome_b,
                count(*)::BIGINT AS n, sum(occurrences)::BIGINT AS occurrences,
                count(*) FILTER (WHERE change_type NOT IN ('same', 'same_entities_other_fields',
                    'missing_result'))::BIGINT AS n_entities_differ
                FROM joined GROUP BY ALL""",
            "change_types": """SELECT resource, library, change_type,
                count(*)::BIGINT AS n, sum(occurrences)::BIGINT AS occurrences
                FROM joined GROUP BY ALL""",
        }
        for table, sql in statements.items():
            if created:
                con.execute(f"INSERT INTO {table} {sql}")
            else:
                con.execute(f"CREATE TABLE {table} AS {sql}")
        created = True
    if not created:
        con.close()
        db_path.unlink(missing_ok=True)
        raise RuntimeError(f"Nothing to compare; skipped: {skipped or 'no common resources'}")
    copy = "(FORMAT PARQUET, COMPRESSION ZSTD)"
    con.execute(
        f"COPY (SELECT * FROM differences ORDER BY resource, change_type, input_id) "
        f"TO '{_sql_path(output / 'differences.parquet')}' {copy}"
    )
    con.execute(
        f"COPY (SELECT * FROM crosstab ORDER BY resource, library, outcome_a, outcome_b) "
        f"TO '{_sql_path(output / 'crosstab.parquet')}' {copy}"
    )
    con.execute(
        f"COPY (SELECT * FROM change_types ORDER BY resource, library, change_type) "
        f"TO '{_sql_path(output / 'change_types.parquet')}' {copy}"
    )
    summary = _summary(con, examples, selected, skipped)
    summary["sample_per_resource_and_library"] = sample or None
    con.close()
    db_path.unlink(missing_ok=True)
    (output / "summary.md").write_text(
        render_markdown(summary, str(a_dir), str(b_dir)), encoding="utf-8"
    )
    (output / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    return summary


def _summary(con, examples: int, resources: list[str], skipped: dict[str, str]) -> dict[str, Any]:
    change_types = con.execute(
        "SELECT resource, library, change_type, n, occurrences FROM change_types "
        "ORDER BY resource, library, change_type"
    ).fetchall()
    crosstab = con.execute(
        "SELECT resource, library, outcome_a, outcome_b, n, occurrences, n_entities_differ "
        "FROM crosstab ORDER BY resource, library, n DESC, outcome_a, outcome_b"
    ).fetchall()
    sample = con.execute(
        f"""SELECT resource, library, change_type, input_id, namespace, identifier, taxon,
                   outcome_a, outcome_b, entities_a, entities_b, votes
            FROM (SELECT *, row_number() OVER (
                      PARTITION BY change_type ORDER BY occurrences DESC, resource, input_id) AS r
                  FROM differences WHERE change_type <> 'same_entities_other_fields')
            WHERE r <= {int(examples)} ORDER BY change_type, r"""
    ).fetchall()
    return dict(
        resources=resources,
        skipped=skipped,
        change_types=[
            dict(resource=r, library=lib, change_type=c, n=n, occurrences=o)
            for r, lib, c, n, o in change_types
        ],
        crosstab=[
            dict(
                resource=r,
                library=lib,
                outcome_a=a,
                outcome_b=b,
                n=n,
                occurrences=o,
                n_entities_differ=d,
            )
            for r, lib, a, b, n, o, d in crosstab
        ],
        examples=[
            dict(
                resource=r,
                library=lib,
                change_type=c,
                input_id=i,
                namespace=ns,
                identifier=ident,
                taxon=t,
                outcome_a=oa,
                outcome_b=ob,
                entities_a=ea,
                entities_b=eb,
                votes=v,
            )
            for r, lib, c, i, ns, ident, t, oa, ob, ea, eb, v in sample
        ],
    )


def _table(header: list[str], rows: list[list[Any]]) -> list[str]:
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    out += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return out


def render_markdown(summary: dict[str, Any], a: str, b: str) -> str:
    lines = ["# Resolution regression diff", "", f"- A: `{a}`", f"- B: `{b}`", ""]
    if summary.get("sample_per_resource_and_library"):
        lines += [
            f"Sample: {summary['sample_per_resource_and_library']} observations per resource "
            "and library (deterministic, md5 of the observation fingerprint).",
            "",
        ]
    if summary["skipped"]:
        lines += ["Skipped: " + "; ".join(f"{k} ({v})" for k, v in summary["skipped"].items()), ""]

    by_group: dict[tuple[str, str], dict[str, int]] = defaultdict(dict)
    for row in summary["change_types"]:
        by_group[(row["resource"], row["library"])][row["change_type"]] = row["n"]
    present = [c for c in CHANGE_TYPES if any(c in g for g in by_group.values())]
    lines += ["## Distinct observations by change type", ""]
    rows = []
    totals: dict[str, int] = defaultdict(int)
    for (resource, library), counts in sorted(by_group.items()):
        total = sum(counts.values())
        rows.append([resource, library, total] + [counts.get(c, 0) for c in present])
        for c in present:
            totals[c] += counts.get(c, 0)
    rows.append(["**all**", "", sum(totals.values())] + [totals[c] for c in present])
    lines += _table(["resource", "library", "observations", *present], rows)
    lines += [
        "",
        "`same`: same accepted entities, protein and gene fields. `same_entities_other_fields`: "
        "same entities but different protein/gene fields. Entity changes are named from A to B; "
        "`b_adds_entities` means A's entities are a strict subset of B's.",
        "",
    ]

    lines += ["## Outcome cross-tab (rows: A, columns: B), all resources", ""]
    for library in sorted({r["library"] for r in summary["crosstab"]}):
        cells: dict[tuple[str, str], int] = defaultdict(int)
        for r in summary["crosstab"]:
            if r["library"] == library:
                cells[(r["outcome_a"], r["outcome_b"])] += r["n"]
        outcomes_a = sorted({a_ for a_, _ in cells})
        outcomes_b = sorted({b_ for _, b_ in cells})
        lines += [f"### {library}", ""]
        lines += _table(
            ["A \\ B", *outcomes_b],
            [[oa, *[cells.get((oa, ob), 0) for ob in outcomes_b]] for oa in outcomes_a],
        )
        lines.append("")

    lines += ["## Largest outcome transitions per resource and library", ""]
    per_group: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for r in summary["crosstab"]:
        if r["outcome_a"] != r["outcome_b"] or r["n_entities_differ"]:
            per_group[(r["resource"], r["library"])].append(r)
    rows = []
    for (resource, library), items in sorted(per_group.items()):
        for r in sorted(items, key=lambda r: -r["n"])[:5]:
            rows.append(
                [resource, library, r["outcome_a"], r["outcome_b"], r["n"], r["n_entities_differ"]]
            )
    lines += (
        _table(["resource", "library", "A", "B", "observations", "entities differ"], rows)
        if rows
        else ["No outcome or entity differences."]
    )
    lines.append("")

    lines += ["## Examples (most frequent first)", ""]
    for ex in summary["examples"]:
        lines.append(
            f"- **{ex['change_type']}** `{ex['resource']}` {ex['library']} "
            f"`{ex['namespace']}:{ex['identifier']}` taxon={ex['taxon'] or '-'} "
            f"{ex['outcome_a']} {ex['entities_a']} -> {ex['outcome_b']} {ex['entities_b']}; "
            f"votes: {'; '.join(ex['votes'][:8])}"
        )
    return "\n".join(lines) + "\n"
