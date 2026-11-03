"""Metabolic-domain classification for ``chemical_entity`` entities.

Coarse, flat, single-valued grouping that primarily subdivides metabolites,
decided by the sources in ``metabolic_domain.yaml`` in priority order: the
first source that yields a bucket wins. Four kinds of source:

- a catalog (``data_source`` + ``domain``) gives its domain to every chemical
  that resource has;
- ``names: goslin`` gives its domain to every chemical with a name in lipid
  shorthand that Goslin parses at species level or finer;
- a ChEBI root (``chebi_root`` + ``domain``) gives its domain to every chemical
  that descends from it in the ChEBI ``subclass_of`` hierarchy;
- a class namespace (``namespace``) is read through the labels of the linked
  classes and the curated keywords: the bucket that most of the chemical's
  labels match wins, ties broken by bucket precedence.

Runs during ``derive`` (no reload), on the canonical tables.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import os
from pathlib import Path

from psycopg2 import sql
import psycopg2.extensions
from psycopg2.extras import execute_values
import yaml


CHEMICAL_ENTITY_TYPE = "chemical_entity"
_RULES_PATH = Path(__file__).with_name("metabolic_domain.yaml")
# Names worth handing to Goslin: shorthand carries a chain such as 16:0. Trivial
# names stay out ("Ethanol" parses as the fatty alcohol FOH 2:0).
_SHORTHAND_NAME = "[0-9]:[0-9]"
# Goslin takes milliseconds on a long name, and a build has some 500,000: the
# names are parsed in chunks by a pool of processes.
_PARSE_CHUNK = 5000
_PARSE_WORKERS = 16


@dataclass(frozen=True)
class MetabolicDomainStats:
    classified: int = 0
    by_default: int = 0


def _load_rules() -> dict:
    with _RULES_PATH.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def classify_metabolic_domain(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = "public",
) -> MetabolicDomainStats:
    """Seed ``vocab_metabolic_domain`` and populate ``entity.metabolic_domain_id``."""

    rules = _load_rules()
    schema_id = sql.Identifier(schema)
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL(
                "SELECT entity_type_id FROM {}.vocab_entity_type WHERE name IN ('chemical_entity','small_molecule')"
            ).format(schema_id),
        )
        rows = cur.fetchall()
        if not rows:
            return MetabolicDomainStats()
        chemical_type_id = [row[0] for row in rows]

        # Seed the controlled vocabulary from the rule file.
        domain_rows = [
            (entry["precedence"], entry["name"], entry["precedence"]) for entry in rules["domains"]
        ]
        cur.executemany(
            sql.SQL(
                """
                INSERT INTO {}.vocab_metabolic_domain
                  (metabolic_domain_id, name, precedence)
                VALUES (%s, %s, %s)
                ON CONFLICT (metabolic_domain_id) DO UPDATE
                SET name = EXCLUDED.name, precedence = EXCLUDED.precedence
                """
            )
            .format(schema_id)
            .as_string(cur.connection),
            domain_rows,
        )

        # Materialise the class namespaces read through their labels.
        cur.execute(
            "CREATE TEMP TABLE _metabolic_domain_source "
            "(namespace text NOT NULL, priority integer NOT NULL) ON COMMIT DROP"
        )
        cur.executemany(
            "INSERT INTO _metabolic_domain_source (namespace, priority) VALUES (%s, %s)",
            [
                (entry["namespace"], entry["priority"])
                for entry in rules["sources"]
                if entry.get("namespace")
            ],
        )

        _assert_domains(cur, schema_id, chemical_type_id, rules["sources"])

        # Materialise the domain -> keyword map (keywords lowercased once).
        cur.execute(
            "CREATE TEMP TABLE _metabolic_domain_keyword "
            "(domain_name text NOT NULL, keyword text NOT NULL) ON COMMIT DROP"
        )
        keyword_rows = [
            (domain_name, keyword.lower())
            for domain_name, keywords in (rules.get("keywords") or {}).items()
            for keyword in keywords
        ]
        if keyword_rows:
            cur.executemany(
                "INSERT INTO _metabolic_domain_keyword (domain_name, keyword) VALUES (%s, %s)",
                keyword_rows,
            )

        # Reset, then assign the first-source / most-specific matching bucket.
        cur.execute(
            sql.SQL(
                "UPDATE {}.entity SET metabolic_domain_id = NULL WHERE entity_type_id = ANY(%s)"
            ).format(schema_id),
            [chemical_type_id],
        )
        cur.execute(
            sql.SQL(
                """
                WITH chem AS (
                  SELECT entity_id FROM {schema}.entity
                  WHERE entity_type_id = ANY(%(chem)s)
                ),
                -- Either side: a class link is `associated_with`, which has
                -- no direction, so the chemical is its subject or its object.
                linked AS (
                  SELECT r.subject_entity_id AS entity_id, r.object_entity_id AS other_id
                  FROM {schema}.relation r
                  JOIN chem ON chem.entity_id = r.subject_entity_id
                  UNION
                  SELECT r.object_entity_id, r.subject_entity_id
                  FROM {schema}.relation r
                  JOIN chem ON chem.entity_id = r.object_entity_id
                ),
                ann AS (
                  SELECT
                    linked.entity_id,
                    min(src.priority) AS term_priority,
                    lower(o.label) AS val
                  FROM linked
                  JOIN {schema}.entity o ON o.entity_id = linked.other_id
                  JOIN {schema}.vocab_identifier_type it
                    ON it.identifier_type_id = o.canonical_identifier_type_id
                  JOIN _metabolic_domain_source src ON src.namespace = it.name
                  WHERE o.label IS NOT NULL
                  GROUP BY linked.entity_id, lower(o.label)
                ),
                matched AS (
                  SELECT
                    ann.entity_id,
                    d.metabolic_domain_id,
                    ann.term_priority,
                    d.precedence AS bucket_precedence,
                    count(DISTINCT ann.val) AS n
                  FROM ann
                  JOIN _metabolic_domain_keyword k
                    ON ann.val LIKE '%%' || k.keyword || '%%'
                  JOIN {schema}.vocab_metabolic_domain d
                    ON d.name = k.domain_name
                  GROUP BY ann.entity_id, d.metabolic_domain_id, ann.term_priority, d.precedence
                  UNION ALL
                  -- Catalogs, lipid names and ChEBI roots: the domain is asserted.
                  SELECT DISTINCT
                    a.entity_id,
                    d.metabolic_domain_id,
                    a.priority,
                    d.precedence,
                    1
                  FROM _metabolic_domain_asserted a
                  JOIN {schema}.vocab_metabolic_domain d ON d.name = a.domain
                ),
                ranked AS (
                  SELECT
                    entity_id,
                    metabolic_domain_id,
                    row_number() OVER (
                      PARTITION BY entity_id
                      ORDER BY term_priority ASC, n DESC, bucket_precedence ASC
                    ) AS rn
                  FROM matched
                )
                UPDATE {schema}.entity e
                SET metabolic_domain_id = ranked.metabolic_domain_id
                FROM ranked
                WHERE ranked.entity_id = e.entity_id
                  AND ranked.rn = 1
                """
            ).format(schema=schema_id),
            {"chem": chemical_type_id},
        )
        classified = int(cur.rowcount)

        # Fallback: any remaining chemical -> default bucket (no NULLs).
        cur.execute(
            sql.SQL(
                "SELECT metabolic_domain_id FROM {}.vocab_metabolic_domain WHERE name = %s"
            ).format(schema_id),
            [rules["default"]],
        )
        default_id = cur.fetchone()[0]
        cur.execute(
            sql.SQL(
                "UPDATE {}.entity SET metabolic_domain_id = %s "
                "WHERE entity_type_id = ANY(%s) AND metabolic_domain_id IS NULL"
            ).format(schema_id),
            [default_id, chemical_type_id],
        )
        by_default = int(cur.rowcount)
    conn.commit()
    return MetabolicDomainStats(classified=classified, by_default=by_default)


def _assert_domains(cur, schema_id, chemical_type_id, sources) -> None:
    """Stage the chemicals whose domain a catalog, a lipid name or ChEBI asserts.

    Catalogs and ChEBI first: a chemical they already place at a priority no
    lower than the names' needs no name parsed.
    """
    cur.execute(
        sql.SQL(
            "CREATE TEMP TABLE _metabolic_domain_asserted ON COMMIT DROP AS "
            "SELECT entity_id, ''::text AS domain, 0 AS priority FROM {}.entity WITH NO DATA"
        ).format(schema_id)
    )
    for entry in sources:
        if entry.get("data_source"):
            cur.execute(
                sql.SQL(
                    """
                    INSERT INTO _metabolic_domain_asserted
                    SELECT DISTINCT r.entity_id, %(domain)s, %(priority)s
                    FROM {schema}.entity_evidence_resolution r
                    JOIN {schema}.data_source ds ON ds.source_id = r.source_id
                    JOIN {schema}.entity e ON e.entity_id = r.entity_id
                    WHERE ds.name = %(source)s AND e.entity_type_id = ANY(%(chem)s)
                    """
                ).format(schema=schema_id),
                {
                    "domain": entry["domain"],
                    "priority": entry["priority"],
                    "source": entry["data_source"],
                    "chem": chemical_type_id,
                },
            )
    chebi_roots = [entry for entry in sources if entry.get("chebi_root")]
    if chebi_roots:
        _assert_chebi_descendants(cur, schema_id, chemical_type_id, chebi_roots)
    cur.execute("ANALYZE _metabolic_domain_asserted")
    for entry in sources:
        if entry.get("names") == "goslin":
            _assert_lipid_names(cur, schema_id, chemical_type_id, entry)


def _assert_chebi_descendants(cur, schema_id, chemical_type_id, roots) -> None:
    """Chemicals that descend from a ChEBI root through ``subclass_of``."""
    cur.execute(
        "CREATE TEMP TABLE _metabolic_domain_chebi_root "
        "(identifier text NOT NULL, domain text NOT NULL, priority integer NOT NULL) "
        "ON COMMIT DROP"
    )
    cur.executemany(
        "INSERT INTO _metabolic_domain_chebi_root VALUES (%s, %s, %s)",
        [(entry["chebi_root"], entry["domain"], entry["priority"]) for entry in roots],
    )
    cur.execute(
        sql.SQL(
            """
            INSERT INTO _metabolic_domain_asserted
            WITH RECURSIVE edge AS MATERIALIZED (
              SELECT DISTINCT r.subject_entity_id AS child, r.object_entity_id AS parent
              FROM {schema}.entity_ontology_relation r
              JOIN {schema}.vocab_relation_predicate p
                ON p.relation_predicate_id = r.predicate_id
              WHERE r.ontology_id = 'chebi' AND p.name = 'subclass_of'
            ),
            descendant (entity_id, domain, priority) AS (
              SELECT e.entity_id, root.domain, root.priority
              FROM _metabolic_domain_chebi_root root
              JOIN {schema}.entity e ON e.canonical_identifier = root.identifier
              UNION
              SELECT edge.child, descendant.domain, descendant.priority
              FROM descendant
              JOIN edge ON edge.parent = descendant.entity_id
            )
            SELECT DISTINCT d.entity_id, d.domain, d.priority
            FROM descendant d
            JOIN {schema}.entity e ON e.entity_id = d.entity_id
            WHERE e.entity_type_id = ANY(%(chem)s)
            """
        ).format(schema=schema_id),
        {"chem": chemical_type_id},
    )


def _assert_lipid_names(cur, schema_id, chemical_type_id, entry) -> None:
    """Chemicals with a name Goslin parses as a lipid at species level or finer."""
    cur.execute(
        sql.SQL(
            """
            SELECT DISTINCT n.entity_id, n.name
            FROM {schema}.entity_name n
            JOIN {schema}.entity e ON e.entity_id = n.entity_id
            WHERE e.entity_type_id = ANY(%(chem)s)
              AND n.name ~ %(shorthand)s
              AND NOT EXISTS (
                SELECT 1 FROM _metabolic_domain_asserted a
                WHERE a.entity_id = n.entity_id AND a.priority <= %(priority)s
              )
            """
        ).format(schema=schema_id),
        {"chem": chemical_type_id, "shorthand": _SHORTHAND_NAME, "priority": entry["priority"]},
    )
    rows = cur.fetchall()
    names = sorted({name for _, name in rows})
    chunks = [names[i : i + _PARSE_CHUNK] for i in range(0, len(names), _PARSE_CHUNK)]
    if len(chunks) > 1:
        workers = min(_PARSE_WORKERS, len(chunks), os.cpu_count() or 1)
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(_parse_lipid_names, chunks))
    else:
        results = [_parse_lipid_names(chunk) for chunk in chunks]
    lipid_names = {
        name for chunk, ok in zip(chunks, results) for name, hit in zip(chunk, ok) if hit
    }
    lipids = {entity_id for entity_id, name in rows if name in lipid_names}
    execute_values(
        cur,
        "INSERT INTO _metabolic_domain_asserted (entity_id, domain, priority) VALUES %s",
        [(entity_id, entry["domain"], entry["priority"]) for entity_id in lipids],
        page_size=10000,
    )


def _parse_lipid_names(names: list[str]) -> list[bool]:
    """Whether Goslin parses each name as a lipid at species level or finer."""
    from omnipath_resolver.goslin import normalize

    return [normalize(name)["status"] == "parsed" for name in names]
