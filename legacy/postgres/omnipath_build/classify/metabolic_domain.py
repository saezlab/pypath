"""Metabolic-domain classification for ``Chemical:OM:0037`` entities.

Coarse, flat, single-valued grouping that primarily subdivides metabolites. A
chemical's domain is read from its taxonomy/pathway annotations (HMDB ClassyFire
super/class, KEGG pathway, Recon3D subsystem) via the curated keyword rules in
``metabolic_domain.yaml``: the first source (in priority order) that yields a
matching bucket wins, ties broken by bucket precedence. Runs during ``derive``
(no reload). Annotations reach the canonical entity through
``entity_evidence_annotation`` -> ``entity_evidence_resolution``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from psycopg2 import sql
import psycopg2.extensions
import yaml


CHEMICAL_ENTITY_TYPE = 'Chemical:OM:0037'
_RULES_PATH = Path(__file__).with_name('metabolic_domain.yaml')


@dataclass(frozen=True)
class MetabolicDomainStats:
    classified: int = 0
    by_default: int = 0


def _load_rules() -> dict:
    with _RULES_PATH.open('r', encoding='utf-8') as handle:
        return yaml.safe_load(handle)


def classify_metabolic_domain(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
) -> MetabolicDomainStats:
    """Seed ``vocab_metabolic_domain`` and populate ``entity.metabolic_domain_id``."""

    rules = _load_rules()
    schema_id = sql.Identifier(schema)
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL(
                'SELECT entity_type_id FROM {}.vocab_entity_type WHERE name = %s'
            ).format(schema_id),
            [CHEMICAL_ENTITY_TYPE],
        )
        row = cur.fetchone()
        if row is None:
            return MetabolicDomainStats()
        chemical_type_id = row[0]

        # Seed the controlled vocabulary from the rule file.
        domain_rows = [
            (entry['precedence'], entry['name'], entry['precedence'])
            for entry in rules['domains']
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
            ).format(schema_id).as_string(cur.connection),
            domain_rows,
        )

        # Materialise the consulted annotation terms (with priority).
        cur.execute(
            'CREATE TEMP TABLE _metabolic_domain_term '
            '(term text PRIMARY KEY, priority integer NOT NULL) ON COMMIT DROP'
        )
        cur.executemany(
            'INSERT INTO _metabolic_domain_term (term, priority) VALUES (%s, %s) '
            'ON CONFLICT (term) DO NOTHING',
            [(entry['term'], entry['priority']) for entry in rules['source_terms']],
        )

        # Materialise the domain -> keyword map (keywords lowercased once).
        cur.execute(
            'CREATE TEMP TABLE _metabolic_domain_keyword '
            '(domain_name text NOT NULL, keyword text NOT NULL) ON COMMIT DROP'
        )
        keyword_rows = [
            (domain_name, keyword.lower())
            for domain_name, keywords in (rules.get('keywords') or {}).items()
            for keyword in keywords
        ]
        if keyword_rows:
            cur.executemany(
                'INSERT INTO _metabolic_domain_keyword (domain_name, keyword) '
                'VALUES (%s, %s)',
                keyword_rows,
            )

        # Reset, then assign the first-source / most-specific matching bucket.
        cur.execute(
            sql.SQL(
                'UPDATE {}.entity SET metabolic_domain_id = NULL '
                'WHERE entity_type_id = %s'
            ).format(schema_id),
            [chemical_type_id],
        )
        cur.execute(
            sql.SQL(
                """
                WITH chem AS (
                  SELECT entity_id FROM {schema}.entity
                  WHERE entity_type_id = %(chem)s
                ),
                ann AS (
                  SELECT
                    eer.entity_id,
                    t.priority AS term_priority,
                    lower(a.value) AS val
                  FROM {schema}.entity_evidence_annotation eea
                  JOIN {schema}.annotation a
                    ON a.annotation_key = eea.annotation_key
                  JOIN {schema}.entity_evidence_resolution eer
                    ON eer.entity_evidence_id = eea.entity_evidence_id
                  JOIN chem ON chem.entity_id = eer.entity_id
                  JOIN _metabolic_domain_term t ON t.term = a.term
                  WHERE a.value IS NOT NULL
                ),
                matched AS (
                  SELECT
                    ann.entity_id,
                    d.metabolic_domain_id,
                    ann.term_priority,
                    d.precedence AS bucket_precedence
                  FROM ann
                  JOIN _metabolic_domain_keyword k
                    ON ann.val LIKE '%%' || k.keyword || '%%'
                  JOIN {schema}.vocab_metabolic_domain d
                    ON d.name = k.domain_name
                ),
                ranked AS (
                  SELECT
                    entity_id,
                    metabolic_domain_id,
                    row_number() OVER (
                      PARTITION BY entity_id
                      ORDER BY term_priority ASC, bucket_precedence ASC
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
            {'chem': chemical_type_id},
        )
        classified = int(cur.rowcount)

        # Fallback: any remaining chemical -> default bucket (no NULLs).
        cur.execute(
            sql.SQL(
                'SELECT metabolic_domain_id FROM {}.vocab_metabolic_domain '
                'WHERE name = %s'
            ).format(schema_id),
            [rules['default']],
        )
        default_id = cur.fetchone()[0]
        cur.execute(
            sql.SQL(
                'UPDATE {}.entity SET metabolic_domain_id = %s '
                'WHERE entity_type_id = %s AND metabolic_domain_id IS NULL'
            ).format(schema_id),
            [default_id, chemical_type_id],
        )
        by_default = int(cur.rowcount)
    conn.commit()
    return MetabolicDomainStats(classified=classified, by_default=by_default)
