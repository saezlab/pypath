"""Chemical-class classification for ``Chemical:OM:0037`` entities.

Flat, single-valued: a chemical's class is the most specific (lowest
``precedence``) class implied by the resources that contributed it, per the
curated rules in ``chemical_class.yaml``. Runs during ``derive`` against the
canonical graph + ``entity_evidence_resolution`` (every contributing source,
regardless of resolution status), so a chemical is classed by what actually
produced it rather than only its resolved evidence (no reload).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from psycopg2 import sql
import psycopg2.extensions
import yaml


CHEMICAL_ENTITY_TYPE = 'Chemical:OM:0037'
_RULES_PATH = Path(__file__).with_name('chemical_class.yaml')


@dataclass(frozen=True)
class ChemicalClassStats:
    classified: int = 0
    by_default: int = 0


def _load_rules() -> dict:
    with _RULES_PATH.open('r', encoding='utf-8') as handle:
        return yaml.safe_load(handle)


def classify_chemical_class(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
) -> ChemicalClassStats:
    """Seed ``vocab_chemical_class`` and populate ``entity.chemical_class_id``."""

    rules = _load_rules()
    schema_id = sql.Identifier(schema)
    with conn.cursor() as cur:
        # Resolve the chemical entity type; if absent, nothing to do.
        cur.execute(
            sql.SQL(
                'SELECT entity_type_id FROM {}.vocab_entity_type WHERE name = %s'
            ).format(schema_id),
            [CHEMICAL_ENTITY_TYPE],
        )
        row = cur.fetchone()
        if row is None:
            return ChemicalClassStats()
        chemical_type_id = row[0]

        # Seed the controlled vocabulary from the rule file.
        class_rows = [
            (entry['precedence'], entry['name'], entry['precedence'])
            for entry in rules['classes']
        ]
        cur.executemany(
            sql.SQL(
                """
                INSERT INTO {}.vocab_chemical_class
                  (chemical_class_id, name, precedence)
                VALUES (%s, %s, %s)
                ON CONFLICT (chemical_class_id) DO UPDATE
                SET name = EXCLUDED.name, precedence = EXCLUDED.precedence
                """
            ).format(schema_id).as_string(cur.connection),
            class_rows,
        )

        # Materialise the source -> class map from the rules.
        cur.execute(
            'CREATE TEMP TABLE _source_chemical_class '
            '(source_name text PRIMARY KEY, class_name text NOT NULL) '
            'ON COMMIT DROP'
        )
        source_rows = [
            (source_name, class_name)
            for class_name, sources in (rules.get('sources') or {}).items()
            for source_name in sources
        ]
        if source_rows:
            cur.executemany(
                'INSERT INTO _source_chemical_class (source_name, class_name) '
                'VALUES (%s, %s) ON CONFLICT (source_name) DO NOTHING',
                source_rows,
            )

        # Reset, then assign the most-specific class implied by an entity's
        # contributing sources.
        cur.execute(
            sql.SQL(
                'UPDATE {}.entity SET chemical_class_id = NULL '
                'WHERE entity_type_id = %s'
            ).format(schema_id),
            [chemical_type_id],
        )
        cur.execute(
            sql.SQL(
                """
                WITH ranked AS (
                  SELECT
                    er.entity_id,
                    vcc.chemical_class_id,
                    row_number() OVER (
                      PARTITION BY er.entity_id ORDER BY vcc.precedence ASC
                    ) AS rn
                  FROM {schema}.entity_evidence_resolution er
                  JOIN {schema}.entity ce
                    ON ce.entity_id = er.entity_id
                   AND ce.entity_type_id = %(chem)s
                  JOIN {schema}.data_source ds
                    ON ds.source_id = er.source_id
                  JOIN _source_chemical_class m
                    ON m.source_name = ds.name
                  JOIN {schema}.vocab_chemical_class vcc
                    ON vcc.name = m.class_name
                  WHERE er.entity_id IS NOT NULL
                )
                UPDATE {schema}.entity e
                SET chemical_class_id = ranked.chemical_class_id
                FROM ranked
                WHERE ranked.entity_id = e.entity_id
                  AND ranked.rn = 1
                """
            ).format(schema=schema_id),
            {'chem': chemical_type_id},
        )
        classified = int(cur.rowcount)

        # Fallback: any remaining chemical (sources present but unmapped, or no
        # source_count row) gets the default class.
        cur.execute(
            sql.SQL(
                'SELECT chemical_class_id FROM {}.vocab_chemical_class '
                'WHERE name = %s'
            ).format(schema_id),
            [rules['default']],
        )
        default_id = cur.fetchone()[0]
        cur.execute(
            sql.SQL(
                'UPDATE {}.entity SET chemical_class_id = %s '
                'WHERE entity_type_id = %s AND chemical_class_id IS NULL'
            ).format(schema_id),
            [default_id, chemical_type_id],
        )
        by_default = int(cur.rowcount)
    conn.commit()
    return ChemicalClassStats(classified=classified, by_default=by_default)
