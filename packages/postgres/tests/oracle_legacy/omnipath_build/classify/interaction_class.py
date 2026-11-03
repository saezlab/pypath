"""Interaction-class mapping for ``vocab_relation_predicate``.

Every relation predicate maps to exactly one coarse interaction class, per the
curated rules in ``interaction_class.yaml``. Runs during ``derive`` after load
(the predicate vocabulary is data-driven), populating
``vocab_relation_predicate.interaction_class_id`` with no reload.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from psycopg2 import sql
import psycopg2.extensions
import yaml


_RULES_PATH = Path(__file__).with_name('interaction_class.yaml')


@dataclass(frozen=True)
class InteractionClassStats:
    mapped: int = 0
    by_default: int = 0
    default_predicates: tuple[str, ...] = field(default_factory=tuple)


def _load_rules() -> dict:
    with _RULES_PATH.open('r', encoding='utf-8') as handle:
        return yaml.safe_load(handle)


def classify_interaction_class(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
) -> InteractionClassStats:
    """Seed ``vocab_interaction_class`` and populate predicate FK column."""

    rules = _load_rules()
    schema_id = sql.Identifier(schema)
    with conn.cursor() as cur:
        # Seed the controlled vocabulary (id = declaration order). `name` is the
        # snake_case slug everything filters by, `label` its capitalised display
        # form, `controlled_vocabulary_mapping` a nullable ontology CURIE.
        class_rows = [
            (
                index,
                entry['name'],
                entry['label'],
                entry.get('controlled_vocabulary_mapping'),
            )
            for index, entry in enumerate(rules['classes'], start=1)
        ]
        cur.executemany(
            sql.SQL(
                """
                INSERT INTO {}.vocab_interaction_class
                  (interaction_class_id, name, label,
                   controlled_vocabulary_mapping)
                VALUES (%s, %s, %s, %s)
                ON CONFLICT (interaction_class_id) DO UPDATE
                SET name = EXCLUDED.name,
                    label = EXCLUDED.label,
                    controlled_vocabulary_mapping =
                      EXCLUDED.controlled_vocabulary_mapping
                """
            ).format(schema_id).as_string(cur.connection),
            class_rows,
        )

        # Materialise the predicate -> class map from the rules.
        cur.execute(
            'CREATE TEMP TABLE _predicate_interaction_class '
            '(predicate_name text PRIMARY KEY, class_name text NOT NULL) '
            'ON COMMIT DROP'
        )
        predicate_rows = [
            (predicate, class_name)
            for predicate, class_name in (rules.get('predicates') or {}).items()
        ]
        if predicate_rows:
            cur.executemany(
                'INSERT INTO _predicate_interaction_class '
                '(predicate_name, class_name) VALUES (%s, %s) '
                'ON CONFLICT (predicate_name) DO NOTHING',
                predicate_rows,
            )

        # Reset, then assign the mapped class to each predicate.
        cur.execute(
            sql.SQL(
                'UPDATE {}.vocab_relation_predicate '
                'SET interaction_class_id = NULL'
            ).format(schema_id)
        )
        cur.execute(
            sql.SQL(
                """
                UPDATE {schema}.vocab_relation_predicate p
                SET interaction_class_id = vic.interaction_class_id
                FROM _predicate_interaction_class m
                JOIN {schema}.vocab_interaction_class vic
                  ON vic.name = m.class_name
                WHERE m.predicate_name = p.name
                """
            ).format(schema=schema_id)
        )
        mapped = int(cur.rowcount)

        # Fallback: every uncovered predicate -> default class (no NULLs).
        cur.execute(
            sql.SQL(
                'SELECT interaction_class_id FROM {}.vocab_interaction_class '
                'WHERE name = %s'
            ).format(schema_id),
            [rules['default']],
        )
        default_id = cur.fetchone()[0]
        cur.execute(
            sql.SQL(
                'UPDATE {}.vocab_relation_predicate SET interaction_class_id = %s '
                'WHERE interaction_class_id IS NULL'
            ).format(schema_id),
            [default_id],
        )
        by_default = int(cur.rowcount)

        # The predicates that fell to default (surface new/uncurated ones).
        cur.execute(
            sql.SQL(
                'SELECT name FROM {}.vocab_relation_predicate '
                'WHERE interaction_class_id = %s ORDER BY name'
            ).format(schema_id),
            [default_id],
        )
        default_predicates = tuple(row[0] for row in cur.fetchall())
    conn.commit()
    return InteractionClassStats(
        mapped=mapped,
        by_default=by_default,
        default_predicates=default_predicates,
    )
