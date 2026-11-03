"""Stored entity labels — gene symbol + universal fallback.

Runs during ``derive`` against the canonical graph. Gene-based entities
(``Gene:MI:0250``) are labelled by their **primary gene symbol** — the
``Gene Name Primary`` identifier the resolver attached (``GENE_NAME_PRIMARY``,
distinct from synonyms). A gene carrying several primary-symbol identifiers
(across sources) takes the most-attested one (then shortest, then alphabetical)
for determinism. Any entity still without a label falls back to its canonical
identifier so every entity has a non-empty ``label``; the chemical and lipid
cascades overwrite their types' labels when they land.
"""

from __future__ import annotations

from dataclasses import dataclass

from psycopg2 import sql
import psycopg2.extensions


# CV label:accession strings (hardcoded per the classify/* convention — these
# are stable vocabulary names, also produced by cv_term_label_accession()).
GENE_ENTITY_TYPE = 'Gene:MI:0250'
GENE_NAME_PRIMARY_TYPE = 'Gene Name Primary:OM:0200'

GENE_SYMBOL_RULE = 'gene_symbol'
IDENTIFIER_FALLBACK_RULE = 'identifier_fallback'


@dataclass(frozen=True)
class EntityLabelStats:
    gene_symbol: int = 0
    identifier_fallback: int = 0
    without_label: int = 0


def _scalar(cur: psycopg2.extensions.cursor, query: sql.SQL, params: list) -> int:
    cur.execute(query, params)
    row = cur.fetchone()
    return int(row[0]) if row and row[0] is not None else 0


def populate_entity_labels(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
) -> EntityLabelStats:
    """Populate ``entity.label`` / ``entity.label_rule`` (gene symbol + fallback)."""

    schema_id = sql.Identifier(schema)
    with conn.cursor() as cur:
        # Resolve the vocabulary ids; skip gracefully if the build predates them.
        gene_type_id = _scalar(
            cur,
            sql.SQL(
                'SELECT entity_type_id FROM {}.vocab_entity_type WHERE name = %s'
            ).format(schema_id),
            [GENE_ENTITY_TYPE],
        )
        gene_name_primary_id = _scalar(
            cur,
            sql.SQL(
                'SELECT identifier_type_id FROM {}.vocab_identifier_type '
                'WHERE name = %s'
            ).format(schema_id),
            [GENE_NAME_PRIMARY_TYPE],
        )

        gene_symbol = 0
        if gene_type_id and gene_name_primary_id:
            # Gene symbol = the most-attested primary-symbol identifier per gene
            # (tie-break shortest, then alphabetical) for a deterministic label.
            cur.execute(
                sql.SQL(
                    """
                    WITH gene_symbol AS (
                      SELECT
                        ei.entity_id,
                        ie.value AS symbol,
                        count(*) AS attestations
                      FROM {schema}.entity_identifier ei
                      JOIN {schema}.entity e
                        ON e.entity_id = ei.entity_id
                      JOIN {schema}.identifier_evidence ie
                        ON ie.identifier_id = ei.identifier_id
                      WHERE e.entity_type_id = %(gene_type)s
                        AND ie.identifier_type_id = %(symbol_type)s
                        AND ie.value IS NOT NULL
                        AND ie.value <> ''
                      GROUP BY ei.entity_id, ie.value
                    ),
                    ranked AS (
                      SELECT
                        entity_id,
                        symbol,
                        row_number() OVER (
                          PARTITION BY entity_id
                          ORDER BY attestations DESC, length(symbol), symbol
                        ) AS rk
                      FROM gene_symbol
                    )
                    UPDATE {schema}.entity e
                    SET label = ranked.symbol,
                        label_rule = %(rule)s
                    FROM ranked
                    WHERE ranked.entity_id = e.entity_id
                      AND ranked.rk = 1
                    """
                ).format(schema=schema_id),
                {
                    'gene_type': gene_type_id,
                    'symbol_type': gene_name_primary_id,
                    'rule': GENE_SYMBOL_RULE,
                },
            )
            gene_symbol = cur.rowcount

        # Universal fallback: any entity still without a label takes its
        # canonical identifier, so every entity has a non-empty label.
        cur.execute(
            sql.SQL(
                """
                UPDATE {schema}.entity
                SET label = canonical_identifier,
                    label_rule = %(rule)s
                WHERE (label IS NULL OR label = '')
                  AND canonical_identifier IS NOT NULL
                  AND canonical_identifier <> ''
                """
            ).format(schema=schema_id),
            {'rule': IDENTIFIER_FALLBACK_RULE},
        )
        identifier_fallback = cur.rowcount

        without_label = _scalar(
            cur,
            sql.SQL(
                "SELECT count(*) FROM {}.entity "
                "WHERE label IS NULL OR label = ''"
            ).format(schema_id),
            [],
        )

    conn.commit()
    return EntityLabelStats(
        gene_symbol=gene_symbol,
        identifier_fallback=identifier_fallback,
        without_label=without_label,
    )
