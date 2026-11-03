"""Populate ``entity_name`` (spec 011 data-model.md section 4).

Runs during ``derive``, alongside :func:`chemical_labels.populate_chemical_labels`
-- entity ids only stabilize after resolution, so this reads the same
resolved ``entity_identifier`` data that module's own label cascade reads,
rather than attaching names during per-source ingestion. One row per
(source, entity, normalized name): the explicit ``is_preferred`` flag
replaces the implicit split between the ``Name`` and ``Synonym`` identifier
types (FR-016) other consumers -- the label cascade, the translation
service's name index -- read instead of re-deriving the split themselves.

``name_normalized`` folds case and whitespace the same way
``omnipath_utils.mapping._id_types.normalize_name`` does on the translation
side (spec 011 T071) -- kept in sync by hand, the same cross-repo mirroring
precedent as ``resolver/chemical_normalization.py`` (T040).
"""

from __future__ import annotations

from dataclasses import dataclass

from psycopg2 import sql
import psycopg2.extensions

from omnipath_build.labels.chemical_labels import (
    CHEMICAL_ENTITY_TYPE,
    INCHIKEY_RE,
    _scalar,
    _values_clause,
)

#: Name-bearing identifier types (the same set chemical_labels.py's own
#: cascade reads), classified by preference and kind. A name is preferred
#: when it is the minting resource's own recommended name (FR-016) -- the
#: ``Name`` type specifically, not a systematic IUPAC form or a synonym.
_NAME_KIND: tuple[tuple[str, bool, str], ...] = (
    ('Name:OM:0202', True, 'trivial'),
    ('Inn:OM:0120', False, 'trivial'),
    ('Abbreviated Name:OM:0208', False, 'abbreviation'),
    ('Iupac Traditional Name:OM:0211', False, 'systematic'),
    ('Synonym:OM:0203', False, 'trivial'),
    ('Iupac Name:OM:0210', False, 'systematic'),
)


@dataclass(frozen=True)
class EntityNameStats:
    rows_written: int = 0


def populate_entity_name(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
) -> EntityNameStats:
    """Write one ``entity_name`` row per (source, entity, normalized name)."""

    schema_id = sql.Identifier(schema)

    with conn.cursor() as cur:
        chem_type_id = _scalar(
            cur,
            sql.SQL(
                'SELECT entity_type_id FROM {}.vocab_entity_type WHERE name = %s'
            ).format(schema_id),
            [CHEMICAL_ENTITY_TYPE],
        )
        if not chem_type_id:
            return EntityNameStats()

        cur.execute(
            sql.SQL(
                """
                WITH name_type(type_name, is_preferred, name_kind)
                  AS (VALUES {name_values}),
                candidate AS (
                  -- One source name may pack several names with a ``|``
                  -- delimiter and a trailing ``(PTF#####)`` source-id suffix
                  -- (FooDB/PTFI); split and strip so each clean name is its
                  -- own candidate, matching chemical_labels.py's own pass.
                  SELECT
                    ei.entity_id,
                    ei.source_id,
                    nt.is_preferred,
                    nt.name_kind,
                    lower(c.cleaned) AS norm,
                    c.cleaned        AS val
                  FROM {schema}.entity e
                  JOIN {schema}.entity_identifier ei
                    ON ei.entity_id = e.entity_id
                  JOIN {schema}.identifier_evidence ie
                    ON ie.identifier_id = ei.identifier_id
                  JOIN {schema}.vocab_identifier_type it
                    ON it.identifier_type_id = ie.identifier_type_id
                  JOIN name_type nt ON nt.type_name = it.name
                  CROSS JOIN LATERAL (
                    SELECT btrim(
                      regexp_replace(part, '\\s*\\(PTF[0-9]+\\)\\s*$', '')
                    ) AS cleaned
                    FROM unnest(string_to_array(ie.value, '|')) AS part
                  ) c
                  WHERE e.entity_type_id = %(chem)s
                    AND c.cleaned <> ''
                    AND length(c.cleaned) BETWEEN 2 AND 120
                    AND c.cleaned !~ '^[0-9]+$'
                    AND c.cleaned !~ %(inchikey)s
                    AND c.cleaned IS DISTINCT FROM e.canonical_identifier
                ),
                -- The same normalized name can reach one entity/source
                -- through more than one identifier type (e.g. a ChEBI
                -- primary name that another CV term also lists as a
                -- synonym); keep the most preferred classification.
                ranked AS (
                  SELECT DISTINCT ON (entity_id, source_id, norm)
                    entity_id, source_id, val, norm, is_preferred, name_kind
                  FROM candidate
                  ORDER BY entity_id, source_id, norm,
                           is_preferred DESC, length(val), val
                )
                INSERT INTO {schema}.entity_name
                  (source_id, entity_id, name, name_normalized,
                   is_preferred, name_kind)
                SELECT source_id, entity_id, val, norm, is_preferred, name_kind
                FROM ranked
                ON CONFLICT (source_id, entity_id, name_normalized) DO NOTHING
                """
            ).format(
                schema=schema_id,
                name_values=_values_clause(_NAME_KIND),
            ),
            dict(chem=chem_type_id, inchikey=INCHIKEY_RE),
        )
        rows_written = cur.rowcount

    conn.commit()
    return EntityNameStats(rows_written=rows_written)
