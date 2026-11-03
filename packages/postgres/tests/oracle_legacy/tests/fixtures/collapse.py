"""The collapse of the interaction record, for the tests that assert its shape.

As amended, the collapse of ``interaction_fact_resource`` for a resource scope
is what a **query** produces at request time, and it is not a table. The build
stops when the record lands, so ``interaction_fact_combined`` — the
all-resources scope materialised — is gone, and with it
``collapse_interaction_scope``, the build routine that filled it.

The assertions that table carried are not gone. The three-valued flags, the
provenance over every contributor, the ordered key and the class precedence are
all statements about **the collapse**, and they hold or fail whether the
collapse is stored or computed. So this module is the oracle those tests read:
the fold, written once, in SQL, on the test side.

**It is a test fixture and not a second implementation.** The production fold
has exactly one caller — the api-service — which is the whole of what removing
the materialisation bought. What lives here is the specification of that fold
in executable form: the same aggregate on the same column, so an engine that
folds differently fails a test in the repository that owns the record rather
than agreeing with itself in the repository that owns the query.

The scope is the resource set, exactly as it was on the removed routine.
``None`` means every resource — the all-resources scope, which is the shape the
removed table held — and an **empty** sequence is an empty scope, which
collapses to no rows. That is a different answer from "no restriction", and the
difference is the point of stating the scope at all.
"""

from __future__ import annotations

from collections.abc import Sequence

import psycopg2.extensions

#: The record the fold reads.
RECORD_TABLE = 'interaction_fact_resource'

#: The view :func:`create_collapse_view` leaves behind, holding the fold over
#: every resource. A test that used to read the materialised table reads this,
#: which keeps the assertion about the collapse and drops the assumption that
#: the collapse is stored.
COLLAPSE_VIEW = 'interaction_collapse_all_resources'

#: The fold's output columns, in the order the statement selects them.
COLLAPSE_COLUMNS = (
    'subject_entity_id',
    'object_entity_id',
    'interaction_class_id',
    'subject_organism',
    'object_organism',
    'is_directed',
    'is_stimulation',
    'is_inhibition',
    'sign_source_count',
    'direction_source_count',
    'affinity',
    'pchembl',
    'score',
    'sources',
    'source_count',
    'dataset_tags',
    'curation_flags',
    'reference_pubmed_ids',
    'reference_dois',
    'reference_count',
    'attributes',
    'interaction_id',
)

_COLLAPSE_SQL = """
WITH scoped_provenance AS (
  SELECT
    r.subject_entity_id,
    r.object_entity_id,
    r.interaction_class_id,
    array_agg(DISTINCT contribution.value)
      FILTER (WHERE contribution.kind = 'pubmed') AS reference_pubmed_ids,
    array_agg(DISTINCT contribution.value)
      FILTER (WHERE contribution.kind = 'doi') AS reference_dois,
    array_agg(DISTINCT contribution.value)
      FILTER (WHERE contribution.kind = 'curation') AS curation_flags
  FROM {record} r
  CROSS JOIN LATERAL (
    SELECT 'pubmed'::text AS kind, pubmed.value
      FROM unnest(r.reference_pubmed_ids) AS pubmed(value)
    UNION ALL
    SELECT 'doi'::text, doi.value
      FROM unnest(r.reference_dois) AS doi(value)
    UNION ALL
    SELECT 'curation'::text, flag.value
      FROM unnest(r.curation_flags) AS flag(value)
  ) AS contribution
  WHERE {scope}
  GROUP BY 1, 2, 3
)
SELECT
  r.subject_entity_id,
  r.object_entity_id,
  r.interaction_class_id,
  min(r.subject_organism) AS subject_organism,
  min(r.object_organism) AS object_organism,
  bool_or(r.is_directed) AS is_directed,
  bool_or(r.is_stimulation) AS is_stimulation,
  bool_or(r.is_inhibition) AS is_inhibition,
  (count(DISTINCT r.source_id) FILTER (
     WHERE r.is_stimulation IS NOT NULL OR r.is_inhibition IS NOT NULL
   ))::smallint AS sign_source_count,
  (count(DISTINCT r.source_id) FILTER (
     WHERE r.is_directed IS NOT NULL
   ))::smallint AS direction_source_count,
  min(r.affinity) AS affinity,
  max(r.pchembl) AS pchembl,
  max(r.score) AS score,
  array_agg(DISTINCT contributor.name) AS sources,
  count(DISTINCT r.source_id)::int AS source_count,
  NULL::text[] AS dataset_tags,
  provenance.curation_flags,
  provenance.reference_pubmed_ids,
  provenance.reference_dois,
  (coalesce(cardinality(provenance.reference_pubmed_ids), 0)
   + coalesce(cardinality(provenance.reference_dois), 0))::int
    AS reference_count,
  NULL::jsonb AS attributes,
  -- Postgres has no min(uuid); the text form orders the same way for the
  -- canonical rendering, and every record row of a group carries the same
  -- header anyway, so any one of them is the answer.
  min(r.interaction_id::text)::uuid AS interaction_id
FROM {record} r
JOIN {data_source} contributor ON contributor.source_id = r.source_id
LEFT JOIN scoped_provenance provenance
  ON provenance.subject_entity_id = r.subject_entity_id
 AND provenance.object_entity_id = r.object_entity_id
 AND provenance.interaction_class_id = r.interaction_class_id
WHERE {scope}
GROUP BY
  r.subject_entity_id,
  r.object_entity_id,
  r.interaction_class_id,
  provenance.reference_pubmed_ids,
  provenance.reference_dois,
  provenance.curation_flags
"""


def _quoted(identifier: str) -> str:
    """A double-quoted SQL identifier, for a statement built as text."""
    return '"{}"'.format(identifier.replace('"', '""'))


def collapse_sql(
    *,
    schema: str,
    sources: Sequence[str] | None = None,
) -> tuple[str, tuple[object, ...]]:
    """``(statement, parameters)`` for the collapse over ``sources``.

    The statement is a bare ``SELECT`` producing :data:`COLLAPSE_COLUMNS` in
    order, so a caller embeds it as a subquery. ``sources`` names resources by
    ``data_source.name``; ``None`` is the all-resources scope.

    The scope predicate appears twice — once in the provenance CTE and once in
    the outer aggregate — so the parameter is passed twice. A fold that applied
    it to only one of them would report the scope's sources beside every
    contributor's references, which is the scope-rule defect with the halves
    the other way round.
    """
    schema_sql = _quoted(schema)
    record = f'{schema_sql}.{_quoted(RECORD_TABLE)}'
    data_source = f'{schema_sql}.{_quoted("data_source")}'
    if sources is None:
        scope = 'true'
        parameters: tuple[object, ...] = ()
    else:
        scope = (
            'r.source_id IN (SELECT source_id FROM '
            f'{data_source} WHERE name = ANY(%s::text[]))'
        )
        parameters = (list(sources), list(sources))
    statement = _COLLAPSE_SQL.format(
        record=record,
        data_source=data_source,
        scope=scope,
    )
    return statement, parameters


def create_collapse_view(
    conn: psycopg2.extensions.connection,
    schema: str,
) -> str:
    """Create :data:`COLLAPSE_VIEW` over the record, and return its name.

    A view rather than a table on purpose: it holds no rows and costs no build
    step, so a test reading it is reading the fold as a query computes it. The
    tests that used to name ``interaction_fact_combined`` name this instead,
    and what they assert is unchanged — which is the evidence that those
    assertions were about the collapse and never about its storage.
    """
    statement, parameters = collapse_sql(schema=schema)
    assert not parameters, 'the all-resources scope takes no parameters'
    with conn.cursor() as cur:
        cur.execute(
            f'CREATE OR REPLACE VIEW {_quoted(schema)}.{_quoted(COLLAPSE_VIEW)}'
            f' AS {statement}'
        )
    conn.commit()
    return COLLAPSE_VIEW
