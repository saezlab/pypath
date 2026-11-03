"""Chemical resolution-level grouping.

Resolved chemical entities are grouped into **selectable structural-specificity
levels** derived purely from the InChIKey block structure — *without RDKit*
(Constitution II keeps RDKit out of the core build). A relation expressed on
different sub-level structures by different resources then projects to **one**
edge at a coarser level, carrying per-resource provenance.

InChIKey layout (standard key, 27 chars)::

    AAAAAAAAAAAAAA - BBBBBBBBFV - P
    `------ 14 ----'   `-- 10 --'   `1
       block 1            block 2    final
    skeleton /         stereo /      protonation /
    connectivity       isotope /     charge
                       tautomer (+ standard flag)

Levels (coarse -> fine):

* ``connectivity``             -> block 1 (first 14 chars). Collapses L/D/racemic
  and other stereo/charge variants of one carbon skeleton.
* ``stereo_isotope_tautomer``  -> blocks 1+2 (first 25 chars, through the 2nd
  block, i.e. *up to* the 2nd dash). Distinguishes stereoisomers; still ignores
  protonation/charge.
* ``full``                     -> the full InChIKey. Effectively one group per
  distinct structure (InChIKey-canonicalised chemicals already merge here).

The three alanine stereo/charge variants (L/D/DL) share block 1
(``QNAYBMKLOCPYGJ``) so they collapse at ``connectivity`` only — at
``stereo_isotope_tautomer`` their block-2 hashes differ. beta-alanine (a
positional isomer) and N-acetyl-L-alanine (the peptide-bond / residue-context
form) have a *different* block 1 and never collapse with alanine at any level.

The materialised tables (built post-load over the full Postgres graph, like the
other derived summaries) are:

* ``chemical_resolution_level``        — the level seed (this module is the
  single source of truth; :func:`ensure_chemical_resolution_schema`).
* ``chemical_resolution_group`` / ``chemical_resolution_group_member`` —
  per level, the InChIKey-prefix group and its member chemical entities.
* ``chemical_resolution_relation`` — per level, each relation touching a chemical
  re-pointed to its group representative and collapsed, with union provenance
  (``source_ids``). No request-time GROUP BY: the API serves a level straight
  from this table.
"""

from __future__ import annotations

from dataclasses import dataclass
import re
import time

from psycopg2 import sql
import psycopg2.extensions

from omnipath_build.cv_terms import CHEMICAL_ENTITY_TYPE

#: Identifier type name of a standard InChIKey (matches ``chemical_fallback``).
STANDARD_INCHI_KEY_TYPE = 'Standard Inchi Key:MI:1101'

#: Trivial-name identifier types (a name borne by >1 structure is ambiguous).
CHEMICAL_NAME_TYPES = ('Name:OM:0202', 'Synonym:OM:0203')

#: A well-formed standard InChIKey (14-10-1, all upper-case A-Z).
INCHIKEY_REGEX = r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$'
_INCHIKEY_RE = re.compile(INCHIKEY_REGEX)


@dataclass(frozen=True)
class ChemicalResolutionLevel:
    """One structural-specificity level keyed by an InChIKey prefix length."""

    level_id: int
    name: str
    inchikey_prefix_length: int
    specificity_rank: int  # 1 = coarsest .. higher = finer
    description: str


#: Single source of truth for the level seed (coarse -> fine).
LEVELS: tuple[ChemicalResolutionLevel, ...] = (
    ChemicalResolutionLevel(
        1,
        'connectivity',
        14,
        1,
        'InChIKey block 1 (first 14 chars) — skeleton/connectivity only; '
        'ignores stereochemistry, isotopes, tautomers and protonation. '
        'Collapses L/D/racemic and charge variants of one skeleton.',
    ),
    ChemicalResolutionLevel(
        2,
        'stereo_isotope_tautomer',
        25,
        2,
        'InChIKey blocks 1+2 (first 25 chars, through the 2nd block) — adds the '
        'stereo/isotope/tautomer layer; still ignores protonation/charge.',
    ),
    ChemicalResolutionLevel(
        3,
        'full',
        27,
        3,
        'Full InChIKey — adds the protonation/charge layer; the most specific '
        'level (one group per distinct structure).',
    ),
)

LEVELS_BY_NAME: dict[str, ChemicalResolutionLevel] = {
    level.name: level for level in LEVELS
}


def resolution_level_key(
    inchikey: str,
    level: str | ChemicalResolutionLevel,
) -> str | None:
    """The InChIKey prefix that is the group key for ``level``.

    Returns ``None`` when ``inchikey`` is not a well-formed standard InChIKey
    (so a junk value never seeds a spurious group). This is the pure-Python
    mirror of the ``left(inchikey, prefix_length)`` used in the SQL build.
    """

    if not inchikey or not _INCHIKEY_RE.match(inchikey):
        return None
    resolved = (
        level
        if isinstance(level, ChemicalResolutionLevel)
        else LEVELS_BY_NAME[level]
    )
    return inchikey[: resolved.inchikey_prefix_length]


# ---------------------------------------------------------------------------
# Schema (DDL + seed) — called from db.schema._ensure_resolution_schema
# ---------------------------------------------------------------------------


def ensure_chemical_resolution_schema(
    cur: psycopg2.extensions.cursor,
    schema: str = 'public',
) -> None:
    """Create the resolution-level tables and (re)seed ``*_level``."""

    schema_id = sql.Identifier(schema)
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {0}.chemical_resolution_level (
              level_id smallint PRIMARY KEY,
              name text NOT NULL UNIQUE,
              inchikey_prefix_length smallint NOT NULL,
              specificity_rank smallint NOT NULL,
              description text NOT NULL
            )
            """
        ).format(schema_id)
    )
    cur.executemany(
        sql.SQL(
            """
            INSERT INTO {0}.chemical_resolution_level
              (level_id, name, inchikey_prefix_length, specificity_rank,
               description)
            VALUES (%s, %s, %s, %s, %s)
            ON CONFLICT (level_id) DO UPDATE SET
              name = EXCLUDED.name,
              inchikey_prefix_length = EXCLUDED.inchikey_prefix_length,
              specificity_rank = EXCLUDED.specificity_rank,
              description = EXCLUDED.description
            """
        ).format(schema_id),
        [
            (
                level.level_id,
                level.name,
                level.inchikey_prefix_length,
                level.specificity_rank,
                level.description,
            )
            for level in LEVELS
        ],
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {0}.chemical_resolution_group (
              level_id smallint NOT NULL
                REFERENCES {0}.chemical_resolution_level(level_id),
              group_key text NOT NULL,
              representative_entity_id uuid NOT NULL,
              representative_inchikey text NOT NULL,
              member_count integer NOT NULL,
              PRIMARY KEY (level_id, group_key)
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {0}.chemical_resolution_group_member (
              level_id smallint NOT NULL,
              group_key text NOT NULL,
              entity_id uuid NOT NULL,
              inchikey text NOT NULL,
              PRIMARY KEY (level_id, entity_id)
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {0}.chemical_ambiguous_name_candidate (
              ambiguous_name text NOT NULL,
              candidate_entity_id uuid NOT NULL,
              candidate_inchikey text NOT NULL,
              resolution_mechanism text NOT NULL
                DEFAULT 'ambiguous_name_class',
              PRIMARY KEY (ambiguous_name, candidate_entity_id)
            )
            """
        ).format(schema_id)
    )
    cur.execute(
        sql.SQL(
            """
            CREATE TABLE IF NOT EXISTS {0}.chemical_resolution_relation (
              level_id smallint NOT NULL,
              subject_entity_id uuid NOT NULL,
              predicate_id bigint NOT NULL,
              object_entity_id uuid NOT NULL,
              relation_category_id bigint,
              member_relation_count integer NOT NULL,
              source_ids bigint[] NOT NULL,
              source_count integer NOT NULL,
              PRIMARY KEY
                (level_id, subject_entity_id, predicate_id, object_entity_id)
            )
            """
        ).format(schema_id)
    )


# ---------------------------------------------------------------------------
# Build (post-load, over the full Postgres graph) — called from cli derive
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ChemicalResolutionStats:
    """Summary counts from the resolution-level build (cost notes / manifest)."""

    chemical_entities_with_inchikey: int = 0
    group_members: int = 0
    groups: int = 0
    relations: int = 0


def _log(progress: bool, step: str, event: str, **fields: object) -> None:
    if not progress:
        return
    details = ' '.join(f'{key}={value}' for key, value in fields.items())
    print(
        f'[chem-resolution-level] step={step} event={event}'
        + (f' {details}' if details else ''),
        flush=True,
    )


def rebuild_chemical_resolution_levels(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
    progress: bool = False,
) -> ChemicalResolutionStats:
    """Materialise the per-level chemical group + collapsed-relation tables.

    Runs AFTER the canonical graph and ``entity_identifier_lookup`` exist
    (i.e. after :func:`rebuild_derived_tables`), reflecting the *full* graph.
    Idempotent: truncates and rebuilds from ``entity`` / ``relation``.
    """

    schema_id = sql.Identifier(schema)
    chem_type = sql.Literal(CHEMICAL_ENTITY_TYPE)
    ik_type = sql.Literal(STANDARD_INCHI_KEY_TYPE)
    ik_regex = sql.Literal(INCHIKEY_REGEX)
    started = time.perf_counter()

    with conn.cursor() as cur:
        ensure_chemical_resolution_schema(cur, schema)
        for table in (
            'chemical_resolution_relation',
            'chemical_resolution_group',
            'chemical_resolution_group_member',
        ):
            cur.execute(
                sql.SQL('TRUNCATE {0}.{1}').format(
                    schema_id, sql.Identifier(table)
                )
            )

        # 1) one primary InChIKey per chemical entity (canonical preferred,
        #    else the lexicographically-min attached InChIKey).
        _log(progress, 'primary_inchikey', 'start')
        cur.execute(
            sql.SQL(
                """
                CREATE TEMP TABLE _chem_entity_inchikey ON COMMIT DROP AS
                WITH chem AS (
                  SELECT
                    e.entity_id,
                    e.canonical_identifier_type_id,
                    e.canonical_identifier
                  FROM {0}.entity e
                  JOIN {0}.vocab_entity_type vet
                    ON vet.entity_type_id = e.entity_type_id
                  WHERE vet.name = {1}
                ),
                ik_type AS (
                  SELECT identifier_type_id
                  FROM {0}.vocab_identifier_type
                  WHERE name = {2}
                ),
                ik_per_entity AS (
                  SELECT c.entity_id, c.canonical_identifier AS inchikey, 0 AS pref
                  FROM chem c
                  WHERE c.canonical_identifier_type_id
                          = (SELECT identifier_type_id FROM ik_type)
                    AND c.canonical_identifier ~ {3}
                  UNION ALL
                  SELECT c.entity_id, ie.value AS inchikey, 1 AS pref
                  FROM chem c
                  JOIN {0}.entity_identifier_lookup eil
                    ON eil.entity_id = c.entity_id
                  JOIN {0}.identifier_evidence ie
                    ON ie.identifier_id = eil.identifier_id
                  WHERE ie.identifier_type_id
                          = (SELECT identifier_type_id FROM ik_type)
                    AND ie.value ~ {3}
                )
                SELECT DISTINCT ON (entity_id) entity_id, inchikey
                FROM ik_per_entity
                ORDER BY entity_id, pref, inchikey
                """
            ).format(schema_id, chem_type, ik_type, ik_regex)
        )
        cur.execute('SELECT count(*) FROM _chem_entity_inchikey')
        n_chem = int(cur.fetchone()[0])

        # 2) members (one row per entity per level) + groups.
        _log(progress, 'groups', 'start')
        cur.execute(
            sql.SQL(
                """
                INSERT INTO {0}.chemical_resolution_group_member
                  (level_id, group_key, entity_id, inchikey)
                SELECT
                  lvl.level_id,
                  left(p.inchikey, lvl.inchikey_prefix_length),
                  p.entity_id,
                  p.inchikey
                FROM _chem_entity_inchikey p
                CROSS JOIN {0}.chemical_resolution_level lvl
                """
            ).format(schema_id)
        )
        n_members = int(cur.rowcount)
        cur.execute(
            sql.SQL(
                """
                INSERT INTO {0}.chemical_resolution_group
                  (level_id, group_key, representative_entity_id,
                   representative_inchikey, member_count)
                WITH rep AS (
                  SELECT DISTINCT ON (level_id, group_key)
                    level_id, group_key,
                    entity_id AS representative_entity_id,
                    inchikey AS representative_inchikey
                  FROM {0}.chemical_resolution_group_member
                  ORDER BY level_id, group_key, entity_id
                ),
                cnt AS (
                  SELECT level_id, group_key, count(*) AS member_count
                  FROM {0}.chemical_resolution_group_member
                  GROUP BY level_id, group_key
                )
                SELECT
                  rep.level_id, rep.group_key,
                  rep.representative_entity_id, rep.representative_inchikey,
                  cnt.member_count
                FROM rep
                JOIN cnt
                  ON cnt.level_id = rep.level_id
                 AND cnt.group_key = rep.group_key
                """
            ).format(schema_id)
        )
        n_groups = int(cur.rowcount)

        # 3) collapsed relations per level, with union provenance.
        #    Each chemical endpoint is re-pointed to its group representative;
        #    non-chemical endpoints (genes) stay as themselves. Only relations
        #    touching at least one grouped chemical are materialised.
        _log(progress, 'relations', 'start')
        cur.execute(
            sql.SQL(
                """
                INSERT INTO {0}.chemical_resolution_relation
                  (level_id, subject_entity_id, predicate_id, object_entity_id,
                   relation_category_id, member_relation_count,
                   source_ids, source_count)
                WITH node_map AS (
                  SELECT
                    m.level_id, m.entity_id,
                    g.representative_entity_id AS projected_entity_id
                  FROM {0}.chemical_resolution_group_member m
                  JOIN {0}.chemical_resolution_group g
                    ON g.level_id = m.level_id
                   AND g.group_key = m.group_key
                ),
                projected AS (
                  SELECT
                    lvl.level_id,
                    r.relation_id,
                    coalesce(sm.projected_entity_id, r.subject_entity_id)
                      AS subject_entity_id,
                    r.predicate_id,
                    coalesce(om.projected_entity_id, r.object_entity_id)
                      AS object_entity_id,
                    r.relation_category_id
                  FROM {0}.relation r
                  CROSS JOIN {0}.chemical_resolution_level lvl
                  LEFT JOIN node_map sm
                    ON sm.level_id = lvl.level_id
                   AND sm.entity_id = r.subject_entity_id
                  LEFT JOIN node_map om
                    ON om.level_id = lvl.level_id
                   AND om.entity_id = r.object_entity_id
                  WHERE sm.projected_entity_id IS NOT NULL
                     OR om.projected_entity_id IS NOT NULL
                )
                SELECT
                  pr.level_id,
                  pr.subject_entity_id,
                  pr.predicate_id,
                  pr.object_entity_id,
                  min(pr.relation_category_id) AS relation_category_id,
                  count(DISTINCT pr.relation_id) AS member_relation_count,
                  coalesce(
                    array_agg(DISTINCT rer.source_id)
                      FILTER (WHERE rer.source_id IS NOT NULL),
                    '{{}}'::bigint[]
                  ) AS source_ids,
                  count(DISTINCT rer.source_id) AS source_count
                FROM projected pr
                LEFT JOIN {0}.relation_evidence_relation rer
                  ON rer.relation_id = pr.relation_id
                GROUP BY
                  pr.level_id, pr.subject_entity_id,
                  pr.predicate_id, pr.object_entity_id
                """
            ).format(schema_id)
        )
        n_relations = int(cur.rowcount)

        _create_chemical_resolution_indexes(cur, schema)

    conn.commit()
    _log(
        progress,
        'all',
        'done',
        chemical_entities=n_chem,
        members=n_members,
        groups=n_groups,
        relations=n_relations,
        seconds=f'{time.perf_counter() - started:.3f}',
    )
    return ChemicalResolutionStats(
        chemical_entities_with_inchikey=n_chem,
        group_members=n_members,
        groups=n_groups,
        relations=n_relations,
    )


@dataclass(frozen=True)
class ChemicalAmbiguousNameStats:
    """Summary counts from the ambiguous-name→candidate build."""

    ambiguous_names: int = 0
    candidate_links: int = 0


def rebuild_chemical_ambiguous_name_candidates(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
    progress: bool = False,
) -> ChemicalAmbiguousNameStats:
    """Attach ambiguous trivial names to their candidate structures.

    A trivial Name/Synonym borne by **>1 distinct chemical structure** (full
    InChIKey) is ambiguous (e.g. ``alanine`` → L-/D-/racemic). Rather than
    collapsing distinct compounds under the name or dropping it, link the name to
    **each** candidate structure entity (``resolution_mechanism =
    'ambiguous_name_class'``) — the destination for ambiguous/variable chemicals
    now that the SMILES tier is gone and the fallback abstains on ambiguity.
    Broad ChEBI **class** nodes are retained separately as ontology nodes
    (``entity_ontology_relation``, class → member edges) — this builds only the
    name→candidate attachment.

    Reads ``chemical_resolution_group_member`` (full level) for the structure
    entities + their InChIKeys, so it must run AFTER
    :func:`rebuild_chemical_resolution_levels`. The name scan is restricted to
    those structure entities (a candidate must be a resolved structure).
    """

    schema_id = sql.Identifier(schema)
    chem_type = sql.Literal(CHEMICAL_ENTITY_TYPE)
    ik_type = sql.Literal(STANDARD_INCHI_KEY_TYPE)
    ik_regex = sql.Literal(INCHIKEY_REGEX)
    name_types = sql.SQL(', ').join(
        sql.Literal(name) for name in CHEMICAL_NAME_TYPES
    )
    started = time.perf_counter()

    with conn.cursor() as cur:
        ensure_chemical_resolution_schema(cur, schema)
        cur.execute(
            sql.SQL(
                'TRUNCATE {0}.chemical_ambiguous_name_candidate'
            ).format(schema_id)
        )
        cur.execute(
            sql.SQL(
                """
                INSERT INTO {0}.chemical_ambiguous_name_candidate
                  (ambiguous_name, candidate_entity_id, candidate_inchikey,
                   resolution_mechanism)
                WITH name_type AS (
                  SELECT identifier_type_id FROM {0}.vocab_identifier_type
                  WHERE name IN ({4})
                ),
                -- structure entities + their InChIKey from the already-built
                -- full-level groups (one row per chemical structure entity).
                struct AS (
                  SELECT m.entity_id, m.inchikey
                  FROM {0}.chemical_resolution_group_member m
                  JOIN {0}.chemical_resolution_level l
                    ON l.level_id = m.level_id AND l.name = 'full'
                ),
                -- (structure entity, trivial name, inchikey), junk-guarded.
                struct_name AS (
                  SELECT DISTINCT
                    s.entity_id,
                    s.inchikey,
                    lower(trim(nm.value)) AS ambiguous_name
                  FROM struct s
                  JOIN {0}.entity_identifier_lookup eil
                    ON eil.entity_id = s.entity_id
                  JOIN {0}.identifier_evidence nm
                    ON nm.identifier_id = eil.identifier_id
                   AND nm.identifier_type_id IN (SELECT identifier_type_id FROM name_type)
                  WHERE length(trim(nm.value)) >= 2
                    AND trim(nm.value) !~ '^[0-9]+$'
                    AND trim(nm.value) !~ {3}
                ),
                -- ambiguous = a name spanning >1 distinct structure (InChIKey).
                ambiguous AS (
                  SELECT ambiguous_name
                  FROM struct_name
                  GROUP BY ambiguous_name
                  HAVING count(DISTINCT inchikey) > 1
                )
                SELECT
                  sn.ambiguous_name,
                  sn.entity_id,
                  min(sn.inchikey) AS candidate_inchikey,
                  'ambiguous_name_class' AS resolution_mechanism
                FROM struct_name sn
                JOIN ambiguous a USING (ambiguous_name)
                GROUP BY sn.ambiguous_name, sn.entity_id
                """
            ).format(schema_id, chem_type, ik_type, ik_regex, name_types)
        )
        links = int(cur.rowcount)
        cur.execute(
            sql.SQL(
                'SELECT count(DISTINCT ambiguous_name) '
                'FROM {0}.chemical_ambiguous_name_candidate'
            ).format(schema_id)
        )
        names = int(cur.fetchone()[0])
        cur.execute(
            sql.SQL(
                """
                CREATE INDEX IF NOT EXISTS
                  chemical_ambiguous_name_candidate_entity_idx
                ON {0}.chemical_ambiguous_name_candidate (candidate_entity_id)
                """
            ).format(schema_id)
        )

    conn.commit()
    _log(
        progress,
        'ambiguous_name_candidates',
        'done',
        ambiguous_names=names,
        candidate_links=links,
        seconds=f'{time.perf_counter() - started:.3f}',
    )
    return ChemicalAmbiguousNameStats(
        ambiguous_names=names,
        candidate_links=links,
    )


def _create_chemical_resolution_indexes(
    cur: psycopg2.extensions.cursor,
    schema: str,
) -> None:
    schema_id = sql.Identifier(schema)
    statements = (
        """
        CREATE INDEX IF NOT EXISTS chemical_resolution_group_member_entity_idx
        ON {0}.chemical_resolution_group_member (entity_id)
        """,
        """
        CREATE INDEX IF NOT EXISTS chemical_resolution_group_member_group_idx
        ON {0}.chemical_resolution_group_member (level_id, group_key)
        """,
        """
        CREATE INDEX IF NOT EXISTS chemical_resolution_relation_subject_idx
        ON {0}.chemical_resolution_relation (level_id, subject_entity_id)
        """,
        """
        CREATE INDEX IF NOT EXISTS chemical_resolution_relation_object_idx
        ON {0}.chemical_resolution_relation (level_id, object_entity_id)
        """,
    )
    for statement in statements:
        cur.execute(sql.SQL(statement).format(schema_id))
