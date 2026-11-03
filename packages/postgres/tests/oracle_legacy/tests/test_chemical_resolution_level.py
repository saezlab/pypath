"""Chemical resolution-level grouping tests (spec-003 Phase 6).

Two layers:

* **Pure unit** (no DB): the InChIKey-prefix level-key logic on the alanine
  benchmark structures — L/D/racemic collapse at ``connectivity``; beta-alanine
  and the N-acetyl (residue-context) form stay distinct at every level.

* **Integration** (needs ``DATABASE_URL``): a synthetic schema is built in the
  connected Postgres, the real :func:`rebuild_chemical_resolution_levels` runs
  over it, and the materialised group/member/relation tables are asserted —
  including that a relation expressed on different sub-level structures by
  different resources collapses to **one** projected edge carrying per-resource
  provenance. Needs only a Postgres connection, not a full build.
"""

from __future__ import annotations

import os

import pytest

from omnipath_build.chemical_resolution_level import (
    LEVELS,
    LEVELS_BY_NAME,
    resolution_level_key,
)
from tests.fixtures.resolution_benchmarks import (
    ALANINE_CONNECTIVITY_BLOCK1,
    ALANINE_CONNECTIVITY_GROUP,
    ALANINE_DISTINCT_STRUCTURES,
    BETA_ALANINE,
    DL_ALANINE,
    D_ALANINE,
    L_ALANINE,
    N_ACETYL_L_ALANINE,
)


# ---------------------------------------------------------------------------
# Pure unit — level-key logic
# ---------------------------------------------------------------------------


def test_levels_seed_shape():
    names = [level.name for level in LEVELS]
    assert names == ['connectivity', 'stereo_isotope_tautomer', 'full']
    assert [level.inchikey_prefix_length for level in LEVELS] == [14, 25, 27]
    # ranks are strictly increasing (coarse -> fine).
    ranks = [level.specificity_rank for level in LEVELS]
    assert ranks == sorted(ranks) and len(set(ranks)) == len(ranks)


def test_connectivity_collapses_stereo_and_charge_variants():
    keys = {
        resolution_level_key(s.inchikey, 'connectivity')
        for s in ALANINE_CONNECTIVITY_GROUP
    }
    # L/D/racemic alanine all share one connectivity group key (block 1).
    assert keys == {ALANINE_CONNECTIVITY_BLOCK1}
    # beta-alanine and the residue-context form are a *different* skeleton.
    distinct = {
        resolution_level_key(s.inchikey, 'connectivity')
        for s in ALANINE_DISTINCT_STRUCTURES
    }
    assert ALANINE_CONNECTIVITY_BLOCK1 not in distinct
    assert len(distinct) == len(ALANINE_DISTINCT_STRUCTURES)


def test_stereo_level_keeps_stereoisomers_distinct():
    keys = {
        resolution_level_key(s.inchikey, 'stereo_isotope_tautomer')
        for s in ALANINE_CONNECTIVITY_GROUP
    }
    # L, D and DL differ in block 2 → three distinct stereo-level groups.
    assert len(keys) == 3
    # each is the 25-char prefix (through the 2nd block, up to the 2nd dash).
    for key in keys:
        assert len(key) == 25 and key.startswith(ALANINE_CONNECTIVITY_BLOCK1)


def test_full_level_is_the_whole_key():
    for s in ALANINE_CONNECTIVITY_GROUP + ALANINE_DISTINCT_STRUCTURES:
        assert resolution_level_key(s.inchikey, 'full') == s.inchikey


def test_malformed_inchikey_yields_no_key():
    for bad in ('', 'not-an-inchikey', 'QNAYBMKLOCPYGJ', '12345', None):
        for level in LEVELS_BY_NAME:
            assert resolution_level_key(bad, level) is None


# ---------------------------------------------------------------------------
# Integration — materialised tables over a synthetic Postgres schema
# ---------------------------------------------------------------------------

DATABASE_URL = os.environ.get('DATABASE_URL')

pg = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; resolution-level build test needs Postgres',
)

# Deterministic ids; the connectivity representative is min(entity_id) = E_L.
E_L = '11111111-1111-1111-1111-111111111111'
E_D = '22222222-2222-2222-2222-222222222222'
E_DL = '33333333-3333-3333-3333-333333333333'
E_B = '44444444-4444-4444-4444-444444444444'
E_N = '55555555-5555-5555-5555-555555555555'
E_T = '99999999-9999-9999-9999-999999999999'  # protein target (not grouped)

CHEM_TYPE_ID = 1
GENE_TYPE_ID = 2
INCHIKEY_TYPE_ID = 10
UNIPROT_TYPE_ID = 11
NAME_TYPE_ID = 16
PREDICATE_ID = 7
SRC_A, SRC_B, SRC_C = 101, 102, 103

# (struct, trivial name) — L/D/DL share the ambiguous name 'Alanine'.
_CHEMS = {
    E_L: (L_ALANINE, 'Alanine'),
    E_D: (D_ALANINE, 'Alanine'),
    E_DL: (DL_ALANINE, 'Alanine'),
    E_B: (BETA_ALANINE, 'beta-Alanine'),
    E_N: (N_ACETYL_L_ALANINE, 'N-acetyl-L-alanine'),
}


@pytest.fixture()
def synthetic_schema():
    import psycopg2

    schema = 'crl_test_phase6'
    conn = psycopg2.connect(DATABASE_URL)
    conn.autocommit = False
    with conn.cursor() as cur:
        cur.execute(f'DROP SCHEMA IF EXISTS {schema} CASCADE')
        cur.execute(f'CREATE SCHEMA {schema}')
        cur.execute(
            f"""
            CREATE TABLE {schema}.vocab_entity_type (
              entity_type_id bigint PRIMARY KEY, name text);
            CREATE TABLE {schema}.vocab_identifier_type (
              identifier_type_id bigint PRIMARY KEY, name text);
            CREATE TABLE {schema}.entity (
              entity_id uuid PRIMARY KEY,
              entity_type_id bigint,
              taxonomy_id bigint,
              canonical_identifier_type_id bigint,
              canonical_identifier text);
            CREATE TABLE {schema}.identifier_evidence (
              identifier_id uuid PRIMARY KEY,
              identifier_type_id bigint, value text);
            CREATE TABLE {schema}.entity_identifier_lookup (
              entity_id uuid, identifier_id uuid);
            CREATE TABLE {schema}.relation (
              relation_id uuid PRIMARY KEY,
              subject_entity_id uuid, predicate_id bigint,
              object_entity_id uuid, relation_category_id bigint);
            CREATE TABLE {schema}.relation_evidence_relation (
              source_id bigint, relation_evidence_id uuid, relation_id uuid);
            """
        )
        cur.execute(
            f'INSERT INTO {schema}.vocab_entity_type VALUES (%s,%s),(%s,%s)',
            (CHEM_TYPE_ID, 'Chemical:OM:0037', GENE_TYPE_ID, 'Gene:MI:0250'),
        )
        cur.execute(
            f'INSERT INTO {schema}.vocab_identifier_type '
            f'VALUES (%s,%s),(%s,%s),(%s,%s)',
            (
                INCHIKEY_TYPE_ID, 'Standard Inchi Key:MI:1101',
                UNIPROT_TYPE_ID, 'Uniprot:MI:1097',
                NAME_TYPE_ID, 'Name:OM:0202',
            ),
        )
        for n, (eid, (struct, name)) in enumerate(_CHEMS.items()):
            cur.execute(
                f'INSERT INTO {schema}.entity VALUES (%s,%s,%s,%s,%s)',
                (eid, CHEM_TYPE_ID, None, INCHIKEY_TYPE_ID, struct.inchikey),
            )
            # attach the InChIKey and the trivial name as identifiers.
            ik_id = f'cccccccc-0000-0000-0000-0000000000{n:02d}'
            nm_id = f'dddddddd-0000-0000-0000-0000000000{n:02d}'
            cur.execute(
                f'INSERT INTO {schema}.identifier_evidence VALUES '
                f'(%s,%s,%s),(%s,%s,%s)',
                (ik_id, INCHIKEY_TYPE_ID, struct.inchikey,
                 nm_id, NAME_TYPE_ID, name),
            )
            cur.execute(
                f'INSERT INTO {schema}.entity_identifier_lookup VALUES '
                f'(%s,%s),(%s,%s)',
                (eid, ik_id, eid, nm_id),
            )
        # the protein target — a non-chemical endpoint, never grouped.
        cur.execute(
            f'INSERT INTO {schema}.entity VALUES (%s,%s,%s,%s,%s)',
            (E_T, GENE_TYPE_ID, 9606, UNIPROT_TYPE_ID, 'P00533'),
        )
        # relations: L/D/DL/beta-alanine each -> the same target, varied source.
        rels = [
            ('aaaaaaaa-0000-0000-0000-000000000001', E_L, SRC_A),
            ('aaaaaaaa-0000-0000-0000-000000000002', E_D, SRC_B),
            ('aaaaaaaa-0000-0000-0000-000000000003', E_DL, SRC_A),
            ('aaaaaaaa-0000-0000-0000-000000000004', E_B, SRC_A),
            ('aaaaaaaa-0000-0000-0000-000000000005', E_N, SRC_C),
        ]
        for i, (rid, subj, src) in enumerate(rels):
            cur.execute(
                f'INSERT INTO {schema}.relation VALUES (%s,%s,%s,%s,%s)',
                (rid, subj, PREDICATE_ID, E_T, None),
            )
            cur.execute(
                f'INSERT INTO {schema}.relation_evidence_relation '
                f'VALUES (%s,%s,%s)',
                (src, f'bbbbbbbb-0000-0000-0000-00000000000{i + 1}', rid),
            )
    conn.commit()
    try:
        yield conn, schema
    finally:
        with conn.cursor() as cur:
            cur.execute(f'DROP SCHEMA IF EXISTS {schema} CASCADE')
        conn.commit()
        conn.close()


def _rows(conn, query):
    with conn.cursor() as cur:
        cur.execute(query)
        return cur.fetchall()


@pg
def test_build_groups_collapse_and_stay_distinct(synthetic_schema):
    from omnipath_build.chemical_resolution_level import (
        rebuild_chemical_resolution_levels,
    )

    conn, schema = synthetic_schema
    rebuild_chemical_resolution_levels(conn, schema=schema)

    # connectivity: one group of 3 (alanine) + beta + N-acetyl singletons.
    conn_groups = {
        key: count
        for key, count in _rows(
            conn,
            f"""
            SELECT group_key, member_count
            FROM {schema}.chemical_resolution_group g
            JOIN {schema}.chemical_resolution_level l USING (level_id)
            WHERE l.name = 'connectivity'
            """,
        )
    }
    assert conn_groups.get(ALANINE_CONNECTIVITY_BLOCK1) == 3
    assert conn_groups.get(BETA_ALANINE.block1) == 1
    assert conn_groups.get(N_ACETYL_L_ALANINE.block1) == 1

    # members of the alanine connectivity group are recoverable.
    members = {
        eid
        for (eid,) in _rows(
            conn,
            f"""
            SELECT entity_id
            FROM {schema}.chemical_resolution_group_member m
            JOIN {schema}.chemical_resolution_level l USING (level_id)
            WHERE l.name = 'connectivity'
              AND m.group_key = '{ALANINE_CONNECTIVITY_BLOCK1}'
            """,
        )
    }
    assert members == {E_L, E_D, E_DL}

    # stereo level: the three alanines are three distinct groups.
    stereo_alanine = _rows(
        conn,
        f"""
        SELECT count(*)
        FROM {schema}.chemical_resolution_group g
        JOIN {schema}.chemical_resolution_level l USING (level_id)
        WHERE l.name = 'stereo_isotope_tautomer'
          AND g.group_key LIKE '{ALANINE_CONNECTIVITY_BLOCK1}%'
        """,
    )[0][0]
    assert stereo_alanine == 3

    # full level: one group per distinct structure (5 total).
    full_total = _rows(
        conn,
        f"""
        SELECT count(*)
        FROM {schema}.chemical_resolution_group g
        JOIN {schema}.chemical_resolution_level l USING (level_id)
        WHERE l.name = 'full'
        """,
    )[0][0]
    assert full_total == len(_CHEMS)


@pg
def test_relation_collapses_with_union_provenance(synthetic_schema):
    from omnipath_build.chemical_resolution_level import (
        rebuild_chemical_resolution_levels,
    )

    conn, schema = synthetic_schema
    rebuild_chemical_resolution_levels(conn, schema=schema)

    rep = _rows(
        conn,
        f"""
        SELECT representative_entity_id
        FROM {schema}.chemical_resolution_group g
        JOIN {schema}.chemical_resolution_level l USING (level_id)
        WHERE l.name = 'connectivity'
          AND g.group_key = '{ALANINE_CONNECTIVITY_BLOCK1}'
        """,
    )[0][0]
    assert str(rep) == E_L  # deterministic min(entity_id)

    # connectivity: the 3 alanine relations collapse to ONE projected edge
    # (rep -> target) carrying both source resources (A and B; r3 is A again).
    collapsed = _rows(
        conn,
        f"""
        SELECT member_relation_count, source_count, source_ids
        FROM {schema}.chemical_resolution_relation r
        JOIN {schema}.chemical_resolution_level l USING (level_id)
        WHERE l.name = 'connectivity'
          AND r.subject_entity_id = '{E_L}'
          AND r.object_entity_id = '{E_T}'
        """,
    )
    assert len(collapsed) == 1
    member_count, source_count, source_ids = collapsed[0]
    assert member_count == 3
    assert source_count == 2
    assert sorted(source_ids) == [SRC_A, SRC_B]

    # beta-alanine is a separate connectivity edge (own subject, single source).
    beta = _rows(
        conn,
        f"""
        SELECT member_relation_count, source_count
        FROM {schema}.chemical_resolution_relation r
        JOIN {schema}.chemical_resolution_level l USING (level_id)
        WHERE l.name = 'connectivity'
          AND r.subject_entity_id = '{E_B}'
        """,
    )
    assert beta == [(1, 1)]

    # full level: the alanine relations stay distinct (no collapse).
    full_alanine_edges = _rows(
        conn,
        f"""
        SELECT count(*)
        FROM {schema}.chemical_resolution_relation r
        JOIN {schema}.chemical_resolution_level l USING (level_id)
        WHERE l.name = 'full'
          AND r.subject_entity_id IN ('{E_L}', '{E_D}', '{E_DL}')
          AND r.object_entity_id = '{E_T}'
        """,
    )[0][0]
    assert full_alanine_edges == 3


@pg
def test_ambiguous_name_attaches_to_candidates(synthetic_schema):
    # The trivial name 'alanine' is borne by L/D/DL (3 distinct
    # structures) → it attaches to all three candidates; an unambiguous name
    # (one structure) produces no candidate links.
    from omnipath_build.chemical_resolution_level import (
        rebuild_chemical_ambiguous_name_candidates,
        rebuild_chemical_resolution_levels,
    )

    conn, schema = synthetic_schema
    # the candidate build reads the full-level groups for structure InChIKeys.
    rebuild_chemical_resolution_levels(conn, schema=schema)
    rebuild_chemical_ambiguous_name_candidates(conn, schema=schema)

    candidates = {
        eid
        for (eid,) in _rows(
            conn,
            f"""
            SELECT candidate_entity_id
            FROM {schema}.chemical_ambiguous_name_candidate
            WHERE ambiguous_name = 'alanine'
            """,
        )
    }
    assert candidates == {E_L, E_D, E_DL}

    mechanisms = {
        m
        for (m,) in _rows(
            conn,
            f'SELECT DISTINCT resolution_mechanism '
            f'FROM {schema}.chemical_ambiguous_name_candidate',
        )
    }
    assert mechanisms == {'ambiguous_name_class'}

    # beta-alanine / N-acetyl names map to a single structure → not ambiguous.
    unambiguous = _rows(
        conn,
        f"""
        SELECT count(*)
        FROM {schema}.chemical_ambiguous_name_candidate
        WHERE ambiguous_name IN ('beta-alanine', 'n-acetyl-l-alanine')
        """,
    )[0][0]
    assert unambiguous == 0
