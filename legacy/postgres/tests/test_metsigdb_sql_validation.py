"""Direct SQL validation of the published MetSigDB substrate (cycle 010).

These are the V-queries of `specs/010-metsigdb-subset/sql-validation.md`, run as
tests. They run after every resource load, not only at the end. Run against a
built instance::

    DATABASE_URL=postgresql://omnipath:omnipath@localhost:55435/omnipath \
        uv run --with pytest --with psycopg2-binary \
        pytest tests/test_metsigdb_sql_validation.py -v
"""

from __future__ import annotations

import os

import pytest

DATABASE_URL = os.environ.get('DATABASE_URL')

pytestmark = pytest.mark.skipif(
    not DATABASE_URL,
    reason='DATABASE_URL not set; the MetSigDB substrate needs a built database',
)

TABLE = 'metsigdb_membership'

# Measured on 2026-08-25 and frozen into contracts/row-contract.md. A build that
# departs from these needs an explanation before it is accepted.
EXPECTED = {
    'Reactome': {'rows': 24235, 'sets': 2239, 'metabolites': 2191},
    'WikiPathways': {'rows': 5431, 'sets': 818, 'metabolites': 2789},
    'KEGG': {'rows': 4969, 'sets': 176, 'metabolites': 1799},
    'MACdb': {'rows': 20291, 'sets': 269, 'metabolites': 5389},
    'ClassyFire': {'rows': 3453875, 'sets': 2863, 'metabolites': 145937},
}

TOLERANCE = 0.02


@pytest.fixture(scope='module')
def conn():
    import psycopg2

    connection = psycopg2.connect(DATABASE_URL)
    connection.autocommit = True
    try:
        yield connection
    finally:
        connection.close()


def _rows(conn, query, params=None):
    with conn.cursor() as cur:
        cur.execute(query, params)
        return cur.fetchall()


def _loaded_resources(conn):
    return {name for (name,) in _rows(conn, f'SELECT DISTINCT resource FROM {TABLE}')}


def test_v1_every_row_has_a_canonical_metabolite(conn):
    rows = _rows(
        conn,
        f"""
        SELECT count(*) FROM {TABLE} m
        LEFT JOIN entity e
          ON e.entity_id = m.metabolite_entity_id AND e.entity_type_id = 5
        WHERE e.entity_id IS NULL
        """,
    )
    assert rows[0][0] == 0


def test_v2_no_duplicate_memberships(conn):
    rows = _rows(
        conn,
        f"""
        SELECT resource, set_source_id, metabolite_entity_id, count(*)
        FROM {TABLE} GROUP BY 1, 2, 3 HAVING count(*) > 1
        """,
    )
    assert rows == []


def test_v3_mandatory_fields_are_present(conn):
    rows = _rows(
        conn,
        f"""
        SELECT count(*) FROM {TABLE}
        WHERE metabolite_label IS NULL
           OR metabolite_entity_type IS NULL
           OR set_type IS NULL
           OR set_size IS NULL
           OR provenance_source IS NULL
           OR build_id IS NULL
        """,
    )
    assert rows[0][0] == 0


def test_v4_set_type_matches_the_resource_rule(conn):
    rows = _rows(
        conn,
        f"""
        SELECT resource, set_type, count(*) FROM {TABLE}
        WHERE (resource, set_type) NOT IN (
          ('KEGG', 'pathway'), ('Reactome', 'pathway'),
          ('WikiPathways', 'pathway'), ('MACdb', 'disease'),
          ('ClassyFire', 'chemical_class'))
        GROUP BY 1, 2
        """,
    )
    assert rows == []


def test_v5_set_size_matches_the_published_population(conn):
    rows = _rows(
        conn,
        f"""
        SELECT resource, set_source_id, min(set_size), count(*)
        FROM {TABLE} GROUP BY 1, 2
        HAVING min(set_size) <> count(*) OR min(set_size) <> max(set_size)
        """,
    )
    assert rows == []


def test_v6_only_in_scope_resources_appear(conn):
    assert _loaded_resources(conn) <= set(EXPECTED) | {
        'KEGG',
        'WikiPathways',
        'MACdb',
        'ClassyFire',
    }


def test_v7_organism_comes_from_the_source(conn):
    """Reactome is human; WikiPathways spans many species; the rest are null.

    An earlier version derived the organism from an `R-HSA-` prefix, which
    published null for every WikiPathways pathway even though the source
    records a species for most of them.
    """
    rows = _rows(
        conn,
        f"""
        SELECT resource, count(organism), count(DISTINCT organism)
        FROM {TABLE} GROUP BY 1
        """,
    )
    seen = {resource: (n, taxa) for resource, n, taxa in rows}
    if 'Reactome' in seen:
        assert seen['Reactome'][1] == 1
        assert _rows(
            conn, f"SELECT DISTINCT organism FROM {TABLE} WHERE resource='Reactome'"
        ) == [(9606,)]
    if 'WikiPathways' in seen:
        assert seen['WikiPathways'][0] > 0
        assert seen['WikiPathways'][1] > 10, 'WikiPathways is many-species'
    for resource in ('KEGG', 'MACdb', 'ClassyFire'):
        if resource in seen:
            assert seen[resource][0] == 0, resource


def test_v8_row_counts_match_the_frozen_population(conn):
    loaded = _loaded_resources(conn)
    measured = {
        resource: {'rows': rows, 'sets': sets, 'metabolites': metabolites}
        for resource, rows, sets, metabolites in _rows(
            conn,
            f"""
            SELECT resource, count(*), count(DISTINCT set_source_id),
                   count(DISTINCT metabolite_entity_id)
            FROM {TABLE} GROUP BY 1
            """,
        )
    }
    for resource, expected in EXPECTED.items():
        if resource not in loaded:
            pytest.skip(f'{resource} is not loaded yet')
        for key, want in expected.items():
            got = measured[resource][key]
            assert abs(got - want) <= want * TOLERANCE, (
                f'{resource}.{key}: {got} against a frozen {want}'
            )


def test_v9_every_row_carries_one_build_stamp(conn):
    if not _loaded_resources(conn):
        pytest.skip('the substrate is empty')
    stamps = _rows(conn, f'SELECT DISTINCT build_id FROM {TABLE}')
    assert len(stamps) == 1
    manifest = _rows(conn, 'SELECT build_id FROM build_manifest')
    assert stamps[0][0] == manifest[0][0]


def test_set_labels_come_from_a_term_or_from_the_name(conn):
    """A set is named by its ontology term, or failing that by its `Name`.

    Three resources are ontology terms and carry a label. WikiPathways
    pathways are not terms at all, and the source publishes their title as a
    `Name` identifier instead; cycle 012 reads that as the fallback, which took
    WikiPathways from no named set to 745 of 818.

    Two populations are still unnamed and neither is a naming defect:

    - **73 WikiPathways sets** have no name anywhere in the database and no
      duplicate entity that carries one. The name never arrived.
    - **All 176 KEGG sets** are anchored on an entity that carries no name,
      while a duplicate of the same accession does. 2,605 KEGG pathway
      accessions map to more than one entity. That is the entity-duplication
      defect cycle 011 owns, not something a label fallback can reach.

    Asserted as counts rather than as "every set is named", because the honest
    state is partial and a test that demanded completeness would have to be
    skipped instead of read.
    """
    rows = _rows(
        conn,
        f'SELECT resource, count(DISTINCT set_source_id), '
        f'count(DISTINCT set_source_id) FILTER (WHERE set_label IS NOT NULL) '
        f'FROM {TABLE} GROUP BY 1',
    )
    seen = {resource: (total, named) for resource, total, named in rows}

    for resource in ('Reactome', 'MACdb', 'ClassyFire'):
        if resource in seen:
            total, named = seen[resource]
            assert named == total, f'{resource}: {named} of {total} named'

    if 'WikiPathways' in seen:
        total, named = seen['WikiPathways']
        assert named > 0, 'the Name fallback stopped working'
        assert total - named <= 73, (
            f'WikiPathways: {total - named} unnamed, was 73. More sets lost '
            'their name than the recorded gap.'
        )

    if 'KEGG' in seen:
        # Zero until the duplicate pathway entities are merged upstream.
        assert seen['KEGG'][1] == 0, (
            'KEGG sets gained names: the entity duplication may be fixed, and '
            'this expectation should be raised rather than left passing.'
        )


def test_a_named_set_reads_as_a_name(conn):
    """A spot check, because a non-null column can still hold the identifier."""
    rows = _rows(
        conn,
        f"""
        SELECT DISTINCT set_source_id, set_label FROM {TABLE}
        WHERE set_source_id IN ('R-HSA-1059683', 'CHEMONTID:0000118', '46')
          AND set_label IS NOT NULL
        """,
    )
    labels = dict(rows)
    assert labels.get('R-HSA-1059683') == 'Interleukin-6 signaling'
    assert labels.get('CHEMONTID:0000118') == 'Ketones'
    assert labels.get('46') == 'Colorectal Cancer (CRC)'


def test_the_sub_type_comes_from_the_source(conn):
    """MACdb records a trait type, and KEGG's overview maps are a known list.

    MACdb calls all 269 of its traits a disease, and 116 of them are
    interventions, phenotypes, genotypes or gene abnormalities. The sub-type is
    where that difference lives, and it is read rather than inferred from a
    label.
    """
    seen = {
        (resource, sub_type)
        for resource, sub_type in _rows(
            conn,
            f'SELECT DISTINCT resource, set_sub_type FROM {TABLE} '
            f'WHERE set_sub_type IS NOT NULL',
        )
    }
    macdb = {s for r, s in seen if r == 'MACdb'}
    assert macdb == {
        'cancer', 'phenotype', 'medical intervention',
        'gene abnormality', 'genotype',
    }
    assert {s for r, s in seen if r == 'KEGG'} == {'overview_map', 'metabolic_map'}
    # No other resource publishes a structured sub-type in this build.
    assert {r for r, _ in seen} == {'MACdb', 'KEGG'}


def test_kegg_overview_maps_are_the_declared_ones(conn):
    """The eleven whole-metabolism maps, not a guess from an identifier shape."""
    rows = _rows(
        conn,
        f"SELECT DISTINCT set_source_id FROM {TABLE} "
        f"WHERE resource = 'KEGG' AND set_sub_type = 'overview_map'",
    )
    from omnipath_build.metsigdb.mapping import KEGG_OVERVIEW_MAPS

    assert {set_id for (set_id,) in rows} <= set(KEGG_OVERVIEW_MAPS)
