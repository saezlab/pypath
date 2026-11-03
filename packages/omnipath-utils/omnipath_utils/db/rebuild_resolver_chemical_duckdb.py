"""Rebuild ``resolver_chemical`` via DuckDB instead of native Postgres SQL.

Standalone maintenance script, not part of the normal build path
(``create_resolver_views`` in ``_build.py`` still owns the SQL-only path and
is what a routine reload calls). Reach for this only when that native
rebuild is unworkably slow.

Why this exists: on this project's sandbox database, running
``resolver_chemical.sql`` directly in Postgres ran for 14+ hours without
finishing. That query is a plain ``UNION`` plus per-row regex validation
over the ~130M rows in ``id_mapping`` that carry a chemical structure key.
No disk spilling occurred, and no error was raised. The query plan's own
cost estimate implied minutes, not hours. Nobody diagnosed the bottleneck,
and it did not reproduce with DuckDB.

The equivalent computation, run here via DuckDB's ``postgres`` extension
(``ATTACH ... TYPE POSTGRES``) reading the same tables, completed in under
two minutes end to end: read, compute, write back, on the same database.
DuckDB is a vectorized, multi-threaded analytical engine with its own
memory management. It never touches Postgres's ``/dev/shm``-based
parallel-worker mechanism, which independently crashed twice on this same
host at even two workers.

Mirrors ``resolver_chemical.sql`` exactly: the same four contributions
(direct, pubchem-hub, chebi-hub, pubchem-substance-hub), the same
per-namespace ``source_pattern`` validation, the same InChIKey syntax
check. Keep the two in sync. This is not a replacement for that file, only
an alternate way to run the same logic when Postgres itself cannot.

Usage::

    uv run --with duckdb python -m omnipath_utils.db.rebuild_resolver_chemical_duckdb \
        postgresql://user:pass@host:5432/omnipath_utils

Leaves behind ``omnipath_utils.resolver_chemical_staging``: the
materialized view's ``CREATE ... AS SELECT * FROM`` reads that table, so
Postgres keeps a catalog dependency on it. The data itself is fully
copied. This is only a definition-time reference. Do not drop the staging
table without ``CASCADE``, and know that doing so drops the materialized
view too.
"""

from __future__ import annotations

import sys
import time

_SOURCE_PATTERNS = (
    ('chebi', r'^CHEBI:\d+$'),
    ('hmdb', r'^HMDB\d{7}$'),
    ('swisslipids', r'^SLM:\d+$'),
    ('lipidmaps', r'^LM[A-Z0-9]+$'),
    ('kegg', r'^C\d{5}$'),
    ('pubchem', r'^([1-9]\d*|0)$'),
    ('pubchem_substance', r'^([1-9]\d*|0)$'),
    ('chembl', r'^CHEMBL\d+$'),
)

_INCHIKEY_PATTERN = r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$'


def rebuild(pg_url: str) -> int:
    """Compute and write ``resolver_chemical`` via DuckDB. Returns the row
    count written."""

    import duckdb

    con = duckdb.connect()
    con.execute('INSTALL postgres; LOAD postgres;')
    con.execute(f"ATTACH '{pg_url}' AS pg (TYPE POSTGRES)")

    started = time.time()
    con.execute('CREATE TEMP TABLE source_pattern (source_type VARCHAR, pattern VARCHAR)')
    con.executemany(
        'INSERT INTO source_pattern VALUES (?, ?)', list(_SOURCE_PATTERNS),
    )
    con.execute(
        'CREATE TEMP TABLE id_type_local AS SELECT id, name FROM pg.omnipath_utils.id_type'
    )
    con.execute("""
        CREATE TEMP TABLE pubchem_ik AS
        SELECT m.source_id AS hub, m.target_id AS inchikey
        FROM pg.omnipath_utils.id_mapping m
        JOIN id_type_local pc ON pc.id = m.source_type_id AND pc.name = 'pubchem'
        JOIN id_type_local ik ON ik.id = m.target_type_id AND ik.name = 'inchikey'
    """)
    con.execute("""
        CREATE TEMP TABLE chebi_ik AS
        SELECT m.source_id AS hub, m.target_id AS inchikey
        FROM pg.omnipath_utils.id_mapping m
        JOIN id_type_local ce ON ce.id = m.source_type_id AND ce.name = 'chebi'
        JOIN id_type_local ik ON ik.id = m.target_type_id AND ik.name = 'inchikey'
    """)
    con.execute(f"""
        CREATE TEMP TABLE result_all AS
        SELECT st.name AS source_type, m.source_id, m.target_id AS inchikey
        FROM pg.omnipath_utils.id_mapping m
        JOIN id_type_local st ON st.id = m.source_type_id
        JOIN id_type_local ik ON ik.id = m.target_type_id AND ik.name = 'inchikey'
        LEFT JOIN source_pattern sp ON sp.source_type = st.name
        WHERE m.source_id IS NOT NULL
          AND regexp_matches(m.target_id, '{_INCHIKEY_PATTERN}')
          AND (sp.pattern IS NULL OR regexp_matches(m.source_id, sp.pattern))

        UNION

        SELECT st.name AS source_type, m.source_id, h.inchikey
        FROM pg.omnipath_utils.id_mapping m
        JOIN id_type_local st ON st.id = m.source_type_id
        JOIN id_type_local pc ON pc.id = m.target_type_id AND pc.name = 'pubchem'
        JOIN pubchem_ik h ON h.hub = m.target_id
        LEFT JOIN source_pattern sp ON sp.source_type = st.name
        WHERE st.name IN ('chembl', 'hmdb', 'chebi', 'drugbank')
          AND m.source_id IS NOT NULL
          AND regexp_matches(h.inchikey, '{_INCHIKEY_PATTERN}')
          AND (sp.pattern IS NULL OR regexp_matches(m.source_id, sp.pattern))

        UNION

        SELECT st.name AS source_type, m.source_id, h.inchikey
        FROM pg.omnipath_utils.id_mapping m
        JOIN id_type_local st ON st.id = m.source_type_id
        JOIN id_type_local ce ON ce.id = m.target_type_id AND ce.name = 'chebi'
        JOIN chebi_ik h ON h.hub = m.target_id
        LEFT JOIN source_pattern sp ON sp.source_type = st.name
        WHERE st.name IN ('kegg', 'hmdb', 'reactome')
          AND m.source_id IS NOT NULL
          AND regexp_matches(h.inchikey, '{_INCHIKEY_PATTERN}')
          AND (sp.pattern IS NULL OR regexp_matches(m.source_id, sp.pattern))

        UNION

        SELECT st.name AS source_type, m.source_id, h.inchikey
        FROM pg.omnipath_utils.id_mapping m
        JOIN id_type_local st ON st.id = m.source_type_id AND st.name = 'pubchem_substance'
        JOIN id_type_local pc ON pc.id = m.target_type_id AND pc.name = 'pubchem'
        JOIN pubchem_ik h ON h.hub = m.target_id
        LEFT JOIN source_pattern sp ON sp.source_type = st.name
        WHERE m.source_id IS NOT NULL
          AND regexp_matches(h.inchikey, '{_INCHIKEY_PATTERN}')
          AND (sp.pattern IS NULL OR regexp_matches(m.source_id, sp.pattern))
    """)

    n = con.execute('SELECT count(*) FROM result_all').fetchone()[0]
    print(f'computed {n} rows in {time.time() - started:.1f}s', file=sys.stderr)

    con.execute('DROP TABLE IF EXISTS pg.omnipath_utils.resolver_chemical_staging')
    con.execute("""
        CREATE TABLE pg.omnipath_utils.resolver_chemical_staging (
            source_type VARCHAR, source_id VARCHAR, inchikey VARCHAR
        )
    """)
    con.execute(
        'INSERT INTO pg.omnipath_utils.resolver_chemical_staging SELECT * FROM result_all'
    )
    print(f'written to Postgres staging table in {time.time() - started:.1f}s', file=sys.stderr)
    print(
        'Next: run resolver_chemical.sql\'s own DROP-and-CREATE block, but '
        'with the CREATE MATERIALIZED VIEW query replaced by '
        "'SELECT source_type, source_id, inchikey FROM "
        "omnipath_utils.resolver_chemical_staging', then create the index.",
        file=sys.stderr,
    )
    return n


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f'usage: {sys.argv[0]} <postgres-url>', file=sys.stderr)
        raise SystemExit(2)
    rebuild(sys.argv[1])
