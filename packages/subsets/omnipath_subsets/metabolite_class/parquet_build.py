"""
Build the metabolite class tables as Parquet from a published release.

The product reads the pinned resources' Parquet through DuckDB. It resolves no identifier,
calls no service and downloads nothing, so the same release and the same rule file always
give the same rows.

The build runs in four steps: each source of the rule file claims molecules, the claims
collapse to one class per molecule, the connectivity groups spread the best class to every
member, and the readable columns and the report are added.
"""

from __future__ import annotations

import logging
from pathlib import Path

import duckdb

from omnipath_core.subset_layout import BuildContext
from omnipath_subsets.duckdb_script import run_script
from omnipath_subsets.exclusion_metabolites.version import list_version
from omnipath_subsets.metabolite_class.rules import Rules, Source, read_rules
from omnipath_subsets.release_source import CHEMICAL_ENTITY_TYPES, open_release

_SQL_DIR = Path(__file__).with_name('duckdb_sql')

logger = logging.getLogger(__name__)

TABLE = 'metabolite_class'
REPORT_TABLE = 'metabolite_class_report'

_CLAIMS = {
    'resource': 'claim_resource',
    'chebi_role': 'claim_chebi_role',
    'network': 'claim_network',
}


def _sql(name: str) -> str:
    """
    The text of the SQL file ``name``.
    """

    return (_SQL_DIR / f'{name}.sql').read_text(encoding = 'utf-8')


def _literals(values) -> str:
    """
    The values as a list of SQL string literals.
    """

    return ', '.join("'" + value.replace("'", "''") + "'" for value in values)


def _claim(con: duckdb.DuckDBPyConnection, source: Source) -> int:
    """
    Record the claims of one source.

    Args:
        con: The connection with the release views and the working tables.
        source: The source of the rule file.

    Returns:
        The number of molecules the source claims.
    """

    text = _sql(_CLAIMS[source.signal])

    if source.signal == 'network':

        text = text.replace('{resources}', _literals(source.resources))

    run_script(
        con,
        text,
        {
            'metabolite_class': source.metabolite_class,
            'priority': source.priority,
            'entity_types': list(CHEMICAL_ENTITY_TYPES),
            'resource': source.resource,
            'root': source.chebi_role_root,
            'predicates': list(source.predicates),
        },
    )
    claimed = con.execute(
        'SELECT count(DISTINCT entity_key) FROM mc_claim WHERE priority = ?',
        [source.priority],
    ).fetchone()[0]

    if claimed == 0:

        logger.warning(
            'metabolite_class: source of priority %s claims nothing',
            source.priority,
        )

    return claimed


def build_tables(con: duckdb.DuckDBPyConnection, *, build_id: str, rules: Rules) -> dict:
    """
    Fill ``metabolite_class`` and its report in ``con``.

    Args:
        con: A connection with the release views.
        build_id: The digest of the release manifest, stamped on every row.
        rules: The rule file.

    Returns:
        The counts of the build, for the subset manifest.
    """

    run_script(con, _sql('tables'))
    con.executemany(
        'INSERT INTO mc_class VALUES (?, ?)',
        [(c.name, c.precedence) for c in rules.classes],
    )
    claims = {source.priority: _claim(con, source) for source in rules.sources}

    run_script(con, _sql('assign'), {'entity_types': list(CHEMICAL_ENTITY_TYPES)})
    run_script(con, _sql('publish'), {'build_id': build_id})
    run_script(
        con,
        _sql('report'),
        {'build_id': build_id, 'classes': list(rules.claimable)},
    )

    sizes = dict(
        con.execute(f'SELECT metabolite_class, count(*) FROM {TABLE} GROUP BY 1').fetchall()
    )
    variants = con.execute(
        f"SELECT count(*) FROM {TABLE} WHERE match_kind = 'structure_key'"
    ).fetchone()[0]

    return {
        'build_id': build_id,
        'rows': sum(sizes.values()),
        'classes': {name: sizes.get(name, 0) for name in rules.claimable},
        'structure_key_rows': variants,
        'sources': [
            {
                'priority': s.priority,
                'class': s.metabolite_class,
                'signal': s.signal,
                'claimed': claims[s.priority],
            }
            for s in rules.sources
        ],
    }


def write_metabolite_class(
    directory: Path,
    context: BuildContext,
    *,
    seed_directory = None,
    rules_path: Path | str | None = None,
) -> dict:
    """
    Write both tables of the product for ``context.release`` into ``directory``.

    Args:
        directory: The working directory of the product build.
        context: The release the product is built from.
        seed_directory: The exclusion product's seed files; not read here.
        rules_path: The rule file; the packaged one by default.

    Returns:
        The counts of the build and the version of the rule file.
    """

    rules = read_rules(rules_path)
    version = list_version([rules.path])
    con = open_release(context.data_root, context.release)

    try:

        stats = build_tables(con, build_id = context.release_manifest_sha256, rules = rules)

        for table, order in (
            (TABLE, 'entity_key'),
            (REPORT_TABLE, 'kind, metabolite_class NULLS FIRST, detail'),
        ):

            target = str(directory / f'{table}.parquet').replace("'", "''")
            con.execute(
                f"COPY (SELECT * FROM {table} ORDER BY {order}) TO '{target}' "
                '(FORMAT parquet, COMPRESSION zstd)'
            )

    finally:

        con.close()

    stats['rules_version'] = version.value
    stats['uncommitted_rules'] = version.uncommitted
    logger.info('metabolite_class: %s', stats)

    return stats
