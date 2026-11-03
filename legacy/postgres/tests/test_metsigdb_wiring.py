"""MetSigDB is a step of the build, not something a person remembers to run.

Cycle 010 shipped the substrate as a Python entry point with no caller, so a
fresh build produced a database with no ``metsigdb_membership`` at all and the
API answered 503. This covers the wiring that fixed it, in the shape cycle 008
established for the network views:

- the ``derive`` runs it by default, and ``--no-metsigdb`` opts out;
- it runs **after** the build manifest, whose stamp every published row carries;
- a failure is loud but does not abort the build, exactly like the network
  views and unlike the interaction projection;
- ``metsigdb`` is also a subcommand of its own, for a rebuild without a derive;
- ``reset-content`` truncates the substrate, so a reload cannot leave rows
  keyed on entity ids that no longer exist.

No database is touched: the build step is stubbed.

    uv run --with pytest pytest tests/test_metsigdb_wiring.py -v
"""

from __future__ import annotations

import inspect
import logging
import re

import pytest

import omnipath_build  # noqa: F401  (imports `_session`, configuring logging)
from omnipath_build import cli
from omnipath_build.db.schema import CONTENT_TABLES
from omnipath_build.metsigdb import BuildStats, ResourceLoadStats

STATS = BuildStats(
    build_id='0ed5afcb26c6',
    partial=False,
    rows=3_508_801,
    seconds=99.4,
    resources=(
        ResourceLoadStats(
            resource='KEGG',
            rows=4969,
            sets=176,
            metabolites=1799,
            removed=0,
            seconds=1.5,
        ),
    ),
)


class _FakeConn:
    """The little the CLI asks of a connection before the step runs."""

    def __init__(self) -> None:
        self.rolled_back = False

    def __enter__(self) -> '_FakeConn':
        return self

    def __exit__(self, *exc_info: object) -> bool:
        return False

    def rollback(self) -> None:
        self.rolled_back = True


@pytest.fixture
def stubbed(monkeypatch):
    """A CLI whose database and session tuning are inert."""
    conn = _FakeConn()
    monkeypatch.setattr(cli.psycopg2, 'connect', lambda url: conn)
    monkeypatch.setattr(cli, '_apply_build_session_tuning', lambda conn: None)
    return conn


def test_the_subcommand_runs_the_step(stubbed, monkeypatch, capsys):
    """`omnipath_build.cli metsigdb` publishes the substrate and says what it wrote."""
    seen = {}

    def _build(conn, *, max_records=None, schema='public', **kwargs):
        seen['schema'] = schema
        seen['max_records'] = max_records
        return STATS

    monkeypatch.setattr(cli, 'build_metsigdb', _build)

    assert cli.main(['--database-url', 'postgresql:///x', 'metsigdb']) == 0

    assert seen['schema'] == 'public'
    out = capsys.readouterr().out
    assert 'build=0ed5afcb26c6' in out
    assert 'rows=3508801' in out
    assert 'KEGG=4969' in out


def test_the_subcommand_applies_a_record_cap(stubbed, monkeypatch):
    """A capped run is for the dev loop, and the step is told the cap."""
    seen = {}
    monkeypatch.setattr(
        cli,
        'build_metsigdb',
        lambda conn, *, max_records=None, schema='public': (
            seen.update(max_records=max_records) or STATS
        ),
    )

    cli.main(
        ['--database-url', 'postgresql:///x', 'metsigdb', '--max-records', '1000']
    )

    assert seen['max_records'] == 1000


@pytest.mark.parametrize(
    ('value', 'expected'),
    [(None, None), ('', None), ('0', None), ('1000', 1000), ('nonsense', None)],
)
def test_an_unusable_cap_means_no_cap(value, expected):
    """An unset or unparseable MAX_RECORDS is the authoritative case, not an error."""
    assert cli._record_cap(value) == expected


def test_the_derive_runs_it_by_default_and_can_be_told_not_to():
    """The flag exists, defaults to on, and is a BooleanOptionalAction."""
    source = inspect.getsource(cli.main)
    flag = source.index("'--metsigdb'")
    declaration = source[flag:flag + 260]
    assert 'BooleanOptionalAction' in declaration, (
        '--no-metsigdb has to be spellable, like --no-network-views'
    )
    assert 'default=True' in declaration, (
        'the substrate is part of the build; opting out is the explicit act'
    )


def _derive_step() -> str:
    """The derive's MetSigDB block, from its start line to the end of its finally.

    Anchored on the log event rather than on how the call is spelled, so a
    change to the call site does not silently stop testing it.
    """
    source = inspect.getsource(cli.main)
    start = source.index("_derive_log('metsigdb_start')")
    end = source.index('metsigdb_conn.close()', start)
    return source[start:end + len('metsigdb_conn.close()')]


def test_it_runs_after_the_manifest_it_stamps_rows_with():
    """`build_id()` reads build_manifest, so the order is a requirement."""
    source = inspect.getsource(cli.main)
    step = source.index("_derive_log('metsigdb_start')")
    assert source.index("'build_manifest_done'") < step


def test_the_derive_gives_the_step_its_own_connection():
    """The step stages temp tables and must set `temp_buffers` to size them.

    Postgres refuses that once a session has touched a temp table, and by the
    time the step runs the derive's session has touched many. Lending it that
    connection leaves the step on the 8 MB default, where it loads Reactome and
    WikiPathways and then dies on KEGG with "no empty local buffer available".
    """
    step = _derive_step()
    assert 'psycopg2.connect(args.database_url)' in step, (
        'the derive must open a connection for the step, not lend its own'
    )
    assert 'metsigdb_conn' in step, 'the step has to receive that connection'
    assert 'conn.commit()' in step, (
        'the derive owes a commit before a second session reads what it wrote'
    )
    assert 'metsigdb_conn.close()' in step, "the step's connection is leaked"


def test_the_derive_does_not_reapply_the_load_cap():
    """The substrate is a projection of what the load wrote, capped already."""
    assert 'max_records' not in _derive_step(), (
        'the derive passes a record cap it should not re-apply'
    )


def test_a_failure_is_loud_and_does_not_abort_the_build():
    """Like the network views: reported at error level, and the build goes on."""
    step = _derive_step()
    # The failure goes through the shared helper, which prints as well as logs
    # and records the step for the derive's exit code. See
    # tests/test_derive_supplementary_steps.py.
    assert '_supplementary_step_failed(' in step
    assert "'metsigdb'," in step
    assert 'conn.rollback()' in step
    assert not re.search(r'^\s*raise\b', step, re.MULTILINE), (
        'a published dataset must not abort the build that produced its content'
    )


def test_reset_content_truncates_the_substrate():
    """A reload cannot leave rows keyed on entity ids that no longer exist."""
    assert 'metsigdb_membership' in CONTENT_TABLES


def test_reset_content_includes_interaction_fact_resource():
    """interaction_fact_resource carries a foreign key to interaction, so a
    reset that truncates interaction without it fails with
    psycopg2.errors.FeatureNotSupported (caught live during a spec 011 T047
    rebuild) -- both must always be truncated together."""
    assert 'interaction_fact_resource' in CONTENT_TABLES


def test_reset_content_includes_every_data_source_dependent():
    """Five more tables carry a foreign key to data_source and were missing
    from CONTENT_TABLES for the same reason as interaction_fact_resource
    (caught live, same rebuild): each fails TRUNCATE the moment it has any
    rows, since Postgres refuses regardless of the referencing table's
    row count, not just when one is populated."""
    for table in (
        'data_source_license',
        'entity_annotation_relation_default',
        'identifier_authority',
        'identifier_role',
        'resource_overlap_summary',
    ):
        assert table in CONTENT_TABLES
