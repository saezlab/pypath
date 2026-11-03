"""A supplementary step that fails cannot leave the build looking successful.

MetSigDB and the network views are caught on purpose: a published dataset over
the core content must not abort the build that produced the content it reads.
The gap was that nobody could tell it had happened. The error went to the
package logger, `make` exports `config/pkg_infra_quiet.yaml`, and a derive that
published nothing printed its usual summary and exited 0. A real build did
exactly that on 2026-09-02: the substrate was never written and the only trace
was `metsigdb_membership.last_vacuum` still reading a week old.

Two things close it. The failure prints on the channel the steps' own progress
already uses, which no logging config suppresses. And the derive exits 1, so a
caller that reads exit 0 as "the database is ready" is not misled.

No database is touched.

    uv run --with pytest pytest tests/test_derive_supplementary_steps.py -v
"""

from __future__ import annotations

import inspect
import logging

import omnipath_build  # noqa: F401  (imports `_session`, configuring logging)
from omnipath_build import cli


def test_a_failure_prints_where_logging_cannot_swallow_it(capsys, caplog):
    """The whole point: `print`, not only the logger `make` silences."""
    caplog.set_level(logging.INFO, logger='omnipath_build.cli')
    failed: list[str] = []

    cli._supplementary_step_failed(
        failed, 'metsigdb', RuntimeError('no empty local buffer available'), 7.2
    )

    out = capsys.readouterr().out
    assert 'STEP FAILED' in out
    assert 'step=metsigdb' in out
    assert 'no empty local buffer available' in out
    assert failed == ['metsigdb']
    # And still in the structured shape, for the sinks that do read it.
    errors = [r for r in caplog.records if r.levelno == logging.ERROR]
    assert any('event=metsigdb_failed' in r.getMessage() for r in errors)


def test_a_complete_derive_still_exits_zero(capsys):
    assert cli._derive_exit_code([]) == 0
    assert capsys.readouterr().out == ''


def test_a_derive_that_skipped_a_step_exits_non_zero(capsys):
    """The build completed and committed what ran; it was not the whole build."""
    assert cli._derive_exit_code(['metsigdb', 'network_views']) == 1
    out = capsys.readouterr().out
    assert '2 step(s) failed' in out
    assert 'metsigdb' in out and 'network_views' in out


def test_every_supplementary_step_reports_through_the_helper():
    """None of them may go back to logging the failure and moving on."""
    source = inspect.getsource(cli.main)
    for step in ('network_views', 'metsigdb', 'cosmos'):
        anchor = source.index(f"_derive_log('{step}_start')")
        tail = source[anchor:anchor + 2600]
        assert '_supplementary_step_failed(' in tail, step
        assert f"'{step}'," in tail, step


def test_the_derive_returns_the_computed_code():
    """A bare `return 0` inside the derive branch would undo all of the above.

    `main` keeps a final `return 0` fallthrough at its own indentation; this
    looks only inside the derive branch, where the code has to be computed.
    """
    source = inspect.getsource(cli.main)
    derive = source[source.index("if args.command == 'derive':"):]
    assert 'return _derive_exit_code(failed_steps)' in derive
    assert '\n            return 0\n' not in derive
