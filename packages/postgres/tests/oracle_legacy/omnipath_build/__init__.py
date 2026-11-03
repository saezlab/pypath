"""DuckDB/PostgreSQL OmniPath build pipeline."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

# Configures the logging tree once, for every module in the package. The
# import is the part that matters — without it `get_session` never runs and the
# build's structured progress output is silent — so it must not be removed.
from omnipath_build._session import session as _session

__all__ = ['TMPDIR_ENV', 'configure_build_tmpdir']

#: Environment variable that sets the default scratch/temp base for every
#: OmniPath builder. Overridable per-invocation with the ``--tmpdir`` flag.
TMPDIR_ENV = 'OMNIPATH_BUILD_TMPDIR'


def configure_build_tmpdir(tmpdir: str | os.PathLike | None = None) -> str | None:
    """Route all builder temporary files to a single base directory.

    Resolves the temp base from ``tmpdir`` (e.g. the ``--tmpdir`` CLI flag) or,
    when that is not given, from the ``OMNIPATH_BUILD_TMPDIR`` environment
    variable. When a value is found the directory is created and made the
    process-wide default for both Python's :mod:`tempfile` and child processes
    (via ``TMPDIR``), so every ``tempfile.TemporaryDirectory()`` and DuckDB
    file-backed spill lands there instead of the root filesystem's ``/tmp``.

    A full build can spill several hundred GB of DuckDB temp storage; pointing
    this at large scratch storage keeps that off ``/`` (which otherwise fills
    and can take down co-located services). Returns the resolved absolute path,
    or ``None`` when neither source is set (callers keep the system default).
    """

    tmpdir = tmpdir or os.environ.get(TMPDIR_ENV)
    if not tmpdir:
        return None
    path = Path(tmpdir).expanduser()
    path.mkdir(parents=True, exist_ok=True)
    resolved = str(path)
    tempfile.tempdir = resolved       # Python's tempfile.* default
    os.environ['TMPDIR'] = resolved   # inherited by child processes (DuckDB, ...)
    return resolved
