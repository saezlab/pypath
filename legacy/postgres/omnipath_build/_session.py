"""Session and logging setup via ``pkg_infra``.

Importing this module configures the whole logging tree once, so every other
module in the package takes a plain ``logging.getLogger(__name__)`` and
inherits the ``pkg_infra`` handlers and format. The package ``__init__``
imports it — without that import the configuration never runs, which is
exactly the bug the omnipath-utils reference implementation carries.
"""

from __future__ import annotations

import logging

from pkg_infra.session import get_session

session = get_session(workspace = '.')

# `pkg_infra` installs the root handlers but leaves the root logger at WARNING,
# which would silence every INFO line the build emits — including the derive
# steps' structured `step=… event=…` progress output that the sign-conflict
# summary and the manifest's derive-cost figures are read from. Set the level
# the package already expected on the package logger, where it is not at the
# mercy of a host that never called `basicConfig` (and where a `basicConfig`
# call would be a no-op anyway, root now carrying handlers).
logging.getLogger(__package__).setLevel(logging.INFO)
