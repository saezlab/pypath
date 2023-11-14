import functools as _ft

from pypath_common import Logger
from pypath_common import log as _log, session as _session

session = _ft.partial(_session, 'pypath')
log = _ft.partial(_log, 'pypath')
