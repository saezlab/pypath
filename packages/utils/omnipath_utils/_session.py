"""Session and logging setup via pkg_infra."""

from pkg_infra.session import get_session

session = get_session(workspace='.')
