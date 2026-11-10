"""OmniPath from Python: the hosted API (the default), or resources installed locally and
queried with DuckDB in the database the API serves them from."""

from .client import Client, connect
from .remote import ApiClient
from .errors import ClientError
from .installation import available, data_directory, init

__all__ = ["ApiClient", "Client", "ClientError", "available", "connect", "data_directory", "init"]
__version__ = "0.1.0"
