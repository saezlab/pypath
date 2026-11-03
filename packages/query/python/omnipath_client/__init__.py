"""Install versioned OmniPath Parquet resources and query them locally with DuckDB."""

from .client import Client, connect
from .errors import ClientError
from .installation import available, data_directory, init

__all__ = ["Client", "ClientError", "available", "connect", "data_directory", "init"]
__version__ = "0.1.0"
