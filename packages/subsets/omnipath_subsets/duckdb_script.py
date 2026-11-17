"""Run a file of DuckDB statements, giving each statement only the parameters it names."""

from __future__ import annotations

import duckdb


def run_script(con: duckdb.DuckDBPyConnection, text: str, parameters: dict | None = None) -> None:
    parameters = parameters or {}
    for statement in con.extract_statements(text):
        names = statement.named_parameters
        con.execute(statement.query, {name: parameters[name] for name in names} if names else None)
