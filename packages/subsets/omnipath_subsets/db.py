"""Small PostgreSQL helpers shared by release loading and product builds."""

import re


def validate_schema(schema: str) -> None:
    if (
        not isinstance(schema, str)
        or not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]{0,62}", schema)
        or schema.lower() in {"public", "information_schema"}
        or schema.lower().startswith("pg_")
    ):
        raise ValueError("Destination must be a new, non-system PostgreSQL schema identifier")
