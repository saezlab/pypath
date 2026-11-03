"""Shared extension namespace for independently versioned release schemas."""

from psycopg2 import sql


def ensure_public_extension(cursor, name: str) -> None:
    """Install in public; require explicit setup for an existing wrong namespace.

    Never relocate, drop or recreate an existing extension as part of a build:
    other release schemas can own dependent data. Public is already included in
    the main derivation search path and preserves the main installed catalogue.
    """
    if name not in {"pg_trgm", "roaringbitmap"}:
        raise ValueError("Unsupported main PostgreSQL extension")
    cursor.execute(
        sql.SQL("CREATE EXTENSION IF NOT EXISTS {} WITH SCHEMA public").format(sql.Identifier(name))
    )
    cursor.execute(
        "SELECT namespace.nspname FROM pg_extension extension "
        "JOIN pg_namespace namespace ON namespace.oid=extension.extnamespace "
        "WHERE extension.extname=%s",
        (name,),
    )
    row = cursor.fetchone()
    if row is None or row[0] != "public":
        actual = row[0] if row else "missing"
        raise ValueError(
            f"Extension {name} is installed in schema {actual}; main-compatible "
            "release builds require public. Relocate it explicitly before building."
        )
