"""Main identifier authority publication using pinned metadata-only declarations.

The main algorithm is retained. No resource discovery, input parser, resolver
or download is called, and no authority is inferred when mints is unavailable.
"""

from __future__ import annotations


def populate_identifier_authority(
    cur,
    schema: str,
    configs: dict[str, object],
    *,
    progress: bool = False,
) -> int:
    """Write the ``mints`` declarations from ``discover_resources`` into
    ``identifier_authority``: one row per namespace a resource is the
    authority for.

    It looks up against this database's live ``vocab_identifier_type`` and
    ``data_source`` tables, not the resolver's static namespace list. Each
    build assigns its own identifier type IDs, and that narrower,
    resolver-only list does not cover every namespace a resource mints.

    The walk skips a namespace with no matching ``vocab_identifier_type``
    row, or a source with no ``data_source`` row, rather than raising an
    error. The resource walk runs before ingest populates those tables.

    Returns the number of rows written.
    """
    from psycopg2 import sql

    from pypath.internals.cv_terms import cv_term_label_accession

    # spec 011 T058: a mints declaration for one of the resolver's own
    # structure-bearing chemical namespaces (the same set resolver_chemical
    # supplies a lookup for) makes that resource the *structure* authority
    # for it, not just the identifier authority -- what
    # resolver_candidate_authority in duckdb_load.py checks to arbitrate a
    # genuine cross-skeleton disagreement.
    from .identifier_types import (  # noqa: PLC0415
        RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE,
    )

    structure_bearing_type_names = frozenset(RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE.values())

    schema_id = sql.Identifier(schema)
    written = 0

    for source_slug, config in configs.items():
        if not config.mints:
            continue

        cur.execute(
            sql.SQL("SELECT source_id FROM {}.data_source WHERE name = %s").format(schema_id),
            [source_slug],
        )
        source_row = cur.fetchone()
        if source_row is None:
            if progress:
                print(
                    f"[identifier_authority] no data_source row for {source_slug!r}, skipping",
                    flush=True,
                )
            continue
        source_id = source_row[0]

        for namespace in config.mints:
            type_name = cv_term_label_accession(namespace)
            cur.execute(
                sql.SQL(
                    "SELECT identifier_type_id FROM {}.vocab_identifier_type WHERE name = %s"
                ).format(schema_id),
                [type_name],
            )
            type_row = cur.fetchone()
            if type_row is None:
                if progress:
                    print(
                        f"[identifier_authority] no vocab_identifier_type "
                        f"row for {type_name!r} ({source_slug}), skipping",
                        flush=True,
                    )
                continue

            is_structure_authority = type_name in structure_bearing_type_names
            cur.execute(
                sql.SQL(
                    """
                    INSERT INTO {}.identifier_authority
                        (identifier_type_id, source_id, is_structure_authority)
                    VALUES (%s, %s, %s)
                    ON CONFLICT (identifier_type_id) DO UPDATE
                    SET source_id = EXCLUDED.source_id,
                        is_structure_authority = EXCLUDED.is_structure_authority
                    """
                ).format(schema_id),
                [type_row[0], source_id, is_structure_authority],
            )
            written += 1

    if progress:
        print(
            f"[identifier_authority] wrote {written} row(s) from {len(configs)} resource config(s)",
            flush=True,
        )
    return written
