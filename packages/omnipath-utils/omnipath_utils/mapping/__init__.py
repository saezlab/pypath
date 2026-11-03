"""ID translation -- the core service of omnipath-utils.

Example::

    from omnipath_utils.mapping import map_name, map_names

    # Translate gene symbol to UniProt
    map_name('TP53', 'genesymbol', 'uniprot')
    # {'P04637'}

    # Translate multiple
    map_names(['TP53', 'EGFR'], 'genesymbol', 'uniprot')
    # {'P04637', 'P00533'}

    # Translate a DataFrame column (pandas, polars, or pyarrow)
    import pandas as pd
    from omnipath_utils.mapping import translate_column

    df = pd.DataFrame({'protein': ['P04637', 'P00533']})
    translate_column(df, 'protein', 'uniprot', 'genesymbol')
"""

from __future__ import annotations

from typing import Any
from collections.abc import Iterable

import narwhals as nw

from omnipath_utils.mapping._mapper import Mapper


def map_name(
    name: str,
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int | None = None,
    strict: bool = False,
    raw: bool = False,
    backend: str | None = None,
    full_uniprot: str = 'fallback',
) -> set[str]:
    """Translate a single identifier.

    Args:
        name: The identifier to translate.
        id_type: Source ID type (e.g. 'genesymbol').
        target_id_type: Target ID type (e.g. 'uniprot').
        ncbi_tax_id: Organism (default: 9606).
        strict: Skip fuzzy fallbacks (gene symbol case, append "1", etc.).
        raw: Skip special-case handling. Direct table lookup only.
        backend: Force a specific backend (e.g. 'uniprot', 'biomart').
        full_uniprot: Use of the comprehensive full-UniProt table (DB mode):
            ``'fallback'`` (default), ``'never'``, ``'both'``, ``'only'``.

    Returns:
        Set of target identifiers.
    """

    from omnipath_utils.mapping._translate import translate_core

    result = translate_core(
        [name],
        id_type,
        target_id_type,
        ncbi_tax_id,
        strict=strict,
        raw=raw,
        backend=backend,
        full_uniprot=full_uniprot,
    )
    return result.get(name, set())


def map_names(
    names: object,
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int | None = None,
    raw: bool = False,
    backend: str | None = None,
    full_uniprot: str = 'fallback',
) -> set[str]:
    """Translate multiple identifiers, return union of results.

    Args:
        names: Iterable of identifiers.
        id_type: Source ID type.
        target_id_type: Target ID type.
        ncbi_tax_id: Organism (default: 9606).
        raw: Skip special-case handling. Direct table lookup only.
        backend: Force a specific backend.
        full_uniprot: Use of the comprehensive full-UniProt table (DB mode):
            ``'fallback'`` (default), ``'never'``, ``'both'``, ``'only'``.

    Returns:
        Union of all target identifiers.
    """

    from omnipath_utils.mapping._translate import translate_core

    result = translate_core(
        list(names),
        id_type,
        target_id_type,
        ncbi_tax_id,
        raw=raw,
        backend=backend,
        full_uniprot=full_uniprot,
    )
    merged = set()
    for targets in result.values():
        merged.update(targets)
    return merged


def map_name0(
    name: str,
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int | None = None,
    raw: bool = False,
    backend: str | None = None,
    full_uniprot: str = 'fallback',
) -> str | None:
    """Translate, returning a single result.

    Args:
        name: The identifier to translate.
        id_type: Source ID type.
        target_id_type: Target ID type.
        ncbi_tax_id: Organism (default: 9606).
        raw: Skip special-case handling.
        backend: Force a specific backend.
        full_uniprot: Use of the comprehensive full-UniProt table (DB mode):
            ``'fallback'`` (default), ``'never'``, ``'both'``, ``'only'``.

    Returns:
        A single target identifier, or None.
    """

    result = map_name(
        name,
        id_type,
        target_id_type,
        ncbi_tax_id,
        raw=raw,
        backend=backend,
        full_uniprot=full_uniprot,
    )
    return next(iter(result)) if result else None


def translate(
    identifiers: Iterable[str],
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int | None = None,
    raw: bool = False,
    backend: str | None = None,
    full_uniprot: str = 'fallback',
) -> dict[str, set[str]]:
    """Batch translate. Returns dict mapping source -> set of targets.

    Uses vectorized table lookup when available. Falls back to per-ID
    map_name (with full special-case handling) when the table is not
    available.

    Args:
        identifiers: Iterable of source identifiers.
        id_type: Source ID type.
        target_id_type: Target ID type.
        ncbi_tax_id: Organism (default: 9606).
        raw: Skip special-case handling. Direct table lookup only.
        backend: Force a specific backend.
        full_uniprot: Use of the comprehensive full-UniProt table (DB mode):
            ``'fallback'`` (default), ``'never'``, ``'both'``, ``'only'``.

    Returns:
        Dict mapping each source ID to a set of target IDs.
    """

    from omnipath_utils.mapping._translate import translate_core

    return translate_core(
        list(identifiers),
        id_type,
        target_id_type,
        ncbi_tax_id,
        raw=raw,
        full_uniprot=full_uniprot,
        backend=backend,
    )


def translation_table(
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int | None = None,
) -> dict[str, set[str]]:
    """Get the full translation table.

    Uses PostgreSQL if ``OMNIPATH_UTILS_DB_URL`` is set,
    otherwise falls back to in-memory backends.
    """

    from omnipath_utils.mapping._translate import _get_db_session

    session = _get_db_session()

    if session is not None:
        try:
            from omnipath_utils.db._query import get_full_table
            from omnipath_utils.mapping._id_types import IdTypeRegistry

            reg = IdTypeRegistry.get()
            ncbi_tax_id = ncbi_tax_id if ncbi_tax_id is not None else 9606
            src = reg.resolve(id_type) or id_type
            tgt = reg.resolve(target_id_type) or target_id_type
            return get_full_table(session, src, tgt, ncbi_tax_id)
        except Exception:
            pass
        finally:
            session.close()

    return Mapper.get().translation_table(
        id_type,
        target_id_type,
        ncbi_tax_id,
    )


def id_types() -> list[str]:
    """List all known ID types."""

    return Mapper.get().id_types()


def translate_column(
    df: Any,
    column: str,
    id_type: str,
    target_id_type: str,
    ncbi_tax_id: int | None = None,
    new_column: str | None = None,
    keep_untranslated: bool = True,
    expand: bool = True,
    raw: bool = False,
    backend: str | None = None,
    full_uniprot: str = 'fallback',
) -> Any:
    """Translate a column of identifiers in a DataFrame.

    Performs vectorized lookup using the full mapping table,
    not per-row translation. Accepts any DataFrame backend
    supported by narwhals (pandas, polars, or pyarrow).

    Args:
        df: DataFrame (pandas, polars, or pyarrow).
        column: Column name with source IDs.
        id_type: Source ID type.
        target_id_type: Target ID type.
        ncbi_tax_id: Organism (default: 9606).
        new_column: Output column name (default: target_id_type).
        keep_untranslated: Keep rows with no translation (default: True).
        expand: Expand one-to-many to multiple rows (default: True).
        raw: Skip special-case handling. Direct table lookup only.
        backend: Force a specific backend.

    Returns:
        DataFrame with translated column added, in the same
        backend format as the input.
    """

    from omnipath_utils.mapping._translate import translate_core

    mapper = Mapper.get()
    ncbi_tax_id = ncbi_tax_id or mapper.ncbi_tax_id
    new_col = new_column or target_id_type

    nw_df = nw.from_native(df, eager_only=True)
    native_ns = nw.get_native_namespace(nw_df)

    # Get unique IDs for batch translation
    unique_ids = nw_df[column].drop_nulls().unique().to_list()
    unique_ids = [str(v) for v in unique_ids]

    trans = translate_core(
        unique_ids,
        id_type,
        target_id_type,
        ncbi_tax_id,
        raw=raw,
        backend=backend,
        full_uniprot=full_uniprot,
    )

    if expand:
        # Build mapping rows: (source, target) for each pair
        src_col_vals = []
        tgt_col_vals = []

        for src, targets in trans.items():
            if targets:
                for tgt in sorted(targets):
                    src_col_vals.append(src)
                    tgt_col_vals.append(tgt)
            elif keep_untranslated:
                src_col_vals.append(src)
                tgt_col_vals.append(None)

        # Handle source values not present in trans (not queried)
        all_sources = set(unique_ids)
        for src in all_sources - set(trans.keys()):
            if keep_untranslated:
                src_col_vals.append(src)
                tgt_col_vals.append(None)

        tmp_col = nw.generate_temporary_column_name(
            n_bytes=8, columns=nw_df.columns,
        )

        mapping_native = native_ns.DataFrame(
            {tmp_col: src_col_vals, new_col: tgt_col_vals},
        )
        mapping_nw = nw.from_native(mapping_native, eager_only=True)

        # Cast column to string for joining
        result = nw_df.with_columns(
            nw.col(column).cast(nw.String).alias(tmp_col),
        )
        result = result.join(
            mapping_nw,
            on=tmp_col,
            how='left' if keep_untranslated else 'inner',
        )
        result = result.drop(tmp_col)

        if not keep_untranslated:
            result = result.filter(~nw.col(new_col).is_null())
    else:
        # No expansion: pick one result per source
        mapping = {}
        for src, targets in trans.items():
            if targets:
                mapping[src] = next(iter(sorted(targets)))
            else:
                mapping[src] = None

        values = []
        for v in nw_df[column].to_list():
            sv = str(v) if v is not None else None
            values.append(mapping.get(sv) if sv else None)

        new_series = nw.new_series(
            name=new_col,
            values=values,
            backend=native_ns,
        )
        result = nw_df.with_columns(new_series)

        if not keep_untranslated:
            result = result.filter(~nw.col(new_col).is_null())

    return nw.to_native(result)


def translate_columns(
    df: Any,
    *translations: tuple[str, ...],
    ncbi_tax_id: int | None = None,
    keep_untranslated: bool = True,
    expand: bool = True,
    raw: bool = False,
    backend: str | None = None,
    full_uniprot: str = 'fallback',
) -> Any:
    """Translate multiple columns in a DataFrame.

    Each translation is a tuple of (column, id_type, target_id_type)
    or (column, id_type, target_id_type, new_column_name).

    Args:
        df: DataFrame (pandas, polars, or pyarrow).
        *translations: Tuples defining translations.
        ncbi_tax_id: Organism (default: 9606).
        keep_untranslated: Keep rows with no translation.
        expand: Expand one-to-many mappings.
        raw: Skip special-case handling. Direct table lookup only.
        backend: Force a specific backend.

    Returns:
        DataFrame with all translated columns added, in the same
        backend format as the input.

    Example::

        translate_columns(
            df,
            ('uniprot', 'uniprot', 'genesymbol'),
            ('uniprot', 'uniprot', 'entrez', 'entrez_id'),
        )
    """

    for t in translations:
        col, src, tgt = t[0], t[1], t[2]
        new_col = t[3] if len(t) > 3 else None
        df = translate_column(
            df,
            col,
            src,
            tgt,
            ncbi_tax_id=ncbi_tax_id,
            new_column=new_col,
            keep_untranslated=keep_untranslated,
            expand=expand,
            raw=raw,
            backend=backend,
            full_uniprot=full_uniprot,
        )

    return df


def identify(
    identifiers: list[str],
    ncbi_tax_id: int = 9606,
) -> dict[str, list[dict]]:
    """Identify the type of given identifiers.

    Requires database mode (PostgreSQL).
    """
    from sqlalchemy.orm import Session

    from omnipath_utils.db._query import identify_ids
    from omnipath_utils.db._connection import get_engine

    engine = get_engine()
    with Session(engine) as session:
        return identify_ids(session, identifiers, ncbi_tax_id)


def all_mappings(
    identifiers: list[str],
    id_type: str,
    ncbi_tax_id: int = 9606,
) -> dict[str, dict[str, list[str]]]:
    """Get all known mappings for identifiers across all target types.

    Requires database mode (PostgreSQL).
    """
    from sqlalchemy.orm import Session

    from omnipath_utils.db._query import get_all_mappings
    from omnipath_utils.db._connection import get_engine

    engine = get_engine()
    with Session(engine) as session:
        return get_all_mappings(session, identifiers, id_type, ncbi_tax_id)
