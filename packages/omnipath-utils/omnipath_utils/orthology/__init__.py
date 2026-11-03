"""Orthology -- cross-species gene translation.

Translates identifiers between organisms using orthologous gene pairs
from OMA, Ensembl Compara, HomoloGene, and HCOP.

Example::

    from omnipath_utils.orthology import translate, translate_column

    # Human TP53 -> mouse orthologs
    translate(['TP53'], source=9606, target=10090)
    # {'TP53': {'Trp53'}}

    # Translate a DataFrame column (pandas, polars, or pyarrow)
    translate_column(df, 'gene', source=9606, target=10090)
"""

from __future__ import annotations

import narwhals as nw

from omnipath_utils._constants import DEFAULT_ORGANISM


def translate(
    identifiers: list[str] | set[str],
    source: int = 9606,
    target: int = 10090,
    id_type: str = 'genesymbol',
    only_swissprot: bool = True,
    resource: str | None = None,
    min_sources: int = 1,
    rel_type: set[str] | None = None,
    raw: bool = False,
) -> dict[str, set[str]]:
    """Translate identifiers to orthologs in another organism.

    Args:
        identifiers: Source identifiers.
        source: Source organism NCBI Taxonomy ID (default: 9606 human).
        target: Target organism NCBI Taxonomy ID (default: 10090 mouse).
        id_type: Identifier type (default: 'genesymbol').
        only_swissprot: Prefer SwissProt for UniProt IDs.
        resource: Force specific resource ('oma', 'ensembl', 'homologene', 'hcop').
        min_sources: For HCOP, minimum number of supporting databases.
        rel_type: Filter by relationship type (e.g. {'1:1'} for one-to-one).
        raw: Skip post-processing (ID translation, SwissProt preference).

    Returns:
        Dict mapping source IDs to sets of ortholog IDs.
    """
    from omnipath_utils.orthology._manager import OrthologyManager

    mgr = OrthologyManager.get()
    return mgr.translate(
        identifiers=list(identifiers),
        source=source,
        target=target,
        id_type=id_type,
        only_swissprot=only_swissprot,
        resource=resource,
        min_sources=min_sources,
        rel_type=rel_type,
        raw=raw,
    )


def translate_column(
    df,
    column: str,
    source: int = 9606,
    target: int = 10090,
    id_type: str = 'genesymbol',
    new_column: str | None = None,
    keep_untranslated: bool = True,
    expand: bool = True,
    only_swissprot: bool = True,
    resource: str | None = None,
    min_sources: int = 1,
):
    """Translate a DataFrame column to orthologs.

    Accepts any DataFrame backend supported by narwhals
    (pandas, polars, or pyarrow).

    Args:
        df: DataFrame (pandas, polars, or pyarrow).
        column: Column with source identifiers.
        source: Source organism.
        target: Target organism.
        id_type: Identifier type.
        new_column: Output column name (default: '{column}_{target}').
        keep_untranslated: Keep rows with no ortholog.
        expand: Expand one-to-many to multiple rows.
        only_swissprot: Prefer SwissProt for UniProt IDs.
        resource: Force specific resource.
        min_sources: For HCOP, minimum supporting databases.

    Returns:
        DataFrame with ortholog column added, in the same
        backend format as the input.
    """

    new_col = new_column or f'{column}_{target}'

    nw_df = nw.from_native(df, eager_only=True)
    native_ns = nw.get_native_namespace(nw_df)

    # Get unique source IDs
    unique_ids = nw_df[column].drop_nulls().unique().to_list()
    source_ids = [str(v) for v in unique_ids]

    trans = translate(
        source_ids,
        source=source,
        target=target,
        id_type=id_type,
        only_swissprot=only_swissprot,
        resource=resource,
        min_sources=min_sources,
    )

    if expand:
        # Build mapping rows
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

        # Handle source values not in trans
        all_sources = set(source_ids)
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


def get_table(
    source: int = 9606,
    target: int = 10090,
    id_type: str = 'genesymbol',
    resource: str | None = None,
    min_sources: int = 1,
) -> dict[str, set[str]]:
    """Get the full orthology table for a pair of organisms."""
    from omnipath_utils.orthology._manager import OrthologyManager

    return OrthologyManager.get().get_table(
        source=source,
        target=target,
        id_type=id_type,
        resource=resource,
        min_sources=min_sources,
    )
