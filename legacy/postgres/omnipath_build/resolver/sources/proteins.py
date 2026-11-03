"""Build protein identifier resolver mappings from UniProt inputs.

The protein resolver emits rows from evidence identifier namespaces to canonical
primary UniProt accessions. Reference mappings keep taxonomy so canonicalization
can require species agreement; secondary UniProt accessions are accepted only
when the primary accession has a single known taxonomy in the source snapshot.
Ambiguous key/taxonomy pairs are split into a separate audit table instead of
being used as resolver candidates.
"""

from __future__ import annotations

import os
from pathlib import Path
from collections.abc import Iterable

import polars as pl

from pypath.internals.cv_terms import (
    IdentifierNamespaceCv,
    cv_term_label_accession,
)
from omnipath_build.resolver.paths import (
    ensure_proteins_data_dir,
    activate_raw_download_data_dir,
)
from omnipath_build.resolver.identifier_types import (
    IDENTIFIER_TYPE_SCHEMA,
    identifier_type_id,
    identifier_type_rows,
)

PROTEIN_IDENTIFIER_LOOKUP_SCHEMA: dict[str, pl.DataType] = {
    'key_identifier_type_id': pl.UInt32,
    'key_value': pl.Utf8,
    'taxonomy_id': pl.Utf8,
    'canonical_identifier_type_id': pl.UInt32,
    'canonical_identifier': pl.Utf8,
}

UNIPROT_TYPE = cv_term_label_accession(IdentifierNamespaceCv.UNIPROT)
# Gene-anchored model: the canonical identifier is the gene anchor
# = NCBI Gene (Entrez); the protein/isoform is record-level state, not identity.
ENTREZ_TYPE = cv_term_label_accession(IdentifierNamespaceCv.ENTREZ)
# omnipath_utils.resolver_gene source_type slug -> build CV id-type.
RESOLVER_GENE_SLUG_TO_KEY = {
    'genesymbol': cv_term_label_accession(IdentifierNamespaceCv.GENE_NAME_PRIMARY),
    'entrez': ENTREZ_TYPE,
    'ensg': cv_term_label_accession(IdentifierNamespaceCv.ENSEMBL),
    'ensp': cv_term_label_accession(IdentifierNamespaceCv.ENSEMBL),
    'uniprot': UNIPROT_TYPE,
}
# Default in-scope organisms when no --taxonomy-id is given (full long-tail
# coverage = a follow-up materialisation of resolver_gene). Covers the EGFR
# ortholog benchmark + the main model/agricultural organisms.
DEFAULT_ORGANISMS = (
    9606, 10090, 10116, 7955, 7227, 6239, 3702, 4932, 559292, 83333,
    9913, 9823, 9615, 9031, 9598, 9544, 8364, 9796, 9940, 9986,
)
KEY_TYPE_ALIASES = {
    'MI:1097:Uniprot': UNIPROT_TYPE,
    'MI:0476:Ensembl': cv_term_label_accession(IdentifierNamespaceCv.ENSEMBL),
    'MI:0477:Entrez': cv_term_label_accession(IdentifierNamespaceCv.ENTREZ),
    'MI:1095:HGNC': cv_term_label_accession(IdentifierNamespaceCv.HGNC),
    'OM:0200:Gene Name Primary': cv_term_label_accession(
        IdentifierNamespaceCv.GENE_NAME_PRIMARY
    ),
    'OM:0201:Gene Name Synonym': cv_term_label_accession(
        IdentifierNamespaceCv.GENE_NAME_SYNONYM
    ),
    'OM:0221:Uniprot Entry Name': cv_term_label_accession(
        IdentifierNamespaceCv.UNIPROT_ENTRY_NAME
    ),
}
PROTEIN_IDENTIFIER_LOOKUP_OUTPUT_FILENAME = 'protein_identifier_lookup.parquet'
PROTEIN_IDENTIFIER_LOOKUP_AMBIGUOUS_OUTPUT_FILENAME = (
    'protein_identifier_lookup_ambiguous.parquet'
)
IDENTIFIER_TYPE_OUTPUT_FILENAME = 'identifier_type.parquet'
GENE_PROTEIN_REPRESENTATIVE_OUTPUT_FILENAME = (
    'gene_protein_representative.parquet'
)
GENE_PROTEIN_REPRESENTATIVE_SCHEMA: dict[str, pl.DataType] = {
    'taxonomy_id': pl.Utf8,
    'canonical_identifier': pl.Utf8,  # the gene's Entrez id (gene-entity key)
    'representative_uniprot': pl.Utf8,
    'is_reviewed': pl.Boolean,
    'uniprot_all': pl.List(pl.Utf8),
}


def _utils_pg_url() -> str:
    url = os.environ.get('OMNIPATH_BUILD_UTILS_PG_URL')
    if not url:
        raise RuntimeError(
            'OMNIPATH_BUILD_UTILS_PG_URL is not set; the gene-anchored protein '
            'resolver reads omnipath_utils.resolver_gene from the omnipath-utils '
            'Postgres.'
        )
    return url


def _protein_identifier_rows(
    taxonomy_ids: Iterable[int | str] | None = None,
) -> Iterable[dict]:
    """Yield gene-anchored resolver rows from ``omnipath_utils.resolver_gene``.

    Each row maps a source identifier (gene symbol / Entrez / Ensembl / UniProt)
    to its **NCBI Gene (Entrez) anchor** for its organism — the canonical
    collapsing identity. The asserted protein/isoform becomes record-level state
    downstream, not here. Read directly from the utils Postgres, per-taxon so the
    ``resolver_gene`` taxon filter pushes into the covering index.
    ``primary_uniprot`` carries the Entrez gene id (kept name = downstream
    contract; it is the canonical anchor, not a UniProt).
    """
    import psycopg2

    taxa = [int(t) for t in (taxonomy_ids or DEFAULT_ORGANISMS)]
    conn = psycopg2.connect(_utils_pg_url())
    try:
        for tax in taxa:
            with conn.cursor() as cur:
                cur.execute(
                    'SELECT source_type, source_id, entrez '
                    'FROM omnipath_utils.resolver_gene '
                    'WHERE ncbi_tax_id = %s',
                    (tax,),
                )
                rows = cur.fetchall()
            for source_type, source_id, entrez in rows:
                key_type = RESOLVER_GENE_SLUG_TO_KEY.get(source_type)
                if not (key_type and source_id and entrez):
                    continue
                yield {
                    'key_type': key_type,
                    'key_value': str(source_id),
                    'taxonomy_id': str(tax),
                    'primary_uniprot': str(entrez),
                }
    finally:
        conn.close()


def _global_uniprot_resolver_enabled(flag: bool | None) -> bool:
    """Whether to fold in the taxon-agnostic UniProt/Entrez slice.

    Off by default (the ~24M-row global map would slow capped dev loops); enabled
    explicitly for full/showcase builds via ``OMNIPATH_BUILD_GLOBAL_UNIPROT_RESOLVER``
    or the ``materialize_proteins`` argument.
    """
    if flag is not None:
        return flag
    return os.environ.get('OMNIPATH_BUILD_GLOBAL_UNIPROT_RESOLVER', '') not in (
        '',
        '0',
        'false',
        'False',
    )


def _global_protein_lookup_df() -> pl.DataFrame:
    """Taxon-agnostic UniProt/Entrez -> Entrez gene rows.

    Reads ``omnipath_utils.resolver_gene_protein_global`` (the global UniProt->Entrez
    map, ~24M rows over ~50k taxa) so UniProt-referenced proteins that arrive with
    no taxonomy (Rhea/Brenda/TCDB/ChEMBL participants) still gene-anchor — the
    resolved entity inherits the organism derived from the AC. Read via DuckDB
    ATTACH so the 24M rows never materialise as Python objects; the result feeds the
    same ambiguity split as the per-taxon rows (an AC->>1 gene falls to the
    ambiguous lookup / multi-gene split, never a force-merge).
    """
    import duckdb

    url = _utils_pg_url()
    literal = "'" + url.replace("'", "''") + "'"
    con = duckdb.connect()
    try:
        con.execute('INSTALL postgres; LOAD postgres;')
        con.execute(f'ATTACH {literal} AS up (TYPE postgres, READ_ONLY)')
        df = con.execute(
            'SELECT source_type, source_id, '
            'CAST(ncbi_tax_id AS VARCHAR) AS taxonomy_id, '
            'CAST(entrez AS VARCHAR) AS entrez '
            'FROM up.omnipath_utils.resolver_gene_protein_global '
            'WHERE source_id IS NOT NULL AND entrez IS NOT NULL'
        ).pl()
    finally:
        con.close()

    uniprot_type_id = identifier_type_id(UNIPROT_TYPE)
    entrez_type_id = identifier_type_id(ENTREZ_TYPE)
    return df.select(
        pl.when(pl.col('source_type') == 'uniprot')
        .then(pl.lit(uniprot_type_id))
        .otherwise(pl.lit(entrez_type_id))
        .cast(pl.UInt32)
        .alias('key_identifier_type_id'),
        pl.col('source_id').cast(pl.Utf8).alias('key_value'),
        pl.col('taxonomy_id').cast(pl.Utf8).alias('taxonomy_id'),
        pl.lit(entrez_type_id)
        .cast(pl.UInt32)
        .alias('canonical_identifier_type_id'),
        pl.col('entrez').cast(pl.Utf8).alias('canonical_identifier'),
    )


def _single_taxonomy_id(values: set[str] | None) -> str | None:
    if not values or len(values) != 1:
        return None
    return next(iter(values))


def _gene_protein_representative_rows(
    taxonomy_ids: Iterable[int | str] | None = None,
) -> Iterable[dict]:
    """Yield per-gene representative UniProt rows.

    From ``omnipath_utils.resolver_protein`` (Entrez → UniProt, already
    SwissProt-preferred), grouped per gene: ``representative_uniprot`` = a reviewed
    (SwissProt) AC if any, else the chosen (sorted) AC; ``uniprot_all`` = every AC;
    ``is_reviewed`` = whether the representative is SwissProt. omnipath-build joins
    this to the gene entities (by Entrez) to fill ``gene_protein_representative`` —
    the no-state-join gene+UniProt output column.
    """
    import psycopg2

    taxa = [int(t) for t in (taxonomy_ids or DEFAULT_ORGANISMS)]
    conn = psycopg2.connect(_utils_pg_url())
    try:
        for tax in taxa:
            with conn.cursor() as cur:
                cur.execute(
                    'SELECT rp.source_id, rp.uniprot, '
                    '(sw.identifier IS NOT NULL) '
                    'FROM omnipath_utils.resolver_protein rp '
                    "LEFT JOIN omnipath_utils.reflist sw "
                    "  ON sw.list_name = 'swissprot' "
                    '  AND sw.identifier = rp.uniprot '
                    "WHERE rp.source_type = 'entrez' AND rp.ncbi_tax_id = %s",
                    (tax,),
                )
                rows = cur.fetchall()
            per_gene: dict[str, dict[str, set]] = {}
            for entrez, uniprot, is_reviewed in rows:
                if not entrez or not uniprot:
                    continue
                g = per_gene.setdefault(
                    str(entrez), {'all': set(), 'reviewed': set()}
                )
                g['all'].add(str(uniprot))
                if is_reviewed:
                    g['reviewed'].add(str(uniprot))
            for entrez, g in per_gene.items():
                uniprot_all = sorted(g['all'])
                reviewed = sorted(g['reviewed'])
                yield {
                    'taxonomy_id': str(tax),
                    'canonical_identifier': entrez,
                    'representative_uniprot': (
                        reviewed[0] if reviewed else uniprot_all[0]
                    ),
                    'is_reviewed': bool(reviewed),
                    'uniprot_all': uniprot_all,
                }
    finally:
        conn.close()


def build_protein_identifier_lookup(
    taxonomy_ids: Iterable[int | str] | None = None,
) -> pl.DataFrame:
    """Return the non-ambiguous protein identifier lookup as a dataframe."""

    activate_raw_download_data_dir()
    return _split_protein_identifier_lookup(
        _protein_identifier_rows(taxonomy_ids=taxonomy_ids)
    )[0]


def materialize_proteins(
    output_dir: str | Path | None = None,
    taxonomy_ids: Iterable[int | str] | None = None,
    skip_existing: bool = True,
    include_global_uniprot: bool | None = None,
) -> dict[str, int]:
    """Write protein resolver parquet files and return output row counts.

    When ``include_global_uniprot`` (env ``OMNIPATH_BUILD_GLOBAL_UNIPROT_RESOLVER``)
    is on, the per-taxon ``DEFAULT_ORGANISMS`` rows are augmented with the
    taxon-agnostic global UniProt/Entrez slice so UniProt-referenced proteins
    with no taxonomy still gene-anchor.
    """

    output_dir = (
        Path(output_dir)
        if output_dir is not None
        else ensure_proteins_data_dir()
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    lookup_path = output_dir / PROTEIN_IDENTIFIER_LOOKUP_OUTPUT_FILENAME
    ambiguous_path = (
        output_dir / PROTEIN_IDENTIFIER_LOOKUP_AMBIGUOUS_OUTPUT_FILENAME
    )
    identifier_type_path = output_dir / IDENTIFIER_TYPE_OUTPUT_FILENAME
    gpr_path = output_dir / GENE_PROTEIN_REPRESENTATIVE_OUTPUT_FILENAME
    if (
        skip_existing
        and taxonomy_ids is None
        and lookup_path.exists()
        and ambiguous_path.exists()
        and identifier_type_path.exists()
        and gpr_path.exists()
    ):
        print(
            f'[resolver] skip source=uniprot existing_dir={output_dir}',
            flush=True,
        )
        return {
            'protein_identifier_lookup_rows': _parquet_row_count(lookup_path),
            'protein_identifier_lookup_ambiguous_rows': _parquet_row_count(
                ambiguous_path
            ),
            'identifier_type_rows': _parquet_row_count(identifier_type_path),
            'gene_protein_representative_rows': _parquet_row_count(gpr_path),
        }

    activate_raw_download_data_dir()
    extra_lookup_df = None
    if _global_uniprot_resolver_enabled(include_global_uniprot):
        print('[resolver] uniprot global taxon-agnostic slice ON', flush=True)
        extra_lookup_df = _global_protein_lookup_df()
        print(
            f'[resolver] uniprot global slice rows={extra_lookup_df.height}',
            flush=True,
        )
    lookup, ambiguous, identifier_types = _split_protein_identifier_lookup(
        _protein_identifier_rows(taxonomy_ids=taxonomy_ids),
        extra_lookup_df=extra_lookup_df,
    )
    lookup.write_parquet(lookup_path)
    ambiguous.write_parquet(ambiguous_path)
    identifier_types.write_parquet(identifier_type_path)

    gpr = pl.DataFrame(
        list(_gene_protein_representative_rows(taxonomy_ids=taxonomy_ids)),
        schema=GENE_PROTEIN_REPRESENTATIVE_SCHEMA,
    )
    gpr.write_parquet(gpr_path)

    return {
        'protein_identifier_lookup_rows': lookup.height,
        'protein_identifier_lookup_ambiguous_rows': ambiguous.height,
        'identifier_type_rows': identifier_types.height,
        'gene_protein_representative_rows': gpr.height,
    }


def _parquet_row_count(path: Path) -> int:
    return pl.scan_parquet(path).select(pl.len()).collect().item()


def _split_protein_identifier_lookup(
    rows: Iterable[dict],
    extra_lookup_df: pl.DataFrame | None = None,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    normalized_rows: list[dict[str, object]] = []
    type_names = {ENTREZ_TYPE}
    for row in rows:
        key_type = row.get('key_type')
        if key_type is None:
            continue
        key_type = str(key_type)
        type_names.add(key_type)
        normalized_rows.append(
            {
                'key_identifier_type_id': identifier_type_id(key_type),
                'key_value': row.get('key_value'),
                'taxonomy_id': row.get('taxonomy_id'),
                'canonical_identifier_type_id': identifier_type_id(
                    ENTREZ_TYPE
                ),
                'canonical_identifier': row.get('primary_uniprot'),
            }
        )

    if extra_lookup_df is not None and extra_lookup_df.height:
        # the global slice keys on UniProt + Entrez; make sure both id types
        # reach the identifier_type parquet even if no per-taxon rows did.
        type_names.update({UNIPROT_TYPE, ENTREZ_TYPE})

    identifier_types = pl.DataFrame(
        identifier_type_rows(type_names),
        schema=IDENTIFIER_TYPE_SCHEMA,
    )
    base_lookup = pl.DataFrame(
        normalized_rows, schema=PROTEIN_IDENTIFIER_LOOKUP_SCHEMA
    )
    if extra_lookup_df is not None and extra_lookup_df.height:
        base_lookup = pl.concat(
            [base_lookup, extra_lookup_df.select(base_lookup.columns)],
            how='vertical',
        )
    if base_lookup.is_empty():
        empty = pl.DataFrame(schema=PROTEIN_IDENTIFIER_LOOKUP_SCHEMA)
        return empty, empty, identifier_types

    return _split_lookup_frame(base_lookup, type_names)


def _split_lookup_frame(
    base_lookup: pl.DataFrame,
    type_names: set[str],
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    identifier_types = pl.DataFrame(
        identifier_type_rows(type_names),
        schema=IDENTIFIER_TYPE_SCHEMA,
    )
    if base_lookup.is_empty():
        empty = pl.DataFrame(schema=PROTEIN_IDENTIFIER_LOOKUP_SCHEMA)
        return empty, empty, identifier_types

    lookup = (
        base_lookup
        .filter(
            pl.col('key_value').is_not_null()
            & (pl.col('key_value') != '')
            & pl.col('canonical_identifier').is_not_null()
            & (pl.col('canonical_identifier') != '')
        )
        .unique()
    )
    ambiguous_keys = (
        lookup.group_by(
            [
                'key_identifier_type_id',
                'key_value',
                'taxonomy_id',
                'canonical_identifier_type_id',
            ]
        )
        .agg(
            pl.col('canonical_identifier')
            .n_unique()
            .alias('canonical_identifier_count')
        )
        .filter(pl.col('canonical_identifier_count') > 1)
        .select(
            [
                'key_identifier_type_id',
                'key_value',
                'taxonomy_id',
                'canonical_identifier_type_id',
            ]
        )
    )
    if ambiguous_keys.is_empty():
        return (
            lookup,
            pl.DataFrame(schema=PROTEIN_IDENTIFIER_LOOKUP_SCHEMA),
            identifier_types,
        )

    join_keys = [
        'key_identifier_type_id',
        'key_value',
        'taxonomy_id',
        'canonical_identifier_type_id',
    ]
    ambiguous = lookup.join(ambiguous_keys, on=join_keys, how='semi')
    unambiguous = lookup.join(ambiguous_keys, on=join_keys, how='anti')
    return unambiguous, ambiguous, identifier_types
