"""Low-level DuckDB/PostgreSQL load helpers for the COPY pipeline."""

from __future__ import annotations

import os
import logging
from pathlib import Path
import tempfile
from itertools import islice
from collections.abc import Iterable

import duckdb
import pyarrow as pa
import psycopg2
from psycopg2 import sql
from psycopg2.extras import execute_values

from omnipath_build.cv_terms import (
    CV_TERM_ID_TYPE,
    GENE_ENTITY_TYPE,
    CV_TERM_ENTITY_TYPE,
    CHEMICAL_ENTITY_TYPE,
)
from pypath.internals.cv_terms import (
    BiologicalRoleCv,
    EntityTypeCv,
    IdentifierNamespaceCv,
    cv_term_label_accession,
)
from omnipath_build.ingest.common import unwrap_record
from pypath.internals.silver_schema import Entity
from omnipath_build.evidence_projector import (
    ProjectionStats,
    EvidenceProjectorBase,
    _MutableProjectionStats,
)
from omnipath_build.chemical_fallback import (
    build_chemical_anchor_map,
    build_chemical_fallback_resolution,
    chemical_fallback_fires_sql,
)
from omnipath_build.multigene_split import explode_multi_gene_protein_mentions
from omnipath_build.resolver.identifier_types import (
    UNRESOLVED_ID_TYPE,
    IDENTIFIER_TYPE_NAMES,
    COMPLEX_MEMBER_HASH_ID_TYPE,
    REACTION_MEMBER_HASH_ID_TYPE,
    identifier_type_id,
    identifier_type_rows,
)

_log = logging.getLogger(__name__)

PROTEIN_ENTITY_TYPE = 'Protein:MI:0326'
COMPLEX_ENTITY_TYPE = 'Complex:MI:0314'
REACTION_ENTITY_TYPE = 'Reaction:OM:0015'
MIRNA_ENTITY_TYPE = cv_term_label_accession(EntityTypeCv.MIRNA)
PATHWAY_ENTITY_TYPE = cv_term_label_accession(EntityTypeCv.PATHWAY)
REACTANT_ROLE_TERMS = (
    cv_term_label_accession(BiologicalRoleCv.REACTANT),
    str(BiologicalRoleCv.REACTANT),
    cv_term_label_accession(BiologicalRoleCv.SUBSTRATE),
    str(BiologicalRoleCv.SUBSTRATE),
)
PRODUCT_ROLE_TERMS = (
    cv_term_label_accession(BiologicalRoleCv.PRODUCT),
    str(BiologicalRoleCv.PRODUCT),
)
PROTEIN_TAXONOMY_OPTIONAL_IDENTIFIER_TYPES = (
    cv_term_label_accession(IdentifierNamespaceCv.UNIPROT),
    cv_term_label_accession(IdentifierNamespaceCv.ENSEMBL),
    cv_term_label_accession(IdentifierNamespaceCv.ENTREZ),
    cv_term_label_accession(IdentifierNamespaceCv.HGNC),
    cv_term_label_accession(IdentifierNamespaceCv.UNIPROT_ENTRY_NAME),
)
RESOLVER_ALIAS_EXPANSION_EXCLUDED_IDENTIFIER_TYPES = (
    cv_term_label_accession(IdentifierNamespaceCv.ENSEMBL),
)
UNIPROT_TYPE = cv_term_label_accession(IdentifierNamespaceCv.UNIPROT)
GENE_NAME_PRIMARY_TYPE = cv_term_label_accession(
    IdentifierNamespaceCv.GENE_NAME_PRIMARY
)
UNIPROT_ENTRY_NAME_TYPE = cv_term_label_accession(
    IdentifierNamespaceCv.UNIPROT_ENTRY_NAME
)
ENSEMBL_TYPE = cv_term_label_accession(IdentifierNamespaceCv.ENSEMBL)
ENTREZ_TYPE = cv_term_label_accession(IdentifierNamespaceCv.ENTREZ)
RESOLVER_PROTEIN_SLUG_TO_IDENTIFIER_TYPE = {
    'genesymbol': GENE_NAME_PRIMARY_TYPE,
    'entrez': ENTREZ_TYPE,
    'ensg': ENSEMBL_TYPE,
    'ensp': ENSEMBL_TYPE,
    'uniprot': UNIPROT_TYPE,
    'uniprot_entry': UNIPROT_ENTRY_NAME_TYPE,
}
# Deterministic best-identifier priority for an entity that could not be resolved
# to a canonical accession. Such an entity is keyed by its single highest-priority
# raw identifier (rather than an opaque hash of its whole identifier bag), so the
# key is human-readable, stable across builds, and repeated mentions of the same
# identifier collapse into one entity. Ties break deterministically by identifier
# string. Identifier types not listed here fall to a common lowest rank (still
# deterministic via the identifier-type-id and identifier tie-breakers).
UNRESOLVED_KEY_IDENTIFIER_PRIORITY = (
    UNIPROT_TYPE,
    ENSEMBL_TYPE,
    ENTREZ_TYPE,
    GENE_NAME_PRIMARY_TYPE,
    UNIPROT_ENTRY_NAME_TYPE,
)
# The UniProtKB accession format (both 6- and 10-character forms). A protein
# mention that carries a well-formed UniProt accession names a real gene product,
# so even when the identifier-resolution database has no gene for it we can type
# it as a gene keyed by that accession rather than leaving a bare protein node.
# An optional ``-<n>`` isoform suffix is stripped before matching.
UNIPROT_AC_REGEX = (
    '^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$'
)
# protein molecular_type_id (matches the entity_evidence_resolution CASE seed);
# the per-record asserted UniProt form is a protein state.
PROTEIN_MOLECULAR_TYPE_ID = 2
STANDARD_INCHI_KEY_TYPE = cv_term_label_accession(
    IdentifierNamespaceCv.STANDARD_INCHI_KEY
)
CHEBI_TYPE = cv_term_label_accession(IdentifierNamespaceCv.CHEBI)
CHEMBL_COMPOUND_TYPE = cv_term_label_accession(
    IdentifierNamespaceCv.CHEMBL_COMPOUND
)
HMDB_TYPE = cv_term_label_accession(IdentifierNamespaceCv.HMDB)
LIPIDMAPS_TYPE = cv_term_label_accession(IdentifierNamespaceCv.LIPIDMAPS)
PUBCHEM_COMPOUND_TYPE = cv_term_label_accession(
    IdentifierNamespaceCv.PUBCHEM_COMPOUND
)
SWISSLIPIDS_TYPE = cv_term_label_accession(IdentifierNamespaceCv.SWISSLIPIDS)
KEGG_COMPOUND_TYPE = cv_term_label_accession(IdentifierNamespaceCv.KEGG_COMPOUND)
DRUGBANK_TYPE = cv_term_label_accession(IdentifierNamespaceCv.DRUGBANK)
REACTOME_STABLE_ID_TYPE = cv_term_label_accession(
    IdentifierNamespaceCv.REACTOME_STABLE_ID
)
LIPID_NAME_TYPE = cv_term_label_accession(IdentifierNamespaceCv.LIPID_NAME)
# Name-axis identifier types that may carry a lipid shorthand (spec 011 T116) --
# matches omnipath-metabo's own _lipid_layer.py _NAME_TYPES exactly, so both
# layers consider the same candidate identifiers lipid-shaped.
LIPID_NAME_SOURCE_TYPES = (
    'Name:OM:0202',
    'Synonym:OM:0203',
    'Abbreviated Name:OM:0208',
    'Iupac Traditional Name:OM:0211',
)
RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE = {
    'pubchem': PUBCHEM_COMPOUND_TYPE,
    'chembl': CHEMBL_COMPOUND_TYPE,
    'chebi': CHEBI_TYPE,
    'hmdb': HMDB_TYPE,
    'lipidmaps': LIPIDMAPS_TYPE,
    'swisslipids': SWISSLIPIDS_TYPE,
    'kegg': KEGG_COMPOUND_TYPE,
    # T092: every namespace with a structure route through
    # omnipath_utils.resolver_chemical, not just the original seven.
    # DrugBank was already loaded (id_mapping, via UniChem's pubchem
    # bridge -- 15,293 rows, research.md R7) and simply never read here.
    # Reactome is this cycle's new chebi-bridged load (T089). RefMet
    # (T088) is NOT added here: its cross-reference lives in
    # id_mapping_long (the name-axis table), which resolver_chemical does
    # not read at all -- RefMet translates through the name axis
    # (translate_ids), a separate, already-adequate mechanism, not
    # through this structure-route resolver.
    'drugbank': DRUGBANK_TYPE,
    'reactome': REACTOME_STABLE_ID_TYPE,
}
WIKIPATHWAYS_ID_TYPE = cv_term_label_accession(IdentifierNamespaceCv.WIKIPATHWAYS)


def _sql_literal(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def _duckdb_identifier(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'


def _duckdb_pg_table(schema: str, table: str) -> str:
    return f'pg.{_duckdb_identifier(schema)}.{_duckdb_identifier(table)}'


def _resolver_component_dir(
    resolver_dir: Path,
    *,
    component: str,
    lookup_filename: str,
) -> Path:
    """Resolve the directory that contains one resolver component's parquet outputs."""

    candidates = (
        resolver_dir / component,
        resolver_dir,
        resolver_dir.parent / component,
    )
    for candidate in candidates:
        if (candidate / lookup_filename).exists():
            return candidate
    raise FileNotFoundError(
        f'Could not find {lookup_filename!r} for resolver component {component!r} '
        f'under {resolver_dir}.'
    )


def _chemical_resolver_component_dir(resolver_dir: Path) -> Path:
    """Resolve the directory that contains partitioned chemical lookup files."""

    candidates = (
        resolver_dir / 'chemicals',
        resolver_dir,
        resolver_dir.parent / 'chemicals',
    )
    for candidate in candidates:
        if (
            (candidate / 'lookup').is_dir()
            and _has_parquet_files(candidate / 'lookup')
            and (candidate / 'identifier_type.parquet').exists()
        ):
            return candidate
    raise FileNotFoundError(
        'Could not find partitioned chemical resolver outputs under '
        f'{resolver_dir}. Expected chemicals/lookup/*.parquet and '
        'chemicals/identifier_type.parquet.'
    )


def _parquet_glob(path: Path) -> str:
    return str(path).replace("'", "''")


def _has_parquet_files(path: Path) -> bool:
    return path.is_dir() and any(path.glob('*.parquet'))


def _has_nested_parquet_files(path: Path) -> bool:
    return path.is_dir() and any(path.rglob('*.parquet'))


def _combined_resolver_parquet_glob(resolver_dir: Path) -> str | None:
    """Return a glob for the combined gene/protein resolver parquet dataset."""

    configured = os.environ.get('OMNIPATH_BUILD_COMBINED_RESOLVER_PARQUET')
    candidates: list[Path] = []
    if configured:
        candidates.append(Path(configured).expanduser())
    candidates.extend(
        (
            resolver_dir / 'resolver_gene_protein_combined',
            resolver_dir / 'resolver_gene_protein_combined.parquet',
        )
    )
    for candidate in candidates:
        if _has_nested_parquet_files(candidate):
            return _parquet_glob(candidate / '**' / '*.parquet')
    if configured:
        raise FileNotFoundError(
            'OMNIPATH_BUILD_COMBINED_RESOLVER_PARQUET does not contain '
            f'any parquet files: {configured}'
        )
    return None


def _combined_resolver_duckdb_path(resolver_dir: Path) -> Path | None:
    """Return a persistent DuckDB path for the combined gene/protein resolver."""

    configured = os.environ.get('OMNIPATH_BUILD_COMBINED_RESOLVER_DUCKDB')
    candidates: list[Path] = []
    if configured:
        candidates.append(Path(configured).expanduser())
    candidates.append(resolver_dir / 'resolver_gene_protein_combined.duckdb')
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    if configured:
        raise FileNotFoundError(
            'OMNIPATH_BUILD_COMBINED_RESOLVER_DUCKDB does not exist: '
            f'{configured}'
        )
    return None


def _attach_combined_resolver_duckdb(
    con: duckdb.DuckDBPyConnection,
    path: Path,
    *,
    alias: str = 'combined_resolver_db',
) -> None:
    attached = any(
        row[1] == alias for row in con.execute('PRAGMA database_list').fetchall()
    )
    if not attached:
        con.execute(
            f'ATTACH {_sql_literal(str(path))} AS {_duckdb_identifier(alias)} '
            '(READ_ONLY)'
        )


def _duckdb_load_postgres_extension(con: duckdb.DuckDBPyConnection) -> None:
    con.execute('INSTALL postgres; LOAD postgres;')


def _duckdb_attach_utils_postgres(
    con: duckdb.DuckDBPyConnection,
    *,
    alias: str = 'utils_pg',
) -> bool:
    url = os.environ.get('OMNIPATH_BUILD_UTILS_PG_URL')
    if not url:
        return False
    attached = any(
        row[1] == alias for row in con.execute('PRAGMA database_list').fetchall()
    )
    if not attached:
        _duckdb_load_postgres_extension(con)
        con.execute(
            f'ATTACH {_sql_literal(url)} AS {_duckdb_identifier(alias)} '
            '(TYPE postgres, READ_ONLY)'
        )
    return True


def _attached_utils_relation_exists(
    con: duckdb.DuckDBPyConnection,
    relation_name: str,
) -> bool:
    return bool(
        con.execute(
            """
            SELECT count(*) > 0
            FROM utils_pg.pg_catalog.pg_class c
            JOIN utils_pg.pg_catalog.pg_namespace n
              ON n.oid = c.relnamespace
            WHERE n.nspname = 'omnipath_utils'
              AND c.relname = ?
            """,
            [relation_name],
        ).fetchone()[0]
    )


def _create_empty_protein_uniprot_fallback_view(
    con: duckdb.DuckDBPyConnection,
) -> None:
    con.execute(
        f"""
        CREATE VIEW protein_uniprot_fallback_lookup AS
        SELECT
          {_sql_literal(PROTEIN_ENTITY_TYPE)} AS entity_type,
          NULL::BIGINT AS key_identifier_type_id,
          NULL::VARCHAR AS key_value,
          NULL::VARCHAR AS taxonomy_id,
          NULL::BIGINT AS canonical_identifier_type_id,
          NULL::VARCHAR AS canonical_identifier
        WHERE false
        """
    )


def _create_protein_uniprot_fallback_view(
    con: duckdb.DuckDBPyConnection,
) -> None:
    """Expose utils ``resolver_protein`` as the protein-only fallback lookup."""

    if not _duckdb_attach_utils_postgres(con):
        con.execute(
            f"""
            CREATE VIEW protein_uniprot_fallback_lookup AS
            SELECT
              {_sql_literal(PROTEIN_ENTITY_TYPE)} AS entity_type,
              NULL::BIGINT AS key_identifier_type_id,
              NULL::VARCHAR AS key_value,
              NULL::VARCHAR AS taxonomy_id,
              NULL::BIGINT AS canonical_identifier_type_id,
              NULL::VARCHAR AS canonical_identifier
            WHERE false
            """
        )
        return

    key_type_case = '\n'.join(
        f"              WHEN '{source_type}' THEN "
        f'{identifier_type_id(identifier_type)}'
        for source_type, identifier_type
        in RESOLVER_PROTEIN_SLUG_TO_IDENTIFIER_TYPE.items()
    )
    source_types_sql = ', '.join(
        _sql_literal(source_type)
        for source_type in RESOLVER_PROTEIN_SLUG_TO_IDENTIFIER_TYPE
    )
    con.execute(
        f"""
        CREATE VIEW protein_uniprot_fallback_lookup AS
        SELECT
          {_sql_literal(PROTEIN_ENTITY_TYPE)} AS entity_type,
          CASE rp.source_type
{key_type_case}
          END::BIGINT AS key_identifier_type_id,
          rp.source_id::VARCHAR AS key_value,
          rp.ncbi_tax_id::VARCHAR AS taxonomy_id,
          {identifier_type_id(UNIPROT_TYPE)}::BIGINT
            AS canonical_identifier_type_id,
          rp.uniprot::VARCHAR AS canonical_identifier
        FROM utils_pg.omnipath_utils.resolver_protein rp
        WHERE rp.source_type IN ({source_types_sql})
          AND rp.source_id IS NOT NULL
          AND rp.uniprot IS NOT NULL
        UNION ALL
        SELECT
          {_sql_literal(PROTEIN_ENTITY_TYPE)} AS entity_type,
          {identifier_type_id(UNIPROT_TYPE)}::BIGINT AS key_identifier_type_id,
          m.source_id::VARCHAR AS key_value,
          NULLIF(m.ncbi_tax_id, 0)::VARCHAR AS taxonomy_id,
          {identifier_type_id(UNIPROT_TYPE)}::BIGINT
            AS canonical_identifier_type_id,
          m.target_id::VARCHAR AS canonical_identifier
        FROM utils_pg.omnipath_utils.id_mapping m
        JOIN utils_pg.omnipath_utils.id_type st
          ON m.source_type_id = st.id
         AND st.name = 'uniprot-sec'
        JOIN utils_pg.omnipath_utils.id_type tt
          ON m.target_type_id = tt.id
         AND tt.name = 'uniprot-pri'
        WHERE m.source_id IS NOT NULL
          AND m.target_id IS NOT NULL
        """
    )


def _source_type_case_sql(
    mapping: dict[str, str],
    *,
    expression: str,
    indent: str = '          ',
) -> str:
    lines = [f'CASE {expression}']
    lines.extend(
        f"{indent}WHEN {_sql_literal(source_type)} THEN "
        f'{identifier_type_id(identifier_type)}'
        for source_type, identifier_type in mapping.items()
    )
    lines.append(f'{indent}END')
    return '\n'.join(lines)


def _source_type_in_sql(mapping: dict[str, str]) -> str:
    return ', '.join(_sql_literal(source_type) for source_type in mapping)


def _static_identifier_type_select_sql(
    names: set[str] | None = None,
) -> str:
    rows = ',\n'.join(
        f'({row["identifier_type_id"]}, {_sql_literal(row["name"])})'
        for row in identifier_type_rows(names)
    )
    return f"""
        SELECT *
        FROM (VALUES
        {rows}
        ) AS t(identifier_type_id, name)
        """


def _mirna_resolver_sql_parts(
    resolver_dir: Path,
) -> tuple[str, str, str]:
    try:
        mirna_dir = _resolver_component_dir(
            resolver_dir,
            component='mirna',
            lookup_filename='mirna_identifier_lookup.parquet',
        )
    except FileNotFoundError:
        return '', '', ''

    mirna_lookup_path = mirna_dir / 'mirna_identifier_lookup.parquet'
    mirna_type_path = mirna_dir / 'identifier_type.parquet'
    mirna_identifier_type_sql = f"""
        UNION
        SELECT *
        FROM read_parquet({_sql_literal(mirna_type_path)})
        """
    mirna_lookup_sql = f"""
        UNION ALL
        SELECT
          {_sql_literal(MIRNA_ENTITY_TYPE)} AS entity_type,
          key_identifier_type_id,
          key_value,
          NULL::VARCHAR AS taxonomy_id,
          canonical_identifier_type_id,
          canonical_identifier
        FROM read_parquet({_sql_literal(mirna_lookup_path)})
        WHERE key_value IS NOT NULL
          AND canonical_identifier IS NOT NULL
        """
    mirna_canonical_sql = f"""
        UNION ALL
        SELECT
          row_number() OVER (
            ORDER BY canonical_identifier_type_id, canonical_identifier
          )::BIGINT AS resolver_entity_id,
          {_sql_literal(MIRNA_ENTITY_TYPE)} AS entity_type,
          NULL::VARCHAR AS taxonomy_id,
          canonical_identifier_type_id,
          canonical_identifier,
          list_distinct(list(key_identifier_type_id)) AS key_identifier_type_ids,
          count(*)::BIGINT AS lookup_rows
        FROM read_parquet({_sql_literal(mirna_lookup_path)})
        WHERE canonical_identifier IS NOT NULL
        GROUP BY
          canonical_identifier_type_id,
          canonical_identifier
        """
    return mirna_identifier_type_sql, mirna_lookup_sql, mirna_canonical_sql


def _create_live_utils_resolver_views(
    con: duckdb.DuckDBPyConnection,
    *,
    resolver_dir: Path,
) -> None:
    """Expose resolver inputs directly from attached omnipath-utils Postgres."""

    mirna_identifier_type_sql, mirna_lookup_sql, mirna_canonical_sql = (
        _mirna_resolver_sql_parts(resolver_dir)
    )
    gene_key_case = _source_type_case_sql(
        RESOLVER_PROTEIN_SLUG_TO_IDENTIFIER_TYPE,
        expression='rg.source_type',
    )
    gene_source_types = _source_type_in_sql(
        RESOLVER_PROTEIN_SLUG_TO_IDENTIFIER_TYPE
    )
    combined_key_case = _source_type_case_sql(
        RESOLVER_PROTEIN_SLUG_TO_IDENTIFIER_TYPE,
        expression='rgp.source_type',
    )
    chemical_key_case = _source_type_case_sql(
        RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE,
        expression='rc.source_type',
    )
    chemical_source_types = _source_type_in_sql(
        RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE
    )
    con.execute(
        f"""
        CREATE VIEW identifier_type AS
        {_static_identifier_type_select_sql()}
        {mirna_identifier_type_sql}
        """
    )
    _create_duckdb_identifier_type_all_view(con)
    combined_duckdb_path = _combined_resolver_duckdb_path(resolver_dir)
    combined_parquet_glob = _combined_resolver_parquet_glob(resolver_dir)
    if combined_duckdb_path is not None:
        _attach_combined_resolver_duckdb(con, combined_duckdb_path)
        con.execute(
            f"""
            CREATE VIEW resolver_lookup AS
            SELECT
              entity_type,
              key_identifier_type_id,
              key_value,
              taxonomy_id,
              canonical_identifier_type_id,
              canonical_identifier
            FROM combined_resolver_db.resolver_lookup_gene_protein
            WHERE key_value IS NOT NULL
              AND canonical_identifier IS NOT NULL
            UNION ALL
            SELECT
              {_sql_literal(CHEMICAL_ENTITY_TYPE)} AS entity_type,
              ({chemical_key_case})::BIGINT AS key_identifier_type_id,
              rc.source_id::VARCHAR AS key_value,
              NULL::VARCHAR AS taxonomy_id,
              {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}::BIGINT
                AS canonical_identifier_type_id,
              rc.inchikey::VARCHAR AS canonical_identifier
            FROM utils_pg.omnipath_utils.resolver_chemical rc
            WHERE rc.source_type IN ({chemical_source_types})
              AND rc.source_id IS NOT NULL
              AND rc.inchikey IS NOT NULL
            UNION ALL
            SELECT DISTINCT
              {_sql_literal(CHEMICAL_ENTITY_TYPE)} AS entity_type,
              {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}::BIGINT
                AS key_identifier_type_id,
              rc.inchikey::VARCHAR AS key_value,
              NULL::VARCHAR AS taxonomy_id,
              {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}::BIGINT
                AS canonical_identifier_type_id,
              rc.inchikey::VARCHAR AS canonical_identifier
            FROM utils_pg.omnipath_utils.resolver_chemical rc
            WHERE rc.inchikey IS NOT NULL
              AND rc.inchikey <> ''
            {mirna_lookup_sql}
            """
        )
        _create_empty_protein_uniprot_fallback_view(con)
        con.execute(
            f"""
            CREATE VIEW resolver_canonical_entity AS
            SELECT
              row_number() OVER (
                ORDER BY
                  entity_type,
                  taxonomy_id NULLS FIRST,
                  canonical_identifier_type_id,
                  canonical_identifier
              )::BIGINT AS resolver_entity_id,
              entity_type,
              taxonomy_id,
              canonical_identifier_type_id,
              canonical_identifier,
              list_distinct(list(key_identifier_type_id))
                AS key_identifier_type_ids,
              count(*)::BIGINT AS lookup_rows
            FROM combined_resolver_db.resolver_lookup_gene_protein
            WHERE canonical_identifier IS NOT NULL
            GROUP BY
              entity_type,
              taxonomy_id,
              canonical_identifier_type_id,
              canonical_identifier
            UNION ALL
            SELECT
              row_number() OVER (
                ORDER BY {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}, rc.inchikey
              )::BIGINT AS resolver_entity_id,
              {_sql_literal(CHEMICAL_ENTITY_TYPE)} AS entity_type,
              NULL::VARCHAR AS taxonomy_id,
              {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}::BIGINT
                AS canonical_identifier_type_id,
              rc.inchikey::VARCHAR AS canonical_identifier,
              list_distinct(list(({chemical_key_case})::BIGINT))
                AS key_identifier_type_ids,
              count(*)::BIGINT AS lookup_rows
            FROM utils_pg.omnipath_utils.resolver_chemical rc
            WHERE rc.source_type IN ({chemical_source_types})
              AND rc.inchikey IS NOT NULL
            GROUP BY rc.inchikey
            {mirna_canonical_sql}
            """
        )
    elif combined_parquet_glob is not None:
        con.execute(
            f"""
            CREATE VIEW resolver_lookup AS
            SELECT
              CASE rgp.canonical_type
                WHEN 'entrez' THEN {_sql_literal(GENE_ENTITY_TYPE)}
                WHEN 'uniprot' THEN {_sql_literal(PROTEIN_ENTITY_TYPE)}
              END AS entity_type,
              ({combined_key_case})::BIGINT AS key_identifier_type_id,
              rgp.source_id::VARCHAR AS key_value,
              NULLIF(rgp.ncbi_tax_id, 0)::VARCHAR AS taxonomy_id,
              CASE rgp.canonical_type
                WHEN 'entrez' THEN {identifier_type_id(ENTREZ_TYPE)}
                WHEN 'uniprot' THEN {identifier_type_id(UNIPROT_TYPE)}
              END::BIGINT AS canonical_identifier_type_id,
              rgp.canonical_id::VARCHAR AS canonical_identifier
            FROM read_parquet(
              '{combined_parquet_glob}',
              hive_partitioning = true
            ) rgp
            WHERE rgp.source_type IN ({gene_source_types})
              AND rgp.source_id IS NOT NULL
              AND rgp.canonical_id IS NOT NULL
              AND rgp.canonical_type IN ('entrez', 'uniprot')
            UNION ALL
            SELECT
              {_sql_literal(CHEMICAL_ENTITY_TYPE)} AS entity_type,
              ({chemical_key_case})::BIGINT AS key_identifier_type_id,
              rc.source_id::VARCHAR AS key_value,
              NULL::VARCHAR AS taxonomy_id,
              {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}::BIGINT
                AS canonical_identifier_type_id,
              rc.inchikey::VARCHAR AS canonical_identifier
            FROM utils_pg.omnipath_utils.resolver_chemical rc
            WHERE rc.source_type IN ({chemical_source_types})
              AND rc.source_id IS NOT NULL
              AND rc.inchikey IS NOT NULL
            UNION ALL
            SELECT DISTINCT
              {_sql_literal(CHEMICAL_ENTITY_TYPE)} AS entity_type,
              {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}::BIGINT
                AS key_identifier_type_id,
              rc.inchikey::VARCHAR AS key_value,
              NULL::VARCHAR AS taxonomy_id,
              {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}::BIGINT
                AS canonical_identifier_type_id,
              rc.inchikey::VARCHAR AS canonical_identifier
            FROM utils_pg.omnipath_utils.resolver_chemical rc
            WHERE rc.inchikey IS NOT NULL
              AND rc.inchikey <> ''
            {mirna_lookup_sql}
            """
        )
        _create_empty_protein_uniprot_fallback_view(con)
        con.execute(
            f"""
            CREATE VIEW resolver_canonical_entity AS
            WITH gene_product AS (
              SELECT
                CASE rgp.canonical_type
                  WHEN 'entrez' THEN {_sql_literal(GENE_ENTITY_TYPE)}
                  WHEN 'uniprot' THEN {_sql_literal(PROTEIN_ENTITY_TYPE)}
                END AS entity_type,
                NULLIF(rgp.ncbi_tax_id, 0)::VARCHAR AS taxonomy_id,
                CASE rgp.canonical_type
                  WHEN 'entrez' THEN {identifier_type_id(ENTREZ_TYPE)}
                  WHEN 'uniprot' THEN {identifier_type_id(UNIPROT_TYPE)}
                END::BIGINT AS canonical_identifier_type_id,
                rgp.canonical_id::VARCHAR AS canonical_identifier,
                ({combined_key_case})::BIGINT AS key_identifier_type_id
              FROM read_parquet(
                '{combined_parquet_glob}',
                hive_partitioning = true
              ) rgp
              WHERE rgp.source_type IN ({gene_source_types})
                AND rgp.canonical_id IS NOT NULL
                AND rgp.canonical_type IN ('entrez', 'uniprot')
            )
            SELECT
              row_number() OVER (
                ORDER BY
                  entity_type,
                  taxonomy_id NULLS FIRST,
                  canonical_identifier_type_id,
                  canonical_identifier
              )::BIGINT AS resolver_entity_id,
              entity_type,
              taxonomy_id,
              canonical_identifier_type_id,
              canonical_identifier,
              list_distinct(list(key_identifier_type_id))
                AS key_identifier_type_ids,
              count(*)::BIGINT AS lookup_rows
            FROM gene_product
            GROUP BY
              entity_type,
              taxonomy_id,
              canonical_identifier_type_id,
              canonical_identifier
            UNION ALL
            SELECT
              row_number() OVER (
                ORDER BY {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}, rc.inchikey
              )::BIGINT AS resolver_entity_id,
              {_sql_literal(CHEMICAL_ENTITY_TYPE)} AS entity_type,
              NULL::VARCHAR AS taxonomy_id,
              {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}::BIGINT
                AS canonical_identifier_type_id,
              rc.inchikey::VARCHAR AS canonical_identifier,
              list_distinct(list(({chemical_key_case})::BIGINT))
                AS key_identifier_type_ids,
              count(*)::BIGINT AS lookup_rows
            FROM utils_pg.omnipath_utils.resolver_chemical rc
            WHERE rc.source_type IN ({chemical_source_types})
              AND rc.inchikey IS NOT NULL
            GROUP BY rc.inchikey
            {mirna_canonical_sql}
            """
        )
    elif _attached_utils_relation_exists(
        con,
        'resolver_gene_protein_combined',
    ):
        con.execute(
            f"""
            CREATE VIEW resolver_lookup AS
            SELECT
              CASE rgp.canonical_type
                WHEN 'entrez' THEN {_sql_literal(GENE_ENTITY_TYPE)}
                WHEN 'uniprot' THEN {_sql_literal(PROTEIN_ENTITY_TYPE)}
              END AS entity_type,
              ({combined_key_case})::BIGINT AS key_identifier_type_id,
              rgp.source_id::VARCHAR AS key_value,
              NULLIF(rgp.ncbi_tax_id, 0)::VARCHAR AS taxonomy_id,
              CASE rgp.canonical_type
                WHEN 'entrez' THEN {identifier_type_id(ENTREZ_TYPE)}
                WHEN 'uniprot' THEN {identifier_type_id(UNIPROT_TYPE)}
              END::BIGINT AS canonical_identifier_type_id,
              rgp.canonical_id::VARCHAR AS canonical_identifier
            FROM utils_pg.omnipath_utils.resolver_gene_protein_combined rgp
            WHERE rgp.source_type IN ({gene_source_types})
              AND rgp.source_id IS NOT NULL
              AND rgp.canonical_id IS NOT NULL
              AND rgp.canonical_type IN ('entrez', 'uniprot')
            UNION ALL
            SELECT
              {_sql_literal(CHEMICAL_ENTITY_TYPE)} AS entity_type,
              ({chemical_key_case})::BIGINT AS key_identifier_type_id,
              rc.source_id::VARCHAR AS key_value,
              NULL::VARCHAR AS taxonomy_id,
              {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}::BIGINT
                AS canonical_identifier_type_id,
              rc.inchikey::VARCHAR AS canonical_identifier
            FROM utils_pg.omnipath_utils.resolver_chemical rc
            WHERE rc.source_type IN ({chemical_source_types})
              AND rc.source_id IS NOT NULL
              AND rc.inchikey IS NOT NULL
            UNION ALL
            SELECT DISTINCT
              {_sql_literal(CHEMICAL_ENTITY_TYPE)} AS entity_type,
              {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}::BIGINT
                AS key_identifier_type_id,
              rc.inchikey::VARCHAR AS key_value,
              NULL::VARCHAR AS taxonomy_id,
              {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}::BIGINT
                AS canonical_identifier_type_id,
              rc.inchikey::VARCHAR AS canonical_identifier
            FROM utils_pg.omnipath_utils.resolver_chemical rc
            WHERE rc.inchikey IS NOT NULL
              AND rc.inchikey <> ''
            {mirna_lookup_sql}
            """
        )
        _create_empty_protein_uniprot_fallback_view(con)
        con.execute(
            f"""
            CREATE VIEW resolver_canonical_entity AS
            WITH gene_product AS (
              SELECT
                CASE rgp.canonical_type
                  WHEN 'entrez' THEN {_sql_literal(GENE_ENTITY_TYPE)}
                  WHEN 'uniprot' THEN {_sql_literal(PROTEIN_ENTITY_TYPE)}
                END AS entity_type,
                NULLIF(rgp.ncbi_tax_id, 0)::VARCHAR AS taxonomy_id,
                CASE rgp.canonical_type
                  WHEN 'entrez' THEN {identifier_type_id(ENTREZ_TYPE)}
                  WHEN 'uniprot' THEN {identifier_type_id(UNIPROT_TYPE)}
                END::BIGINT AS canonical_identifier_type_id,
                rgp.canonical_id::VARCHAR AS canonical_identifier,
                ({combined_key_case})::BIGINT AS key_identifier_type_id
              FROM utils_pg.omnipath_utils.resolver_gene_protein_combined rgp
              WHERE rgp.source_type IN ({gene_source_types})
                AND rgp.canonical_id IS NOT NULL
                AND rgp.canonical_type IN ('entrez', 'uniprot')
            )
            SELECT
              row_number() OVER (
                ORDER BY
                  entity_type,
                  taxonomy_id NULLS FIRST,
                  canonical_identifier_type_id,
                  canonical_identifier
              )::BIGINT AS resolver_entity_id,
              entity_type,
              taxonomy_id,
              canonical_identifier_type_id,
              canonical_identifier,
              list_distinct(list(key_identifier_type_id))
                AS key_identifier_type_ids,
              count(*)::BIGINT AS lookup_rows
            FROM gene_product
            GROUP BY
              entity_type,
              taxonomy_id,
              canonical_identifier_type_id,
              canonical_identifier
            UNION ALL
            SELECT
              row_number() OVER (
                ORDER BY {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}, rc.inchikey
              )::BIGINT AS resolver_entity_id,
              {_sql_literal(CHEMICAL_ENTITY_TYPE)} AS entity_type,
              NULL::VARCHAR AS taxonomy_id,
              {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}::BIGINT
                AS canonical_identifier_type_id,
              rc.inchikey::VARCHAR AS canonical_identifier,
              list_distinct(list(({chemical_key_case})::BIGINT))
                AS key_identifier_type_ids,
              count(*)::BIGINT AS lookup_rows
            FROM utils_pg.omnipath_utils.resolver_chemical rc
            WHERE rc.source_type IN ({chemical_source_types})
              AND rc.inchikey IS NOT NULL
            GROUP BY rc.inchikey
            {mirna_canonical_sql}
            """
        )
    else:
        con.execute(
            f"""
            CREATE VIEW resolver_lookup AS
            SELECT
              {_sql_literal(GENE_ENTITY_TYPE)} AS entity_type,
              ({gene_key_case})::BIGINT AS key_identifier_type_id,
              rg.source_id::VARCHAR AS key_value,
              rg.ncbi_tax_id::VARCHAR AS taxonomy_id,
              {identifier_type_id(ENTREZ_TYPE)}::BIGINT
                AS canonical_identifier_type_id,
              rg.entrez::VARCHAR AS canonical_identifier
            FROM utils_pg.omnipath_utils.resolver_gene rg
            WHERE rg.source_type IN ({gene_source_types})
              AND rg.source_id IS NOT NULL
              AND rg.entrez IS NOT NULL
            UNION ALL
            SELECT
              {_sql_literal(CHEMICAL_ENTITY_TYPE)} AS entity_type,
              ({chemical_key_case})::BIGINT AS key_identifier_type_id,
              rc.source_id::VARCHAR AS key_value,
              NULL::VARCHAR AS taxonomy_id,
              {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}::BIGINT
                AS canonical_identifier_type_id,
              rc.inchikey::VARCHAR AS canonical_identifier
            FROM utils_pg.omnipath_utils.resolver_chemical rc
            WHERE rc.source_type IN ({chemical_source_types})
              AND rc.source_id IS NOT NULL
              AND rc.inchikey IS NOT NULL
            UNION ALL
            SELECT DISTINCT
              {_sql_literal(CHEMICAL_ENTITY_TYPE)} AS entity_type,
              {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}::BIGINT
                AS key_identifier_type_id,
              rc.inchikey::VARCHAR AS key_value,
              NULL::VARCHAR AS taxonomy_id,
              {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}::BIGINT
                AS canonical_identifier_type_id,
              rc.inchikey::VARCHAR AS canonical_identifier
            FROM utils_pg.omnipath_utils.resolver_chemical rc
            WHERE rc.inchikey IS NOT NULL
              AND rc.inchikey <> ''
            {mirna_lookup_sql}
            """
        )
        _create_protein_uniprot_fallback_view(con)
        con.execute(
            f"""
            CREATE VIEW resolver_canonical_entity AS
            SELECT
              row_number() OVER (
                ORDER BY
                  rg.ncbi_tax_id NULLS FIRST,
                  {identifier_type_id(ENTREZ_TYPE)},
                  rg.entrez
              )::BIGINT AS resolver_entity_id,
              {_sql_literal(GENE_ENTITY_TYPE)} AS entity_type,
              rg.ncbi_tax_id::VARCHAR AS taxonomy_id,
              {identifier_type_id(ENTREZ_TYPE)}::BIGINT
                AS canonical_identifier_type_id,
              rg.entrez::VARCHAR AS canonical_identifier,
              list_distinct(list(({gene_key_case})::BIGINT))
                AS key_identifier_type_ids,
              count(*)::BIGINT AS lookup_rows
            FROM utils_pg.omnipath_utils.resolver_gene rg
            WHERE rg.source_type IN ({gene_source_types})
              AND rg.entrez IS NOT NULL
            GROUP BY
              rg.ncbi_tax_id,
              rg.entrez
            UNION ALL
            SELECT
              row_number() OVER (
                ORDER BY {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}, rc.inchikey
              )::BIGINT AS resolver_entity_id,
              {_sql_literal(CHEMICAL_ENTITY_TYPE)} AS entity_type,
              NULL::VARCHAR AS taxonomy_id,
              {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}::BIGINT
                AS canonical_identifier_type_id,
              rc.inchikey::VARCHAR AS canonical_identifier,
              list_distinct(list(({chemical_key_case})::BIGINT))
                AS key_identifier_type_ids,
              count(*)::BIGINT AS lookup_rows
            FROM utils_pg.omnipath_utils.resolver_chemical rc
            WHERE rc.source_type IN ({chemical_source_types})
              AND rc.inchikey IS NOT NULL
            GROUP BY rc.inchikey
            {mirna_canonical_sql}
            """
        )
    con.execute(
        """
        CREATE VIEW gene_protein_representative_src AS
        WITH protein AS (
          SELECT
            rp.ncbi_tax_id::VARCHAR AS taxonomy_id,
            rp.source_id::VARCHAR AS canonical_identifier,
            rp.uniprot::VARCHAR AS uniprot,
            (sw.identifier IS NOT NULL) AS is_reviewed
          FROM utils_pg.omnipath_utils.resolver_protein rp
          LEFT JOIN utils_pg.omnipath_utils.reflist sw
            ON sw.list_name = 'swissprot'
           AND sw.identifier = rp.uniprot
          WHERE rp.source_type = 'entrez'
            AND rp.source_id IS NOT NULL
            AND rp.uniprot IS NOT NULL
        )
        SELECT
          taxonomy_id,
          canonical_identifier,
          coalesce(
            min(uniprot) FILTER (WHERE is_reviewed),
            min(uniprot)
          ) AS representative_uniprot,
          bool_or(is_reviewed) AS is_reviewed,
          list_sort(list_distinct(list(uniprot))) AS uniprot_all
        FROM protein
        GROUP BY taxonomy_id, canonical_identifier
        """
    )


def _create_duckdb_content_uuid_macro(con: duckdb.DuckDBPyConnection) -> None:
    con.execute(
        """
        CREATE OR REPLACE MACRO content_uuid(payload) AS (
          (
            substr(md5(payload), 1, 8) || '-' ||
            substr(md5(payload), 9, 4) || '-' ||
            substr(md5(payload), 13, 4) || '-' ||
            substr(md5(payload), 17, 4) || '-' ||
            substr(md5(payload), 21, 12)
          )::UUID
        )
        """
    )
    con.execute(
        """
        CREATE OR REPLACE MACRO canonical_entity_key(
          entity_type,
          taxonomy_id,
          canonical_identifier_type,
          canonical_identifier
        ) AS (
          to_json(
            list_value(
              entity_type::VARCHAR,
              coalesce(taxonomy_id::VARCHAR, ''),
              canonical_identifier_type::VARCHAR,
              canonical_identifier::VARCHAR
            )
          )
        )
        """
    )
    # The interaction header's identity: a content hash of the **sorted
    # participant multiset** and the interaction class. Sorting is what makes
    # it endpoint-independent — the same participants in any input order name
    # the same interaction, so A→B and B→A are two facts of one interaction —
    # and keeping it a multiset rather than a set keeps a homodimer distinct
    # from a single participant. The participants are cast to lowercase text
    # before sorting so the payload is byte-identical to the one the Postgres
    # derive step builds (`db/derived_tables.py`, `interaction_content_uuid_sql`);
    # the two engines must mint the same uuid for the same content.
    con.execute(
        """
        CREATE OR REPLACE MACRO interaction_content_key(
          participant_ids,
          interaction_class
        ) AS (
          to_json(
            list_prepend(
              interaction_class::VARCHAR,
              list_sort(
                list_transform(participant_ids, x -> lower(x::VARCHAR))
              )
            )
          )
        )
        """
    )
    con.execute(
        """
        CREATE OR REPLACE MACRO interaction_content_uuid(
          participant_ids,
          interaction_class
        ) AS (
          content_uuid(
            interaction_content_key(participant_ids, interaction_class)
          )
        )
        """
    )
    con.execute(
        """
        CREATE OR REPLACE MACRO canonical_entity_uuid(
          entity_type,
          taxonomy_id,
          canonical_identifier_type,
          canonical_identifier
        ) AS (
          content_uuid(
            canonical_entity_key(
              entity_type,
              taxonomy_id,
              canonical_identifier_type,
              canonical_identifier
            )
          )
        )
        """
    )


def _create_duckdb_resolver_views(
    con: duckdb.DuckDBPyConnection,
    *,
    resolver_dir: Path,
) -> None:
    """Expose resolver parquet inputs in the DuckDB shape used by canonicalize."""

    if _duckdb_attach_utils_postgres(con):
        _create_live_utils_resolver_views(con, resolver_dir=resolver_dir)
        return

    protein_dir = _resolver_component_dir(
        resolver_dir,
        component='proteins',
        lookup_filename='protein_identifier_lookup.parquet',
    )
    chemical_dir = _chemical_resolver_component_dir(resolver_dir)
    protein_lookup_path = protein_dir / 'protein_identifier_lookup.parquet'
    protein_type_path = protein_dir / 'identifier_type.parquet'
    gene_protein_representative_path = (
        protein_dir / 'gene_protein_representative.parquet'
    )
    chemical_lookup_dir = chemical_dir / 'lookup'
    chemical_lookup_glob = chemical_lookup_dir / '*.parquet'
    chemical_type_path = chemical_dir / 'identifier_type.parquet'

    # miRNA: organism-agnostic name/accession -> MI#/MIMAT#. The
    # lookup already carries the identity rows (MI#->MI#, MIMAT#->MIMAT#), so it
    # both resolves precursor/mature names and collapses the maturation-stub
    # matures onto their MIMAT#. Guarded so resolver snapshots predating the
    # mirbase source still load.
    try:
        mirna_dir = _resolver_component_dir(
            resolver_dir,
            component='mirna',
            lookup_filename='mirna_identifier_lookup.parquet',
        )
    except FileNotFoundError:
        mirna_dir = None
    if mirna_dir is not None:
        mirna_lookup_path = mirna_dir / 'mirna_identifier_lookup.parquet'
        mirna_type_path = mirna_dir / 'identifier_type.parquet'
        mirna_identifier_type_sql = f"""
        UNION
        SELECT *
        FROM read_parquet({_sql_literal(mirna_type_path)})
        """
        mirna_lookup_sql = f"""
        UNION ALL
        SELECT
          {_sql_literal(MIRNA_ENTITY_TYPE)} AS entity_type,
          key_identifier_type_id,
          key_value,
          NULL::VARCHAR AS taxonomy_id,
          canonical_identifier_type_id,
          canonical_identifier
        FROM read_parquet({_sql_literal(mirna_lookup_path)})
        WHERE key_value IS NOT NULL
          AND canonical_identifier IS NOT NULL
        """
        mirna_canonical_sql = f"""
        UNION ALL
        SELECT
          row_number() OVER (
            ORDER BY canonical_identifier_type_id, canonical_identifier
          )::BIGINT AS resolver_entity_id,
          {_sql_literal(MIRNA_ENTITY_TYPE)} AS entity_type,
          NULL::VARCHAR AS taxonomy_id,
          canonical_identifier_type_id,
          canonical_identifier,
          list_distinct(list(key_identifier_type_id)) AS key_identifier_type_ids,
          count(*)::BIGINT AS lookup_rows
        FROM read_parquet({_sql_literal(mirna_lookup_path)})
        WHERE canonical_identifier IS NOT NULL
        GROUP BY
          canonical_identifier_type_id,
          canonical_identifier
        """
    else:
        mirna_identifier_type_sql = ''
        mirna_lookup_sql = ''
        mirna_canonical_sql = ''
    chemical_identifier_type_sql = f"""
        UNION
        SELECT *
        FROM read_parquet({_sql_literal(chemical_type_path)})
        """
    chemical_lookup_sql = f"""
        UNION ALL
        SELECT
          {_sql_literal(CHEMICAL_ENTITY_TYPE)} AS entity_type,
          key_identifier_type_id,
          key_value,
          NULL::VARCHAR AS taxonomy_id,
          canonical_identifier_type_id,
          canonical_identifier
        FROM read_parquet('{_parquet_glob(chemical_lookup_glob)}')
        WHERE key_value IS NOT NULL
          AND canonical_identifier IS NOT NULL
        UNION ALL
        SELECT DISTINCT
          {_sql_literal(CHEMICAL_ENTITY_TYPE)} AS entity_type,
          {identifier_type_id(STANDARD_INCHI_KEY_TYPE)} AS key_identifier_type_id,
          canonical_identifier AS key_value,
          NULL::VARCHAR AS taxonomy_id,
          canonical_identifier_type_id,
          canonical_identifier
        FROM read_parquet('{_parquet_glob(chemical_lookup_glob)}')
        WHERE canonical_identifier_type_id =
              {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}
          AND canonical_identifier IS NOT NULL
          AND canonical_identifier <> ''
        """
    con.execute(
        f"""
        CREATE VIEW identifier_type AS
        SELECT *
        FROM read_parquet({_sql_literal(protein_type_path)})
        {chemical_identifier_type_sql}
        {mirna_identifier_type_sql}
        """
    )
    _create_duckdb_identifier_type_all_view(con)
    con.execute(
        f"""
        CREATE VIEW resolver_lookup AS
        SELECT
          {_sql_literal(GENE_ENTITY_TYPE)} AS entity_type,
          key_identifier_type_id,
          key_value,
          taxonomy_id,
          canonical_identifier_type_id,
          canonical_identifier
        FROM read_parquet({_sql_literal(protein_lookup_path)})
        WHERE key_value IS NOT NULL
          AND canonical_identifier IS NOT NULL
        {chemical_lookup_sql}
        {mirna_lookup_sql}
        """
    )
    _create_protein_uniprot_fallback_view(con)
    con.execute(
        f"""
        CREATE VIEW resolver_canonical_entity AS
        SELECT
          row_number() OVER (
            ORDER BY
              taxonomy_id NULLS FIRST,
              canonical_identifier_type_id,
              canonical_identifier
          )::BIGINT AS resolver_entity_id,
          {_sql_literal(GENE_ENTITY_TYPE)} AS entity_type,
          taxonomy_id,
          canonical_identifier_type_id,
          canonical_identifier,
          list_distinct(list(key_identifier_type_id)) AS key_identifier_type_ids,
          count(*)::BIGINT AS lookup_rows
        FROM read_parquet({_sql_literal(protein_lookup_path)})
        WHERE canonical_identifier IS NOT NULL
        GROUP BY
          taxonomy_id,
          canonical_identifier_type_id,
          canonical_identifier
        UNION ALL
        SELECT
          row_number() OVER (
            ORDER BY canonical_identifier_type_id, canonical_identifier
          )::BIGINT AS resolver_entity_id,
          {_sql_literal(CHEMICAL_ENTITY_TYPE)} AS entity_type,
          NULL::VARCHAR AS taxonomy_id,
          canonical_identifier_type_id,
          canonical_identifier,
          list_distinct(list(key_identifier_type_id)) AS key_identifier_type_ids,
          count(*)::BIGINT AS lookup_rows
        FROM read_parquet('{_parquet_glob(chemical_lookup_glob)}')
        WHERE canonical_identifier IS NOT NULL
        GROUP BY
          canonical_identifier_type_id,
          canonical_identifier
        {mirna_canonical_sql}
        """
    )
    # Per-gene representative UniProt. Guarded so resolver snapshots predating
    # the gene_protein_representative output still load; the canonicalize step
    # skips the table when the view is absent.
    if gene_protein_representative_path.exists():
        con.execute(
            f"""
            CREATE VIEW gene_protein_representative_src AS
            SELECT
              taxonomy_id,
              canonical_identifier,
              representative_uniprot,
              is_reviewed,
              uniprot_all
            FROM read_parquet({_sql_literal(gene_protein_representative_path)})
            WHERE canonical_identifier IS NOT NULL
            """
        )


def _create_duckdb_identifier_type_all_view(
    con: duckdb.DuckDBPyConnection,
) -> None:
    """Expose resolver identifier types plus local synthetic namespaces."""

    static_rows_sql = ',\n'.join(
        f'({identifier_type_id(name)}, {_sql_literal(name)})'
        for name in IDENTIFIER_TYPE_NAMES
    )
    has_entity_identifier_raw = bool(
        con.execute(
            """
            SELECT COUNT(*)
            FROM information_schema.tables
            WHERE table_name = 'entity_identifier_raw'
            """
        ).fetchone()[0]
    )
    has_annotation_relation_evidence_raw = bool(
        con.execute(
            """
            SELECT COUNT(*)
            FROM information_schema.tables
            WHERE table_name = 'annotation_relation_evidence_raw'
            """
        ).fetchone()[0]
    )
    has_ontology_relation_raw = bool(
        con.execute(
            """
            SELECT COUNT(*)
            FROM information_schema.tables
            WHERE table_name = 'ontology_relation_raw'
            """
        ).fetchone()[0]
    )
    evidence_identifier_type_sql = (
        """
          UNION
          SELECT DISTINCT identifier_type AS name
          FROM entity_identifier_raw
          WHERE identifier_type IS NOT NULL
            AND identifier_type <> ''
        """
        if has_entity_identifier_raw
        else ''
    )
    annotation_relation_identifier_type_sql = (
        """
          UNION
          SELECT DISTINCT object_id_type AS name
          FROM annotation_relation_evidence_raw
          WHERE object_id_type IS NOT NULL
            AND object_id_type <> ''
        """
        if has_annotation_relation_evidence_raw
        else ''
    )
    ontology_relation_identifier_type_sql = (
        """
          UNION
          SELECT DISTINCT subject_identifier_type AS name
          FROM ontology_relation_raw
          WHERE subject_identifier_type IS NOT NULL
            AND subject_identifier_type <> ''
          UNION
          SELECT DISTINCT object_identifier_type AS name
          FROM ontology_relation_raw
          WHERE object_identifier_type IS NOT NULL
            AND object_identifier_type <> ''
        """
        if has_ontology_relation_raw
        else ''
    )
    con.execute(
        f"""
        CREATE OR REPLACE VIEW identifier_type_all AS
        WITH static_identifier_type(identifier_type_id, name) AS (
          VALUES
          {static_rows_sql}
        ),
        base_identifier_type AS (
          SELECT * FROM identifier_type
          UNION
          SELECT * FROM static_identifier_type
        ),
        base AS (
          SELECT coalesce(max(identifier_type_id), 0) AS max_id
          FROM base_identifier_type
        ),
        required_name AS (
          SELECT {_sql_literal(CV_TERM_ID_TYPE)} AS name
          UNION ALL
          SELECT {_sql_literal(UNRESOLVED_ID_TYPE)} AS name
          {evidence_identifier_type_sql}
          {annotation_relation_identifier_type_sql}
          {ontology_relation_identifier_type_sql}
        ),
        missing AS (
          SELECT DISTINCT required_name.name
          FROM required_name
          LEFT JOIN base_identifier_type base_type
            ON base_type.name = required_name.name
          WHERE base_type.identifier_type_id IS NULL
        )
        SELECT * FROM base_identifier_type
        UNION ALL
        SELECT
          base.max_id + row_number() OVER (ORDER BY missing.name)
            AS identifier_type_id,
          missing.name
        FROM missing
        CROSS JOIN base
        """
    )


class DuckDBEvidenceProjector(EvidenceProjectorBase):
    """Flatten silver entity streams into DuckDB evidence tables."""

    def __init__(
        self,
        con: duckdb.DuckDBPyConnection,
        *,
        chunk_size: int = 100_000,
    ) -> None:
        super().__init__(chunk_size=chunk_size)
        self.con = con

    def project_records(
        self,
        records: Iterable[object],
        *,
        source: str,
        dataset: str,
        max_records: int | None = None,
        row_offset: int = 0,
    ) -> ProjectionStats:
        """Project source records into loaded DuckDB evidence tables."""

        if max_records is not None:
            records = islice(records, max_records)

        writers = _DuckDBEvidenceWriters(self.con, chunk_size=self.chunk_size)
        seen_annotations: set[tuple[str, str, str | None, str | None]] = set()
        stats = _MutableProjectionStats()
        try:
            for index, item in enumerate(records, start=row_offset + 1):
                entity, _ = unwrap_record(item)
                if not isinstance(entity, Entity):
                    continue
                self._flatten_entity_tree(
                    entity,
                    source=source,
                    dataset=dataset,
                    row_id=index,
                    occurrence_id=f'{dataset}:{index}:parent',
                    parent_entity_evidence_id=None,
                    entity_role='parent',
                    writers=writers,
                    seen_annotations=seen_annotations,
                    stats=stats,
                )
                stats.source_rows += 1
        finally:
            writers.close()
        return stats.freeze()


class _DuckDBEvidenceWriters:
    def __init__(
        self,
        con: duckdb.DuckDBPyConnection,
        *,
        chunk_size: int,
    ) -> None:
        self.entity = _DuckDBRowWriter(
            con,
            'entity_evidence_raw',
            _ENTITY_EVIDENCE_SCHEMA,
            chunk_size=chunk_size,
        )
        self.identifier = _DuckDBRowWriter(
            con,
            'entity_identifier_raw',
            _ENTITY_IDENTIFIER_SCHEMA,
            chunk_size=chunk_size,
        )
        self.entity_annotation = _DuckDBRowWriter(
            con,
            'entity_annotation_raw',
            _ANNOTATION_REF_SCHEMA,
            chunk_size=chunk_size,
        )
        self.relation_annotation = _DuckDBRowWriter(
            con,
            'relation_annotation_raw',
            _RELATION_ANNOTATION_REF_SCHEMA,
            chunk_size=chunk_size,
        )
        self.annotation = _DuckDBRowWriter(
            con,
            'annotation_value',
            _ANNOTATION_VALUE_SCHEMA,
            chunk_size=chunk_size,
        )
        self.relation = _DuckDBRowWriter(
            con,
            'relation_evidence_raw',
            _RELATION_EVIDENCE_SCHEMA,
            chunk_size=chunk_size,
        )
        self.annotation_relation = _DuckDBRowWriter(
            con,
            'annotation_relation_evidence_raw',
            _ANNOTATION_RELATION_EVIDENCE_SCHEMA,
            chunk_size=chunk_size,
        )
        self.ontology_relation = _DuckDBRowWriter(
            con,
            'ontology_relation_raw',
            _ONTOLOGY_RELATION_SCHEMA,
            chunk_size=chunk_size,
        )

    def close(self) -> None:
        self.entity.close()
        self.identifier.close()
        self.entity_annotation.close()
        self.relation_annotation.close()
        self.annotation.close()
        self.relation.close()
        self.annotation_relation.close()
        self.ontology_relation.close()


class _DuckDBRowWriter:
    def __init__(
        self,
        con: duckdb.DuckDBPyConnection,
        table: str,
        schema: pa.Schema,
        *,
        chunk_size: int,
    ) -> None:
        self.con = con
        self.table = table
        self.schema = schema
        self.chunk_size = chunk_size
        self.rows: list[dict[str, object]] = []

    def write(self, row: dict[str, object]) -> None:
        self.rows.append(row)
        if len(self.rows) >= self.chunk_size:
            self.flush()

    def flush(self) -> None:
        if not self.rows:
            return
        batch = f'_batch_{self.table}'
        table = pa.Table.from_pylist(self.rows, schema=self.schema)
        self.con.register(batch, table)
        try:
            self.con.execute(f'INSERT INTO {self.table} SELECT * FROM {batch}')
        finally:
            self.con.unregister(batch)
        self.rows.clear()

    def close(self) -> None:
        self.flush()


_ENTITY_EVIDENCE_SCHEMA = pa.schema(
    [
        ('source', pa.string()),
        ('dataset', pa.string()),
        ('row_id', pa.int64()),
        ('entity_evidence_id', pa.string()),
        ('parent_entity_evidence_id', pa.string()),
        ('entity_role', pa.string()),
        ('entity_type', pa.string()),
        ('taxonomy_id', pa.string()),
    ]
)


_ENTITY_IDENTIFIER_SCHEMA = pa.schema(
    [
        ('source', pa.string()),
        ('entity_evidence_id', pa.string()),
        ('identifier_id', pa.string()),
        ('identifier_type', pa.string()),
        ('identifier', pa.string()),
        # WP1's normalized key, chemical identifier types only (spec 011
        # T040). NULL for every other identifier type -- the resolver join
        # falls back to `identifier` itself when this is unset.
        ('identifier_normalized', pa.string()),
    ]
)


_ANNOTATION_REF_SCHEMA = pa.schema(
    [
        ('source', pa.string()),
        ('evidence_id', pa.string()),
        ('annotation_key', pa.string()),
        ('term', pa.string()),
        ('value', pa.string()),
        ('unit', pa.string()),
    ]
)


_RELATION_ANNOTATION_REF_SCHEMA = pa.schema(
    [
        ('source', pa.string()),
        ('evidence_id', pa.string()),
        ('annotation_key', pa.string()),
        ('annotation_scope', pa.string()),
        ('term', pa.string()),
        ('value', pa.string()),
        ('unit', pa.string()),
    ]
)


_ANNOTATION_VALUE_SCHEMA = pa.schema(
    [
        ('annotation_key', pa.string()),
        ('term', pa.string()),
        ('value', pa.string()),
        ('unit', pa.string()),
    ]
)


_RELATION_EVIDENCE_SCHEMA = pa.schema(
    [
        ('source', pa.string()),
        ('dataset', pa.string()),
        ('row_id', pa.int64()),
        ('relation_evidence_id', pa.string()),
        ('subject_entity_evidence_id', pa.string()),
        ('predicate', pa.string()),
        ('object_entity_evidence_id', pa.string()),
        ('relation_category', pa.string()),
    ]
)


_ANNOTATION_RELATION_EVIDENCE_SCHEMA = pa.schema(
    [
        ('relation_evidence_id', pa.string()),
        ('source', pa.string()),
        ('dataset', pa.string()),
        ('row_id', pa.int64()),
        ('subject_entity_evidence_id', pa.string()),
        ('predicate', pa.string()),
        ('object_entity_type', pa.string()),
        ('object_id_type', pa.string()),
        ('object_id', pa.string()),
        ('relation_category', pa.string()),
    ]
)


_ONTOLOGY_RELATION_SCHEMA = pa.schema(
    [
        ('source', pa.string()),
        ('dataset', pa.string()),
        ('subject_entity_evidence_id', pa.string()),
        ('ontology_id', pa.string()),
        ('subject_entity_type', pa.string()),
        ('subject_identifier_type', pa.string()),
        ('subject_identifier', pa.string()),
        ('predicate', pa.string()),
        ('object_entity_type', pa.string()),
        ('object_identifier_type', pa.string()),
        ('object_identifier', pa.string()),
    ]
)


def _create_duckdb_evidence_tables(con: duckdb.DuckDBPyConnection) -> None:
    con.execute(
        """
        CREATE OR REPLACE TABLE entity_evidence_raw (
          source VARCHAR,
          dataset VARCHAR,
          row_id BIGINT,
          entity_evidence_id VARCHAR,
          parent_entity_evidence_id VARCHAR,
          entity_role VARCHAR,
          entity_type VARCHAR,
          taxonomy_id VARCHAR
        )
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TABLE ontology_terms_raw (
          source VARCHAR,
          dataset VARCHAR,
          term_id VARCHAR,
          term_entity_type VARCHAR,
          term_identifier_type VARCHAR,
          term_identifier VARCHAR,
          ontology_prefix VARCHAR,
          label VARCHAR,
          definition VARCHAR,
          ontology_id VARCHAR,
          synonyms VARCHAR[],
          synonyms_text VARCHAR,
          sources VARCHAR[]
        )
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TABLE ontology_relation_raw (
          source VARCHAR,
          dataset VARCHAR,
          subject_entity_evidence_id VARCHAR,
          ontology_id VARCHAR,
          subject_entity_type VARCHAR,
          subject_identifier_type VARCHAR,
          subject_identifier VARCHAR,
          predicate VARCHAR,
          object_entity_type VARCHAR,
          object_identifier_type VARCHAR,
          object_identifier VARCHAR
        )
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TABLE entity_identifier_raw (
          source VARCHAR,
          entity_evidence_id VARCHAR,
          identifier_id VARCHAR,
          identifier_type VARCHAR,
          identifier VARCHAR,
          identifier_normalized VARCHAR
        )
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TABLE entity_annotation_raw (
          source VARCHAR,
          evidence_id VARCHAR,
          annotation_key VARCHAR,
          term VARCHAR,
          value VARCHAR,
          unit VARCHAR
        )
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TABLE relation_annotation_raw (
          source VARCHAR,
          evidence_id VARCHAR,
          annotation_key VARCHAR,
          annotation_scope VARCHAR,
          term VARCHAR,
          value VARCHAR,
          unit VARCHAR
        )
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TABLE annotation_value (
          annotation_key VARCHAR,
          term VARCHAR,
          value VARCHAR,
          unit VARCHAR
        )
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TABLE relation_evidence_raw (
          source VARCHAR,
          dataset VARCHAR,
          row_id BIGINT,
          relation_evidence_id VARCHAR,
          subject_entity_evidence_id VARCHAR,
          predicate VARCHAR,
          object_entity_evidence_id VARCHAR,
          relation_category VARCHAR
        )
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TABLE annotation_relation_evidence_raw (
          relation_evidence_id VARCHAR,
          source VARCHAR,
          dataset VARCHAR,
          row_id BIGINT,
          subject_entity_evidence_id VARCHAR,
          predicate VARCHAR,
          object_entity_type VARCHAR,
          object_id_type VARCHAR,
          object_id VARCHAR,
          relation_category VARCHAR
        )
        """
    )


def _ensure_duckdb_canonical_caches(con: duckdb.DuckDBPyConnection) -> None:
    con.execute(
        """
        CREATE TABLE IF NOT EXISTS entity_key_cache (
          entity_id UUID,
          entity_key VARCHAR PRIMARY KEY,
          entity_type VARCHAR,
          taxonomy_id VARCHAR,
          canonical_identifier_type VARCHAR,
          canonical_identifier_type_id BIGINT,
          canonical_identifier VARCHAR,
          sources VARCHAR,
          first_seen_at TIMESTAMP,
          last_seen_at TIMESTAMP
        )
        """
    )
    con.execute(
        'ALTER TABLE entity_key_cache DROP COLUMN IF EXISTS identifiers_json'
    )
    con.execute(
        """
        CREATE TABLE IF NOT EXISTS relation_key_cache (
          relation_id UUID,
          relation_key VARCHAR PRIMARY KEY,
          subject_entity_id UUID,
          predicate VARCHAR,
          object_entity_id UUID,
          sources VARCHAR,
          first_seen_at TIMESTAMP,
          last_seen_at TIMESTAMP
        )
        """
    )


def _drop_duckdb_batch_tables(con: duckdb.DuckDBPyConnection) -> None:
    for table in (
        'evidence_identifier_key',
        'resolver_entity_type_match',
        'resolver_evidence_identifier_key',
        'ontology_endpoint_identifier_key',
        'taxonomy_optional_resolver_key_type',
        'taxonomy_optional_unambiguous_key',
        'needed_resolver_lookup',
        'needed_ontology_endpoint_resolver_lookup',
        'protein_uniprot_fallback_taxonomy_optional_unambiguous_key',
        'needed_protein_uniprot_fallback_lookup',
        'entity_identifier_group',
        'cv_term_evidence_resolution',
        'pathway_identifier_evidence_resolution',
        'pathway_stable_id_evidence_resolution',
        'standard_inchi_key_evidence_resolution',
        'entity_resolution_base',
        'complex_member_signature_base',
        'complex_member_signature',
        'reaction_member_signature',
        'entity_resolution',
        'ontology_term_resolution',
        'batch_entity_candidate',
        'new_entity',
        'canonical_entity',
        'canonical_entity_identifier',
        'gene_protein_representative',
        'entity_evidence_resolution',
        'evidence_state_link',
        'state',
        'state_component',
        'evidence_state',
        'entity_ontology_relation',
        'relation_candidate_evidence',
        'batch_relation_candidate',
        'new_relation',
        'relation',
        'relation_evidence_relation',
    ):
        con.execute(f'DROP TABLE IF EXISTS {table}')


def _live_utils_attached(con: duckdb.DuckDBPyConnection) -> bool:
    """True when the utils Postgres is ATTACHed as ``utils_pg`` in this DuckDB."""
    return any(
        row[1] == 'utils_pg'
        for row in con.execute('PRAGMA database_list').fetchall()
    )


def _bulk_insert_resolver_lookup_rows(
    con: duckdb.DuckDBPyConnection,
    table: str,
    rows: list[tuple[str, int, str, str | None, int, str]],
) -> None:
    """Bulk-insert fetched resolver rows into a local DuckDB ``resolver_lookup``.

    Row-by-row ``executemany`` is pathologically slow in DuckDB for large result
    sets — observed making per-shard canonicalize balloon to ~30 min on id-dense
    sources (ChEMBL), both Postgres instances idle (so the cost was the DuckDB
    insert, not the keyed fetch). Insert via an Arrow table (vectorised, one
    statement) instead. Column types are pinned so an all-null column (e.g. the
    organism-agnostic ``taxonomy_id`` for chemicals) still matches the schema.
    """
    if not rows:
        return
    columns = list(zip(*rows))
    arrow_table = pa.table({
        'entity_type': pa.array(columns[0], type=pa.string()),
        'key_identifier_type_id': pa.array(columns[1], type=pa.int64()),
        'key_value': pa.array(columns[2], type=pa.string()),
        'taxonomy_id': pa.array(columns[3], type=pa.string()),
        'canonical_identifier_type_id': pa.array(columns[4], type=pa.int64()),
        'canonical_identifier': pa.array(columns[5], type=pa.string()),
    })
    con.register('_resolver_lookup_rows', arrow_table)
    try:
        con.execute(
            f'INSERT INTO {table} SELECT * FROM _resolver_lookup_rows'
        )
    finally:
        con.unregister('_resolver_lookup_rows')


def _fetch_live_utils_rows_for_keys(
    *,
    url: str,
    table: str,
    mapping: dict[str, str],
    entity_type: str,
    canonical_identifier_type_id: int,
    canonical_column: str,
    key_rows: list[tuple[str, str, int | None]],
    has_taxonomy: bool = True,
) -> list[tuple[str, int, str, str | None, int, str]]:
    """Fetch resolver rows for the current shard's keys from utils Postgres.

    Keyed lookup (2026-07-04, adopted+finished from Jonathan's keyed-lookup WIP):
    only this shard's ids are shipped to Postgres (a temp table) and joined against
    the MATERIALIZED, indexed resolver, so index probes return just the needed
    rows -- never the full-table scan into DuckDB that the full-view path did (the
    ~500M-row-per-shard remote scan + multi-GB blow-up). ALL resolution and
    normalization now lives in the utils resolver (ADR 0006 + the authoritative
    gene2ensembl/BioMart gene-space); this only READS it -- so Jonathan's original
    inline re-derivations from ``id_mapping``/``id_mapping_ftp`` (one special-case
    per resolver) are dropped in favour of this single generic keyed join.
    ``has_taxonomy`` is False for the organism-agnostic ``resolver_chemical``.
    """
    if not key_rows:
        return []

    key_type_case = ' '.join(
        f"WHEN {source_type!r} THEN {identifier_type_id(identifier_type)}"
        for source_type, identifier_type in mapping.items()
    )
    tax_join = (
        'AND (k.ncbi_tax_id IS NULL OR k.ncbi_tax_id = r.ncbi_tax_id)'
        if has_taxonomy else ''
    )
    tax_select = 'NULLIF(r.ncbi_tax_id, 0)::text' if has_taxonomy else 'NULL::text'
    query = f"""
        SELECT
          %s AS entity_type,
          (CASE r.source_type {key_type_case} END)::bigint
            AS key_identifier_type_id,
          r.source_id::text AS key_value,
          {tax_select} AS taxonomy_id,
          %s::bigint AS canonical_identifier_type_id,
          r.{canonical_column}::text AS canonical_identifier
        FROM omnipath_utils.{table} r
        JOIN tmp_resolver_key k
          ON k.source_type = r.source_type
         AND k.source_id = r.source_id
         {tax_join}
        WHERE r.{canonical_column} IS NOT NULL
    """
    with psycopg2.connect(url) as pg:
        with pg.cursor() as cur:
            # source_id is TEXT (not varchar(64)): evidence can carry ids longer
            # than the resolver's varchar(64) source_id (BRENDA/MONDO), which used
            # to raise StringDataRightTruncation on insert. A >64-char key simply
            # never matches a (<=64-char) resolver row, so TEXT is both safe and
            # correct here.
            cur.execute(
                'CREATE TEMP TABLE tmp_resolver_key '
                '(source_type text, source_id text, ncbi_tax_id integer) '
                'ON COMMIT DROP'
            )
            # None-safe ordering: ncbi_tax_id is None for no-taxon evidence, so a
            # plain tuple sort raised "'<' not supported between int and NoneType"
            # whenever the same (source_type, source_id) appeared both with and
            # without a taxon (guidetopharma/bindingdb/drugcentral/wikipathways).
            execute_values(
                cur,
                'INSERT INTO tmp_resolver_key '
                '(source_type, source_id, ncbi_tax_id) VALUES %s',
                sorted(
                    set(key_rows),
                    key=lambda r: (r[0], r[1], -1 if r[2] is None else r[2]),
                ),
                page_size=10000,
            )
            cur.execute(query, (entity_type, canonical_identifier_type_id))
            return cur.fetchall()


def _live_utils_key_rows(
    con: duckdb.DuckDBPyConnection,
    *,
    mapping: dict[str, str],
    evidence_entity_types: tuple[str, ...],
) -> list[tuple[str, str, int | None]]:
    """The shard's distinct (source_type, source_id, ncbi_tax_id) resolver keys."""
    id_to_source_types: dict[int, list[str]] = {}
    for source_type, identifier_type in mapping.items():
        id_to_source_types.setdefault(
            identifier_type_id(identifier_type), []
        ).append(source_type)

    entity_placeholders = ', '.join('?' for _ in evidence_entity_types)
    id_placeholders = ', '.join('?' for _ in id_to_source_types)
    rows = con.execute(
        f"""
        SELECT DISTINCT
          key_identifier_type_id,
          key_value,
          try_cast(taxonomy_id AS INTEGER) AS ncbi_tax_id
        FROM evidence_identifier_key
        WHERE entity_type IN ({entity_placeholders})
          AND key_identifier_type_id IN ({id_placeholders})
          AND key_value IS NOT NULL
          AND key_value <> ''
        """,
        [*evidence_entity_types, *id_to_source_types.keys()],
    ).fetchall()

    keys: list[tuple[str, str, int | None]] = []
    for key_identifier_type_id, key_value, ncbi_tax_id in rows:
        for source_type in id_to_source_types[int(key_identifier_type_id)]:
            keys.append((source_type, str(key_value), ncbi_tax_id))
    return keys


# The identifier-resolution (utils) relations the live build reads, and the
# evidence entity types each one can resolve. An entity type is only left
# unresolvable *for lack of reference data* when EVERY relation that could
# resolve it is absent or empty.
REQUIRED_TRANSLATION_TABLES: dict[str, tuple[str, ...]] = {
    'resolver_gene': (GENE_ENTITY_TYPE, PROTEIN_ENTITY_TYPE),
    'resolver_gene_protein_global': (GENE_ENTITY_TYPE, PROTEIN_ENTITY_TYPE),
    'resolver_protein': (PROTEIN_ENTITY_TYPE,),
    'resolver_chemical': (CHEMICAL_ENTITY_TYPE,),
}


def _attached_utils_relation_nonempty(
    con: duckdb.DuckDBPyConnection,
    relation_name: str,
) -> bool:
    """Whether an attached utils relation holds at least one row (cheap probe)."""
    return bool(
        con.execute(
            'SELECT EXISTS (SELECT 1 FROM utils_pg.omnipath_utils.'
            f'{_duckdb_identifier(relation_name)})'
        ).fetchone()[0]
    )


def preflight_translation_tables(
    con: duckdb.DuckDBPyConnection,
) -> list[dict[str, object]]:
    """Report whether each required identifier-resolution table is usable.

    Before resolving evidence against the connected resolver database, check that
    every relation the build depends on is present and non-empty, and log a
    warning for any that is missing or empty — so a silently degraded resolver
    database shows up in the build log instead of quietly producing unresolved
    records. Records the evidence entity types that no present relation can
    resolve, so canonicalisation can mark them as blocked by missing reference
    data rather than genuinely unresolvable. Returns the per-relation report
    (also surfaced in the build manifest). A no-op returning an empty report when
    no resolver database is attached.
    """
    con.execute(
        'CREATE TABLE IF NOT EXISTS missing_translation_entity_type '
        '(entity_type VARCHAR)'
    )
    if not _live_utils_attached(con):
        return []

    report: list[dict[str, object]] = []
    resolvable_entity_type: dict[str, bool] = {}
    for relation, entity_types in REQUIRED_TRANSLATION_TABLES.items():
        present = _attached_utils_relation_exists(con, relation)
        usable = present and _attached_utils_relation_nonempty(con, relation)
        status = 'present' if usable else ('empty' if present else 'absent')
        report.append({'relation': relation, 'status': status})
        if not usable:
            _log.warning(
                'identifier-resolution table %r is %s in the connected resolver '
                'database; records that rely on it cannot be resolved',
                relation,
                status,
            )
        for entity_type in entity_types:
            resolvable_entity_type[entity_type] = (
                resolvable_entity_type.get(entity_type, False) or usable
            )

    blocked = sorted(
        entity_type
        for entity_type, resolvable in resolvable_entity_type.items()
        if not resolvable
    )
    con.execute('DELETE FROM missing_translation_entity_type')
    if blocked:
        con.executemany(
            'INSERT INTO missing_translation_entity_type VALUES (?)',
            [(entity_type,) for entity_type in blocked],
        )
    _log.info(
        'identifier-resolution pre-flight: %s',
        ', '.join(f'{r["relation"]}={r["status"]}' for r in report),
    )
    return report


def _recover_no_taxon_protein_from_uniprot(
    con: duckdb.DuckDBPyConnection,
) -> None:
    """Recover the organism of a protein mention that arrived without one.

    Some sources cite a UniProt accession but state no organism (reaction and
    binding-affinity resources in particular). The accession itself names the
    organism, so look it up in the identifier-resolution database's
    accession -> organism map and record it on the mention *before* resolution.
    A recovered organism both lets the mention reach its gene through the
    per-organism gene map and, when no gene is known, lets it be typed as a gene
    keyed by the accession (which requires a known organism). No-op when the
    resolution database is not attached or carries no accession -> organism map.
    """
    _fetch_uniprot_taxon_lookup(con)
    _apply_uniprot_taxon_recovery(con)


def _fetch_uniprot_taxon_lookup(
    con: duckdb.DuckDBPyConnection,
) -> None:
    """Fetch accession -> organism for this shard's organism-less proteins.

    Ships the shard's well-formed accessions (protein mentions with no organism)
    to the identifier-resolution Postgres and reads back their organism from the
    ``resolver_uniprot_taxon`` map, into a local ``uniprot_taxon_lookup`` table.
    No-op (and leaves any existing lookup untouched) when the resolver is not
    attached or the map is absent.
    """
    url = os.environ.get('OMNIPATH_BUILD_UTILS_PG_URL')
    if not url or not _live_utils_attached(con):
        return
    accessions = [
        row[0]
        for row in con.execute(
            f"""
            SELECT DISTINCT regexp_replace(ei.identifier, '-[0-9]+$', '') AS ac
            FROM entity_evidence_raw ee
            JOIN entity_identifier_raw ei
              ON ei.source = ee.source
             AND ei.entity_evidence_id = ee.entity_evidence_id
            WHERE ee.entity_type = {_sql_literal(PROTEIN_ENTITY_TYPE)}
              AND ee.taxonomy_id IS NULL
              AND ei.identifier_type = {_sql_literal(UNIPROT_TYPE)}
              AND ei.identifier IS NOT NULL
              AND ei.identifier <> ''
              AND regexp_matches(
                regexp_replace(ei.identifier, '-[0-9]+$', ''),
                {_sql_literal(UNIPROT_AC_REGEX)}
              )
            """
        ).fetchall()
    ]
    con.execute('DROP TABLE IF EXISTS uniprot_taxon_lookup')
    con.execute(
        'CREATE TABLE uniprot_taxon_lookup (uniprot VARCHAR, taxonomy_id VARCHAR)'
    )
    if not accessions:
        return
    with psycopg2.connect(url) as pg:
        with pg.cursor() as cur:
            cur.execute(
                "SELECT to_regclass('omnipath_utils.resolver_uniprot_taxon')"
                ' IS NOT NULL'
            )
            if not cur.fetchone()[0]:
                return
            cur.execute(
                'CREATE TEMP TABLE tmp_uniprot_key (uniprot text) ON COMMIT DROP'
            )
            execute_values(
                cur,
                'INSERT INTO tmp_uniprot_key (uniprot) VALUES %s',
                [(a,) for a in sorted(set(accessions))],
                page_size=10000,
            )
            cur.execute(
                """
                SELECT r.uniprot, NULLIF(r.ncbi_tax_id, 0)::text
                FROM omnipath_utils.resolver_uniprot_taxon r
                JOIN tmp_uniprot_key k ON k.uniprot = r.uniprot
                WHERE r.ncbi_tax_id IS NOT NULL
                """
            )
            rows = cur.fetchall()
    if not rows:
        return
    columns = list(zip(*rows, strict=True))
    arrow_table = pa.table({
        'uniprot': pa.array(columns[0], type=pa.string()),
        'taxonomy_id': pa.array(columns[1], type=pa.string()),
    })
    con.register('_uniprot_taxon_rows', arrow_table)
    try:
        con.execute(
            'INSERT INTO uniprot_taxon_lookup SELECT * FROM _uniprot_taxon_rows'
        )
    finally:
        con.unregister('_uniprot_taxon_rows')
    _log.info(
        'organism recovery: %d of %d organism-less protein accessions '
        'matched an organism',
        len(rows),
        len(accessions),
    )


def _fetch_identifier_authority(
    con: duckdb.DuckDBPyConnection,
) -> None:
    """Fetch the structure-authority declarations into a local table (spec
    011 T058, WP2).

    ``identifier_authority`` lives in the *build's own* Postgres (populated
    by the resource walk before ``load`` runs), not the identifier-resolution
    Postgres this module otherwise reads via ``OMNIPATH_BUILD_UTILS_PG_URL``
    -- so this reads ``DATABASE_URL`` instead, the same env var the CLI's own
    ``--database-url`` default already reads. Small, static table (one row
    per identifier type), so no keyed/sharded fetch is needed -- the whole
    thing is read every time. Degrades to an empty local table (arbitration
    then finds no authority match, so a genuine disagreement stays
    unresolved, never silently guessed) when the env var is unset or the
    table does not exist yet -- a build's very first run, before the
    resource walk has populated it.
    """
    url = os.environ.get('DATABASE_URL')
    if not url:
        con.execute(
            'CREATE TABLE IF NOT EXISTS identifier_authority '
            '(identifier_type_id BIGINT, source_id BIGINT, '
            'is_structure_authority BOOLEAN)'
        )
        return
    with psycopg2.connect(url) as pg:
        with pg.cursor() as cur:
            cur.execute(
                "SELECT to_regclass('public.identifier_authority') IS NOT NULL"
            )
            if not cur.fetchone()[0]:
                con.execute(
                    'CREATE TABLE IF NOT EXISTS identifier_authority '
                    '(identifier_type_id BIGINT, source_id BIGINT, '
                    'is_structure_authority BOOLEAN)'
                )
                return
            cur.execute(
                'SELECT identifier_type_id, source_id, is_structure_authority '
                'FROM identifier_authority WHERE is_structure_authority'
            )
            rows = cur.fetchall()
    con.execute(
        'CREATE TABLE IF NOT EXISTS identifier_authority '
        '(identifier_type_id BIGINT, source_id BIGINT, '
        'is_structure_authority BOOLEAN)'
    )
    if not rows:
        return
    columns = list(zip(*rows, strict=True))
    arrow_table = pa.table({
        'identifier_type_id': pa.array(columns[0], type=pa.int64()),
        'source_id': pa.array(columns[1], type=pa.int64()),
        'is_structure_authority': pa.array(columns[2], type=pa.bool_()),
    })
    con.register('_identifier_authority_rows', arrow_table)
    try:
        con.execute(
            'INSERT INTO identifier_authority '
            'SELECT * FROM _identifier_authority_rows'
        )
    finally:
        con.unregister('_identifier_authority_rows')


def _apply_uniprot_taxon_recovery(
    con: duckdb.DuckDBPyConnection,
) -> None:
    """Write the recovered organism onto organism-less protein mentions.

    Reads the local ``uniprot_taxon_lookup`` (built by
    :func:`_fetch_uniprot_taxon_lookup`, or supplied directly in tests) and sets
    ``taxonomy_id`` on every protein mention that had none, from the organism of
    its UniProt accession. When one mention carries accessions of different
    organisms the lowest taxon id is chosen deterministically. No-op when the
    lookup is absent or empty.
    """
    has_lookup = bool(
        con.execute(
            """
            SELECT COUNT(*)
            FROM information_schema.tables
            WHERE table_name = 'uniprot_taxon_lookup'
            """
        ).fetchone()[0]
    )
    if not has_lookup:
        return
    con.execute(
        f"""
        UPDATE entity_evidence_raw AS ee
        SET taxonomy_id = sub.taxonomy_id
        FROM (
          SELECT ee2.source, ee2.entity_evidence_id,
                 min(utl.taxonomy_id) AS taxonomy_id
          FROM entity_evidence_raw ee2
          JOIN entity_identifier_raw ei
            ON ei.source = ee2.source
           AND ei.entity_evidence_id = ee2.entity_evidence_id
           AND ei.identifier_type = {_sql_literal(UNIPROT_TYPE)}
          JOIN uniprot_taxon_lookup utl
            ON utl.uniprot = regexp_replace(ei.identifier, '-[0-9]+$', '')
          WHERE ee2.entity_type = {_sql_literal(PROTEIN_ENTITY_TYPE)}
            AND ee2.taxonomy_id IS NULL
          GROUP BY ee2.source, ee2.entity_evidence_id
        ) AS sub
        WHERE ee.source = sub.source
          AND ee.entity_evidence_id = sub.entity_evidence_id
          AND ee.taxonomy_id IS NULL
        """
    )


def _canonicalize_taxon_to_species(
    con: duckdb.DuckDBPyConnection,
) -> None:
    """Canonicalize each mention's organism to species level (strain -> species).

    NCBI files some genes under a strain taxon while an Entrez GeneID belongs to
    one species; left unreconciled, the same gene surfaces under several NCBI
    strain taxa and splits into duplicate taxon-specific entities. Map every
    mention's taxon to its species-level taxon via the resolver's ``taxon_species``
    map (built from the NCBI taxdump) so mentions and the resolver agree on one
    taxon per gene. Ships this shard's distinct taxa to the resolver Postgres and
    reads back only the strain taxa that actually differ from their species. No-op
    when the resolver is not attached or the map is absent/empty.
    """
    url = os.environ.get('OMNIPATH_BUILD_UTILS_PG_URL')
    if not url or not _live_utils_attached(con):
        return
    taxa = [
        row[0]
        for row in con.execute(
            """
            SELECT DISTINCT taxonomy_id
            FROM entity_evidence_raw
            WHERE taxonomy_id IS NOT NULL AND taxonomy_id <> ''
            """
        ).fetchall()
    ]
    con.execute('DROP TABLE IF EXISTS taxon_species_lookup')
    con.execute(
        'CREATE TABLE taxon_species_lookup (tax_id VARCHAR, species_tax_id VARCHAR)'
    )
    numeric = sorted({t for t in taxa if t and t.isdigit()})
    if not numeric:
        return
    with psycopg2.connect(url) as pg:
        with pg.cursor() as cur:
            cur.execute(
                "SELECT to_regclass('omnipath_utils.taxon_species') IS NOT NULL"
            )
            if not cur.fetchone()[0]:
                return
            cur.execute(
                'CREATE TEMP TABLE tmp_taxon_key (tax_id integer) ON COMMIT DROP'
            )
            execute_values(
                cur,
                'INSERT INTO tmp_taxon_key (tax_id) VALUES %s',
                [(int(t),) for t in numeric],
                page_size=10000,
            )
            cur.execute(
                """
                SELECT s.tax_id::text, s.species_tax_id::text
                FROM omnipath_utils.taxon_species s
                JOIN tmp_taxon_key k ON k.tax_id = s.tax_id
                WHERE s.species_tax_id IS NOT NULL
                  AND s.species_tax_id <> s.tax_id
                """
            )
            rows = cur.fetchall()
    if not rows:
        return
    columns = list(zip(*rows, strict=True))
    arrow_table = pa.table({
        'tax_id': pa.array(columns[0], type=pa.string()),
        'species_tax_id': pa.array(columns[1], type=pa.string()),
    })
    con.register('_taxon_species_rows', arrow_table)
    try:
        con.execute(
            'INSERT INTO taxon_species_lookup SELECT * FROM _taxon_species_rows'
        )
    finally:
        con.unregister('_taxon_species_rows')
    con.execute(
        """
        UPDATE entity_evidence_raw AS ee
        SET taxonomy_id = tsl.species_tax_id
        FROM taxon_species_lookup tsl
        WHERE ee.taxonomy_id = tsl.tax_id
        """
    )
    _log.info(
        'taxon canonicalization: %d strain taxa mapped to species level',
        len(rows),
    )


def _build_resolver_lookup(
    con: duckdb.DuckDBPyConnection,
) -> None:
    """Replace the full live-utils resolver VIEWS with per-shard keyed tables.

    No-op unless the utils Postgres is attached. Reads the materialized resolvers
    (resolver_gene -> entrez per taxon; resolver_gene_protein_global -> entrez
    taxon-agnostic for no-taxon proteins; resolver_protein -> uniprot fallback;
    resolver_chemical -> inchikey) by keyed lookup and rebuilds the small local
    ``resolver_lookup`` / ``protein_uniprot_fallback_lookup`` tables that
    canonicalize joins against.
    """
    preflight_translation_tables(con)
    url = os.environ.get('OMNIPATH_BUILD_UTILS_PG_URL')
    if not url or not _live_utils_attached(con):
        return

    protein_gene_keys = _live_utils_key_rows(
        con,
        mapping=RESOLVER_PROTEIN_SLUG_TO_IDENTIFIER_TYPE,
        evidence_entity_types=(PROTEIN_ENTITY_TYPE, GENE_ENTITY_TYPE),
    )
    resolver_rows: list[tuple[str, int, str, str | None, int, str]] = []
    resolver_rows.extend(
        _fetch_live_utils_rows_for_keys(
            url=url,
            table='resolver_gene',
            mapping=RESOLVER_PROTEIN_SLUG_TO_IDENTIFIER_TYPE,
            entity_type=GENE_ENTITY_TYPE,
            canonical_identifier_type_id=identifier_type_id(ENTREZ_TYPE),
            canonical_column='entrez',
            key_rows=protein_gene_keys,
        )
    )
    # Global gene anchor for proteins that arrive with NO taxon
    # (Rhea/Brenda/TCDB/ChEMBL reference UniProt without an organism): the
    # per-taxon resolver_gene only materializes a subset of organisms, so a
    # UniProt from an unmaterialized organism has no resolver_gene row at all.
    # resolver_gene_protein_global covers uniprot->entrez across ALL taxa. It
    # DOES carry the gene's real ncbi_tax_id (entrez->taxon is 1:1), so fetch
    # WITH taxonomy (has_taxonomy=True): the no-taxon evidence still matches
    # (the keyed join keeps rows where the evidence taxon IS NULL), and the
    # resolved gene gets its concrete taxon. Fetching it taxon-agnostically
    # (NULL) instead created a second, NULL-taxon copy of every gene, colliding
    # with the real- taxon copy of the same entrez and breaking the
    # entity_evidence_resolution primary key at derive time.
    resolver_rows.extend(
        _fetch_live_utils_rows_for_keys(
            url=url,
            table='resolver_gene_protein_global',
            mapping=RESOLVER_PROTEIN_SLUG_TO_IDENTIFIER_TYPE,
            entity_type=GENE_ENTITY_TYPE,
            canonical_identifier_type_id=identifier_type_id(ENTREZ_TYPE),
            canonical_column='entrez',
            has_taxonomy=True,
            key_rows=protein_gene_keys,
        )
    )
    resolver_rows.extend(
        _fetch_live_utils_rows_for_keys(
            url=url,
            table='resolver_chemical',
            mapping=RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE,
            entity_type=CHEMICAL_ENTITY_TYPE,
            canonical_identifier_type_id=identifier_type_id(
                STANDARD_INCHI_KEY_TYPE
            ),
            canonical_column='inchikey',
            has_taxonomy=False,
            key_rows=_live_utils_key_rows(
                con,
                mapping=RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE,
                evidence_entity_types=(CHEMICAL_ENTITY_TYPE,),
            ),
        )
    )
    fallback_rows = _fetch_live_utils_rows_for_keys(
        url=url,
        table='resolver_protein',
        mapping=RESOLVER_PROTEIN_SLUG_TO_IDENTIFIER_TYPE,
        entity_type=PROTEIN_ENTITY_TYPE,
        canonical_identifier_type_id=identifier_type_id(UNIPROT_TYPE),
        canonical_column='uniprot',
        key_rows=_live_utils_key_rows(
            con,
            mapping=RESOLVER_PROTEIN_SLUG_TO_IDENTIFIER_TYPE,
            evidence_entity_types=(PROTEIN_ENTITY_TYPE,),
        ),
    )

    for name in ('resolver_lookup', 'protein_uniprot_fallback_lookup'):
        con.execute(f'DROP VIEW IF EXISTS {name}')
        con.execute(f'DROP TABLE IF EXISTS {name}')
    schema_sql = """
        (
          entity_type VARCHAR,
          key_identifier_type_id BIGINT,
          key_value VARCHAR,
          taxonomy_id VARCHAR,
          canonical_identifier_type_id BIGINT,
          canonical_identifier VARCHAR
        )
    """
    con.execute(f'CREATE TABLE resolver_lookup {schema_sql}')
    con.execute(f'CREATE TABLE protein_uniprot_fallback_lookup {schema_sql}')
    _bulk_insert_resolver_lookup_rows(con, 'resolver_lookup', resolver_rows)
    _bulk_insert_resolver_lookup_rows(
        con, 'protein_uniprot_fallback_lookup', fallback_rows
    )


def _append_live_utils_ontology_endpoint_resolver_rows(
    con: duckdb.DuckDBPyConnection,
) -> None:
    """Add ontology-endpoint keys to the keyed ``resolver_lookup``.

    Ontology relation endpoints (``ontology_endpoint_identifier_key``) can
    reference ids with NO ``entity_evidence`` row (``subject_entity_evidence_id
    IS NULL``), so they are absent from ``evidence_identifier_key`` and thus from
    the per-shard keyed ``resolver_lookup`` built by
    ``_build_resolver_lookup``. The ontology prefilter
    (``needed_ontology_endpoint_resolver_lookup``) reads ``resolver_lookup``
    ``WHERE taxonomy_id IS NULL`` — i.e. only the organism-agnostic **chemical**
    resolver rows — so fetch the chemical resolver rows for those endpoint keys
    and append them, keeping ontology resolution at full coverage under the keyed
    read (the full-view path saw every row for free). No-op unless utils attached.
    """
    url = os.environ.get('OMNIPATH_BUILD_UTILS_PG_URL')
    if not url or not _live_utils_attached(con):
        return

    chem_id_to_slugs: dict[int, list[str]] = {}
    for slug, idt in RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE.items():
        chem_id_to_slugs.setdefault(identifier_type_id(idt), []).append(slug)
    placeholders = ', '.join('?' for _ in chem_id_to_slugs)
    endpoint_rows = con.execute(
        f"""
        SELECT DISTINCT key_identifier_type_id, key_value
        FROM ontology_endpoint_identifier_key
        WHERE key_identifier_type_id IN ({placeholders})
          AND key_value IS NOT NULL
          AND key_value <> ''
        """,
        list(chem_id_to_slugs.keys()),
    ).fetchall()

    key_rows: list[tuple[str, str, int | None]] = []
    for key_identifier_type_id, key_value in endpoint_rows:
        for slug in chem_id_to_slugs[int(key_identifier_type_id)]:
            key_rows.append((slug, str(key_value), None))

    new_rows = _fetch_live_utils_rows_for_keys(
        url=url,
        table='resolver_chemical',
        mapping=RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE,
        entity_type=CHEMICAL_ENTITY_TYPE,
        canonical_identifier_type_id=identifier_type_id(STANDARD_INCHI_KEY_TYPE),
        canonical_column='inchikey',
        has_taxonomy=False,
        key_rows=key_rows,
    )
    _bulk_insert_resolver_lookup_rows(con, 'resolver_lookup', new_rows)


def _append_live_utils_annotation_object_resolver_rows(
    con: duckdb.DuckDBPyConnection,
) -> None:
    """Add annotation-relation chemical-object keys to the keyed ``resolver_lookup``.

    An annotation relation's object (``annotation_relation_evidence_raw``) is a
    bare external identifier, never an ``entity_evidence`` mention -- so, like
    ontology relation endpoints just above, it can be absent from the per-shard
    keyed ``resolver_lookup`` (``_build_resolver_lookup``, which is keyed off
    ``entity_evidence``/``entity_identifier_raw``). Fetch the chemical resolver
    rows for those object keys and append them.

    This is the fix for the previously-known gap where a chemical annotation
    object (e.g. a reaction annotated with a bare ChEBI id, not a full
    ``entity_evidence`` mention) self-typed by its raw id and minted its own
    entity, permanently separate from whatever entity the same real molecule
    resolves to elsewhere via the normal mention pipeline. Scoped to chemical
    objects only (the documented scope of the gap) -- no-op unless utils
    attached.
    """
    url = os.environ.get('OMNIPATH_BUILD_UTILS_PG_URL')
    if not url or not _live_utils_attached(con):
        return

    chem_id_to_slugs: dict[int, list[str]] = {}
    for slug, idt in RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE.items():
        chem_id_to_slugs.setdefault(identifier_type_id(idt), []).append(slug)
    placeholders = ', '.join('?' for _ in chem_id_to_slugs)
    object_rows = con.execute(
        f"""
        SELECT DISTINCT object_type.identifier_type_id, ar.object_id
        FROM annotation_relation_evidence_raw ar
        JOIN identifier_type_all object_type
          ON object_type.name = ar.object_id_type
        WHERE ar.object_entity_type = ?
          AND object_type.identifier_type_id IN ({placeholders})
          AND ar.object_id IS NOT NULL
          AND ar.object_id <> ''
        """,
        [CHEMICAL_ENTITY_TYPE] + list(chem_id_to_slugs.keys()),
    ).fetchall()

    key_rows: list[tuple[str, str, int | None]] = []
    for key_identifier_type_id, key_value in object_rows:
        for slug in chem_id_to_slugs[int(key_identifier_type_id)]:
            key_rows.append((slug, str(key_value), None))

    new_rows = _fetch_live_utils_rows_for_keys(
        url=url,
        table='resolver_chemical',
        mapping=RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE,
        entity_type=CHEMICAL_ENTITY_TYPE,
        canonical_identifier_type_id=identifier_type_id(STANDARD_INCHI_KEY_TYPE),
        canonical_column='inchikey',
        has_taxonomy=False,
        key_rows=key_rows,
    )
    _bulk_insert_resolver_lookup_rows(con, 'resolver_lookup', new_rows)


def _build_annotation_object_resolution(con: duckdb.DuckDBPyConnection) -> None:
    """Resolve every distinct annotation-relation object key once, as a
    persisted table -- shared by ``annotation_object_entity`` (which entity
    it mints/joins), the entity ``sources`` attribution, and
    ``annotation_projected`` (which entity a relation's object edge points
    at), so all three agree on the same promoted identity instead of drifting
    if only one were fixed.

    A key promotes to its resolver hit when ``resolver_lookup`` has one
    (chemical keys only, see ``_append_live_utils_annotation_object_resolver_rows``
    just above); anything else -- non-chemical object types, or a chemical id
    with no resolver coverage -- falls through to the previous raw self-typed
    behaviour via ``coalesce``, unchanged.
    """
    con.execute(
        """
        CREATE TABLE annotation_object_resolution AS
        SELECT DISTINCT
          ar.object_entity_type AS entity_type,
          object_type.identifier_type_id AS raw_identifier_type_id,
          ar.object_id AS raw_identifier,
          coalesce(rl.canonical_identifier_type_id, object_type.identifier_type_id)
            AS canonical_identifier_type_id,
          coalesce(rl.canonical_identifier, ar.object_id) AS canonical_identifier
        FROM annotation_relation_evidence_raw ar
        JOIN identifier_type_all object_type
          ON object_type.name = ar.object_id_type
        LEFT JOIN resolver_lookup rl
          ON rl.entity_type = ar.object_entity_type
         AND rl.key_identifier_type_id = object_type.identifier_type_id
         AND rl.key_value = ar.object_id
         AND rl.taxonomy_id IS NULL
        WHERE ar.object_id IS NOT NULL AND ar.object_id <> ''
        QUALIFY row_number() OVER (
          PARTITION BY
            ar.object_entity_type, object_type.identifier_type_id, ar.object_id
          ORDER BY
            rl.canonical_identifier IS NULL,
            rl.canonical_identifier_type_id,
            rl.canonical_identifier
        ) = 1
        """
    )


def _canonicalize_loaded_duckdb(
    con: duckdb.DuckDBPyConnection,
) -> tuple[int, int, int]:
    """Canonicalize already-loaded DuckDB evidence and resolver tables."""

    _create_duckdb_identifier_type_all_view(con)
    _ensure_duckdb_canonical_caches(con)
    _drop_duckdb_batch_tables(con)
    # Fill in the organism of protein mentions that cited a UniProt accession but
    # no species, from the accession itself — before the resolver keys are built,
    # so the recovered organism drives both gene resolution and the accession-keyed
    # gene mint (which requires a known organism).
    _recover_no_taxon_protein_from_uniprot(con)
    # Canonicalize each mention's organism to species level (strain -> species),
    # after recovery so a recovered strain taxon is canonicalized too. Keeps
    # mentions and the resolver on one taxon per gene, so a gene is not split
    # across NCBI strain-vs-species taxa.
    _canonicalize_taxon_to_species(con)
    # WP2 skeleton-collapse/authority-arbitration (T056-T058) needs to know
    # which identifier types are declared structure authorities; fetch it
    # once per shard, before entity_resolution_base is built.
    _fetch_identifier_authority(con)
    resolver_alias_expansion_excluded_type_ids = ', '.join(
        str(identifier_type_id(name))
        for name in RESOLVER_ALIAS_EXPANSION_EXCLUDED_IDENTIFIER_TYPES
    )
    con.execute(
        """
        CREATE TABLE evidence_identifier_key AS
        SELECT DISTINCT
          ee.entity_type,
          ee.taxonomy_id,
          kit.identifier_type_id AS key_identifier_type_id,
          COALESCE(ei.identifier_normalized, ei.identifier) AS key_value
        FROM entity_evidence_raw ee
        JOIN entity_identifier_raw ei
          ON ei.source = ee.source
         AND ei.entity_evidence_id = ee.entity_evidence_id
        JOIN identifier_type_all kit
          ON kit.name = ei.identifier_type
        WHERE ei.identifier IS NOT NULL
          AND ei.identifier <> ''
        """
    )
    # Keyed lookup: now that the shard's evidence keys exist, replace the full
    # live-utils resolver VIEWS with small per-shard tables fetched by key (only
    # this shard's ids scan the remote resolver, not the whole 500M-row tables).
    _build_resolver_lookup(con)
    con.execute(
        """
        CREATE TABLE resolver_entity_type_match AS
        SELECT ? AS evidence_entity_type, ? AS resolver_entity_type
        UNION ALL
        SELECT ?, ?
        UNION ALL
        SELECT ?, ?
        UNION ALL
        SELECT ?, ?
        """,
        # Protein evidence first tries the gene-product resolver's gene anchor
        # rows, but can also resolve to protein rows when Entrez is unavailable.
        [
            PROTEIN_ENTITY_TYPE,
            GENE_ENTITY_TYPE,
            PROTEIN_ENTITY_TYPE,
            PROTEIN_ENTITY_TYPE,
            GENE_ENTITY_TYPE,
            GENE_ENTITY_TYPE,
            CHEMICAL_ENTITY_TYPE,
            CHEMICAL_ENTITY_TYPE,
        ],
    )
    con.execute(
        """
        CREATE TABLE resolver_evidence_identifier_key AS
        SELECT DISTINCT
          k.entity_type AS evidence_entity_type,
          etm.resolver_entity_type,
          k.taxonomy_id,
          k.key_identifier_type_id,
          k.key_value
        FROM evidence_identifier_key k
        JOIN resolver_entity_type_match etm
          ON etm.evidence_entity_type = k.entity_type
        """
    )
    con.execute(
        """
        CREATE TABLE taxonomy_optional_resolver_key_type AS
        SELECT identifier_type_id
        FROM identifier_type_all
        WHERE name IN ({})
        """.format(
            ', '.join(
                '?' for _ in PROTEIN_TAXONOMY_OPTIONAL_IDENTIFIER_TYPES
            )
        ),
        list(PROTEIN_TAXONOMY_OPTIONAL_IDENTIFIER_TYPES),
    )
    has_resolver_keys = bool(
        con.execute(
            'SELECT count(*) > 0 FROM resolver_evidence_identifier_key'
        ).fetchone()[0]
    )
    if has_resolver_keys:
        con.execute(
            """
            CREATE TABLE taxonomy_optional_unambiguous_key AS
            SELECT
              rl.key_identifier_type_id,
              rl.key_value,
              rl.canonical_identifier_type_id
            FROM resolver_lookup rl
            JOIN taxonomy_optional_resolver_key_type opt
              ON opt.identifier_type_id = rl.key_identifier_type_id
            JOIN resolver_evidence_identifier_key k
              ON k.evidence_entity_type = ?
             AND k.resolver_entity_type = rl.entity_type
             AND k.key_identifier_type_id = rl.key_identifier_type_id
             AND k.key_value = rl.key_value
            WHERE rl.entity_type = ?
            GROUP BY
              rl.key_identifier_type_id,
              rl.key_value,
              rl.canonical_identifier_type_id
            HAVING count(DISTINCT rl.canonical_identifier) = 1
            """,
            [PROTEIN_ENTITY_TYPE, GENE_ENTITY_TYPE],
        )
        con.execute(
            """
            CREATE TABLE needed_resolver_lookup AS
            SELECT DISTINCT
              k.evidence_entity_type,
              rl.*,
              opt.key_identifier_type_id IS NOT NULL AS taxonomy_optional_match
            FROM resolver_lookup rl
            JOIN resolver_evidence_identifier_key k
              ON k.resolver_entity_type = rl.entity_type
             AND k.key_identifier_type_id = rl.key_identifier_type_id
             AND k.key_value = rl.key_value
            LEFT JOIN taxonomy_optional_unambiguous_key opt
              ON opt.key_identifier_type_id = rl.key_identifier_type_id
             AND opt.key_value = rl.key_value
             AND opt.canonical_identifier_type_id =
                 rl.canonical_identifier_type_id
            WHERE
              rl.taxonomy_id = k.taxonomy_id
              OR rl.taxonomy_id IS NULL
              OR opt.key_identifier_type_id IS NOT NULL
            """
        )
    else:
        con.execute(
            """
            CREATE TABLE taxonomy_optional_unambiguous_key AS
            SELECT
              NULL::BIGINT AS key_identifier_type_id,
              NULL::VARCHAR AS key_value,
              NULL::BIGINT AS canonical_identifier_type_id
            WHERE false
            """
        )
        con.execute(
            """
            CREATE TABLE needed_resolver_lookup AS
            SELECT
              NULL::VARCHAR AS evidence_entity_type,
              NULL::VARCHAR AS entity_type,
              NULL::BIGINT AS key_identifier_type_id,
              NULL::VARCHAR AS key_value,
              NULL::VARCHAR AS taxonomy_id,
              NULL::BIGINT AS canonical_identifier_type_id,
              NULL::VARCHAR AS canonical_identifier,
              NULL::BOOLEAN AS taxonomy_optional_match
            WHERE false
            """
        )
    con.execute(
        """
        CREATE TABLE ontology_endpoint_identifier_key AS
        SELECT DISTINCT
          endpoint.entity_type AS evidence_entity_type,
          etm.resolver_entity_type,
          kit.identifier_type_id AS key_identifier_type_id,
          endpoint.identifier AS key_value
        FROM (
          SELECT
            subject_entity_type AS entity_type,
            subject_identifier_type AS identifier_type,
            subject_identifier AS identifier
          FROM ontology_relation_raw
          WHERE subject_entity_evidence_id IS NULL
            AND subject_identifier IS NOT NULL
            AND subject_identifier <> ''
          UNION
          SELECT
            object_entity_type AS entity_type,
            object_identifier_type AS identifier_type,
            object_identifier AS identifier
          FROM ontology_relation_raw
          WHERE object_identifier IS NOT NULL
            AND object_identifier <> ''
        ) endpoint
        JOIN resolver_entity_type_match etm
          ON etm.evidence_entity_type = endpoint.entity_type
        JOIN identifier_type_all kit
          ON kit.name = endpoint.identifier_type
        WHERE endpoint.entity_type IS NOT NULL
        """
    )
    # Keyed read only fetched entity-evidence keys; ontology endpoints may carry
    # ids absent from that set — append their (chemical) resolver rows so the
    # prefilter below keeps full coverage.
    _append_live_utils_ontology_endpoint_resolver_rows(con)
    has_ontology_endpoint_resolver_keys = bool(
        con.execute(
            'SELECT count(*) > 0 FROM ontology_endpoint_identifier_key'
        ).fetchone()[0]
    )
    if has_ontology_endpoint_resolver_keys:
        con.execute(
            """
            CREATE TABLE needed_ontology_endpoint_resolver_lookup AS
            SELECT DISTINCT
              k.evidence_entity_type,
              rl.*
            FROM resolver_lookup rl
            JOIN ontology_endpoint_identifier_key k
              ON k.resolver_entity_type = rl.entity_type
             AND k.key_identifier_type_id = rl.key_identifier_type_id
             AND k.key_value = rl.key_value
            WHERE rl.taxonomy_id IS NULL
            """
        )
    else:
        con.execute(
            """
            CREATE TABLE needed_ontology_endpoint_resolver_lookup AS
            SELECT
              NULL::VARCHAR AS evidence_entity_type,
              NULL::VARCHAR AS entity_type,
              NULL::BIGINT AS key_identifier_type_id,
              NULL::VARCHAR AS key_value,
              NULL::VARCHAR AS taxonomy_id,
              NULL::BIGINT AS canonical_identifier_type_id,
              NULL::VARCHAR AS canonical_identifier
            WHERE false
            """
        )
    # Annotation-relation chemical objects promote through the same keyed
    # resolver_lookup, for the same reason ontology endpoints just above do --
    # built once here so batch_entity_candidate (below) and the later
    # annotation_projected relation-building step both read one agreed answer.
    _append_live_utils_annotation_object_resolver_rows(con)
    _build_annotation_object_resolution(con)
    con.execute(
        """
        CREATE TABLE protein_uniprot_fallback_taxonomy_optional_unambiguous_key AS
        SELECT
          puf.key_identifier_type_id,
          puf.key_value,
          puf.canonical_identifier_type_id
        FROM protein_uniprot_fallback_lookup puf
        JOIN taxonomy_optional_resolver_key_type opt
          ON opt.identifier_type_id = puf.key_identifier_type_id
        JOIN evidence_identifier_key k
          ON k.entity_type = ?
         AND k.key_identifier_type_id = puf.key_identifier_type_id
         AND k.key_value = puf.key_value
        GROUP BY
          puf.key_identifier_type_id,
          puf.key_value,
          puf.canonical_identifier_type_id
        HAVING count(
          DISTINCT coalesce(puf.taxonomy_id, '') || chr(31) ||
          puf.canonical_identifier
        ) = 1
        """,
        [PROTEIN_ENTITY_TYPE],
    )
    con.execute(
        """
        CREATE TABLE needed_protein_uniprot_fallback_lookup AS
        SELECT DISTINCT
          puf.*,
          opt.key_identifier_type_id IS NOT NULL AS taxonomy_optional_match
        FROM protein_uniprot_fallback_lookup puf
        LEFT JOIN protein_uniprot_fallback_taxonomy_optional_unambiguous_key opt
          ON opt.key_identifier_type_id = puf.key_identifier_type_id
         AND opt.key_value = puf.key_value
         AND opt.canonical_identifier_type_id =
             puf.canonical_identifier_type_id
        JOIN evidence_identifier_key k
          ON k.entity_type = ?
         AND k.key_identifier_type_id = puf.key_identifier_type_id
         AND k.key_value = puf.key_value
         AND (
           puf.taxonomy_id = k.taxonomy_id
           OR puf.taxonomy_id IS NULL
           OR opt.key_identifier_type_id IS NOT NULL
         )
        """,
        [PROTEIN_ENTITY_TYPE],
    )
    # Multi-gene protein split: a UniProt mapping to >1 gene is duplicated 1:1
    # per gene *before* resolution, so the existing 1:1 machinery resolves each
    # copy to its gene. Rewrites the raw tables in place + leaves
    # `multigene_resolution` for entity_resolution_base's direct-resolution arm.
    explode_multi_gene_protein_mentions(con)
    # Non-lipid chemical fallback: best non-structure id per chemical mention,
    # consumed by entity_resolution_base's unresolved branch below. Stage 2
    # anchor map (1:1 id->structure/ChEBI) feeds the fallback's translated
    # candidates, so it must be built first.
    build_chemical_anchor_map(con)
    build_chemical_fallback_resolution(con)
    con.execute(
        """
        CREATE TABLE entity_identifier_group AS
        SELECT
          ee.source,
          ee.entity_evidence_id,
          coalesce(
            string_agg(
              kit.identifier_type_id::VARCHAR || '=' || ei.identifier,
              '|'
              ORDER BY kit.identifier_type_id, ei.identifier
            ),
            ee.source || ':' || ee.entity_evidence_id
          ) AS unresolved_identifier_key
        FROM entity_evidence_raw ee
        LEFT JOIN entity_identifier_raw ei
          ON ei.source = ee.source
         AND ei.entity_evidence_id = ee.entity_evidence_id
         AND ei.identifier IS NOT NULL
         AND ei.identifier <> ''
        LEFT JOIN identifier_type_all kit
          ON kit.name = ei.identifier_type
        GROUP BY ee.source, ee.entity_evidence_id
        """
    )
    # The single highest-priority raw identifier per evidence record, used as the
    # canonical identifier for records that survive every resolution path (in
    # place of a hash of the whole identifier bag). Deterministic: rank by the
    # best-identifier priority, then by identifier-type-id, then by the identifier
    # string. Records with no usable identifier get no row here and fall back to
    # the per-record composite key carried by entity_identifier_group.
    unresolved_priority_case = '\n'.join(
        [
            '            CASE',
            *(
                f'              WHEN ei.identifier_type = '
                f'{_sql_literal(identifier_type)} THEN {rank}'
                for rank, identifier_type in enumerate(
                    UNRESOLVED_KEY_IDENTIFIER_PRIORITY
                )
            ),
            f'              ELSE {len(UNRESOLVED_KEY_IDENTIFIER_PRIORITY)}',
            '            END',
        ]
    )
    con.execute(
        f"""
        CREATE TABLE entity_best_identifier AS
        SELECT
          source,
          entity_evidence_id,
          best_identifier
        FROM (
          SELECT
            ee.source,
            ee.entity_evidence_id,
            ei.identifier AS best_identifier,
{unresolved_priority_case} AS priority,
            coalesce(kit.identifier_type_id, 2147483647) AS identifier_type_id
          FROM entity_evidence_raw ee
          JOIN entity_identifier_raw ei
            ON ei.source = ee.source
           AND ei.entity_evidence_id = ee.entity_evidence_id
           AND ei.identifier IS NOT NULL
           AND ei.identifier <> ''
          LEFT JOIN identifier_type_all kit
            ON kit.name = ei.identifier_type
        )
        QUALIFY row_number() OVER (
          PARTITION BY source, entity_evidence_id
          ORDER BY priority, identifier_type_id, best_identifier
        ) = 1
        """
    )
    con.execute(
        """
        CREATE TABLE cv_term_evidence_resolution AS
        SELECT
          ee.source,
          ee.entity_evidence_id,
          cv_type.identifier_type_id AS canonical_identifier_type_id,
          min(ei.identifier) AS canonical_identifier
        FROM entity_evidence_raw ee
        JOIN entity_identifier_raw ei
          ON ei.source = ee.source
         AND ei.entity_evidence_id = ee.entity_evidence_id
        CROSS JOIN (
          SELECT identifier_type_id
          FROM identifier_type_all
          WHERE name = ?
        ) cv_type
        WHERE ee.entity_type = ?
          AND ei.identifier_type = ?
          AND ei.identifier IS NOT NULL
          AND ei.identifier <> ''
        GROUP BY
          ee.source,
          ee.entity_evidence_id,
          cv_type.identifier_type_id
        """,
        [CV_TERM_ID_TYPE, CV_TERM_ENTITY_TYPE, CV_TERM_ID_TYPE],
    )
    con.execute(
        """
        CREATE TABLE pathway_identifier_evidence_resolution AS
        WITH candidate AS (
          SELECT
            ee.source,
            ee.entity_evidence_id,
            kit.identifier_type_id AS canonical_identifier_type_id,
            min(ei.identifier) AS canonical_identifier,
            CASE
              WHEN ei.identifier_type = ? THEN 1
              WHEN ei.identifier_type = ? THEN 2
              ELSE 100
            END AS priority
          FROM entity_evidence_raw ee
          JOIN entity_identifier_raw ei
            ON ei.source = ee.source
           AND ei.entity_evidence_id = ee.entity_evidence_id
          JOIN identifier_type_all kit
            ON kit.name = ei.identifier_type
          WHERE ee.entity_type = ?
            AND ei.identifier_type IN (?, ?)
            AND ei.identifier IS NOT NULL
            AND ei.identifier <> ''
          GROUP BY
            ee.source,
            ee.entity_evidence_id,
            kit.identifier_type_id,
            ei.identifier_type
          HAVING count(DISTINCT ei.identifier) = 1
        )
        SELECT
          source,
          entity_evidence_id,
          canonical_identifier_type_id,
          canonical_identifier
        FROM candidate
        QUALIFY row_number() OVER (
          PARTITION BY source, entity_evidence_id
          ORDER BY priority, canonical_identifier
        ) = 1
        """,
        [
            REACTOME_STABLE_ID_TYPE,
            WIKIPATHWAYS_ID_TYPE,
            PATHWAY_ENTITY_TYPE,
            REACTOME_STABLE_ID_TYPE,
            WIKIPATHWAYS_ID_TYPE,
        ],
    )
    # A lipid mention canonicalizes on its standardized nomenclature name,
    # never on structure (research R9): "any InChIKey, ChEBI, LIPID MAPS or
    # SwissLipids identifier is attached as an ordinary identifier" -- so
    # this must be resolved and threaded into direct_resolution's coalesce
    # BEFORE std_inchi_key below, not after. Degrades cleanly (T112/R14):
    # an empty table, not a failure, when the utils Postgres isn't attached
    # or its lipid_name table doesn't exist (an older utils build, or one
    # without the `lipid` extra).
    lipid_name_available = _duckdb_attach_utils_postgres(con) and (
        _attached_utils_relation_exists(con, 'lipid_name')
    )
    if lipid_name_available:
        con.execute(
            f"""
            CREATE TABLE lipid_name_evidence_resolution AS
            WITH candidate AS (
              SELECT
                ee.source,
                ee.entity_evidence_id,
                ln.lipid_name,
                ln.lipid_level,
                ln.chains_listed,
                ln.chains_possible,
                -- a mention may carry more than one name-axis identifier
                -- (a primary name and a synonym, say); the finer-grained
                -- match wins, mirroring omnipath-metabo's own
                -- "highest-granularity parse wins" rule.
                row_number() OVER (
                  PARTITION BY ee.source, ee.entity_evidence_id
                  ORDER BY
                    CASE ln.lipid_level
                      WHEN 'complete_structure' THEN 8
                      WHEN 'full_structure' THEN 7
                      WHEN 'structure_defined' THEN 6
                      WHEN 'sn_position' THEN 5
                      WHEN 'molecular_species' THEN 4
                      WHEN 'partially_specified' THEN 3
                      WHEN 'species' THEN 2
                      ELSE 1
                    END DESC
                ) AS rk
              FROM entity_evidence_raw ee
              JOIN entity_identifier_raw ei
                ON ei.source = ee.source
               AND ei.entity_evidence_id = ee.entity_evidence_id
              JOIN utils_pg.omnipath_utils.lipid_name ln
                ON ln.raw_name = ei.identifier
              WHERE ee.entity_type = ?
                AND ei.identifier_type IN ({
                    ','.join('?' for _ in LIPID_NAME_SOURCE_TYPES)
                })
            )
            SELECT
              source,
              entity_evidence_id,
              lname_type.identifier_type_id AS canonical_identifier_type_id,
              -- the composite key (lipid_name, lipid_level, chains_listed,
              -- chains_possible) encoded as one string: lipid_name alone is
              -- not unique across levels (e.g. a one-chain class's
              -- full_structure and complete_structure renderings can be
              -- textually identical, verified live on this build's own
              -- corpus -- 'CAR 11:0' is one of 2,828 such collisions).
              lipid_name || '|' || lipid_level || '|' || chains_listed || '|'
                || chains_possible AS canonical_identifier
            FROM candidate
            CROSS JOIN (
              SELECT identifier_type_id FROM identifier_type_all
              WHERE name = ?
            ) lname_type
            WHERE rk = 1
            """,
            [CHEMICAL_ENTITY_TYPE, *LIPID_NAME_SOURCE_TYPES, LIPID_NAME_TYPE],
        )
    else:
        con.execute(
            """
            CREATE TABLE lipid_name_evidence_resolution AS
            SELECT
              NULL::VARCHAR AS source,
              NULL::VARCHAR AS entity_evidence_id,
              NULL::BIGINT AS canonical_identifier_type_id,
              NULL::VARCHAR AS canonical_identifier
            WHERE false
            """
        )
    con.execute(
        """
        CREATE TABLE standard_inchi_key_evidence_resolution AS
        SELECT
          ee.source,
          ee.entity_evidence_id,
          std_type.identifier_type_id AS canonical_identifier_type_id,
          min(ei.identifier) AS canonical_identifier
        FROM entity_evidence_raw ee
        JOIN entity_identifier_raw ei
          ON ei.source = ee.source
         AND ei.entity_evidence_id = ee.entity_evidence_id
        CROSS JOIN (
          SELECT identifier_type_id
          FROM identifier_type_all
          WHERE name = ?
        ) std_type
        WHERE ee.entity_type = ?
          AND ei.identifier_type = ?
          AND regexp_matches(
            ei.identifier,
            '^[A-Z]{14}-[A-Z]{10}-[A-Z]$'
          )
        GROUP BY
          ee.source,
          ee.entity_evidence_id,
          std_type.identifier_type_id
        HAVING count(DISTINCT ei.identifier) = 1
        """,
        [
            STANDARD_INCHI_KEY_TYPE,
            CHEMICAL_ENTITY_TYPE,
            STANDARD_INCHI_KEY_TYPE,
        ],
    )
    cf_fires = chemical_fallback_fires_sql()
    protein_fallback_fires = (
        f"ee.entity_type = {_sql_literal(PROTEIN_ENTITY_TYPE)} "
        "AND (rcs.candidate_count IS NULL OR rcs.candidate_count = 0) "
        "AND puf.candidate_count = 1"
    )
    # A protein mention with a well-formed UniProt accession (and a real taxon),
    # for which the resolver produced neither a gene nor a canonicalisable primary
    # UniProt, still names a gene product — type it as a gene keyed by that
    # accession. Requires exactly one distinct accession so the key is unambiguous;
    # mutually exclusive with the resolver and fallback arms by construction.
    protein_uniprot_selftype_fires = (
        f'ee.entity_type = {_sql_literal(PROTEIN_ENTITY_TYPE)} '
        'AND (rcs.candidate_count IS NULL OR rcs.candidate_count = 0) '
        'AND (puf.candidate_count IS NULL OR puf.candidate_count <> 1) '
        'AND st.accession_count = 1'
    )
    con.execute(
        f"""
        CREATE TABLE entity_resolution_base AS
        WITH direct_resolution AS (
          SELECT
            ee.source,
            ee.dataset,
            ee.row_id,
            ee.entity_evidence_id,
            ee.entity_type,
            ee.entity_type AS molecular_entity_type,
            CASE
              WHEN cv_term.canonical_identifier IS NOT NULL THEN NULL
              WHEN pathway_identifier.canonical_identifier IS NOT NULL THEN NULL
              ELSE ee.taxonomy_id
            END AS taxonomy_id,
            coalesce(
              cv_term.canonical_identifier_type_id,
              pathway_identifier.canonical_identifier_type_id,
              lipid_name.canonical_identifier_type_id,
              std_inchi_key.canonical_identifier_type_id
            ) AS canonical_identifier_type_id,
            coalesce(
              cv_term.canonical_identifier,
              pathway_identifier.canonical_identifier,
              lipid_name.canonical_identifier,
              std_inchi_key.canonical_identifier
            ) AS canonical_identifier,
            'resolved' AS status,
            CASE
              WHEN cv_term.canonical_identifier IS NOT NULL THEN 'cv_term'
              WHEN pathway_identifier.canonical_identifier IS NOT NULL
                THEN 'pathway'
              WHEN lipid_name.canonical_identifier IS NOT NULL THEN 'lipid_name'
              ELSE 'inchikey'
            END AS resolution_mechanism
          FROM entity_evidence_raw ee
          LEFT JOIN cv_term_evidence_resolution cv_term
            ON cv_term.source = ee.source
           AND cv_term.entity_evidence_id = ee.entity_evidence_id
          LEFT JOIN pathway_identifier_evidence_resolution pathway_identifier
            ON pathway_identifier.source = ee.source
           AND pathway_identifier.entity_evidence_id = ee.entity_evidence_id
          LEFT JOIN lipid_name_evidence_resolution lipid_name
            ON lipid_name.source = ee.source
           AND lipid_name.entity_evidence_id = ee.entity_evidence_id
          LEFT JOIN standard_inchi_key_evidence_resolution std_inchi_key
            ON std_inchi_key.source = ee.source
           AND std_inchi_key.entity_evidence_id = ee.entity_evidence_id
          WHERE coalesce(
            cv_term.canonical_identifier,
            pathway_identifier.canonical_identifier,
            lipid_name.canonical_identifier,
            std_inchi_key.canonical_identifier
          ) IS NOT NULL
        ),
        remaining_entity AS (
          SELECT ee.*
          FROM entity_evidence_raw ee
          LEFT JOIN direct_resolution direct
            ON direct.source = ee.source
           AND direct.entity_evidence_id = ee.entity_evidence_id
          WHERE direct.entity_evidence_id IS NULL
            -- multi-gene protein copies resolve directly to their assigned
            -- gene below, not through the candidate/unresolved path.
            AND NOT EXISTS (
              SELECT 1 FROM multigene_resolution mr
              WHERE mr.source = ee.source
                AND mr.entity_evidence_id = ee.entity_evidence_id
            )
        ),
        resolver_candidate AS (
          SELECT DISTINCT
            ee.source,
            ee.entity_evidence_id,
            kit.identifier_type_id AS key_identifier_type_id,
            ei.identifier_type AS key_identifier_type_name,
            ei.identifier_key AS key_value,
            rl.entity_type AS resolver_entity_type,
            coalesce(rl.taxonomy_id, ee.taxonomy_id) AS taxonomy_id,
            rl.canonical_identifier_type_id,
            rl.canonical_identifier,
            -- WP2 (spec 011 T056): the connectivity skeleton -- InChIKey
            -- block 1, the first 14 chars -- for a chemical structure
            -- candidate; the full identifier unchanged for anything else
            -- (gene/protein candidates, or a chemical candidate that for
            -- some reason isn't InChIKey-typed). Grouping candidates by this
            -- key, not the full identifier, is what tells a real
            -- cross-skeleton disagreement apart from the same molecule in a
            -- different protonation/tautomer/stereo state.
            CASE
              WHEN rl.canonical_identifier_type_id = {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}
                AND rl.canonical_identifier ~ '^[A-Z]{{14}}-[A-Z]{{10}}-[A-Z]$'
              THEN left(rl.canonical_identifier, 14)
              ELSE rl.canonical_identifier
            END AS skeleton_group_key
          FROM remaining_entity ee
          JOIN (
            SELECT
              source,
              entity_evidence_id,
              identifier_type,
              COALESCE(identifier_normalized, identifier) AS identifier_key
            FROM entity_identifier_raw
          ) ei
            ON ei.source = ee.source
           AND ei.entity_evidence_id = ee.entity_evidence_id
          JOIN identifier_type_all kit
            ON kit.name = ei.identifier_type
          JOIN needed_resolver_lookup rl
            ON rl.key_identifier_type_id = kit.identifier_type_id
           AND rl.key_value = ei.identifier_key
           AND rl.evidence_entity_type = ee.entity_type
           AND (
             rl.taxonomy_id = ee.taxonomy_id
             OR rl.taxonomy_id IS NULL
             OR rl.taxonomy_optional_match
           )
        ),
        resolver_candidate_group_size AS (
          -- how many raw candidates fed each skeleton group -- >1 means real
          -- collapsing happened (T057's resolution_mechanism distinguishes
          -- 'resolver', one candidate all along, from 'skeleton_collapsed').
          SELECT
            source,
            entity_evidence_id,
            skeleton_group_key,
            count(*) AS raw_candidate_count
          FROM resolver_candidate
          GROUP BY source, entity_evidence_id, skeleton_group_key
        ),
        resolver_candidate_group AS (
          -- one row per (mention, skeleton group): the group's own
          -- representative candidate, chosen per T057's fixed, documented
          -- rule -- neutral parent first (the structure key's last
          -- character, 'N'), then the least specified stereo layer (the
          -- InChIKey's own "no stereo/isotope/tautomer defined" hash,
          -- 'UHFFFAOYSA', over one asserting a specific layer some other
          -- candidate in the group disagrees with), then the lowest
          -- identifier as the final deterministic tie-break. Both
          -- conditions are no-ops for a non-InChIKey candidate (gene/
          -- protein, or an ungrouped chemical), which sorts by identifier
          -- alone -- unchanged single-candidate behaviour.
          SELECT
            rc.source,
            rc.entity_evidence_id,
            rc.skeleton_group_key,
            rc.key_identifier_type_id,
            rc.key_identifier_type_name,
            rc.key_value,
            rc.resolver_entity_type,
            rc.taxonomy_id,
            rc.canonical_identifier_type_id,
            rc.canonical_identifier,
            g.raw_candidate_count
          FROM resolver_candidate rc
          JOIN resolver_candidate_group_size g
            ON g.source = rc.source
           AND g.entity_evidence_id = rc.entity_evidence_id
           AND g.skeleton_group_key = rc.skeleton_group_key
          QUALIFY row_number() OVER (
            PARTITION BY rc.source, rc.entity_evidence_id, rc.skeleton_group_key
            ORDER BY
              CASE
                WHEN rc.canonical_identifier_type_id = {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}
                THEN (right(rc.canonical_identifier, 1) <> 'N')
                ELSE false
              END,
              CASE
                WHEN rc.canonical_identifier_type_id = {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}
                THEN (substr(rc.canonical_identifier, 16, 10) <> 'UHFFFAOYSA')
                ELSE false
              END,
              rc.canonical_identifier
          ) = 1
        ),
        resolver_candidate_summary AS (
          -- candidate_count now counts distinct SKELETON GROUPS, not
          -- distinct full identifiers -- several protonation/stereo variants
          -- of one molecule collapse to candidate_count = 1 here.
          SELECT
            source,
            entity_evidence_id,
            count(*) AS candidate_count,
            min(resolver_entity_type) AS resolver_entity_type,
            min(taxonomy_id) AS taxonomy_id,
            min(canonical_identifier_type_id) AS canonical_identifier_type_id,
            min(canonical_identifier) AS canonical_identifier,
            max(raw_candidate_count) AS max_raw_candidate_count
          FROM resolver_candidate_group
          GROUP BY source, entity_evidence_id
        ),
        resolver_candidate_authority AS (
          -- WP2 (spec 011 T058): among a mention's skeleton-group
          -- candidates, the ones whose OWN citing identifier type is a
          -- declared structure authority. Never named a resource in code --
          -- the authority table alone decides.
          SELECT g.*
          FROM resolver_candidate_group g
          JOIN identifier_authority ia
            ON ia.identifier_type_id = g.key_identifier_type_id
           AND ia.is_structure_authority
        ),
        resolver_authority_summary AS (
          -- one row per mention: how many DISTINCT authority-backed skeleton
          -- groups it has. Arbitration only succeeds at exactly one -- two
          -- different authorities backing two different skeletons is still
          -- a genuine, unarbitrated disagreement.
          SELECT
            source,
            entity_evidence_id,
            count(DISTINCT skeleton_group_key) AS authority_candidate_count,
            min(resolver_entity_type) AS resolver_entity_type,
            min(taxonomy_id) AS taxonomy_id,
            min(canonical_identifier_type_id) AS canonical_identifier_type_id,
            min(canonical_identifier) AS canonical_identifier
          FROM resolver_candidate_authority
          GROUP BY source, entity_evidence_id
        ),
        protein_uniprot_fallback_candidate AS (
          SELECT DISTINCT
            ee.source,
            ee.entity_evidence_id,
            coalesce(puf.taxonomy_id, ee.taxonomy_id) AS taxonomy_id,
            puf.canonical_identifier_type_id,
            puf.canonical_identifier
          FROM remaining_entity ee
          JOIN entity_identifier_raw ei
            ON ei.source = ee.source
           AND ei.entity_evidence_id = ee.entity_evidence_id
          JOIN identifier_type_all kit
            ON kit.name = ei.identifier_type
          JOIN needed_protein_uniprot_fallback_lookup puf
            ON puf.key_identifier_type_id = kit.identifier_type_id
           AND puf.key_value = ei.identifier
           AND (
             puf.taxonomy_id = ee.taxonomy_id
             OR puf.taxonomy_id IS NULL
             OR puf.taxonomy_optional_match
           )
          WHERE ee.entity_type = {_sql_literal(PROTEIN_ENTITY_TYPE)}
        ),
        protein_uniprot_fallback_summary AS (
          SELECT
            source,
            entity_evidence_id,
            count(
              DISTINCT coalesce(taxonomy_id, '') || chr(31) ||
              canonical_identifier_type_id::VARCHAR || chr(31) ||
              canonical_identifier
            ) AS candidate_count,
            min(taxonomy_id) AS taxonomy_id,
            min(canonical_identifier_type_id) AS canonical_identifier_type_id,
            min(canonical_identifier) AS canonical_identifier
          FROM protein_uniprot_fallback_candidate
          GROUP BY source, entity_evidence_id
        ),
        protein_uniprot_selftype AS (
          SELECT
            ee.source,
            ee.entity_evidence_id,
            min(ee.taxonomy_id) AS taxonomy_id,
            count(
              DISTINCT regexp_replace(ei.identifier, '-[0-9]+$', '')
            ) AS accession_count,
            min(regexp_replace(ei.identifier, '-[0-9]+$', '')) AS uniprot_ac
          FROM remaining_entity ee
          JOIN entity_identifier_raw ei
            ON ei.source = ee.source
           AND ei.entity_evidence_id = ee.entity_evidence_id
           AND ei.identifier_type = {_sql_literal(UNIPROT_TYPE)}
           AND ei.identifier IS NOT NULL
           AND ei.identifier <> ''
          WHERE ee.entity_type = {_sql_literal(PROTEIN_ENTITY_TYPE)}
            AND ee.taxonomy_id IS NOT NULL
            AND regexp_matches(
              regexp_replace(ei.identifier, '-[0-9]+$', ''),
              {_sql_literal(UNIPROT_AC_REGEX)}
            )
          GROUP BY ee.source, ee.entity_evidence_id
        )
        SELECT * FROM direct_resolution
        UNION ALL
        SELECT
          ee.source,
          ee.dataset,
          ee.row_id,
          ee.entity_evidence_id,
          CASE
            WHEN rcs.candidate_count = 1
            THEN rcs.resolver_entity_type
            -- A gene product the resolver knows only by a UniProt accession, with
            -- no gene, is a gene we cannot name — type it Gene, keyed by that
            -- accession (deduped by it), so the base graph is gene-typed;
            -- molecular_entity_type below stays Protein, so the accession still
            -- becomes a protein `state`, and a later real resolution can merge
            -- this into the known gene. The fallback arm keys by the resolver's
            -- primary UniProt; the self-type arm keys by the accession the source
            -- gave, for proteins whose accession the resolver does not carry.
            WHEN {protein_fallback_fires}
            THEN {_sql_literal(GENE_ENTITY_TYPE)}
            WHEN {protein_uniprot_selftype_fires}
            THEN {_sql_literal(GENE_ENTITY_TYPE)}
            ELSE ee.entity_type
          END AS entity_type,
          ee.entity_type AS molecular_entity_type,
          CASE
            WHEN rcs.candidate_count = 1 THEN rcs.taxonomy_id
            WHEN rcs.candidate_count > 1 AND ras.authority_candidate_count = 1
              THEN ras.taxonomy_id
            WHEN {protein_fallback_fires} THEN puf.taxonomy_id
            WHEN {protein_uniprot_selftype_fires} THEN st.taxonomy_id
            ELSE ee.taxonomy_id
          END AS taxonomy_id,
          -- Chemical fallback: when the structure + resolver-candidate paths
          -- leave a chemical unresolved, take its best non-structure id by
          -- priority (cf) instead of the md5 hash. Gated: cf fires ONLY when
          -- the resolver produced no candidates ({cf_fires}); an ambiguous
          -- candidate_count > 1 that no authority arbitrates (spec 011 T058,
          -- T060) stays unresolved — never a fallback pick over several
          -- distinct structures.
          CASE
            WHEN rcs.candidate_count = 1 THEN rcs.canonical_identifier_type_id
            WHEN rcs.candidate_count > 1 AND ras.authority_candidate_count = 1
              THEN ras.canonical_identifier_type_id
            WHEN {protein_fallback_fires}
              THEN puf.canonical_identifier_type_id
            WHEN {protein_uniprot_selftype_fires}
              THEN {identifier_type_id(UNIPROT_TYPE)}
            WHEN {cf_fires}
              THEN cf.canonical_identifier_type_id
            ELSE unresolved_type.identifier_type_id
          END AS canonical_identifier_type_id,
          CASE
            WHEN rcs.candidate_count = 1 THEN rcs.canonical_identifier
            WHEN rcs.candidate_count > 1 AND ras.authority_candidate_count = 1
              THEN ras.canonical_identifier
            WHEN {protein_fallback_fires} THEN puf.canonical_identifier
            WHEN {protein_uniprot_selftype_fires} THEN st.uniprot_ac
            WHEN {cf_fires} THEN cf.canonical_identifier
            -- An entity that could not be resolved is keyed by its best raw
            -- identifier (deterministic, human-readable, reproducible), so
            -- repeated mentions of the same identifier collapse into one entity.
            -- Records with no usable identifier fall back to the per-record
            -- composite (source:entity_evidence_id).
            ELSE coalesce(ebi.best_identifier, eig.unresolved_identifier_key)
          END AS canonical_identifier,
          CASE
            WHEN rcs.candidate_count = 1 THEN 'resolved'
            WHEN rcs.candidate_count > 1 AND ras.authority_candidate_count = 1
              THEN 'resolved'
            WHEN {protein_fallback_fires} THEN 'resolved'
            WHEN {protein_uniprot_selftype_fires} THEN 'resolved'
            WHEN {cf_fires} THEN 'resolved'
            ELSE 'unresolved'
          END AS status,
          CASE
            -- spec 011 T057: a skeleton that collapsed more than one raw
            -- candidate records that it did, not the plain 'resolver' a
            -- single-candidate mention gets — same status, different
            -- provenance.
            WHEN rcs.candidate_count = 1 AND rcs.max_raw_candidate_count > 1
              THEN 'skeleton_collapsed'
            WHEN rcs.candidate_count = 1 THEN 'resolver'
            WHEN rcs.candidate_count > 1 AND ras.authority_candidate_count = 1
              THEN 'authority_arbitrated'
            WHEN {protein_fallback_fires} THEN 'unknown_gene'
            WHEN {protein_uniprot_selftype_fires} THEN 'unknown_gene'
            WHEN {cf_fires} THEN cf.mechanism
            ELSE 'unresolved'
          END AS resolution_mechanism
        FROM remaining_entity ee
        JOIN entity_identifier_group eig
          ON eig.source = ee.source
         AND eig.entity_evidence_id = ee.entity_evidence_id
        LEFT JOIN entity_best_identifier ebi
          ON ebi.source = ee.source
         AND ebi.entity_evidence_id = ee.entity_evidence_id
        CROSS JOIN (
          SELECT identifier_type_id
          FROM identifier_type_all
          WHERE name = ?
        ) unresolved_type
        LEFT JOIN resolver_candidate_summary rcs
          ON rcs.source = ee.source
         AND rcs.entity_evidence_id = ee.entity_evidence_id
        LEFT JOIN resolver_authority_summary ras
          ON ras.source = ee.source
         AND ras.entity_evidence_id = ee.entity_evidence_id
        LEFT JOIN protein_uniprot_fallback_summary puf
          ON puf.source = ee.source
         AND puf.entity_evidence_id = ee.entity_evidence_id
        LEFT JOIN protein_uniprot_selftype st
          ON st.source = ee.source
         AND st.entity_evidence_id = ee.entity_evidence_id
        LEFT JOIN chemical_fallback_resolution cf
          ON cf.source = ee.source
         AND cf.entity_evidence_id = ee.entity_evidence_id
        UNION ALL
        -- Multi-gene protein copies: resolve each duplicated mention directly
        -- to its assigned gene (entity_type already GENE in
        -- multigene_resolution); molecular_entity_type stays the protein
        -- evidence type so the state layer records the UniProt per gene.
        SELECT
          ee.source,
          ee.dataset,
          ee.row_id,
          ee.entity_evidence_id,
          mr.entity_type AS entity_type,
          ee.entity_type AS molecular_entity_type,
          mr.taxonomy_id,
          mr.canonical_identifier_type_id,
          mr.canonical_identifier,
          'resolved' AS status,
          'gene_anchor' AS resolution_mechanism
        FROM entity_evidence_raw ee
        JOIN multigene_resolution mr
          ON mr.source = ee.source
         AND mr.entity_evidence_id = ee.entity_evidence_id
        """,
        [UNRESOLVED_ID_TYPE],
    )
    # spec 011 T059 (WP2): every unarbitrated chemical disagreement gets a
    # queryable conflict record naming the sources, the identifiers and the
    # competing structures -- scenario 3 requires it never be silently
    # dropped. Scoped to entity_resolution_base's own unresolved chemical
    # rows (a small set), re-deriving each one's skeleton-group candidates
    # directly rather than reusing entity_resolution_base's WITH-scoped CTEs
    # (which do not outlive that one statement) -- same grouping and
    # tie-break rule as resolver_candidate_group in entity_resolution_base
    # (T056/T057); keep the two in sync if either changes.
    con.execute(
        f"""
        CREATE TABLE resolution_conflict AS
        WITH conflicted_mention AS (
          SELECT source, entity_evidence_id
          FROM entity_resolution_base
          WHERE status = 'unresolved'
            AND entity_type = {_sql_literal(CHEMICAL_ENTITY_TYPE)}
        ),
        conflict_candidate AS (
          SELECT
            ee.source,
            ee.entity_evidence_id,
            kit.identifier_type_id AS key_identifier_type_id,
            ei.identifier_type AS key_identifier_type_name,
            COALESCE(ei.identifier_normalized, ei.identifier) AS key_value,
            rl.canonical_identifier_type_id,
            rl.canonical_identifier,
            CASE
              WHEN rl.canonical_identifier_type_id = {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}
                AND rl.canonical_identifier ~ '^[A-Z]{{14}}-[A-Z]{{10}}-[A-Z]$'
              THEN left(rl.canonical_identifier, 14)
              ELSE rl.canonical_identifier
            END AS skeleton_group_key
          FROM conflicted_mention cm
          JOIN entity_evidence_raw ee
            ON ee.source = cm.source
           AND ee.entity_evidence_id = cm.entity_evidence_id
          JOIN entity_identifier_raw ei
            ON ei.source = ee.source
           AND ei.entity_evidence_id = ee.entity_evidence_id
          JOIN identifier_type_all kit
            ON kit.name = ei.identifier_type
          JOIN needed_resolver_lookup rl
            ON rl.key_identifier_type_id = kit.identifier_type_id
           AND rl.key_value = COALESCE(ei.identifier_normalized, ei.identifier)
           AND rl.evidence_entity_type = ee.entity_type
           AND (
             rl.taxonomy_id = ee.taxonomy_id
             OR rl.taxonomy_id IS NULL
             OR rl.taxonomy_optional_match
           )
        ),
        conflict_group AS (
          SELECT *
          FROM conflict_candidate
          QUALIFY row_number() OVER (
            PARTITION BY source, entity_evidence_id, skeleton_group_key
            ORDER BY
              CASE
                WHEN canonical_identifier_type_id = {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}
                THEN (right(canonical_identifier, 1) <> 'N')
                ELSE false
              END,
              CASE
                WHEN canonical_identifier_type_id = {identifier_type_id(STANDARD_INCHI_KEY_TYPE)}
                THEN (substr(canonical_identifier, 16, 10) <> 'UHFFFAOYSA')
                ELSE false
              END,
              canonical_identifier
          ) = 1
        )
        SELECT
          cg.source,
          cg.entity_evidence_id,
          cg.key_identifier_type_name AS identifier_type,
          cg.key_value AS value_normalized,
          cg.canonical_identifier AS candidate_structure,
          -- the candidate came from an identifier attached to this same
          -- mention, so the mention's own source is who supplied it.
          cg.source AS candidate_source,
          CASE
            WHEN ia.identifier_type_id IS NOT NULL THEN 'authoritative'
            ELSE 'cross_reference'
          END AS candidate_role,
          cg.skeleton_group_key AS skeleton
        FROM conflict_group cg
        LEFT JOIN identifier_authority ia
          ON ia.identifier_type_id = cg.key_identifier_type_id
         AND ia.is_structure_authority
        """
    )
    con.execute(
        """
        CREATE TABLE complex_member_signature_base AS
        WITH complex_member AS (
          SELECT DISTINCT
            parent.source,
            parent.entity_evidence_id,
            child_resolution.entity_type,
            child_resolution.taxonomy_id,
            child_resolution.canonical_identifier_type_id,
            child_resolution.canonical_identifier,
            child_resolution.status
          FROM entity_evidence_raw parent
          JOIN entity_evidence_raw child
            ON child.source = parent.source
           AND child.parent_entity_evidence_id = parent.entity_evidence_id
          JOIN entity_resolution_base child_resolution
            ON child_resolution.source = child.source
           AND child_resolution.entity_evidence_id = child.entity_evidence_id
          WHERE parent.entity_type = ?
        )
        SELECT
          complex_member.source,
          complex_member.entity_evidence_id,
          complex_hash_type.identifier_type_id AS canonical_identifier_type_id,
          sha256(
            to_json(
              list(
                struct_pack(
                  entity_type := complex_member.entity_type,
                  taxonomy_id := complex_member.taxonomy_id,
                  canonical_identifier_type_id := complex_member.canonical_identifier_type_id,
                  canonical_identifier := complex_member.canonical_identifier
                )
                ORDER BY
                  complex_member.entity_type,
                  complex_member.taxonomy_id,
                  complex_member.canonical_identifier_type_id,
                  complex_member.canonical_identifier
              )
            )
          ) AS canonical_identifier,
          CASE
            WHEN bool_and(complex_member.status = 'resolved') THEN 'resolved'
            ELSE 'unresolved'
          END AS status
        FROM complex_member
        CROSS JOIN (
          SELECT identifier_type_id
          FROM identifier_type_all
          WHERE name = ?
        ) complex_hash_type
        GROUP BY
          complex_member.source,
          complex_member.entity_evidence_id,
          complex_hash_type.identifier_type_id
        """,
        [COMPLEX_ENTITY_TYPE, COMPLEX_MEMBER_HASH_ID_TYPE],
    )
    con.execute(
        """
        CREATE TABLE complex_member_signature AS
        WITH complex_member AS (
          SELECT DISTINCT
            parent.source,
            parent.entity_evidence_id,
            child_resolution.entity_type,
            child_resolution.taxonomy_id,
            coalesce(
              child_complex.canonical_identifier_type_id,
              child_resolution.canonical_identifier_type_id
            ) AS canonical_identifier_type_id,
            coalesce(
              child_complex.canonical_identifier,
              child_resolution.canonical_identifier
            ) AS canonical_identifier,
            coalesce(child_complex.status, child_resolution.status) AS status
          FROM entity_evidence_raw parent
          JOIN entity_evidence_raw child
            ON child.source = parent.source
           AND child.parent_entity_evidence_id = parent.entity_evidence_id
          JOIN entity_resolution_base child_resolution
            ON child_resolution.source = child.source
           AND child_resolution.entity_evidence_id = child.entity_evidence_id
          LEFT JOIN complex_member_signature_base child_complex
            ON child_complex.source = child.source
           AND child_complex.entity_evidence_id = child.entity_evidence_id
          WHERE parent.entity_type = ?
        )
        SELECT
          complex_member.source,
          complex_member.entity_evidence_id,
          complex_hash_type.identifier_type_id AS canonical_identifier_type_id,
          sha256(
            to_json(
              list(
                struct_pack(
                  entity_type := complex_member.entity_type,
                  taxonomy_id := complex_member.taxonomy_id,
                  canonical_identifier_type_id := complex_member.canonical_identifier_type_id,
                  canonical_identifier := complex_member.canonical_identifier
                )
                ORDER BY
                  complex_member.entity_type,
                  complex_member.taxonomy_id,
                  complex_member.canonical_identifier_type_id,
                  complex_member.canonical_identifier
              )
            )
          ) AS canonical_identifier,
          CASE
            WHEN bool_and(complex_member.status = 'resolved') THEN 'resolved'
            ELSE 'unresolved'
          END AS status
        FROM complex_member
        CROSS JOIN (
          SELECT identifier_type_id
          FROM identifier_type_all
          WHERE name = ?
        ) complex_hash_type
        GROUP BY
          complex_member.source,
          complex_member.entity_evidence_id,
          complex_hash_type.identifier_type_id
        """,
        [COMPLEX_ENTITY_TYPE, COMPLEX_MEMBER_HASH_ID_TYPE],
    )
    con.execute(
        """
        CREATE TABLE reaction_member_signature AS
        WITH reaction_member AS (
          SELECT DISTINCT
            parent.source,
            parent.entity_evidence_id,
            CASE
              WHEN role_annotation.term IN (?, ?, ?, ?) THEN 'reactant'
              WHEN role_annotation.term IN (?, ?) THEN 'product'
            END AS participant_role,
            child_resolution.entity_type,
            child_resolution.taxonomy_id,
            coalesce(
              child_complex.canonical_identifier_type_id,
              child_resolution.canonical_identifier_type_id
            ) AS canonical_identifier_type_id,
            coalesce(
              child_complex.canonical_identifier,
              child_resolution.canonical_identifier
            ) AS canonical_identifier,
            coalesce(child_complex.status, child_resolution.status) AS status
          FROM entity_evidence_raw parent
          JOIN relation_evidence_raw relation
            ON relation.source = parent.source
           AND relation.subject_entity_evidence_id = parent.entity_evidence_id
           AND relation.predicate = 'has_participant'
          JOIN relation_annotation_raw role_annotation
            ON role_annotation.source = relation.source
           AND role_annotation.evidence_id = relation.relation_evidence_id
          JOIN entity_resolution_base child_resolution
            ON child_resolution.source = relation.source
           AND child_resolution.entity_evidence_id =
               relation.object_entity_evidence_id
          LEFT JOIN complex_member_signature child_complex
            ON child_complex.source = child_resolution.source
           AND child_complex.entity_evidence_id =
               child_resolution.entity_evidence_id
          WHERE parent.entity_type = ?
            AND role_annotation.term IN (?, ?, ?, ?, ?, ?)
        )
        SELECT
          reaction_member.source,
          reaction_member.entity_evidence_id,
          reaction_hash_type.identifier_type_id AS canonical_identifier_type_id,
          sha256(
            to_json(
              list(
                struct_pack(
                  participant_role := reaction_member.participant_role,
                  entity_type := reaction_member.entity_type,
                  taxonomy_id := reaction_member.taxonomy_id,
                  canonical_identifier_type_id :=
                    reaction_member.canonical_identifier_type_id,
                  canonical_identifier := reaction_member.canonical_identifier
                )
                ORDER BY
                  reaction_member.participant_role,
                  reaction_member.entity_type,
                  reaction_member.taxonomy_id,
                  reaction_member.canonical_identifier_type_id,
                  reaction_member.canonical_identifier
              )
            )
          ) AS canonical_identifier,
          CASE
            WHEN bool_and(reaction_member.status = 'resolved') THEN 'resolved'
            ELSE 'unresolved'
          END AS status
        FROM reaction_member
        CROSS JOIN (
          SELECT identifier_type_id
          FROM identifier_type_all
          WHERE name = ?
        ) reaction_hash_type
        GROUP BY
          reaction_member.source,
          reaction_member.entity_evidence_id,
          reaction_hash_type.identifier_type_id
        HAVING bool_or(reaction_member.participant_role = 'reactant')
           AND bool_or(reaction_member.participant_role = 'product')
        """,
        [
            *REACTANT_ROLE_TERMS,
            *PRODUCT_ROLE_TERMS,
            REACTION_ENTITY_TYPE,
            *REACTANT_ROLE_TERMS,
            *PRODUCT_ROLE_TERMS,
            REACTION_MEMBER_HASH_ID_TYPE,
        ],
    )
    con.execute(
        """
        CREATE TABLE entity_resolution AS
        SELECT
          base.source,
          base.dataset,
          base.row_id,
          base.entity_evidence_id,
          base.entity_type,
          base.molecular_entity_type,
          base.taxonomy_id,
          coalesce(
            reaction_member.canonical_identifier_type_id,
            complex_member.canonical_identifier_type_id,
            base.canonical_identifier_type_id
          ) AS canonical_identifier_type_id,
          coalesce(
            reaction_member.canonical_identifier,
            complex_member.canonical_identifier,
            base.canonical_identifier
          ) AS canonical_identifier,
          coalesce(
            reaction_member.status,
            complex_member.status,
            base.status
          ) AS status,
          CASE
            WHEN reaction_member.canonical_identifier IS NOT NULL
              THEN 'reaction'
            WHEN complex_member.canonical_identifier IS NOT NULL
              THEN 'complex'
            ELSE base.resolution_mechanism
          END AS resolution_mechanism
        FROM entity_resolution_base base
        LEFT JOIN complex_member_signature complex_member
          ON complex_member.source = base.source
         AND complex_member.entity_evidence_id = base.entity_evidence_id
        LEFT JOIN reaction_member_signature reaction_member
          ON reaction_member.source = base.source
         AND reaction_member.entity_evidence_id = base.entity_evidence_id
        """
    )
    con.execute(
        """
        CREATE TABLE ontology_term_resolution AS
        WITH term_key AS (
          SELECT DISTINCT
            ot.source,
            ot.ontology_id,
            ot.term_id,
            ot.term_entity_type AS entity_type,
            NULL::VARCHAR AS taxonomy_id,
            kit.identifier_type_id AS term_identifier_type_id,
            ot.term_identifier
          FROM ontology_terms_raw ot
          JOIN identifier_type_all kit
            ON kit.name = ot.term_identifier_type
          WHERE ot.term_identifier IS NOT NULL
            AND ot.term_identifier <> ''
            AND ot.term_entity_type IS NOT NULL
        )
        SELECT
          tk.source,
          tk.ontology_id,
          tk.term_id,
          tk.entity_type,
          coalesce(rl.taxonomy_id, tk.taxonomy_id) AS taxonomy_id,
          tk.term_identifier_type_id,
          tk.term_identifier,
          coalesce(
            rl.canonical_identifier_type_id,
            tk.term_identifier_type_id
          ) AS canonical_identifier_type_id,
          coalesce(
            rl.canonical_identifier,
            tk.term_identifier
          ) AS canonical_identifier
        FROM term_key tk
        LEFT JOIN resolver_entity_type_match etm
          ON etm.evidence_entity_type = tk.entity_type
        LEFT JOIN resolver_lookup rl
          ON rl.entity_type = etm.resolver_entity_type
         AND rl.key_identifier_type_id = tk.term_identifier_type_id
         AND rl.key_value = tk.term_identifier
         AND (rl.taxonomy_id = tk.taxonomy_id OR rl.taxonomy_id IS NULL)
        QUALIFY row_number() OVER (
          PARTITION BY
            tk.source,
            tk.ontology_id,
            tk.entity_type,
            tk.term_identifier_type_id,
            tk.term_identifier
          ORDER BY
            rl.canonical_identifier IS NULL,
            rl.canonical_identifier_type_id,
            rl.canonical_identifier
        ) = 1
        """
    )
    con.execute(
        """
        CREATE TABLE batch_entity_candidate AS
        WITH complex_hash_type AS (
          SELECT identifier_type_id
          FROM identifier_type_all
          WHERE name IN (?, ?)
        ),
        needed_resolved_key AS (
          SELECT DISTINCT
            er.entity_type,
            er.taxonomy_id,
            er.canonical_identifier_type_id,
            er.canonical_identifier
          FROM entity_resolution er
          WHERE er.status = 'resolved'
            AND er.canonical_identifier_type_id NOT IN (
              SELECT identifier_type_id
              FROM complex_hash_type
            )
        ),
        needed_resolved_entity AS (
          SELECT DISTINCT
            entity_type,
            taxonomy_id,
            canonical_identifier_type_id,
            canonical_identifier
          FROM needed_resolved_key
        ),
        complex_member_entity AS (
          SELECT
            er.entity_type,
            er.taxonomy_id,
            er.canonical_identifier_type_id,
            er.canonical_identifier,
            er.status AS resolution_status
          FROM entity_resolution er
          WHERE er.canonical_identifier_type_id IN (
            SELECT identifier_type_id
            FROM complex_hash_type
          )
          GROUP BY
            er.entity_type,
            er.taxonomy_id,
            er.canonical_identifier_type_id,
            er.canonical_identifier,
            er.status
        ),
        cv_term_entity AS (
          SELECT DISTINCT
            ? AS entity_type,
            NULL::VARCHAR AS taxonomy_id,
            cv_type.identifier_type_id AS canonical_identifier_type_id,
            object_id AS canonical_identifier
          FROM annotation_relation_evidence_raw
          CROSS JOIN (
            SELECT identifier_type_id
            FROM identifier_type_all
            WHERE name = ?
          ) cv_type
          WHERE object_entity_type = ?
            AND object_id_type = ?
            AND object_id IS NOT NULL
            AND NOT EXISTS (
              SELECT 1
              FROM needed_resolved_entity existing_term
              WHERE existing_term.entity_type = ?
                AND existing_term.taxonomy_id IS NULL
                AND existing_term.canonical_identifier_type_id =
                    cv_type.identifier_type_id
                AND existing_term.canonical_identifier = object_id
            )
        ),
        annotation_object_entity AS (
          -- Fixed (spec 011, found 2026-09-08 during T047, fixed at the
          -- cycle's final rebuild): promotes through annotation_object_resolution
          -- (built above, in _build_annotation_object_resolution) instead of
          -- self-typing the raw (object_id_type, object_id) verbatim -- so a
          -- chemical object (e.g. a reaction annotated with a bare ChEBI id, not
          -- a full entity_evidence mention) lands on whatever entity the same
          -- real molecule resolves to via the normal entity_resolution_base
          -- mention path, instead of minting a permanent, unpromoted duplicate.
          -- Confirmed reproducibly with CHEBI:107644 (piperonylic acid) across
          -- two full rebuilds before this fix, identical orphaned UUID both
          -- times (entity ids are content-addressed) -- inflated SC-001's
          -- ChEBI-canonical count. Non-chemical object types (no resolver
          -- coverage fetched for them) fall through unchanged to the same raw
          -- self-typing as before -- only CV-term objects stay excluded.
          SELECT DISTINCT
            aor.entity_type,
            NULL::VARCHAR AS taxonomy_id,
            aor.canonical_identifier_type_id,
            aor.canonical_identifier
          FROM annotation_object_resolution aor
          WHERE aor.entity_type <> ?
            AND NOT EXISTS (
              SELECT 1
              FROM needed_resolved_entity existing_entity
              WHERE existing_entity.entity_type = aor.entity_type
                AND existing_entity.taxonomy_id IS NULL
                AND existing_entity.canonical_identifier_type_id =
                    aor.canonical_identifier_type_id
                AND existing_entity.canonical_identifier = aor.canonical_identifier
            )
        ),
        ontology_term_identifier_row AS (
          SELECT DISTINCT
            otr.entity_type,
            otr.taxonomy_id,
            otr.canonical_identifier_type_id,
            otr.canonical_identifier,
            otr.term_identifier_type_id AS identifier_type_id,
            otr.term_identifier AS identifier
          FROM ontology_term_resolution otr
          WHERE otr.canonical_identifier IS NOT NULL
            AND otr.canonical_identifier <> ''
          UNION
          SELECT DISTINCT
            otr.entity_type,
            otr.taxonomy_id,
            otr.canonical_identifier_type_id,
            otr.canonical_identifier,
            otr.canonical_identifier_type_id AS identifier_type_id,
            otr.canonical_identifier AS identifier
          FROM ontology_term_resolution otr
          WHERE otr.canonical_identifier IS NOT NULL
            AND otr.canonical_identifier <> ''
        ),
        ontology_term_entity AS (
          SELECT DISTINCT
            otir.entity_type,
            otir.taxonomy_id,
            otir.canonical_identifier_type_id,
            otir.canonical_identifier
          FROM ontology_term_identifier_row otir
          WHERE NOT EXISTS (
            SELECT 1
            FROM needed_resolved_entity existing_entity
            WHERE existing_entity.entity_type = otir.entity_type
              AND existing_entity.taxonomy_id IS NOT DISTINCT FROM otir.taxonomy_id
              AND existing_entity.canonical_identifier_type_id =
                  otir.canonical_identifier_type_id
              AND existing_entity.canonical_identifier =
                  otir.canonical_identifier
          )
        ),
        ontology_relation_endpoint_key AS (
          SELECT DISTINCT
            endpoint.source,
            endpoint.entity_type,
            NULL::VARCHAR AS taxonomy_id,
            kit.identifier_type_id AS endpoint_identifier_type_id,
            endpoint.identifier AS endpoint_identifier,
            coalesce(rl.taxonomy_id, NULL::VARCHAR) AS resolved_taxonomy_id,
            coalesce(
              rl.canonical_identifier_type_id,
              kit.identifier_type_id
            ) AS canonical_identifier_type_id,
            coalesce(
              rl.canonical_identifier,
              endpoint.identifier
            ) AS canonical_identifier
          FROM (
            SELECT
              source,
              subject_entity_type AS entity_type,
              subject_identifier_type AS identifier_type,
              subject_identifier AS identifier
            FROM ontology_relation_raw
            WHERE subject_entity_evidence_id IS NULL
              AND subject_identifier IS NOT NULL
              AND subject_identifier <> ''
            UNION
            SELECT
              source,
              object_entity_type AS entity_type,
              object_identifier_type AS identifier_type,
              object_identifier AS identifier
            FROM ontology_relation_raw
            WHERE object_identifier IS NOT NULL
              AND object_identifier <> ''
          ) endpoint
          JOIN identifier_type_all kit
            ON kit.name = endpoint.identifier_type
          LEFT JOIN needed_ontology_endpoint_resolver_lookup rl
            ON rl.evidence_entity_type = endpoint.entity_type
           AND rl.key_identifier_type_id = kit.identifier_type_id
           AND rl.key_value = endpoint.identifier
           AND rl.taxonomy_id IS NULL
          WHERE endpoint.entity_type IS NOT NULL
          QUALIFY row_number() OVER (
            PARTITION BY
              endpoint.entity_type,
              kit.identifier_type_id,
              endpoint.identifier
            ORDER BY
              rl.canonical_identifier IS NULL,
              rl.canonical_identifier_type_id,
              rl.canonical_identifier
          ) = 1
        ),
        ontology_relation_identifier_row AS (
          SELECT DISTINCT
            entity_type,
            resolved_taxonomy_id AS taxonomy_id,
            canonical_identifier_type_id,
            canonical_identifier,
            endpoint_identifier_type_id AS identifier_type_id,
            endpoint_identifier AS identifier
          FROM ontology_relation_endpoint_key
          WHERE canonical_identifier IS NOT NULL
            AND canonical_identifier <> ''
          UNION
          SELECT DISTINCT
            entity_type,
            resolved_taxonomy_id AS taxonomy_id,
            canonical_identifier_type_id,
            canonical_identifier,
            canonical_identifier_type_id AS identifier_type_id,
            canonical_identifier AS identifier
          FROM ontology_relation_endpoint_key
          WHERE canonical_identifier IS NOT NULL
            AND canonical_identifier <> ''
        ),
        ontology_relation_endpoint_entity AS (
          SELECT DISTINCT
            orir.entity_type,
            orir.taxonomy_id,
            orir.canonical_identifier_type_id,
            orir.canonical_identifier
          FROM ontology_relation_identifier_row orir
          WHERE NOT EXISTS (
            SELECT 1
            FROM needed_resolved_entity existing_entity
            WHERE existing_entity.entity_type = orir.entity_type
              AND existing_entity.taxonomy_id IS NOT DISTINCT FROM orir.taxonomy_id
              AND existing_entity.canonical_identifier_type_id =
                  orir.canonical_identifier_type_id
              AND existing_entity.canonical_identifier =
                  orir.canonical_identifier
          )
        ),
        all_entity AS (
          SELECT
            entity_type,
            taxonomy_id,
            canonical_identifier_type_id,
            canonical_identifier,
            'resolved' AS resolution_status
          FROM needed_resolved_entity
          UNION ALL
          SELECT
            entity_type,
            taxonomy_id,
            canonical_identifier_type_id,
            canonical_identifier,
            'resolved' AS resolution_status
          FROM cv_term_entity
          UNION ALL
          SELECT
            entity_type,
            taxonomy_id,
            canonical_identifier_type_id,
            canonical_identifier,
            'resolved' AS resolution_status
          FROM annotation_object_entity
          UNION ALL
          SELECT
            entity_type,
            taxonomy_id,
            canonical_identifier_type_id,
            canonical_identifier,
            'resolved' AS resolution_status
          FROM ontology_term_entity
          UNION ALL
          SELECT
            entity_type,
            taxonomy_id,
            canonical_identifier_type_id,
            canonical_identifier,
            'resolved' AS resolution_status
          FROM ontology_relation_endpoint_entity
          UNION ALL
          SELECT * FROM complex_member_entity
          UNION ALL
          SELECT
            er.entity_type,
            er.taxonomy_id,
            er.canonical_identifier_type_id,
            er.canonical_identifier,
            er.status AS resolution_status
          FROM entity_resolution er
          WHERE er.status = 'unresolved'
            AND er.canonical_identifier_type_id NOT IN (
              SELECT identifier_type_id
              FROM complex_hash_type
            )
          GROUP BY
            er.entity_type,
            er.taxonomy_id,
            er.canonical_identifier_type_id,
            er.canonical_identifier,
            er.status
        )
        SELECT
          canonical_entity_uuid(
            all_entity.entity_type,
            all_entity.taxonomy_id,
            it.name,
            all_entity.canonical_identifier
          ) AS entity_id,
          canonical_entity_key(
            all_entity.entity_type,
            all_entity.taxonomy_id,
            it.name,
            all_entity.canonical_identifier
                  ) AS entity_key,
                  all_entity.entity_type,
                  all_entity.taxonomy_id,
                  all_entity.canonical_identifier_type_id,
                  it.name AS canonical_identifier_type,
                  all_entity.canonical_identifier,
                  all_entity.resolution_status,
          erm.resolution_mechanism,
          string_agg(DISTINCT source_rows.source, ',' ORDER BY source_rows.source)
            AS sources
        FROM all_entity
        JOIN identifier_type_all it
          ON it.identifier_type_id = all_entity.canonical_identifier_type_id
        LEFT JOIN (
          SELECT DISTINCT
            source,
            entity_type,
            taxonomy_id,
            canonical_identifier_type_id,
            canonical_identifier
          FROM entity_resolution
          UNION
          -- Matches annotation_object_entity's own resolution above (the
          -- promoted identifier when the resolver has one, the raw one when it
          -- doesn't) -- so this source-attribution join actually finds the
          -- all_entity row a promoted annotation object landed on, not the
          -- (now nonexistent) raw-keyed one.
          SELECT DISTINCT
            ar.source,
            aor.entity_type,
            NULL::VARCHAR AS taxonomy_id,
            aor.canonical_identifier_type_id,
            aor.canonical_identifier
          FROM annotation_relation_evidence_raw ar
          JOIN identifier_type_all object_type
            ON object_type.name = ar.object_id_type
          JOIN annotation_object_resolution aor
            ON aor.entity_type = ar.object_entity_type
           AND aor.raw_identifier_type_id = object_type.identifier_type_id
           AND aor.raw_identifier = ar.object_id
          WHERE ar.object_id IS NOT NULL
          UNION
          SELECT DISTINCT
            source,
            entity_type,
            taxonomy_id,
            canonical_identifier_type_id,
            canonical_identifier
          FROM ontology_term_resolution
          UNION
          SELECT DISTINCT
            source,
            entity_type,
            resolved_taxonomy_id AS taxonomy_id,
            canonical_identifier_type_id,
            canonical_identifier
          FROM ontology_relation_endpoint_key
        ) source_rows
          ON source_rows.entity_type = all_entity.entity_type
         AND source_rows.taxonomy_id IS NOT DISTINCT FROM all_entity.taxonomy_id
         AND source_rows.canonical_identifier_type_id =
             all_entity.canonical_identifier_type_id
         AND source_rows.canonical_identifier = all_entity.canonical_identifier
        LEFT JOIN (
          -- representative resolution_mechanism per canonical entity: when
          -- several mentions resolve to the same entity by different
          -- mechanisms, keep the most authoritative one (lower priority wins).
          SELECT
            entity_type,
            taxonomy_id,
            canonical_identifier_type_id,
            canonical_identifier,
            min_by(
              resolution_mechanism,
              CASE resolution_mechanism
                WHEN 'inchikey' THEN 1
                WHEN 'anchored_structure' THEN 2
                WHEN 'smiles' THEN 3
                WHEN 'chebi' THEN 4
                WHEN 'anchored_chebi' THEN 5
                WHEN 'resolver' THEN 6
                WHEN 'cv_term' THEN 7
                WHEN 'pathway' THEN 8
                WHEN 'gene_anchor' THEN 9
                WHEN 'complex' THEN 10
                WHEN 'reaction' THEN 11
                WHEN 'chembl' THEN 12
                WHEN 'pubchem' THEN 13
                WHEN 'swisslipids' THEN 14
                WHEN 'hmdb' THEN 15
                WHEN 'lipidmaps' THEN 16
                WHEN 'original_id' THEN 17
                WHEN 'name' THEN 18
                WHEN 'unresolved' THEN 19
                ELSE 99
              END
            ) AS resolution_mechanism
          FROM entity_resolution
          WHERE resolution_mechanism IS NOT NULL
          GROUP BY
            entity_type,
            taxonomy_id,
            canonical_identifier_type_id,
            canonical_identifier
        ) erm
          ON erm.entity_type = all_entity.entity_type
         AND erm.taxonomy_id IS NOT DISTINCT FROM all_entity.taxonomy_id
         AND erm.canonical_identifier_type_id =
             all_entity.canonical_identifier_type_id
         AND erm.canonical_identifier = all_entity.canonical_identifier
        GROUP BY
          entity_id,
                  entity_key,
                  all_entity.entity_type,
                  all_entity.taxonomy_id,
                  all_entity.canonical_identifier_type_id,
                  it.name,
                  all_entity.canonical_identifier,
                  all_entity.resolution_status,
                  erm.resolution_mechanism
        """,
        [
            COMPLEX_MEMBER_HASH_ID_TYPE,
            REACTION_MEMBER_HASH_ID_TYPE,
            CV_TERM_ENTITY_TYPE,
            CV_TERM_ID_TYPE,
            CV_TERM_ENTITY_TYPE,
            CV_TERM_ID_TYPE,
            CV_TERM_ENTITY_TYPE,
            CV_TERM_ENTITY_TYPE,
        ],
    )
    con.execute(
        """
        CREATE TABLE new_entity AS
        SELECT c.*
        FROM batch_entity_candidate c
        LEFT JOIN entity_key_cache cache
          ON cache.entity_key = c.entity_key
        WHERE cache.entity_key IS NULL
        """
    )
    con.execute(
        """
        CREATE TABLE canonical_entity AS
        SELECT * FROM batch_entity_candidate
        """
    )
    # gene_protein_representative: 1:1 gene entity -> chosen representative
    # UniProt, joined by Entrez anchor. Guarded on the resolver source view
    # (older resolver snapshots omit it) — empty table otherwise so the
    # bulk-copy path stays uniform.
    gene_protein_representative_src_exists = bool(
        con.execute(
            """
            SELECT COUNT(*)
            FROM information_schema.tables
            WHERE table_name = 'gene_protein_representative_src'
            """
        ).fetchone()[0]
    )
    has_gene_entities = bool(
        con.execute(
            f"""
            SELECT count(*) > 0
            FROM canonical_entity
            WHERE entity_type = {_sql_literal(GENE_ENTITY_TYPE)}
              AND canonical_identifier_type_id = {identifier_type_id(ENTREZ_TYPE)}
            """
        ).fetchone()[0]
    )
    # A gene product known only by a UniProt accession is typed as a gene keyed by
    # that accession (an "unknown gene"). Such a gene is its own representative
    # UniProt. This arm is disjoint from the Entrez-anchored arm below — a
    # different canonical identifier type gives a different entity_id — so the 1:1
    # entity_id primary key holds under UNION ALL. It exists regardless of whether
    # the resolver snapshot ships the gene_protein_representative source view.
    uniprot_gene_arm = f"""
        SELECT
          ce.entity_id,
          ce.canonical_identifier AS representative_uniprot,
          NULL::BOOLEAN AS is_reviewed,
          [ce.canonical_identifier] AS uniprot_all
        FROM canonical_entity ce
        WHERE ce.entity_type = {_sql_literal(GENE_ENTITY_TYPE)}
          AND ce.canonical_identifier_type_id = {identifier_type_id(UNIPROT_TYPE)}
    """
    if gene_protein_representative_src_exists and has_gene_entities:
        con.execute(
            f"""
            CREATE TABLE gene_protein_representative AS
            SELECT
              ce.entity_id,
              src.representative_uniprot,
              src.is_reviewed,
              src.uniprot_all
            FROM gene_protein_representative_src src
            JOIN canonical_entity ce
              ON ce.entity_type = {_sql_literal(GENE_ENTITY_TYPE)}
             AND ce.canonical_identifier = src.canonical_identifier
             AND ce.taxonomy_id IS NOT DISTINCT FROM src.taxonomy_id
            UNION ALL
{uniprot_gene_arm}
            """
        )
    else:
        con.execute(
            f"""
            CREATE TABLE gene_protein_representative AS
{uniprot_gene_arm}
            """
        )
    has_resolver_alias_entities = bool(
        con.execute(
            """
            SELECT count(*) > 0
            FROM canonical_entity ce
            JOIN resolver_entity_type_match etm
              ON etm.evidence_entity_type = ce.entity_type
            """
        ).fetchone()[0]
    )
    resolver_alias_identifier_sql = ''
    if has_resolver_alias_entities:
        resolver_alias_identifier_sql = f"""
          UNION
          SELECT
            se.source,
            se.entity_id,
            rl.key_identifier_type_id AS identifier_type_id,
            rl.key_value AS identifier
          FROM source_entity se
          JOIN resolver_entity_type_match etm
            ON etm.evidence_entity_type = se.entity_type
          JOIN resolver_lookup rl
            ON rl.entity_type = etm.resolver_entity_type
           AND rl.canonical_identifier_type_id =
               se.canonical_identifier_type_id
           AND rl.canonical_identifier = se.canonical_identifier
           AND (
             rl.taxonomy_id = se.taxonomy_id
             OR rl.taxonomy_id IS NULL
           )
          WHERE rl.key_value IS NOT NULL
            AND rl.key_value <> ''
            AND rl.key_identifier_type_id NOT IN (
              {resolver_alias_expansion_excluded_type_ids}
            )
        """
    con.execute(
        f"""
        CREATE TABLE canonical_entity_identifier AS
        WITH source_entity AS (
          SELECT DISTINCT
            source_name AS source,
            ce.entity_id,
            ce.entity_type,
            ce.taxonomy_id,
            ce.canonical_identifier_type_id,
            ce.canonical_identifier
          FROM canonical_entity ce,
            unnest(string_split(ce.sources, ',')) AS source_names(source_name)
          WHERE ce.sources IS NOT NULL
            AND ce.sources <> ''
            AND source_name IS NOT NULL
            AND source_name <> ''
        ),
        identifier_rows AS (
          SELECT
            se.source,
            se.entity_id,
            se.canonical_identifier_type_id AS identifier_type_id,
            se.canonical_identifier AS identifier
          FROM source_entity se
          WHERE se.canonical_identifier IS NOT NULL
            AND se.canonical_identifier <> ''
          {resolver_alias_identifier_sql}
          UNION
          SELECT
            er.source,
            ce.entity_id,
            kit.identifier_type_id,
            ei.identifier
          FROM entity_resolution er
          JOIN canonical_entity ce
            ON ce.entity_type = er.entity_type
           AND ce.taxonomy_id IS NOT DISTINCT FROM er.taxonomy_id
           AND ce.canonical_identifier_type_id =
               er.canonical_identifier_type_id
           AND ce.canonical_identifier = er.canonical_identifier
          JOIN entity_identifier_raw ei
            ON ei.source = er.source
           AND ei.entity_evidence_id = er.entity_evidence_id
          JOIN identifier_type_all kit
            ON kit.name = ei.identifier_type
          WHERE ei.identifier IS NOT NULL
            AND ei.identifier <> ''
          UNION
          SELECT
            otr.source,
            ce.entity_id,
            otr.term_identifier_type_id AS identifier_type_id,
            otr.term_identifier AS identifier
          FROM ontology_term_resolution otr
          JOIN canonical_entity ce
            ON ce.entity_type = otr.entity_type
           AND ce.taxonomy_id IS NOT DISTINCT FROM otr.taxonomy_id
           AND ce.canonical_identifier_type_id =
               otr.canonical_identifier_type_id
           AND ce.canonical_identifier = otr.canonical_identifier
          WHERE otr.term_identifier IS NOT NULL
            AND otr.term_identifier <> ''
        )
        SELECT DISTINCT
          ir.source,
          ir.entity_id,
          it.name AS identifier_type,
          ir.identifier_type_id,
          ir.identifier
        FROM identifier_rows ir
        JOIN identifier_type_all it
          ON it.identifier_type_id = ir.identifier_type_id
        WHERE ir.source IS NOT NULL
          AND ir.identifier IS NOT NULL
          AND ir.identifier <> ''
        """
    )
    con.execute(
        """
        CREATE TABLE entity_ontology_relation AS
        WITH raw_edge AS (
          SELECT
            *,
            source || chr(31) || coalesce(ontology_id, '') || chr(31) ||
              coalesce(
                subject_entity_evidence_id,
                subject_entity_type || chr(30) ||
                  subject_identifier_type || chr(30) ||
                  subject_identifier,
                ''
              ) || chr(31) || predicate || chr(31) ||
              object_entity_type || chr(30) ||
                object_identifier_type || chr(30) ||
                object_identifier AS edge_key
          FROM ontology_relation_raw
          WHERE predicate IS NOT NULL
            AND predicate <> ''
            AND object_identifier IS NOT NULL
            AND object_identifier <> ''
        ),
        endpoint_resolution AS (
          SELECT
            raw.source,
            raw.dataset,
            raw.ontology_id,
            raw.edge_key,
            'subject' AS endpoint_role,
            er.entity_type,
            er.taxonomy_id,
            er.canonical_identifier_type_id,
            er.canonical_identifier
          FROM raw_edge raw
          JOIN entity_resolution er
            ON er.source = raw.source
           AND er.entity_evidence_id = raw.subject_entity_evidence_id
          WHERE raw.subject_entity_evidence_id IS NOT NULL
          UNION ALL
          SELECT
            endpoint.source,
            endpoint.dataset,
            endpoint.ontology_id,
            endpoint.edge_key,
            endpoint.endpoint_role,
            endpoint.entity_type,
            coalesce(rl.taxonomy_id, NULL::VARCHAR) AS taxonomy_id,
            coalesce(
              rl.canonical_identifier_type_id,
              kit.identifier_type_id
            ) AS canonical_identifier_type_id,
            coalesce(
              rl.canonical_identifier,
              endpoint.identifier
            ) AS canonical_identifier
          FROM (
            SELECT
              source,
              dataset,
              ontology_id,
              edge_key,
              subject_entity_type AS entity_type,
              subject_identifier_type AS identifier_type,
              subject_identifier AS identifier,
              'subject' AS endpoint_role
            FROM raw_edge
            WHERE subject_entity_evidence_id IS NULL
              AND subject_identifier IS NOT NULL
              AND subject_identifier <> ''
            UNION ALL
            SELECT
              source,
              dataset,
              ontology_id,
              edge_key,
              object_entity_type AS entity_type,
              object_identifier_type AS identifier_type,
              object_identifier AS identifier,
              'object' AS endpoint_role
            FROM raw_edge
          ) endpoint
          JOIN identifier_type_all kit
            ON kit.name = endpoint.identifier_type
          LEFT JOIN needed_ontology_endpoint_resolver_lookup rl
            ON rl.evidence_entity_type = endpoint.entity_type
           AND rl.key_identifier_type_id = kit.identifier_type_id
           AND rl.key_value = endpoint.identifier
           AND rl.taxonomy_id IS NULL
          QUALIFY row_number() OVER (
            PARTITION BY
              endpoint.edge_key,
              endpoint.endpoint_role
            ORDER BY
              rl.canonical_identifier IS NULL,
              rl.canonical_identifier_type_id,
              rl.canonical_identifier
          ) = 1
        ),
        edge_endpoint AS (
          SELECT
            er.source,
            er.dataset,
            er.ontology_id,
            er.edge_key,
            er.endpoint_role,
            ce.entity_id
          FROM endpoint_resolution er
          JOIN canonical_entity ce
            ON ce.entity_type = er.entity_type
           AND ce.taxonomy_id IS NOT DISTINCT FROM er.taxonomy_id
           AND ce.canonical_identifier_type_id =
               er.canonical_identifier_type_id
           AND ce.canonical_identifier = er.canonical_identifier
        )
        SELECT DISTINCT
          raw.source,
          raw.dataset,
          raw.ontology_id,
          subject.entity_id AS subject_entity_id,
          raw.predicate,
          object.entity_id AS object_entity_id
        FROM raw_edge raw
        JOIN edge_endpoint subject
          ON subject.source = raw.source
         AND subject.ontology_id IS NOT DISTINCT FROM raw.ontology_id
         AND subject.endpoint_role = 'subject'
         AND subject.edge_key = raw.edge_key
        JOIN edge_endpoint object
          ON object.source = raw.source
         AND object.ontology_id IS NOT DISTINCT FROM raw.ontology_id
         AND object.endpoint_role = 'object'
         AND object.edge_key = subject.edge_key
        """
    )
    con.execute(
        f"""
        CREATE TABLE entity_evidence_resolution AS
        SELECT
          er.source,
          er.entity_evidence_id,
          er.status,
          ce.entity_id,
          CASE er.molecular_entity_type
            WHEN {_sql_literal(GENE_ENTITY_TYPE)} THEN 1
            WHEN {_sql_literal(PROTEIN_ENTITY_TYPE)} THEN 2
            WHEN {_sql_literal(MIRNA_ENTITY_TYPE)} THEN 4
            ELSE NULL
          END AS molecular_type_id,
          -- Why an unresolved record stayed unresolved: 'missing_translation_table'
          -- (6) when no resolver relation covering its entity type was available
          -- (pre-flight found them absent/empty), otherwise
          -- 'no_accepted_resolver_candidate' (3) — the resolver had data but no
          -- match. Resolved and ambiguous records carry no reason.
          CASE
            WHEN er.status = 'unresolved' AND mtt.entity_type IS NOT NULL THEN 6
            WHEN er.status = 'unresolved' THEN 3
            ELSE NULL
          END::SMALLINT AS reason_id
        FROM entity_resolution er
        LEFT JOIN canonical_entity ce
          ON ce.entity_type = er.entity_type
         AND ce.taxonomy_id IS NOT DISTINCT FROM er.taxonomy_id
         AND ce.canonical_identifier_type_id = er.canonical_identifier_type_id
         AND ce.canonical_identifier = er.canonical_identifier
        LEFT JOIN missing_translation_entity_type mtt
          ON mtt.entity_type = er.entity_type
        """
    )
    # --- molecular state (the heavy opt-in tier) -----------------------------
    # The per-record asserted UniProt AC/isoform a source gave for a
    # gene-resolved protein mention, captured as a `state` (a bag of
    # components) linked to the evidence via `evidence_state` (one-to-many).
    # Only protein records that carry a specific UniProt AC get a state —
    # bare-symbol mentions already have molecular_type=protein on
    # entity_evidence_resolution (the cheap tier) and need no state row.
    # Isoforms (``P12345-2``) split into a ``uniprot`` base-AC component plus
    # an ``isoform`` component.
    con.execute(
        f"""
        CREATE TABLE evidence_state_link AS
        WITH asserted AS (
          SELECT DISTINCT
            eer.source,
            eer.entity_evidence_id,
            eer.entity_id AS gene_entity_id,
            ei.identifier AS uniprot_value
          FROM entity_evidence_resolution eer
          JOIN entity_identifier_raw ei
            ON ei.source = eer.source
           AND ei.entity_evidence_id = eer.entity_evidence_id
          WHERE eer.entity_id IS NOT NULL
            AND eer.molecular_type_id = {PROTEIN_MOLECULAR_TYPE_ID}
            AND ei.identifier_type = {_sql_literal(UNIPROT_TYPE)}
            AND ei.identifier IS NOT NULL
            AND ei.identifier <> ''
        ),
        components AS (
          SELECT
            source,
            entity_evidence_id,
            gene_entity_id,
            regexp_replace(uniprot_value, '-[0-9]+$', '') AS uniprot_ac,
            CASE
              WHEN regexp_matches(uniprot_value, '-[0-9]+$')
              THEN uniprot_value
              ELSE NULL
            END AS isoform
          FROM asserted
        )
        SELECT
          source,
          entity_evidence_id,
          gene_entity_id,
          uniprot_ac,
          isoform,
          content_uuid(
            'state' || chr(31) ||
            gene_entity_id::VARCHAR || chr(31) ||
            '{PROTEIN_MOLECULAR_TYPE_ID}' || chr(31) ||
            uniprot_ac || chr(31) ||
            coalesce(isoform, '')
          ) AS state_id
        FROM components
        """
    )
    con.execute(
        f"""
        CREATE TABLE state AS
        SELECT DISTINCT
          state_id,
          gene_entity_id,
          {PROTEIN_MOLECULAR_TYPE_ID}::SMALLINT AS molecular_type_id
        FROM evidence_state_link
        """
    )
    con.execute(
        """
        CREATE TABLE state_component AS
        SELECT DISTINCT state_id, 'uniprot' AS component_type, uniprot_ac AS value
        FROM evidence_state_link
        UNION
        SELECT DISTINCT state_id, 'isoform' AS component_type, isoform AS value
        FROM evidence_state_link
        WHERE isoform IS NOT NULL
        """
    )
    con.execute(
        """
        CREATE TABLE evidence_state AS
        SELECT DISTINCT source, entity_evidence_id, state_id
        FROM evidence_state_link
        """
    )
    con.execute(
        """
        CREATE TABLE relation_candidate_evidence AS
        WITH member_projected AS (
          SELECT
            rr.source,
            rr.relation_evidence_id,
            subject.entity_id AS subject_entity_id,
            rr.predicate,
            object.entity_id AS object_entity_id,
            rr.relation_category
          FROM relation_evidence_raw rr
          JOIN entity_evidence_resolution subject
            ON subject.source = rr.source
           AND subject.entity_evidence_id = rr.subject_entity_evidence_id
          JOIN entity_evidence_resolution object
            ON object.source = rr.source
           AND object.entity_evidence_id = rr.object_entity_evidence_id
          WHERE subject.entity_id IS NOT NULL
            AND object.entity_id IS NOT NULL
        ),
        annotation_projected AS (
          -- object repoints through annotation_object_resolution (the same
          -- promotion annotation_object_entity applied when minting/joining
          -- the object's own entity row) -- previously joined canonical_entity
          -- on the raw (object_id_type, object_id), a dead end for any object
          -- that promoted to a different, properly-resolved entity elsewhere.
          SELECT
            ar.source,
            ar.relation_evidence_id,
            object.entity_id AS subject_entity_id,
            ar.predicate,
            subject.entity_id AS object_entity_id,
            ar.relation_category
          FROM annotation_relation_evidence_raw ar
          JOIN entity_evidence_resolution subject
            ON subject.source = ar.source
           AND subject.entity_evidence_id = ar.subject_entity_evidence_id
          JOIN identifier_type_all raw_object_type
            ON raw_object_type.name = ar.object_id_type
          JOIN annotation_object_resolution aor
            ON aor.entity_type = ar.object_entity_type
           AND aor.raw_identifier_type_id = raw_object_type.identifier_type_id
           AND aor.raw_identifier = ar.object_id
          JOIN identifier_type_all resolved_object_type
            ON resolved_object_type.identifier_type_id =
               aor.canonical_identifier_type_id
          JOIN canonical_entity object
            ON object.entity_type = aor.entity_type
           AND object.canonical_identifier_type = resolved_object_type.name
           AND object.canonical_identifier = aor.canonical_identifier
          WHERE subject.entity_id IS NOT NULL
        )
        SELECT * FROM member_projected
        UNION ALL
        SELECT * FROM annotation_projected
        """
    )
    con.execute(
        """
        CREATE TABLE batch_relation_candidate AS
        SELECT
          content_uuid(
            subject_entity_id::VARCHAR || '|' ||
            predicate || '|' ||
            object_entity_id::VARCHAR
          ) AS relation_id,
          subject_entity_id::VARCHAR || '|' ||
            predicate || '|' ||
            object_entity_id::VARCHAR AS relation_key,
          subject_entity_id,
          predicate,
          object_entity_id,
          min(relation_category) AS relation_category,
          string_agg(DISTINCT source, ',' ORDER BY source) AS sources
        FROM relation_candidate_evidence
        GROUP BY subject_entity_id, predicate, object_entity_id
        """
    )
    con.execute(
        """
        CREATE TABLE new_relation AS
        SELECT c.*
        FROM batch_relation_candidate c
        LEFT JOIN relation_key_cache cache
          ON cache.relation_key = c.relation_key
        WHERE cache.relation_key IS NULL
        """
    )
    con.execute(
        """
        CREATE TABLE relation AS
        SELECT * FROM batch_relation_candidate
        """
    )
    con.execute(
        """
        CREATE TABLE relation_evidence_relation AS
        SELECT
          evidence.source,
          evidence.relation_evidence_id,
          r.relation_id
        FROM relation_candidate_evidence evidence
        JOIN relation r
          ON r.subject_entity_id = evidence.subject_entity_id
         AND r.predicate = evidence.predicate
         AND r.object_entity_id = evidence.object_entity_id
        """
    )
    _refresh_duckdb_canonical_caches(con)

    entities = int(
        con.sql('SELECT COUNT(*) FROM canonical_entity').fetchone()[0]
    )
    relations = int(con.sql('SELECT COUNT(*) FROM relation').fetchone()[0])
    links = int(
        con.sql('SELECT COUNT(*) FROM relation_evidence_relation').fetchone()[0]
    )
    return entities, relations, links


def _refresh_duckdb_canonical_caches(con: duckdb.DuckDBPyConnection) -> None:
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE entity_cache_source AS
        SELECT
          entity_key,
          unnest(string_split(sources, ',')) AS source
        FROM entity_key_cache
        WHERE sources IS NOT NULL
          AND sources <> ''
        UNION
        SELECT
          entity_key,
          unnest(string_split(sources, ',')) AS source
        FROM batch_entity_candidate
        WHERE sources IS NOT NULL
          AND sources <> ''
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TABLE entity_key_cache AS
        SELECT
          any_value(entity_id) AS entity_id,
          entity_key,
          any_value(entity_type) AS entity_type,
          any_value(taxonomy_id) AS taxonomy_id,
              any_value(canonical_identifier_type) AS canonical_identifier_type,
              any_value(canonical_identifier_type_id) AS canonical_identifier_type_id,
              any_value(canonical_identifier) AS canonical_identifier,
              coalesce(
            (
              SELECT string_agg(DISTINCT source, ',' ORDER BY source)
              FROM entity_cache_source src
              WHERE src.entity_key = merged.entity_key
            ),
            ''
          ) AS sources,
          min(first_seen_at) AS first_seen_at,
          now() AS last_seen_at
        FROM (
          SELECT
            entity_id,
            entity_key,
            entity_type,
            taxonomy_id,
                canonical_identifier_type,
                canonical_identifier_type_id,
                canonical_identifier,
                sources,
            first_seen_at,
            last_seen_at
          FROM entity_key_cache
          UNION ALL
          SELECT
            entity_id,
            entity_key,
            entity_type,
            taxonomy_id,
                canonical_identifier_type,
                canonical_identifier_type_id,
                canonical_identifier,
                sources,
            now(),
            now()
          FROM batch_entity_candidate
        ) merged
        GROUP BY entity_key
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TEMP TABLE relation_cache_source AS
        SELECT
          relation_key,
          unnest(string_split(sources, ',')) AS source
        FROM relation_key_cache
        WHERE sources IS NOT NULL
          AND sources <> ''
        UNION
        SELECT
          relation_key,
          unnest(string_split(sources, ',')) AS source
        FROM batch_relation_candidate
        WHERE sources IS NOT NULL
          AND sources <> ''
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TABLE relation_key_cache AS
        SELECT
          any_value(relation_id) AS relation_id,
          relation_key,
          any_value(subject_entity_id) AS subject_entity_id,
          any_value(predicate) AS predicate,
          any_value(object_entity_id) AS object_entity_id,
          coalesce(
            (
              SELECT string_agg(DISTINCT source, ',' ORDER BY source)
              FROM relation_cache_source src
              WHERE src.relation_key = merged.relation_key
            ),
            ''
          ) AS sources,
          min(first_seen_at) AS first_seen_at,
          now() AS last_seen_at
        FROM (
          SELECT
            relation_id,
            relation_key,
            subject_entity_id,
            predicate,
            object_entity_id,
            sources,
            first_seen_at,
            last_seen_at
          FROM relation_key_cache
          UNION ALL
          SELECT
            relation_id,
            relation_key,
            subject_entity_id,
            predicate,
            object_entity_id,
            sources,
            now(),
            now()
          FROM batch_relation_candidate
        ) merged
        GROUP BY relation_key
        """
    )


_BULK_COPY_CONTENT_TABLES = (
    'identifier_evidence',
    'annotation',
    'entity_evidence',
    'entity_evidence_identifier',
    'entity_evidence_annotation',
    'relation_evidence',
    'relation_evidence_annotation',
    'entity',
    'entity_identifier',
    'ontology_terms',
    'entity_ontology_relation',
    'entity_evidence_resolution',
    'resolution_conflict',
    'relation',
    'relation_evidence_relation',
    'gene_protein_representative',
    'state',
    'state_component',
    'evidence_state',
)


def _bulk_load_create_views_from_loaded_tables(
    con: duckdb.DuckDBPyConnection,
) -> None:
    views = {
        'pq_entity_evidence': 'entity_evidence_raw',
        'pq_entity_identifier': 'entity_identifier_raw',
        'pq_entity_annotation': 'entity_annotation_raw',
        'pq_relation_evidence': 'relation_evidence_raw',
        'pq_annotation_relation_evidence': 'annotation_relation_evidence_raw',
        'pq_relation_annotation': 'relation_annotation_raw',
        'pq_annotation': 'annotation_value',
        'pq_entity': 'canonical_entity',
        'pq_entity_identifier_resolved': 'canonical_entity_identifier',
        'pq_entity_evidence_resolution': 'entity_evidence_resolution',
        'pq_resolution_conflict': 'resolution_conflict',
        'pq_relation': 'relation',
        'pq_relation_evidence_relation': 'relation_evidence_relation',
        'pq_ontology_terms': 'ontology_terms_raw',
        'pq_entity_ontology_relation': 'entity_ontology_relation',
        'pq_gene_protein_representative': 'gene_protein_representative',
        'pq_state': 'state',
        'pq_state_component': 'state_component',
        'pq_evidence_state': 'evidence_state',
    }
    for view, table in views.items():
        con.execute(f'CREATE VIEW {view} AS SELECT * FROM {table}')
    con.execute(
        """
        CREATE VIEW pq_annotation_relation_evidence_resolved AS
        -- object repoints through annotation_object_resolution, same as
        -- annotation_projected above -- the raw (object_id_type, object_id)
        -- join is a dead end for any object promoted to a different,
        -- properly-resolved entity.
        SELECT
          ar.*,
          object.entity_id AS object_entity_id
        FROM annotation_relation_evidence_raw ar
        JOIN identifier_type_all raw_object_type
          ON raw_object_type.name = ar.object_id_type
        JOIN annotation_object_resolution aor
          ON aor.entity_type = ar.object_entity_type
         AND aor.raw_identifier_type_id = raw_object_type.identifier_type_id
         AND aor.raw_identifier = ar.object_id
        JOIN identifier_type_all resolved_object_type
          ON resolved_object_type.identifier_type_id =
             aor.canonical_identifier_type_id
        JOIN canonical_entity object
          ON object.entity_type = aor.entity_type
         AND object.canonical_identifier_type = resolved_object_type.name
         AND object.canonical_identifier = aor.canonical_identifier
        """
    )


def _bulk_load_assert_empty(
    con: duckdb.DuckDBPyConnection,
    schema: str,
) -> None:
    content_tables = (
        'data_source',
        'dataset',
        'identifier_evidence',
        'annotation',
        'entity',
        'entity_identifier',
        'entity_evidence',
        'entity_evidence_identifier',
        'entity_evidence_annotation',
        'relation',
        'relation_evidence',
        'relation_evidence_annotation',
        'entity_evidence_resolution',
        'resolution_conflict',
        'relation_evidence_relation',
        'entity_annotation_relation',
        'ontology_terms',
        'entity_ontology_relation',
        'gene_protein_representative',
        'state',
        'state_component',
        'evidence_state',
    )
    non_empty = []
    for table in content_tables:
        count = int(
            con.sql(f'SELECT COUNT(*) FROM pg.{schema}.{table}').fetchone()[0]
        )
        if count:
            non_empty.append(f'{table}={count}')
    if non_empty:
        raise RuntimeError(
            'bulk_load_parquet_to_postgres requires empty content tables: '
            + ', '.join(non_empty)
        )


def _bulk_load_small_dimensions(
    con: duckdb.DuckDBPyConnection,
    schema: str,
) -> None:
    con.execute(
        f"""
        INSERT INTO pg.{schema}.vocab_entity_type (name)
        SELECT candidate.name
        FROM (
          SELECT DISTINCT entity_type AS name FROM pq_entity_evidence WHERE entity_type IS NOT NULL
          UNION
          SELECT DISTINCT entity_type AS name FROM pq_entity WHERE entity_type IS NOT NULL
          UNION
          SELECT DISTINCT object_entity_type AS name
          FROM pq_annotation_relation_evidence
          WHERE object_entity_type IS NOT NULL
          UNION
          SELECT DISTINCT term_entity_type AS name
          FROM pq_ontology_terms
          WHERE term_id IS NOT NULL
          UNION
          SELECT DISTINCT subject_entity_type AS name
          FROM ontology_relation_raw
          WHERE subject_entity_type IS NOT NULL
          UNION
          SELECT DISTINCT object_entity_type AS name
          FROM ontology_relation_raw
          WHERE object_entity_type IS NOT NULL
        ) candidate
        LEFT JOIN pg.{schema}.vocab_entity_type existing
          ON existing.name = candidate.name
        WHERE existing.entity_type_id IS NULL
        """
    )
    con.execute(
        f"""
        INSERT INTO pg.{schema}.vocab_entity_role (name)
        SELECT candidate.name
        FROM (
          SELECT DISTINCT entity_role AS name
          FROM pq_entity_evidence
          WHERE entity_role IS NOT NULL
        ) candidate
        LEFT JOIN pg.{schema}.vocab_entity_role existing
          ON existing.name = candidate.name
        WHERE existing.entity_role_id IS NULL
        """
    )
    con.execute(
        f"""
        INSERT INTO pg.{schema}.vocab_relation_predicate (name)
        SELECT candidate.name
        FROM (
          SELECT DISTINCT predicate AS name FROM pq_relation WHERE predicate IS NOT NULL
          UNION
          SELECT DISTINCT predicate AS name FROM pq_relation_evidence WHERE predicate IS NOT NULL
          UNION
          SELECT DISTINCT predicate AS name
          FROM pq_annotation_relation_evidence
          WHERE predicate IS NOT NULL
          UNION
          SELECT DISTINCT predicate AS name
          FROM pq_entity_ontology_relation
          WHERE predicate IS NOT NULL
        ) candidate
        LEFT JOIN pg.{schema}.vocab_relation_predicate existing
          ON existing.name = candidate.name
        WHERE existing.relation_predicate_id IS NULL
        """
    )
    con.execute(
        f"""
        INSERT INTO pg.{schema}.vocab_relation_category (name)
        SELECT candidate.name
        FROM (
          SELECT DISTINCT relation_category AS name FROM pq_relation
          UNION
          SELECT DISTINCT relation_category AS name FROM pq_relation_evidence
          UNION
          SELECT DISTINCT relation_category AS name
          FROM pq_annotation_relation_evidence
        ) candidate
        LEFT JOIN pg.{schema}.vocab_relation_category existing
          ON existing.name = candidate.name
        WHERE candidate.name IS NOT NULL
          AND existing.relation_category_id IS NULL
        """
    )
    con.execute(
        f"""
        WITH missing AS (
          SELECT DISTINCT identifier_type AS name
          FROM pq_entity_identifier
          WHERE identifier_type IS NOT NULL
          UNION
          SELECT DISTINCT canonical_identifier_type AS name
          FROM pq_entity
          WHERE canonical_identifier_type IS NOT NULL
          UNION
          SELECT DISTINCT object_id_type AS name
          FROM pq_annotation_relation_evidence
          WHERE object_id_type IS NOT NULL
          UNION
          SELECT DISTINCT term_identifier_type AS name
          FROM pq_ontology_terms
          WHERE term_id IS NOT NULL
          UNION
          SELECT DISTINCT subject_identifier_type AS name
          FROM ontology_relation_raw
          WHERE subject_identifier_type IS NOT NULL
          UNION
          SELECT DISTINCT object_identifier_type AS name
          FROM ontology_relation_raw
          WHERE object_identifier_type IS NOT NULL
        ),
        missing_new AS (
          SELECT missing.name
          FROM missing
          LEFT JOIN pg.{schema}.vocab_identifier_type it
            ON it.name = missing.name
          WHERE it.identifier_type_id IS NULL
        ),
        base AS (
          SELECT coalesce(max(identifier_type_id), 0) AS max_id
          FROM pg.{schema}.vocab_identifier_type
        )
        INSERT INTO pg.{schema}.vocab_identifier_type (identifier_type_id, name)
        SELECT base.max_id + row_number() OVER (ORDER BY missing_new.name),
               missing_new.name
        FROM missing_new
        CROSS JOIN base
        """
    )
    con.execute(
        f"""
        INSERT INTO pg.{schema}.data_source (name)
        SELECT candidate.source
        FROM (
          SELECT source FROM pq_entity_evidence
          UNION
          SELECT source FROM pq_entity_evidence_resolution
          UNION
          SELECT source FROM pq_relation_evidence
          UNION
          SELECT source FROM pq_annotation_relation_evidence
          UNION
          SELECT source FROM pq_ontology_terms
          UNION
          SELECT source FROM pq_entity_ontology_relation
        ) candidate
        LEFT JOIN pg.{schema}.data_source existing
          ON existing.name = candidate.source
        WHERE candidate.source IS NOT NULL
          AND existing.source_id IS NULL
        """
    )
    con.execute(
        f"""
        INSERT INTO pg.{schema}.dataset (source_id, name)
        SELECT DISTINCT ds.source_id, candidate.dataset
        FROM (
          SELECT source, dataset FROM pq_entity_evidence
          UNION
          SELECT source, dataset FROM pq_relation_evidence
          UNION
          SELECT source, dataset FROM pq_annotation_relation_evidence
          UNION
          SELECT source, dataset FROM pq_ontology_terms
          UNION
          SELECT source, dataset FROM pq_entity_ontology_relation
        ) candidate
        JOIN pg.{schema}.data_source ds
          ON ds.name = candidate.source
        LEFT JOIN pg.{schema}.dataset existing
          ON existing.source_id = ds.source_id
         AND existing.name = candidate.dataset
        WHERE candidate.dataset IS NOT NULL
          AND existing.dataset_id IS NULL
        """
    )


def _bulk_load_materialize_dimensions(
    con: duckdb.DuckDBPyConnection,
    schema: str,
) -> None:
    dimension_tables = (
        (
            'load_data_source',
            f"""
            SELECT source_id, name
            FROM pg.{schema}.data_source
            """,
        ),
        (
            'load_dataset',
            f"""
            SELECT source_id, dataset_id, name
            FROM pg.{schema}.dataset
            """,
        ),
        (
            'load_vocab_identifier_type',
            f"""
            SELECT identifier_type_id, name, value_pattern
            FROM pg.{schema}.vocab_identifier_type
            """,
        ),
        (
            'load_vocab_entity_type',
            f"""
            SELECT entity_type_id, name
            FROM pg.{schema}.vocab_entity_type
            """,
        ),
        (
            'load_vocab_entity_role',
            f"""
            SELECT entity_role_id, name
            FROM pg.{schema}.vocab_entity_role
            """,
        ),
        (
            'load_vocab_relation_predicate',
            f"""
            SELECT relation_predicate_id, name
            FROM pg.{schema}.vocab_relation_predicate
            """,
        ),
        (
            'load_vocab_relation_category',
            f"""
            SELECT relation_category_id, name
            FROM pg.{schema}.vocab_relation_category
            """,
        ),
        (
            'load_vocab_resolution_status',
            f"""
            SELECT resolution_status_id, name
            FROM pg.{schema}.vocab_resolution_status
            """,
        ),
        (
            'load_vocab_annotation_scope',
            f"""
            SELECT annotation_scope_id, name
            FROM pg.{schema}.vocab_annotation_scope
            """,
        ),
    )
    for table, query in dimension_tables:
        con.execute(f'CREATE OR REPLACE TEMP TABLE {table} AS {query}')


def _drop_bulk_load_constraints_and_indexes(
    *,
    database_url: str,
    schema: str,
) -> None:
    with psycopg2.connect(database_url) as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT t.relname, c.conname
                FROM pg_constraint c
                JOIN pg_class t ON t.oid = c.conrelid
                JOIN pg_namespace n ON n.oid = t.relnamespace
                WHERE n.nspname = %s
                  AND c.conislocal
                  AND c.contype <> 'n'
                  AND (
                    t.relname = ANY(%s)
                    OR t.relname = ANY(
                      SELECT child.relname
                      FROM pg_inherits i
                      JOIN pg_class parent ON parent.oid = i.inhparent
                      JOIN pg_namespace parent_ns
                        ON parent_ns.oid = parent.relnamespace
                      JOIN pg_class child ON child.oid = i.inhrelid
                      JOIN pg_namespace child_ns
                        ON child_ns.oid = child.relnamespace
                      WHERE parent_ns.nspname = %s
                        AND child_ns.nspname = %s
                        AND parent.relname = ANY(%s)
                    )
                  )
                ORDER BY c.contype = 'f' DESC, t.relname, c.conname
                """,
                [
                    schema,
                    list(_BULK_COPY_CONTENT_TABLES),
                    schema,
                    schema,
                    list(_BULK_COPY_CONTENT_TABLES),
                ],
            )
            constraints = cur.fetchall()
            for table, constraint in constraints:
                cur.execute(
                    sql.SQL(
                        'ALTER TABLE {}.{} DROP CONSTRAINT IF EXISTS {} CASCADE'
                    ).format(
                        sql.Identifier(schema),
                        sql.Identifier(table),
                        sql.Identifier(constraint),
                    )
                )
            cur.execute(
                """
                WITH target_tables AS (
                  SELECT c.oid
                  FROM pg_class c
                  JOIN pg_namespace n ON n.oid = c.relnamespace
                  WHERE n.nspname = %s
                    AND c.relname = ANY(%s)
                  UNION
                  SELECT child.oid
                  FROM pg_inherits i
                  JOIN pg_class parent ON parent.oid = i.inhparent
                  JOIN pg_namespace parent_ns
                    ON parent_ns.oid = parent.relnamespace
                  JOIN pg_class child ON child.oid = i.inhrelid
                  JOIN pg_namespace child_ns
                    ON child_ns.oid = child.relnamespace
                  WHERE parent_ns.nspname = %s
                    AND child_ns.nspname = %s
                    AND parent.relname = ANY(%s)
                )
                SELECT index_class.relname
                FROM pg_index idx
                JOIN pg_class index_class ON index_class.oid = idx.indexrelid
                JOIN pg_namespace index_ns ON index_ns.oid = index_class.relnamespace
                JOIN target_tables target ON target.oid = idx.indrelid
                WHERE index_ns.nspname = %s
                  AND NOT EXISTS (
                    SELECT 1
                    FROM pg_inherits index_inherits
                    WHERE index_inherits.inhrelid = idx.indexrelid
                  )
                """,
                [
                    schema,
                    list(_BULK_COPY_CONTENT_TABLES),
                    schema,
                    schema,
                    list(_BULK_COPY_CONTENT_TABLES),
                    schema,
                ],
            )
            indexes = [row[0] for row in cur.fetchall()]
            for index in indexes:
                cur.execute(
                    sql.SQL('DROP INDEX IF EXISTS {}.{} CASCADE').format(
                        sql.Identifier(schema),
                        sql.Identifier(index),
                    )
                )


def _copy_duckdb_query_to_postgres(
    con: duckdb.DuckDBPyConnection,
    *,
    database_url: str,
    schema: str,
    table: str,
    columns: tuple[str, ...],
    query: str,
    attach_source_partition: bool = False,
    source_id: int | None = None,
) -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        csv_path = Path(tmpdir) / f'{table}.csv'
        con.execute(
            f"""
            COPY (
              {query}
            )
            TO {_sql_literal(csv_path)}
            (
              FORMAT CSV,
              HEADER false,
              DELIMITER ',',
              QUOTE '"',
              ESCAPE '"',
              NULL '\\N'
            )
            """
        )
        column_sql = sql.SQL(', ').join(
            sql.Identifier(column) for column in columns
        )
        with psycopg2.connect(database_url) as conn:
            with conn.cursor() as cur:
                target_table = sql.Identifier(table)
                if attach_source_partition:
                    if source_id is None:
                        raise ValueError(
                            'source_id is required for source partition copy'
                        )
                    staging_table = f'{table}_source_{source_id}_staging'
                    cur.execute(
                        """
                        SELECT EXISTS (
                          SELECT 1
                          FROM pg_inherits i
                          JOIN pg_class parent ON parent.oid = i.inhparent
                          JOIN pg_namespace parent_ns
                            ON parent_ns.oid = parent.relnamespace
                          JOIN pg_class child ON child.oid = i.inhrelid
                          JOIN pg_namespace child_ns
                            ON child_ns.oid = child.relnamespace
                          WHERE parent_ns.nspname = %s
                            AND parent.relname = %s
                            AND child_ns.nspname = %s
                            AND child.relname = %s
                        )
                        """,
                        [schema, table, schema, staging_table],
                    )
                    partition_attached = bool(cur.fetchone()[0])
                    should_attach = not partition_attached
                    if partition_attached:
                        target_table = sql.Identifier(staging_table)
                    else:
                        cur.execute(
                            sql.SQL(
                                'DROP TABLE IF EXISTS {}.{} CASCADE'
                            ).format(
                                sql.Identifier(schema),
                                sql.Identifier(staging_table),
                            )
                        )
                        cur.execute(
                            sql.SQL(
                                """
                                CREATE TABLE {}.{} (
                                  LIKE {}.{}
                                    INCLUDING DEFAULTS
                                    INCLUDING GENERATED
                                    INCLUDING CONSTRAINTS
                                )
                                """
                            ).format(
                                sql.Identifier(schema),
                                sql.Identifier(staging_table),
                                sql.Identifier(schema),
                                sql.Identifier(table),
                            )
                        )
                        target_table = sql.Identifier(staging_table)
                else:
                    should_attach = False
                with csv_path.open('r', encoding='utf-8') as csv_file:
                    cur.copy_expert(
                        sql.SQL(
                            """
                            COPY {}.{} ({})
                            FROM STDIN WITH (
                              FORMAT CSV,
                              NULL '\\N'
                            )
                            """
                        ).format(
                            sql.Identifier(schema),
                            target_table,
                            column_sql,
                        ),
                        csv_file,
                    )
                if attach_source_partition and should_attach:
                    cur.execute(
                        sql.SQL(
                            """
                            ALTER TABLE {}.{}
                            ATTACH PARTITION {}.{}
                            FOR VALUES IN ({})
                            """
                        ).format(
                            sql.Identifier(schema),
                            sql.Identifier(table),
                            sql.Identifier(schema),
                            sql.Identifier(staging_table),
                            sql.Literal(source_id),
                        )
                    )


def _copy_source_partition(
    con: duckdb.DuckDBPyConnection,
    *,
    database_url: str,
    schema: str,
    table: str,
    columns: tuple[str, ...],
    query: str,
    source_id: int,
) -> None:
    _copy_duckdb_query_to_postgres(
        con,
        database_url=database_url,
        schema=schema,
        table=table,
        columns=columns,
        query=query,
        attach_source_partition=True,
        source_id=source_id,
    )


def _duckdb_source_id(
    con: duckdb.DuckDBPyConnection,
    source: str | None = None,
) -> int:
    if source is None:
        rows = con.execute(
            """
            WITH current_source AS (
              SELECT source FROM pq_entity_evidence WHERE source IS NOT NULL
              UNION
              SELECT source FROM pq_entity_evidence_resolution WHERE source IS NOT NULL
              UNION
              SELECT source FROM pq_entity_identifier_resolved WHERE source IS NOT NULL
              UNION
              SELECT source FROM pq_relation_evidence WHERE source IS NOT NULL
              UNION
              SELECT source FROM pq_annotation_relation_evidence WHERE source IS NOT NULL
              UNION
              SELECT source FROM pq_ontology_terms WHERE source IS NOT NULL
              UNION
              SELECT source FROM pq_entity_ontology_relation WHERE source IS NOT NULL
            )
            SELECT ds.source_id, ds.name
            FROM current_source s
            JOIN load_data_source ds
              ON ds.name = s.source
            ORDER BY ds.source_id
            """
        ).fetchall()
        if len(rows) != 1:
            raise ValueError(
                'Expected exactly one source for source-partition COPY load; '
                f'found {len(rows)}: {rows!r}'
            )
        return int(rows[0][0])
    row = con.execute(
        'SELECT source_id FROM load_data_source WHERE name = ?',
        [source],
    ).fetchone()
    if row is None:
        raise ValueError(f'No source_id found for {source!r}')
    return int(row[0])


def _bulk_copy_evidence(
    con: duckdb.DuckDBPyConnection,
    *,
    schema: str,
    database_url: str,
) -> None:
    source_id = _duckdb_source_id(con)
    existing_identifier_evidence = _duckdb_pg_table(
        schema,
        'identifier_evidence',
    )
    # entity_identifier_raw is created with identifier_normalized (T040) on
    # every path that builds it -- but multigene_split.explode_multi_gene_-
    # protein_mentions used to rebuild the table via an explicit column list
    # that predated T040 and dropped this column whenever a shard had a
    # multi-gene UniProt mention (fixed; see multigene_split._explode_one).
    # Kept as a defensive check rather than assuming the fix covers every
    # future rebuild path: degrade to NULL rather than crash the whole load.
    has_identifier_normalized = bool(
        con.execute(
            """
            SELECT COUNT(*) FROM information_schema.columns
            WHERE table_name = 'entity_identifier_raw'
              AND column_name = 'identifier_normalized'
            """
        ).fetchone()[0]
    )
    identifier_normalized_sql = (
        'i.identifier_normalized' if has_identifier_normalized
        else 'NULL::VARCHAR'
    )
    _copy_duckdb_query_to_postgres(
        con,
        database_url=database_url,
        schema=schema,
        table='identifier_evidence',
        columns=(
            'identifier_id', 'identifier_type_id', 'value', 'value_normalized',
        ),
        query="""
          SELECT DISTINCT
            i.identifier_id::UUID,
            it.identifier_type_id,
            i.identifier,
            {identifier_normalized_sql}
          FROM pq_entity_identifier i
          JOIN load_vocab_identifier_type it
            ON it.name = i.identifier_type
          LEFT JOIN {existing_identifier_evidence} existing
            ON existing.identifier_type_id = it.identifier_type_id
           AND existing.value = i.identifier
          WHERE i.identifier_id IS NOT NULL
            AND i.identifier IS NOT NULL
            AND existing.identifier_id IS NULL
        """.format(
            existing_identifier_evidence=existing_identifier_evidence,
            identifier_normalized_sql=identifier_normalized_sql,
        ),
    )
    existing_annotation = _duckdb_pg_table(schema, 'annotation')
    _copy_duckdb_query_to_postgres(
        con,
        database_url=database_url,
        schema=schema,
        table='annotation',
        columns=('annotation_key', 'term', 'value', 'unit'),
        query="""
          SELECT DISTINCT
            pq_annotation.annotation_key::UUID,
            pq_annotation.term,
            pq_annotation.value,
            pq_annotation.unit
          FROM pq_annotation
          LEFT JOIN {existing_annotation} existing
            ON existing.annotation_key = pq_annotation.annotation_key::UUID
          WHERE pq_annotation.term IS NOT NULL
            AND existing.annotation_key IS NULL
        """.format(existing_annotation=existing_annotation),
    )
    _copy_source_partition(
        con,
        database_url=database_url,
        schema=schema,
        table='entity_evidence',
        columns=(
            'source_id',
            'entity_evidence_id',
            'dataset_id',
            'row_id',
            'parent_entity_evidence_id',
            'entity_role_id',
            'entity_type_id',
            'taxonomy_id',
        ),
        query="""
          SELECT DISTINCT
            ds.source_id,
            e.entity_evidence_id::UUID,
            d.dataset_id,
            e.row_id,
            e.parent_entity_evidence_id::UUID,
            er.entity_role_id,
            et.entity_type_id,
            NULLIF(e.taxonomy_id, '')::BIGINT
          FROM pq_entity_evidence e
          JOIN load_data_source ds
            ON ds.name = e.source
          JOIN load_dataset d
            ON d.source_id = ds.source_id
           AND d.name = e.dataset
          JOIN load_vocab_entity_role er
            ON er.name = e.entity_role
          LEFT JOIN load_vocab_entity_type et
            ON et.name = e.entity_type
        """,
        source_id=source_id,
    )
    _copy_source_partition(
        con,
        database_url=database_url,
        schema=schema,
        table='entity_evidence_identifier',
        columns=('source_id', 'entity_evidence_id', 'identifier_id'),
        query=f"""
          SELECT DISTINCT
            ds.source_id,
            i.entity_evidence_id::UUID,
            coalesce(existing.identifier_id, i.identifier_id::UUID) AS identifier_id
          FROM pq_entity_identifier i
          JOIN load_data_source ds
            ON ds.name = i.source
          JOIN load_vocab_identifier_type it
            ON it.name = i.identifier_type
          LEFT JOIN {existing_identifier_evidence} existing
            ON existing.identifier_type_id = it.identifier_type_id
           AND existing.value = i.identifier
          WHERE i.identifier_id IS NOT NULL
        """,
        source_id=source_id,
    )
    _copy_source_partition(
        con,
        database_url=database_url,
        schema=schema,
        table='entity_evidence_annotation',
        columns=('source_id', 'entity_evidence_id', 'annotation_key'),
        query="""
          SELECT DISTINCT
            ds.source_id,
            a.evidence_id::UUID,
            a.annotation_key::UUID
          FROM pq_entity_annotation a
          JOIN load_data_source ds
            ON ds.name = a.source
        """,
        source_id=source_id,
    )


def _bulk_copy_canonical(
    con: duckdb.DuckDBPyConnection,
    *,
    schema: str,
    database_url: str,
) -> None:
    source_id = _duckdb_source_id(con)
    _copy_duckdb_query_to_postgres(
        con,
        database_url=database_url,
        schema=schema,
        table='entity',
        columns=(
            'entity_id',
            'entity_type_id',
            'taxonomy_id',
            'canonical_identifier_type_id',
            'canonical_identifier',
            'resolution_status_id',
            'resolution_mechanism',
        ),
        query="""
          SELECT
            e.entity_id,
            et.entity_type_id,
            NULLIF(e.taxonomy_id, '')::BIGINT,
            it.identifier_type_id,
            e.canonical_identifier,
            rs.resolution_status_id,
            e.resolution_mechanism
          FROM pq_entity e
          JOIN load_vocab_entity_type et
            ON et.name = e.entity_type
          JOIN load_vocab_identifier_type it
            ON it.name = e.canonical_identifier_type
          JOIN load_vocab_resolution_status rs
            ON rs.name = e.resolution_status
          LEFT JOIN {existing_entity} existing
            ON existing.entity_id = e.entity_id
          WHERE existing.entity_id IS NULL
        """.format(existing_entity=_duckdb_pg_table(schema, 'entity')),
    )
    existing_identifier_evidence = _duckdb_pg_table(
        schema,
        'identifier_evidence',
    )
    _copy_duckdb_query_to_postgres(
        con,
        database_url=database_url,
        schema=schema,
        table='identifier_evidence',
        columns=(
            'identifier_id', 'identifier_type_id', 'value', 'value_normalized',
        ),
        query=f"""
          SELECT DISTINCT
            content_uuid(
              'identifier' || chr(31) ||
              i.identifier_type || chr(31) ||
              i.identifier
            ) AS identifier_id,
            it.identifier_type_id,
            i.identifier,
            -- This table's values already come from the resolver's own key
            -- (chemical source_ids are normalized upstream, T033/T034) or
            -- from the canonical identifier itself, so no transform is
            -- needed here -- only a check that the value is already in its
            -- namespace's canonical form (spec 011 T040/T041).
            CASE WHEN it.value_pattern IS NOT NULL
                  AND i.identifier ~ it.value_pattern
                 THEN i.identifier
            END AS value_normalized
          FROM pq_entity_identifier_resolved i
          JOIN load_vocab_identifier_type it
            ON it.name = i.identifier_type
          LEFT JOIN {existing_identifier_evidence} existing
            ON existing.identifier_type_id = it.identifier_type_id
           AND existing.value = i.identifier
          WHERE i.identifier IS NOT NULL
            AND i.identifier <> ''
            AND existing.identifier_id IS NULL
        """,
    )
    existing_entity_identifier = _duckdb_pg_table(
        schema,
        'entity_identifier',
    )
    _copy_source_partition(
        con,
        database_url=database_url,
        schema=schema,
        table='entity_identifier',
        columns=('source_id', 'entity_id', 'identifier_id'),
        query=f"""
          SELECT candidate.*
          FROM (
            SELECT DISTINCT
              ds.source_id,
              i.entity_id::UUID AS entity_id,
              coalesce(
                ie.identifier_id,
                content_uuid(
                  'identifier' || chr(31) ||
                  i.identifier_type || chr(31) ||
                  i.identifier
                )
              ) AS identifier_id
            FROM pq_entity_identifier_resolved i
            JOIN load_data_source ds
              ON ds.name = i.source
            JOIN load_vocab_identifier_type it
              ON it.name = i.identifier_type
            LEFT JOIN {existing_identifier_evidence} ie
              ON ie.identifier_type_id = it.identifier_type_id
             AND ie.value = i.identifier
            WHERE i.identifier IS NOT NULL
              AND i.identifier <> ''
          ) candidate
          LEFT JOIN {existing_entity_identifier} existing
            ON existing.source_id = candidate.source_id
           AND existing.entity_id = candidate.entity_id
           AND existing.identifier_id = candidate.identifier_id
          WHERE existing.source_id IS NULL
        """,
        source_id=source_id,
    )
    _copy_source_partition(
        con,
        database_url=database_url,
        schema=schema,
        table='ontology_terms',
        columns=(
            'source_id',
            'term_entity_id',
            'term_id',
            'ontology_prefix',
            'label',
            'definition',
            'ontology_id',
            'synonyms',
            'synonyms_text',
            'sources',
        ),
        query="""
          SELECT
            ds.source_id,
            ce.entity_id,
            ot.term_id,
            ot.ontology_prefix,
            ot.label,
            ot.definition,
            ot.ontology_id,
            COALESCE(
              '{{' || array_to_string(
                list_transform(
                  ot.synonyms,
                  x -> chr(34)
                    || replace(
                         replace(x, chr(92), chr(92) || chr(92)),
                         chr(34),
                         chr(92) || chr(34)
                       )
                    || chr(34)
                ),
                ','
              ) || '}}',
              '{{}}'
            ),
            COALESCE(ot.synonyms_text, ''),
            COALESCE(
              '{{' || array_to_string(
                list_transform(
                  ot.sources,
                  x -> chr(34)
                    || replace(
                         replace(x, chr(92), chr(92) || chr(92)),
                         chr(34),
                         chr(92) || chr(34)
                       )
                    || chr(34)
                ),
                ','
              ) || '}}',
              '{{}}'
            )
          FROM pq_ontology_terms ot
          JOIN load_data_source ds
            ON ds.name = ot.source
          JOIN ontology_term_resolution otr
            ON otr.source = ot.source
           AND otr.ontology_id = ot.ontology_id
           AND otr.term_id = ot.term_id
          JOIN canonical_entity ce
            ON ce.entity_type = otr.entity_type
           AND ce.taxonomy_id IS NOT DISTINCT FROM otr.taxonomy_id
           AND ce.canonical_identifier_type_id =
               otr.canonical_identifier_type_id
           AND ce.canonical_identifier = otr.canonical_identifier
          WHERE ot.term_id IS NOT NULL
        """,
        source_id=source_id,
    )
    _copy_source_partition(
        con,
        database_url=database_url,
        schema=schema,
        table='entity_ontology_relation',
        columns=(
            'source_id',
            'subject_entity_id',
            'predicate_id',
            'object_entity_id',
            'ontology_id',
        ),
        query="""
          SELECT candidate.*
          FROM (
            SELECT DISTINCT
              ds.source_id,
              eor.subject_entity_id,
              rp.relation_predicate_id,
              eor.object_entity_id,
              eor.ontology_id
            FROM pq_entity_ontology_relation eor
            JOIN load_data_source ds
              ON ds.name = eor.source
            JOIN load_vocab_relation_predicate rp
              ON rp.name = eor.predicate
          ) candidate
          LEFT JOIN {existing_entity_ontology_relation} existing
            ON existing.source_id = candidate.source_id
           AND existing.subject_entity_id = candidate.subject_entity_id
           AND existing.predicate_id = candidate.relation_predicate_id
           AND existing.object_entity_id = candidate.object_entity_id
           AND existing.ontology_id = candidate.ontology_id
          WHERE existing.source_id IS NULL
        """.format(
            existing_entity_ontology_relation=_duckdb_pg_table(
                schema,
                'entity_ontology_relation',
            ),
        ),
        source_id=source_id,
    )
    _copy_source_partition(
        con,
        database_url=database_url,
        schema=schema,
        table='relation_evidence',
        columns=(
            'source_id',
            'relation_evidence_id',
            'dataset_id',
            'row_id',
            'subject_entity_evidence_id',
            'subject_entity_id',
            'predicate_id',
            'object_entity_evidence_id',
            'object_entity_id',
            'relation_category_id',
        ),
        query="""
          SELECT DISTINCT *
          FROM (
            SELECT
              ds.source_id,
              r.relation_evidence_id::UUID,
              d.dataset_id,
              r.row_id,
              r.subject_entity_evidence_id::UUID,
              NULL::UUID AS subject_entity_id,
              rp.relation_predicate_id,
              r.object_entity_evidence_id::UUID,
              NULL::UUID AS object_entity_id,
              rc.relation_category_id
            FROM pq_relation_evidence r
            JOIN load_data_source ds
              ON ds.name = r.source
            JOIN load_dataset d
              ON d.source_id = ds.source_id
             AND d.name = r.dataset
            JOIN load_vocab_relation_predicate rp
              ON rp.name = r.predicate
            JOIN load_vocab_relation_category rc
              ON rc.name = r.relation_category
            UNION ALL
            SELECT
              ds.source_id,
              ar.relation_evidence_id::UUID,
              d.dataset_id,
              ar.row_id,
              NULL::UUID AS subject_entity_evidence_id,
              ar.object_entity_id AS subject_entity_id,
              rp.relation_predicate_id,
              ar.subject_entity_evidence_id::UUID AS object_entity_evidence_id,
              NULL::UUID AS object_entity_id,
              rc.relation_category_id
            FROM pq_annotation_relation_evidence_resolved ar
            JOIN load_data_source ds
              ON ds.name = ar.source
            JOIN load_dataset d
              ON d.source_id = ds.source_id
             AND d.name = ar.dataset
            JOIN load_vocab_relation_predicate rp
              ON rp.name = ar.predicate
            JOIN load_vocab_relation_category rc
              ON rc.name = ar.relation_category
          )
        """,
        source_id=source_id,
    )
    _copy_source_partition(
        con,
        database_url=database_url,
        schema=schema,
        table='relation_evidence_annotation',
        columns=(
            'source_id',
            'relation_evidence_id',
            'annotation_key',
            'annotation_scope_id',
        ),
        query="""
          SELECT DISTINCT
            ds.source_id,
            a.evidence_id::UUID,
            a.annotation_key::UUID,
            sc.annotation_scope_id
          FROM pq_relation_annotation a
          JOIN load_data_source ds
            ON ds.name = a.source
          JOIN load_vocab_annotation_scope sc
            ON sc.name = coalesce(a.annotation_scope, 'relation')
        """,
        source_id=source_id,
    )
    _copy_source_partition(
        con,
        database_url=database_url,
        schema=schema,
        table='entity_evidence_resolution',
        columns=(
            'source_id',
            'entity_evidence_id',
            'status_id',
            'entity_id',
            'reason_id',
            'resolved_at',
            'molecular_type_id',
        ),
        query="""
          SELECT
            ds.source_id,
            er.entity_evidence_id::UUID,
            rs.resolution_status_id,
            er.entity_id,
            er.reason_id::SMALLINT,
            now(),
            er.molecular_type_id::SMALLINT
          FROM pq_entity_evidence_resolution er
          JOIN load_data_source ds
            ON ds.name = er.source
          JOIN load_vocab_resolution_status rs
            ON rs.name = er.status
          WHERE er.entity_id IS NOT NULL
        """,
        source_id=source_id,
    )
    _copy_source_partition(
        con,
        database_url=database_url,
        schema=schema,
        table='resolution_conflict',
        columns=(
            'source_id',
            'entity_evidence_id',
            'identifier_type_id',
            'value_normalized',
            'candidate_structure',
            'candidate_source_id',
            'candidate_role',
            'skeleton',
        ),
        query="""
          SELECT
            ds.source_id,
            rc.entity_evidence_id::UUID,
            it.identifier_type_id,
            rc.value_normalized,
            rc.candidate_structure,
            cs.source_id AS candidate_source_id,
            rc.candidate_role,
            rc.skeleton
          FROM pq_resolution_conflict rc
          JOIN load_data_source ds
            ON ds.name = rc.source
          JOIN load_vocab_identifier_type it
            ON it.name = rc.identifier_type
          JOIN load_data_source cs
            ON cs.name = rc.candidate_source
        """,
        source_id=source_id,
    )
    _copy_duckdb_query_to_postgres(
        con,
        database_url=database_url,
        schema=schema,
        table='relation',
        columns=(
            'relation_id',
            'subject_entity_id',
            'predicate_id',
            'object_entity_id',
            'relation_category_id',
        ),
        query="""
          SELECT
            r.relation_id,
            r.subject_entity_id,
            rp.relation_predicate_id,
            r.object_entity_id,
            rc.relation_category_id
          FROM pq_relation r
          JOIN load_vocab_relation_predicate rp
            ON rp.name = r.predicate
          LEFT JOIN load_vocab_relation_category rc
            ON rc.name = r.relation_category
          LEFT JOIN {existing_relation} existing
            ON existing.relation_id = r.relation_id
          WHERE existing.relation_id IS NULL
        """.format(existing_relation=_duckdb_pg_table(schema, 'relation')),
    )
    _copy_source_partition(
        con,
        database_url=database_url,
        schema=schema,
        table='relation_evidence_relation',
        columns=('source_id', 'relation_id', 'relation_evidence_id'),
        query="""
          SELECT
            ds.source_id,
            rer.relation_id,
            rer.relation_evidence_id::UUID
          FROM pq_relation_evidence_relation rer
          JOIN load_data_source ds
            ON ds.name = rer.source
        """,
        source_id=source_id,
    )
    # gene_protein_representative: global 1:1 table, copied like `entity`. The
    # staged pipeline canonicalises per source-shard, so a gene entity recurs
    # across shards — dedup against the rows already in Postgres (mirrors the
    # `entity` copy). uniprot_all is rendered as a Postgres array literal so
    # the CSV round-trip parses back into text[] (UniProt ACs need no element
    # quoting).
    existing_gene_protein_representative = _duckdb_pg_table(
        schema,
        'gene_protein_representative',
    )
    _copy_duckdb_query_to_postgres(
        con,
        database_url=database_url,
        schema=schema,
        table='gene_protein_representative',
        columns=(
            'entity_id',
            'representative_uniprot',
            'is_reviewed',
            'uniprot_all',
        ),
        query=f"""
          SELECT
            gpr.entity_id,
            gpr.representative_uniprot,
            gpr.is_reviewed,
            CASE
              WHEN gpr.uniprot_all IS NULL THEN NULL
              ELSE '{{' || array_to_string(gpr.uniprot_all, ',') || '}}'
            END AS uniprot_all
          FROM pq_gene_protein_representative gpr
          LEFT JOIN {existing_gene_protein_representative} existing
            ON existing.entity_id = gpr.entity_id
          WHERE existing.entity_id IS NULL
        """,
    )
    # state / state_component: content-hashed global tables — identical (gene,
    # uniprot, isoform) forms across shards/sources collapse to one state.
    # Dedup against the rows already in Postgres (like `entity`).
    existing_state = _duckdb_pg_table(schema, 'state')
    _copy_duckdb_query_to_postgres(
        con,
        database_url=database_url,
        schema=schema,
        table='state',
        columns=('state_id', 'gene_entity_id', 'molecular_type_id'),
        query=f"""
          SELECT
            s.state_id,
            s.gene_entity_id,
            s.molecular_type_id::SMALLINT
          FROM pq_state s
          LEFT JOIN {existing_state} existing
            ON existing.state_id = s.state_id
          WHERE existing.state_id IS NULL
        """,
    )
    existing_state_component = _duckdb_pg_table(schema, 'state_component')
    _copy_duckdb_query_to_postgres(
        con,
        database_url=database_url,
        schema=schema,
        table='state_component',
        columns=('state_id', 'component_type', 'value'),
        query=f"""
          SELECT
            sc.state_id,
            sc.component_type,
            sc.value
          FROM pq_state_component sc
          LEFT JOIN {existing_state_component} existing
            ON existing.state_id = sc.state_id
           AND existing.component_type = sc.component_type
           AND existing.value = sc.value
          WHERE existing.state_id IS NULL
        """,
    )
    # evidence_state: one-to-many assignment of states to source records —
    # partitioned by source_id like entity_evidence_resolution.
    _copy_source_partition(
        con,
        database_url=database_url,
        schema=schema,
        table='evidence_state',
        columns=('source_id', 'entity_evidence_id', 'state_id'),
        query="""
          SELECT
            ds.source_id,
            es.entity_evidence_id::UUID,
            es.state_id
          FROM pq_evidence_state es
          JOIN load_data_source ds
            ON ds.name = es.source
        """,
        source_id=source_id,
    )


def _reset_postgres_sequences(
    *,
    database_url: str,
    schema: str,
) -> None:
    sequence_tables = (
        ('data_source', 'source_id'),
        ('dataset', 'dataset_id'),
    )
    with psycopg2.connect(database_url) as conn:
        with conn.cursor() as cur:
            for table, column in sequence_tables:
                cur.execute(
                    sql.SQL(
                        """
                        SELECT setval(
                          pg_get_serial_sequence(%s, %s),
                          COALESCE((SELECT MAX({column}) FROM {table}), 1),
                          COALESCE((SELECT MAX({column}) FROM {table}), 0) > 0
                        )
                        """
                    ).format(
                        column=sql.Identifier(column),
                        table=sql.SQL('{}.{}').format(
                            sql.Identifier(schema),
                            sql.Identifier(table),
                        ),
                    ),
                    [f'{schema}.{table}', column],
                )
        conn.commit()


__all__ = [
    '_sql_literal',
    '_create_duckdb_content_uuid_macro',
    '_create_duckdb_resolver_views',
    '_create_duckdb_identifier_type_all_view',
    'DuckDBEvidenceProjector',
    '_ENTITY_EVIDENCE_SCHEMA',
    '_ENTITY_IDENTIFIER_SCHEMA',
    '_ANNOTATION_REF_SCHEMA',
    '_ANNOTATION_VALUE_SCHEMA',
    '_RELATION_EVIDENCE_SCHEMA',
    '_ANNOTATION_RELATION_EVIDENCE_SCHEMA',
    '_create_duckdb_evidence_tables',
    '_ensure_duckdb_canonical_caches',
    '_drop_duckdb_batch_tables',
    '_canonicalize_loaded_duckdb',
    '_bulk_load_create_views_from_loaded_tables',
    '_bulk_load_assert_empty',
    '_bulk_load_small_dimensions',
    '_bulk_load_materialize_dimensions',
    '_drop_bulk_load_constraints_and_indexes',
    '_copy_duckdb_query_to_postgres',
    '_copy_source_partition',
    '_duckdb_source_id',
    '_bulk_copy_evidence',
    '_bulk_copy_canonical',
    '_reset_postgres_sequences',
]
