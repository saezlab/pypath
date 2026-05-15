"""
ChEMBL-specific data parsers.
"""

from __future__ import annotations

import os
from pathlib import Path
import shutil
from typing import Any, Generator

try:
    import duckdb
except ImportError:  # pragma: no cover - optional dependency
    duckdb = None

from pypath.inputs_v2.parsers.base import (
    _extract_sqlite_from_opener,
    iter_parquet,
    iter_sqlite,
)

CHEMBL_ACTIVITIES_PARQUET_CHUNK_SIZE = int(
    os.environ.get('OMNIPATH_CHEMBL_ACTIVITIES_PARQUET_CHUNK_SIZE', '50000')
)


DUCKDB_TABLES: dict[str, str] = {
    'molecule_dictionary': """
        SELECT
            molregno, chembl_id, pref_name, molecule_type, max_phase,
            therapeutic_flag, withdrawn_flag, polymer_flag, inorganic_flag,
            natural_product
        FROM s.molecule_dictionary
    """,
    'compound_structures': """
        SELECT molregno, canonical_smiles, standard_inchi, standard_inchi_key
        FROM s.compound_structures
    """,
    'compound_properties': """
        SELECT molregno, full_mwt, full_molformula, alogp
        FROM s.compound_properties
    """,
    'assays': """
        SELECT
            assay_id, tid, doc_id, chembl_id, description, assay_type,
            assay_tax_id, confidence_score, assay_category,
            assay_subcellular_fraction, assay_tissue, assay_cell_type
        FROM s.assays
    """,
    'assay_parameters': """
        SELECT assay_id, standard_type, standard_value, standard_units
        FROM s.assay_parameters
    """,
    'target_dictionary': """
        SELECT tid, target_type, pref_name, tax_id, organism, chembl_id
        FROM s.target_dictionary
    """,
    'target_components': """
        SELECT tid, component_id
        FROM s.target_components
    """,
    'component_sequences': """
        SELECT component_id, component_type, db_source, accession, description
        FROM s.component_sequences
    """,
    'docs': """
        SELECT
            doc_id, chembl_id, pubmed_id, doi, patent_id, journal, year,
            title, abstract, doc_type
        FROM s.docs
    """,
    'drug_mechanism': """
        SELECT *
        FROM s.drug_mechanism
    """,
    'compound_records': """
        SELECT record_id, doc_id
        FROM s.compound_records
    """,
    'mechanism_refs': """
        SELECT mec_id, ref_url
        FROM s.mechanism_refs
    """,
    'activities': """
        SELECT
            activity_id, assay_id, molregno, doc_id, standard_type,
            standard_relation, standard_value, pchembl_value,
            data_validity_comment, action_type
        FROM s.activities
    """,
    'ligand_eff': """
        SELECT activity_id, bei, le, lle, sei
        FROM s.ligand_eff
    """,
    'action_type': """
        SELECT action_type, description, parent_type
        FROM s.action_type
    """,
}


PARQUET_QUERIES: dict[str, str] = {
    'molecules': """
        SELECT
            md.molregno,
            md.chembl_id,
            md.pref_name,
            md.molecule_type,
            md.max_phase,
            md.therapeutic_flag,
            md.withdrawn_flag,
            md.polymer_flag,
            md.inorganic_flag,
            md.natural_product,
            NULL AS chebi_id,
            cs.canonical_smiles,
            cs.standard_inchi,
            cs.standard_inchi_key,
            cp.full_mwt,
            cp.full_molformula,
            cp.alogp
        FROM molecule_dictionary md
        LEFT JOIN compound_structures cs ON md.molregno = cs.molregno
        LEFT JOIN compound_properties cp ON md.molregno = cp.molregno
    """,
    'assays': """
        WITH assay_params_agg AS (
            SELECT
                assay_id,
                GROUP_CONCAT(
                    COALESCE(standard_type, '') || ':' ||
                    COALESCE(CAST(standard_value AS VARCHAR), '') || ' ' ||
                    COALESCE(standard_units, ''),
                    '; '
                ) AS parameters
            FROM assay_parameters
            GROUP BY assay_id
        )
        SELECT
            a.*,
            t.chembl_id AS target_chembl_id,
            d.chembl_id AS document_chembl_id,
            apa.parameters
        FROM assays a
        LEFT JOIN target_dictionary t ON a.tid = t.tid
        LEFT JOIN docs d ON a.doc_id = d.doc_id
        LEFT JOIN assay_params_agg apa ON a.assay_id = apa.assay_id
    """,
    'documents': """
        SELECT *
        FROM docs
    """,
    'targets': """
        SELECT
            td.tid,
            td.target_type,
            td.pref_name,
            td.tax_id,
            td.organism,
            td.chembl_id,
            GROUP_CONCAT(CASE
                WHEN cs.component_type = 'PROTEIN'
                 AND cs.db_source IN ('SWISS-PROT', 'TREMBL')
                THEN COALESCE(cs.accession, '')
                ELSE ''
            END) AS component_uniprot_accessions,
            GROUP_CONCAT(CASE
                WHEN cs.accession LIKE 'ENS%'
                THEN COALESCE(cs.accession, '')
                ELSE ''
            END) AS component_ensembl_accessions,
            GROUP_CONCAT(COALESCE(cs.component_type, '')) AS component_types,
            GROUP_CONCAT(COALESCE(CAST(tc.component_id AS VARCHAR), '')) AS component_ids,
            GROUP_CONCAT(COALESCE(cs.description, '')) AS component_descriptions
        FROM target_dictionary td
        LEFT JOIN target_components tc ON td.tid = tc.tid
        LEFT JOIN component_sequences cs ON tc.component_id = cs.component_id
        GROUP BY
            td.tid, td.target_type, td.pref_name, td.tax_id,
            td.organism, td.chembl_id
    """,
    'mechanisms': """
        WITH target_components_agg AS (
            SELECT
                tc.tid,
                GROUP_CONCAT(DISTINCT CASE
                    WHEN cseq.component_type = 'PROTEIN'
                     AND cseq.db_source IN ('SWISS-PROT', 'TREMBL')
                    THEN cseq.accession
                END) AS target_component_uniprot_accessions,
                GROUP_CONCAT(DISTINCT CASE
                    WHEN cseq.accession LIKE 'ENS%'
                    THEN cseq.accession
                END) AS target_component_ensembl_accessions
            FROM target_components tc
            LEFT JOIN component_sequences cseq ON tc.component_id = cseq.component_id
            GROUP BY tc.tid
        ),
        mechanism_refs_agg AS (
            SELECT mec_id, GROUP_CONCAT(ref_url, '; ') AS mechanism_refs
            FROM mechanism_refs
            GROUP BY mec_id
        )
        SELECT
            dm.*,
            md.chembl_id AS molecule_chembl_id,
            cs.canonical_smiles,
            cs.standard_inchi,
            cs.standard_inchi_key,
            td.chembl_id AS target_chembl_id,
            td.target_type,
            tca.target_component_uniprot_accessions,
            tca.target_component_ensembl_accessions,
            d.chembl_id AS document_chembl_id,
            d.pubmed_id,
            d.doi,
            d.title AS document_title,
            mra.mechanism_refs
        FROM drug_mechanism dm
        LEFT JOIN molecule_dictionary md ON dm.molregno = md.molregno
        LEFT JOIN compound_structures cs ON dm.molregno = cs.molregno
        LEFT JOIN target_dictionary td ON dm.tid = td.tid
        LEFT JOIN target_components_agg tca ON dm.tid = tca.tid
        LEFT JOIN compound_records cr ON dm.record_id = cr.record_id
        LEFT JOIN docs d ON cr.doc_id = d.doc_id
        LEFT JOIN mechanism_refs_agg mra ON dm.mec_id = mra.mec_id
    """,
    'activities': """
        WITH assay_params_agg AS (
            SELECT
                assay_id,
                GROUP_CONCAT(
                    COALESCE(standard_type, '') || ':' ||
                    COALESCE(CAST(standard_value AS VARCHAR), '') || ' ' ||
                    COALESCE(standard_units, ''),
                    '; '
                ) AS assay_parameters
            FROM assay_parameters
            GROUP BY assay_id
        ),
        target_components_agg AS (
            SELECT
                tc.tid,
                GROUP_CONCAT(DISTINCT CASE
                    WHEN cseq.component_type = 'PROTEIN'
                     AND cseq.db_source IN ('SWISS-PROT', 'TREMBL')
                    THEN cseq.accession
                END) AS target_component_uniprot_accessions,
                GROUP_CONCAT(DISTINCT CASE
                    WHEN cseq.accession LIKE 'ENS%'
                    THEN cseq.accession
                END) AS target_component_ensembl_accessions
            FROM target_components tc
            LEFT JOIN component_sequences cseq ON tc.component_id = cseq.component_id
            GROUP BY tc.tid
        )
        SELECT
            act.activity_id,
            act.assay_id,
            act.molregno,
            act.doc_id,
            act.standard_type,
            act.standard_relation,
            act.standard_value,
            act.pchembl_value,
            act.data_validity_comment,
            act.action_type,
            md.chembl_id AS molecule_chembl_id,
            cs.canonical_smiles,
            cs.standard_inchi,
            cs.standard_inchi_key,
            a.chembl_id AS assay_chembl_id,
            a.assay_type,
            a.assay_tax_id,
            a.confidence_score,
            a.assay_category,
            a.assay_subcellular_fraction,
            a.assay_tissue,
            a.assay_cell_type,
            a.description AS assay_description,
            apa.assay_parameters,
            td.chembl_id AS target_chembl_id,
            td.target_type,
            td.tax_id AS target_tax_id,
            tca.target_component_uniprot_accessions,
            tca.target_component_ensembl_accessions,
            d.chembl_id AS document_chembl_id,
            d.pubmed_id,
            d.doi,
            le.bei,
            le.le,
            le.lle,
            le.sei,
            ax.description AS action_description,
            ax.parent_type AS action_parent_type
        FROM activities act
        LEFT JOIN molecule_dictionary md ON act.molregno = md.molregno
        LEFT JOIN compound_structures cs ON act.molregno = cs.molregno
        LEFT JOIN assays a ON act.assay_id = a.assay_id
        LEFT JOIN assay_params_agg apa ON a.assay_id = apa.assay_id
        LEFT JOIN target_dictionary td ON a.tid = td.tid
        LEFT JOIN target_components_agg tca ON a.tid = tca.tid
        LEFT JOIN docs d ON act.doc_id = d.doc_id
        LEFT JOIN ligand_eff le ON act.activity_id = le.activity_id
        LEFT JOIN action_type ax ON act.action_type = ax.action_type
    """,
}


def _chembl_cache_paths(
    sqlite_path: Path | None,
    duckdb_path: Path | None = None,
    parquet_dir: Path | None = None,
) -> tuple[Path, Path]:
    if sqlite_path is None:
        raise ValueError("sqlite_path is required for ChEMBL DuckDB/Parquet cache.")

    if duckdb_path is None:
        duckdb_path = sqlite_path.with_suffix('.duckdb')
    if parquet_dir is None:
        parquet_dir = sqlite_path.parent / f'{sqlite_path.stem}_parquet'

    return Path(duckdb_path), Path(parquet_dir)


def _ensure_sqlite(opener: Any, sqlite_path: Path | None, db_rel_path: str | None) -> None:
    if sqlite_path is not None and not sqlite_path.exists():
        if not opener or not opener.result:
            raise FileNotFoundError(f"ChEMBL SQLite cache not found: {sqlite_path}")
        _extract_sqlite_from_opener(opener, sqlite_path, db_rel_path)


def _ensure_duckdb_cache(sqlite_path: Path, duckdb_path: Path) -> None:
    if duckdb_path.exists():
        return
    if duckdb is None:
        raise ImportError("duckdb is required to create the ChEMBL cache.")

    tmp_path = duckdb_path.with_suffix(f'{duckdb_path.suffix}.tmp')
    if tmp_path.exists():
        tmp_path.unlink()

    print(f"Creating ChEMBL DuckDB cache: {duckdb_path}")
    duckdb_path.parent.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(tmp_path))
    try:
        con.execute('LOAD sqlite')
        con.execute(f"ATTACH '{sqlite_path}' AS s (TYPE sqlite, READ_ONLY)")
        for table_name, select_sql in DUCKDB_TABLES.items():
            con.execute(f"CREATE TABLE {table_name} AS {select_sql}")
        con.execute('CHECKPOINT')
    finally:
        con.close()

    tmp_path.replace(duckdb_path)


def _ensure_parquet_dataset(dataset: str, duckdb_path: Path, parquet_dir: Path) -> Path:
    if duckdb is None:
        raise ImportError("duckdb is required to create ChEMBL Parquet datasets.")

    parquet_path = parquet_dir / f'{dataset}.parquet'
    if parquet_path.exists():
        return parquet_path

    if dataset == 'activities':
        return _ensure_activities_parquet_dataset(duckdb_path, parquet_path)

    tmp_path = parquet_path.with_suffix('.parquet.tmp')
    if tmp_path.exists():
        tmp_path.unlink()

    print(f"Creating ChEMBL Parquet dataset: {parquet_path}")
    parquet_dir.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(duckdb_path), read_only=True)
    try:
        escaped_path = str(tmp_path).replace("'", "''")
        con.execute(
            f"COPY ({PARQUET_QUERIES[dataset]}) TO '{escaped_path}' "
            "(FORMAT PARQUET, COMPRESSION ZSTD)"
        )
    finally:
        con.close()

    tmp_path.replace(parquet_path)
    return parquet_path


def _ensure_activities_parquet_dataset(duckdb_path: Path, parquet_path: Path) -> Path:
    """Write ChEMBL activities as bounded Parquet parts.

    The activities query joins the largest ChEMBL table with several lookup
    tables. Running it as one COPY can make DuckDB build very large
    intermediates, so process stable activity_id ranges instead.
    """
    tmp_path = parquet_path.with_suffix('.parquet.tmp')
    if tmp_path.exists():
        if tmp_path.is_dir():
            shutil.rmtree(tmp_path)
        else:
            tmp_path.unlink()

    print(f"Creating ChEMBL Parquet dataset: {parquet_path}")
    parquet_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(duckdb_path), read_only=True)
    try:
        min_id, max_id = con.execute(
            'SELECT min(activity_id), max(activity_id) FROM activities'
        ).fetchone()
        if min_id is None or max_id is None:
            _copy_activities_chunk(con, tmp_path / 'part-00000.parquet', None, None)
        else:
            chunk_size = max(1, CHEMBL_ACTIVITIES_PARQUET_CHUNK_SIZE)
            start = int(min_id)
            max_id = int(max_id)
            part = 0
            while start <= max_id:
                stop = min(start + chunk_size, max_id + 1)
                print(
                    'Creating ChEMBL activities Parquet part '
                    f'{part:,}: activity_id [{start:,}, {stop:,})',
                    flush=True,
                )
                _copy_activities_chunk(
                    con,
                    tmp_path / f'part-{part:05d}.parquet',
                    start,
                    stop,
                )
                start = stop
                part += 1
    finally:
        con.close()

    tmp_path.replace(parquet_path)
    return parquet_path


def _copy_activities_chunk(
    con: duckdb.DuckDBPyConnection,
    output_path: Path,
    start: int | None,
    stop: int | None,
) -> None:
    escaped_path = str(output_path).replace("'", "''")
    query = PARQUET_QUERIES['activities']
    if start is not None and stop is not None:
        query = (
            f'{query}\n'
            f'        WHERE act.activity_id >= {start}\n'
            f'          AND act.activity_id < {stop}'
        )
    con.execute(
        f"COPY ({query}) TO '{escaped_path}' "
        "(FORMAT PARQUET, COMPRESSION ZSTD)"
    )


def _iter_chembl_prepared(
    dataset: str,
    opener: Any,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    sqlite_path = kwargs.get('sqlite_path')
    db_rel_path = kwargs.get('db_rel_path')
    duckdb_path = kwargs.get('duckdb_path')
    parquet_dir = kwargs.get('parquet_dir')
    use_parquet_cache = kwargs.get('use_parquet_cache', True)

    if use_parquet_cache and sqlite_path is not None:
        try:
            sqlite_path = Path(sqlite_path)
            duckdb_path, parquet_dir = _chembl_cache_paths(
                sqlite_path,
                duckdb_path=duckdb_path,
                parquet_dir=parquet_dir,
            )
            _ensure_sqlite(opener, sqlite_path, db_rel_path)
            _ensure_duckdb_cache(sqlite_path, duckdb_path)
            parquet_path = _ensure_parquet_dataset(dataset, duckdb_path, parquet_dir)
            yield from iter_parquet(path=parquet_path)
            return
        except Exception as error:
            print(
                f"ChEMBL DuckDB/Parquet cache unavailable for {dataset}; "
                f"falling back to SQLite. Reason: {error}"
            )

    yield from iter_sqlite(opener, query=SQLITE_QUERIES[dataset], **kwargs)


SQLITE_QUERIES: dict[str, str] = {
    'molecules': """
        SELECT
            md.*,
            cs.canonical_smiles,
            cs.standard_inchi,
            cs.standard_inchi_key,
            cp.full_mwt,
            cp.full_molformula,
            cp.alogp
        FROM molecule_dictionary md
        LEFT JOIN compound_structures cs ON md.molregno = cs.molregno
        LEFT JOIN compound_properties cp ON md.molregno = cp.molregno
    """,
    'assays': """
        SELECT
            a.*,
            t.chembl_id AS target_chembl_id,
            d.chembl_id AS document_chembl_id,
            GROUP_CONCAT(
                COALESCE(ap.standard_type, '') || ':' ||
                COALESCE(ap.standard_value, '') || ' ' ||
                COALESCE(ap.standard_units, ''),
                '; '
            ) AS parameters
        FROM assays a
        LEFT JOIN target_dictionary t ON a.tid = t.tid
        LEFT JOIN docs d ON a.doc_id = d.doc_id
        LEFT JOIN assay_parameters ap ON a.assay_id = ap.assay_id
        GROUP BY a.assay_id
    """,
    'documents': "SELECT * FROM docs",
    'mechanisms': """
        SELECT
            dm.*,
            md.chembl_id AS molecule_chembl_id,
            cs.canonical_smiles,
            cs.standard_inchi,
            cs.standard_inchi_key,
            td.chembl_id AS target_chembl_id,
            td.target_type,
            GROUP_CONCAT(DISTINCT CASE
                WHEN cseq.component_type = 'PROTEIN'
                 AND cseq.db_source IN ('SWISS-PROT', 'TREMBL')
                THEN cseq.accession
            END) AS target_component_uniprot_accessions,
            GROUP_CONCAT(DISTINCT CASE
                WHEN cseq.accession LIKE 'ENS%'
                THEN cseq.accession
            END) AS target_component_ensembl_accessions,
            d.chembl_id AS document_chembl_id,
            d.pubmed_id,
            d.doi,
            d.title AS document_title,
            GROUP_CONCAT(mr.ref_url, '; ') AS mechanism_refs
        FROM drug_mechanism dm
        LEFT JOIN molecule_dictionary md ON dm.molregno = md.molregno
        LEFT JOIN compound_structures cs ON dm.molregno = cs.molregno
        LEFT JOIN target_dictionary td ON dm.tid = td.tid
        LEFT JOIN target_components tc ON dm.tid = tc.tid
        LEFT JOIN component_sequences cseq ON tc.component_id = cseq.component_id
        LEFT JOIN compound_records cr ON dm.record_id = cr.record_id
        LEFT JOIN docs d ON cr.doc_id = d.doc_id
        LEFT JOIN mechanism_refs mr ON dm.mec_id = mr.mec_id
        GROUP BY dm.mec_id
    """,
    'targets': """
        SELECT
            td.tid,
            td.target_type,
            td.pref_name,
            td.tax_id,
            td.organism,
            td.chembl_id,
            GROUP_CONCAT(CASE
                WHEN cs.component_type = 'PROTEIN'
                 AND cs.db_source IN ('SWISS-PROT', 'TREMBL')
                THEN COALESCE(cs.accession, '')
                ELSE ''
            END) AS component_uniprot_accessions,
            GROUP_CONCAT(CASE
                WHEN cs.accession LIKE 'ENS%'
                THEN COALESCE(cs.accession, '')
                ELSE ''
            END) AS component_ensembl_accessions,
            GROUP_CONCAT(COALESCE(cs.component_type, '')) AS component_types,
            GROUP_CONCAT(COALESCE(tc.component_id, '')) AS component_ids,
            GROUP_CONCAT(COALESCE(cs.description, '')) AS component_descriptions
        FROM target_dictionary td
        LEFT JOIN target_components tc ON td.tid = tc.tid
        LEFT JOIN component_sequences cs ON tc.component_id = cs.component_id
        GROUP BY td.tid
    """,
    'activities': """
        WITH assay_params_agg AS MATERIALIZED (
            SELECT
                ap.assay_id,
                GROUP_CONCAT(
                    COALESCE(ap.standard_type, '') || ':' ||
                    COALESCE(ap.standard_value, '') || ' ' ||
                    COALESCE(ap.standard_units, ''),
                    '; '
                ) AS assay_parameters
            FROM assay_parameters ap
            GROUP BY ap.assay_id
        ),
        target_components_agg AS MATERIALIZED (
            SELECT
                tc.tid,
                GROUP_CONCAT(DISTINCT CASE
                    WHEN cseq.component_type = 'PROTEIN'
                     AND cseq.db_source IN ('SWISS-PROT', 'TREMBL')
                    THEN cseq.accession
                END) AS target_component_uniprot_accessions,
                GROUP_CONCAT(DISTINCT CASE
                    WHEN cseq.accession LIKE 'ENS%'
                    THEN cseq.accession
                END) AS target_component_ensembl_accessions
            FROM target_components tc
            LEFT JOIN component_sequences cseq ON tc.component_id = cseq.component_id
            GROUP BY tc.tid
        )
        SELECT
            act.activity_id,
            act.assay_id,
            act.molregno,
            act.doc_id,
            act.standard_type,
            act.standard_relation,
            act.standard_value,
            act.pchembl_value,
            act.data_validity_comment,
            act.action_type,
            md.chembl_id AS molecule_chembl_id,
            cs.canonical_smiles,
            cs.standard_inchi,
            cs.standard_inchi_key,
            a.chembl_id AS assay_chembl_id,
            a.assay_type,
            a.assay_tax_id,
            a.confidence_score,
            a.assay_category,
            a.assay_subcellular_fraction,
            a.assay_tissue,
            a.assay_cell_type,
            a.description AS assay_description,
            apa.assay_parameters,
            td.chembl_id AS target_chembl_id,
            td.target_type,
            td.tax_id AS target_tax_id,
            tca.target_component_uniprot_accessions,
            tca.target_component_ensembl_accessions,
            d.chembl_id AS document_chembl_id,
            d.pubmed_id,
            d.doi,
            le.bei,
            le.le,
            le.lle,
            le.sei,
            ax.description AS action_description,
            ax.parent_type AS action_parent_type
        FROM activities act
        LEFT JOIN molecule_dictionary md ON act.molregno = md.molregno
        LEFT JOIN compound_structures cs ON act.molregno = cs.molregno
        LEFT JOIN assays a ON act.assay_id = a.assay_id
        LEFT JOIN assay_params_agg apa ON a.assay_id = apa.assay_id
        LEFT JOIN target_dictionary td ON a.tid = td.tid
        LEFT JOIN target_components_agg tca ON a.tid = tca.tid
        LEFT JOIN docs d ON act.doc_id = d.doc_id
        LEFT JOIN ligand_eff le ON act.activity_id = le.activity_id
        LEFT JOIN action_type ax ON act.action_type = ax.action_type
    """,
}


def molecules_parser(opener: Any, **kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Parses ChEMBL molecules."""
    yield from _iter_chembl_prepared('molecules', opener, **kwargs)


def assays_parser(opener: Any, **kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Parses ChEMBL assays."""
    yield from _iter_chembl_prepared('assays', opener, **kwargs)


def mechanisms_parser(opener: Any, **kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Parses ChEMBL drug mechanisms."""
    yield from _iter_chembl_prepared('mechanisms', opener, **kwargs)


def documents_parser(opener: Any, **kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Parses ChEMBL documents."""
    yield from _iter_chembl_prepared('documents', opener, **kwargs)


def targets_parser(opener: Any, **kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Parses ChEMBL targets."""
    yield from _iter_chembl_prepared('targets', opener, **kwargs)


def activities_parser(opener: Any, **kwargs: Any) -> Generator[dict[str, Any], None, None]:
    """Parses ChEMBL activities."""
    yield from _iter_chembl_prepared('activities', opener, **kwargs)
