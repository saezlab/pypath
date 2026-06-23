"""
ChEMBL-specific data parsers.
"""

from __future__ import annotations

import os
import sqlite3
from typing import Any, Generator

from pypath.inputs_v2.parsers.base import (
    iter_sqlite,
    _extract_sqlite_from_opener,
)

# Optional structure/property columns projected from the 1:1 joins. ChEMBL
# renames/drops columns between releases (e.g. v36 has no ``mw_monoisotopic``/
# ``cx_logp``/``cx_logd``/``molecular_species``), so they are selected only when
# present in the deployed schema and NULL-aliased otherwise (version-resilient).
_CS_COLS = ('canonical_smiles', 'standard_inchi', 'standard_inchi_key')
_CP_COLS = (
    'full_mwt', 'full_molformula', 'mw_monoisotopic', 'alogp',
    'cx_logp', 'cx_logd', 'molecular_species',
)


def _table_columns(sqlite_path: Any, table: str) -> set[str]:
    """Return the column names of a table (empty set if it does not exist)."""
    conn = sqlite3.connect(sqlite_path)
    try:
        return {r[1] for r in conn.execute(f'PRAGMA table_info({table})')}
    finally:
        conn.close()


def molecules_parser(
    opener: Any,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Parses ChEMBL molecules by joining dictionary, structures, and properties.

    Resilient to ChEMBL schema differences across releases: optional structure
    and property columns are selected only when present (NULL otherwise), and
    ``synonyms`` aggregates the 1:many ``molecule_synonyms`` table via a
    correlated subquery (joined on the unit separator -- chemical names contain
    commas -- and split back to a list below) when that table exists.

    Args:
        opener: The opener object for handling database connections.
        **kwargs: Additional arguments passed to iter_sqlite (e.g., sqlite_path, db_rel_path).

    Yields:
        Dictionary representing a joined molecule record.
    """
    sqlite_path = kwargs.get('sqlite_path')
    if sqlite_path and not os.path.exists(sqlite_path):
        _extract_sqlite_from_opener(opener, sqlite_path, kwargs.get('db_rel_path'))

    cs_have = _table_columns(sqlite_path, 'compound_structures')
    cp_have = _table_columns(sqlite_path, 'compound_properties')
    has_synonyms = bool(_table_columns(sqlite_path, 'molecule_synonyms'))

    def _select(prefix: str, cols: tuple, have: set[str]) -> str:
        return ',\n        '.join(
            f'{prefix}.{c}' if c in have else f'NULL AS {c}' for c in cols
        )

    syn_select = (
        ',\n        (SELECT GROUP_CONCAT(ms.synonyms, char(31))'
        ' FROM molecule_synonyms ms WHERE ms.molregno = md.molregno) AS synonyms'
        if has_synonyms else ''
    )

    query = f"""
    SELECT
        md.*,
        {_select('cs', _CS_COLS, cs_have)},
        {_select('cp', _CP_COLS, cp_have)}{syn_select}
    FROM molecule_dictionary md
    LEFT JOIN compound_structures cs ON md.molregno = cs.molregno
    LEFT JOIN compound_properties cp ON md.molregno = cp.molregno
    """

    for row in iter_sqlite(opener, query=query, **kwargs):

        syn = row.get('synonyms')

        if syn:
            row['synonyms'] = [s for s in str(syn).split('\x1f') if s]

        yield row


def assays_parser(
    opener: Any,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Parses ChEMBL assays by joining with targets, documents, and parameters.
    Aggregates parameters into a single string using GROUP_CONCAT.

    Args:
        opener: The opener object for handling database connections.
        **kwargs: Additional arguments passed to iter_sqlite.

    Yields:
        Dictionary representing a joined assay record.
    """
    query = """
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
    """
    yield from iter_sqlite(opener, query=query, **kwargs)


def mechanisms_parser(
    opener: Any,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Parses ChEMBL drug mechanisms by joining with molecules, targets, and documents.

    Args:
        opener: The opener object for handling database connections.
        **kwargs: Additional arguments passed to iter_sqlite.

    Yields:
        Dictionary representing a joined mechanism record.
    """
    query = """
    SELECT
        dm.*,
        md.chembl_id AS molecule_chembl_id,
        td.chembl_id AS target_chembl_id,
        d.chembl_id AS document_chembl_id,
        d.pubmed_id,
        d.doi,
        d.title AS document_title,
        GROUP_CONCAT(mr.ref_url, '; ') AS mechanism_refs
    FROM drug_mechanism dm
    LEFT JOIN molecule_dictionary md ON dm.molregno = md.molregno
    LEFT JOIN target_dictionary td ON dm.tid = td.tid
    LEFT JOIN docs d ON dm.doc_id = d.doc_id
    LEFT JOIN mechanism_refs mr ON dm.mec_id = mr.mec_id
    GROUP BY dm.mec_id
    """
    yield from iter_sqlite(opener, query=query, **kwargs)


def documents_parser(
    opener: Any,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Parses ChEMBL documents from the `docs` table.

    Args:
        opener: The opener object for handling database connections.
        **kwargs: Additional arguments passed to iter_sqlite.

    Yields:
        Dictionary representing a document record.
    """
    query = "SELECT * FROM docs"
    yield from iter_sqlite(opener, query=query, **kwargs)


def targets_parser(
    opener: Any,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Parses ChEMBL targets by joining with components and sequences.

    Args:
        opener: The opener object for handling database connections.
        **kwargs: Additional arguments passed to iter_sqlite.

    Yields:
        Dictionary representing a joined target record.
    """
    query = """
    SELECT
        td.tid,
        td.target_type,
        td.pref_name,
        td.tax_id,
        td.organism,
        td.chembl_id,
        GROUP_CONCAT(COALESCE(cs.accession, '')) AS component_accessions,
        GROUP_CONCAT(COALESCE(cs.component_type, '')) AS component_types,
        GROUP_CONCAT(COALESCE(tc.component_id, '')) AS component_ids,
        GROUP_CONCAT(COALESCE(cs.description, '')) AS component_descriptions
    FROM target_dictionary td
    LEFT JOIN target_components tc ON td.tid = tc.tid
    LEFT JOIN component_sequences cs ON tc.component_id = cs.component_id
    GROUP BY td.tid
    """
    yield from iter_sqlite(opener, query=query, **kwargs)


def activities_parser(
    opener: Any,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Parses ChEMBL activities by joining with molecules, assays, targets, documents,
    ligand efficiency, and action types.

    Args:
        opener: The opener object for handling database connections.
        **kwargs: Additional arguments passed to iter_sqlite.

    Yields:
        Dictionary representing a joined activity record.
    """
    query = """
    SELECT
        act.*,
        md.chembl_id AS molecule_chembl_id,
        a.chembl_id AS assay_chembl_id,
        td.chembl_id AS target_chembl_id,
        td.tax_id AS target_tax_id,
        d.chembl_id AS document_chembl_id,
        d.pubmed_id,
        le.bei,
        le.le,
        le.lle,
        le.sei,
        at.description AS action_description,
        at.parent_type AS action_parent_type
    FROM activities act
    LEFT JOIN molecule_dictionary md ON act.molregno = md.molregno
    LEFT JOIN assays a ON act.assay_id = a.assay_id
    LEFT JOIN target_dictionary td ON a.tid = td.tid
    LEFT JOIN docs d ON act.doc_id = d.doc_id
    LEFT JOIN ligand_eff le ON act.activity_id = le.activity_id
    LEFT JOIN action_type at ON act.action_type = at.action_type
    """
    yield from iter_sqlite(opener, query=query, **kwargs)
