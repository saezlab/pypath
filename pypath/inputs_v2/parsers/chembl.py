"""
ChEMBL-specific data parsers.
"""

from __future__ import annotations

from typing import Any, Generator

from pypath.inputs_v2.parsers.base import iter_sqlite


def molecules_parser(
    opener: Any,
    **kwargs: Any,
) -> Generator[dict[str, Any], None, None]:
    """
    Parses ChEMBL molecules by joining dictionary, structures, and properties.

    Args:
        opener: The opener object for handling database connections.
        **kwargs: Additional arguments passed to iter_sqlite (e.g., sqlite_path, db_rel_path).

    Yields:
        Dictionary representing a joined molecule record.
    """
    query = """
    SELECT
        md.*,
        cs.canonical_smiles,
        cs.standard_inchi,
        cs.standard_inchi_key,
        cp.full_mwt,
        cp.full_molformula,
        cp.mw_monoisotopic,
        cp.alogp,
        cp.cx_logp,
        cp.cx_logd,
        cp.molecular_species
    FROM molecule_dictionary md
    LEFT JOIN compound_structures cs ON md.molregno = cs.molregno
    LEFT JOIN compound_properties cp ON md.molregno = cp.molregno
    """
    yield from iter_sqlite(opener, query=query, **kwargs)


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
