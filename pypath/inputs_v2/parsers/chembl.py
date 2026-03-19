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
