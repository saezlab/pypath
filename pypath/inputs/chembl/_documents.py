from collections.abc import Generator

from ._records import ChemblDocument
from . import _raw

__all__ = [
    "document",
]

def document(max_pages: int | None = None) -> Generator[ChemblDocument]:
    """
    Retrieves the Chembl document information.

    This function is a wrapper around the `chembl_general` function.
    It retrieves the document information from ChEMBL and
    yields the data as named tuples of the type `ChemblDocument`.

    Args:
        max_pages (int): The maximum number of pages to retrieve.

    Yields:
        ChemblDocument: The named tuple of the retrieved data.
    """
    documents = _raw.json_pages(data_type="document", max_pages=max_pages)

    yield from (ChemblDocument
        (
            document_chembl_id = document["document_chembl_id"],
            pubmed_id = document["pubmed_id"],
            patent_id = document["patent_id"],
            doc_type = document["doc_type"],
            journal = document["journal"],
            year = document["year"],
            doi = document["doi"],
        )
        for document in documents
    )
