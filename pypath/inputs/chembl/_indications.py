
from collections.abc import Generator

from ._records import ChemblIndication
from . import _raw

def get_indicatons(max_pages: int | None = None) -> Generator[ChemblIndication]:
    """
    Retrieves drug indications from ChEMBL.

    This function is a wrapper around the `chembl_general` function.
    It retrieves the indication information from ChEMBL and
    yields the data as named tuples of the type `ChemblIndication`.

    Args:
        max_pages (int): The maximum number of pages to retrieve.

    Yields:
        ChemblIndication: The named tuple of the retrieved data.
    """
    indications = _raw.json_pages(
        data_type="drug_indication",
        max_pages=max_pages,
    )

    yield from (ChemblDocument
        (
            indication_chembl_id = indication["document_chembl_id"],
            pubmed_id = indication["pubmed_id"],
            patent_id = indication["patent_id"],
            doc_type = indication["doc_type"],
            journal = indication["journal"],
            year = indication["year"],
            doi = indication["doi"],
        )
        for indication in indications
    )

# TODO: check if drug indications are required
def chembl_drug_indications(
    max_phase_threshold: int = 0,
    ) -> list[tuple]:
    """
    Retrieves drug indications data from ChEMBL.

    Args
        max_phase_threshold:
            The threshold for maximum phase of the drug
            for which the indication is valid.
    Returns
        List of drug indications as namedtuples.
    """

    fields_indication = (
        'efo_id',
        'efo_term',
        'max_phase',
        'mesh_heading',
        'mesh_id',
        'molecule_chembl',
    )

    ChemblIndication = collections.namedtuple(
        'ChemblIndication',
        fields_indication,
        defaults = (None,) * len(fields_indication),
    )

    indication_lst = []
    page_dct = {}

    while True:

        if not page_dct:

            url = (
                f"{urls.urls['chembl']['url']}"
                f"{urls.urls['chembl']['drug_indication']}"
            )

        elif page_dct['page_meta']['next']:
            url = (
                f"{urls.urls['chembl']['url']}"
                f"{page_dct['page_meta']['next']}"
            )

        else:
            break

        c = curl.Curl(url, large=True, silent=False)
        fileobj = open(c.fileobj.name, encoding='utf-8')
        page_dct = json.loads(fileobj.read())

        indication_lst.extend(
            ChemblIndication(
                efo_id = ind['efo_id'],
                efo_term = ind['efo_term'],
                max_phase = float(ind['max_phase_for_ind']),
                mesh_heading = ind['mesh_heading'],
                mesh_id = ind['mesh_id'],
                molecule_chembl = ind['molecule_chembl_id'],
            )
            for ind in page_dct['drug_indications']
            if float(ind['max_phase_for_ind']) > max_phase_threshold and max_phase_threshold != 0 \
                or max_phase_threshold == 0
        )

    return indication_lst
