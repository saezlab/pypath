"""
Parse SLC tables data and emit Entity records.

This module converts SLC tables reference data into Entity records using the
declarative schema pattern.
"""

from bs4 import BeautifulSoup

from pypath.inputs_v2.base import (
    ResourceConfig,
    Download,
    Resource,
    Dataset,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
    MoleculeAnnotationsCv,
)

# =================================== SET-UP ===================================

BASE_URL = 'https://slc.bioparadigms.org/'
HUMAN_TAXON_ID = '9606'

config = ResourceConfig(
    id=ResourceCv.SLC_TABLES,
    name='SLC Tables',
    url='https://slc.bioparadigms.org/',
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='23506860',
    primary_category='transporters',
    description=(
        'The SLC (Solute Carrier) series includes the “classical” transporter '
        'families (ion-coupled transporters, exchangers, passive transporters, '
        'etc.). It contains 65 SLC gene families with 458 different human '
        'transporter genes, thus representing a major portion of the human '
        'transporter-related genes.'
    ),
)

# ================================== DOWNLOAD ==================================


def chunk_this(L, n):
    '''
    For a given list *L*, returns another list of *n*-sized chunks from
    it (in the same order).

    * Arguments:
        - *L* [list]: The list to be sliced into sublists of the
          definded size.
        - *n* [int]: The size of the chunks.

    * Returns:
        - [list]: List of *n*-sized chunks from *L*. **NOTE:** If the
          number of items in *L* is not divisible by *n*, the last
          element returned will have an inferior size.

    * Examples:
        >>> L = range(6)
        >>> chunk_this(L, 2)
        [[0, 1], [2, 3], [4, 5]]
        >>> chunk_this(L, 4)
        [[0, 1, 2, 3], [4, 5]]
    '''
    L = list(L)

    return [L[i:i + n] for i in range(0, len(L), n)]


def parser(opener, **_kwargs):

    # Obtaining plain text from Opener instance as single string
    txt = ''.join(opener.result)

    # Parsing HTML text
    soup = BeautifulSoup(txt, 'html.parser')

    # Extracting list of SLC families and corresponding tables
    fams = [t.text for t in soup.find_all('span', attrs={'class': 'slcname'})]
    raw_tables = soup.find_all('table')
    if not raw_tables:
        return

    headers = ['SLC_family'] + [
        t.text
        for t in raw_tables[0].find_all('td', attrs={'class': 'tbl_head'})
    ]

    for slc, table in zip(fams, raw_tables):
        cells = [
            t.text
            for t in table.find_all('td', attrs={'class': 'tbl_cell'})
        ]
        for chunk in chunk_this(cells, len(headers) - 1):
            if len(chunk) != len(headers) - 1:
                continue
            yield dict(zip(headers, [slc] + chunk))


def _slc_name(row: dict) -> str | None:
    value = str(row.get('SLC name') or '').strip()
    return value or None


def _primary_slc_gene_names(row: dict) -> list[str]:
    value = _slc_name(row)
    if value is None:
        return []
    if '_' in value:
        base = value.split('_', 1)[0].strip()
        return [base] if base else []
    return [value]


def _slc_gene_synonyms(row: dict) -> list[str]:
    value = _slc_name(row)
    synonyms = []
    if value and '_' in value:
        synonyms.append(value)
    return synonyms


download = Download(
    url=BASE_URL,
    filename='slctables.html',
    subfolder='slctables',
    large=True,
    ext='.html',
    default_mode='r',
)

# =================================== SCHEMA ===================================

f = FieldConfig(
    delimiter=', ',
    extract={},
    map={},
    transform={},
)

schema = EntityBuilder(
    entity_type=EntityTypeCv.PROTEIN,
    identifiers=IdentifiersBuilder(
        CV(
            term=IdentifierNamespaceCv.GENE_NAME_PRIMARY,
            value=_primary_slc_gene_names,
        ),
        CV(
            term=IdentifierNamespaceCv.GENE_NAME_SYNONYM,
            value=_slc_gene_synonyms,
        ),
        CV(
            term=IdentifierNamespaceCv.NAME,
            value=f('Protein name', delimiter=', ')
        ),
        CV(
            term=IdentifierNamespaceCv.SYNONYM,
            value=f('Aliases', delimiter=', ')
        ),
    ),
    annotations=AnnotationsBuilder(
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=HUMAN_TAXON_ID),
        CV(term=MoleculeAnnotationsCv.PROTEIN_FAMILY, value=f('SLC_family')),
        CV(term=MoleculeAnnotationsCv.DESCRIPTION, value=f('Transport type*')),
        CV(
            term=MoleculeAnnotationsCv.TRANSPORT_SUBSTRATE,
            value=f('Substrates', delimiter=', ')
        ),
        CV(
            term=MoleculeAnnotationsCv.TISSUE_LOCATION,
            value=f('Tissue and cellular expression', delimiter=', ')
        ),
    ),
)

# ================================= RESOURCE ===================================

resource = Resource(
    config=config,
    slc=Dataset(
        download=download,
        mapper=schema,
        raw_parser=parser,
    )
)

# ================================= REFERENCE ==================================

# X SLC_family                        SLC1
# X SLC name                          SLC1A1
# X Protein name                      EAAC1, EAAT3
# X Aliases                           System X-AG
# X Transport type*                   C / Na+, H+, K+
# X Substrates                        L-Glu, D/L-Asp
# X Tissue and cellular expression    brain (neurons), intestine, kidney, ...

# XXX: Thought about splitting the last columns but they have inconsistent
# separators/formatting across thables :(
