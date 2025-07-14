import collections
import re
from collections.abc import Generator

from pypath.utils import taxonomy

from . import _raw
from ._records import BindingdbInteraction, BindingdbLigand, BindingdbTarget

RE_REGION_MUTATION = re.compile(
    r'(.+?)\s?((?:\[[-\d,A-Z\[\]/]+\])?$)'
)

def interactions(dataset: str = 'All',
                 max_lines: int | None = None) -> Generator[BindingdbInteraction]:

    """
    Yields interactions from BindingDB.

    Args:
        dataset : str, optional
            BindingDB dataset to read. The default is 'All'.
        max_lines : int | None, optional
            Maximum number of lines to read from the file. If None, the entire file
            is read.

    Yields:
        BindingdbInteraction
            Namedtuple with ligand and target information.

    Notes:
        The uniprot ID is mapped from the name extracted from the target field 
        according to the regular expression `RE_REGION_MUTATION`.

    """
    uniprot_mapping = _raw.mapping()

    for record in _raw.table(dataset = dataset, max_lines = max_lines):

        organism = record['Target Source Organism According to Curator or DataSource']

        try:
            name, regions_mutations = RE_REGION_MUTATION.match(record['Target Name']).groups()
        except AttributeError:
            print(record['Target Name'])

        yield BindingdbInteraction(
            ligand = BindingdbLigand(
                name = record['BindingDB Ligand Name'],
                smiles = record['Ligand SMILES'],
                inchi = record['Ligand InChI'],
                inchi_key = record['Ligand InChI Key'],
                pubchem = record['PubChem CID'],
            ),
            target = BindingdbTarget(
                name = name,
                organism = organism,
                ncbi_tax_id = taxonomy.ensure_ncbi_tax_id(organism),
                uniprot = uniprot_mapping.get(record['Target Name'], [None])[0], #TODO : create regex split to access all uniprots
                regions_mutations = [x.strip('[]').split(',') for x in regions_mutations.split('/')],
            ),
        )
