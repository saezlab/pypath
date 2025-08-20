from collections import namedtuple
from collections.abc import Generator
import collections
import urllib
import re

from pypath.share import curl
from pypath.inputs import uniprot
import pypath.resources.urls as urls
import pypath.utils.taxonomy as taxonomy
import pypath.internals.intera as intera


Metabolite_cell = collections.namedtuple(
    'Metabolite_cell',
    ['hmdb_id', 'metabolite', 'interaction', 'cell_type', 'experimental_subject', 'disease', 'effect', 'pubmed_id'],
)

def mrclinksdb_raw(organism: int | str = 'human') -> Generator[list, None, None]:
    latin_name = taxonomy.ensure_latin_name(organism)
    url = urls.urls["mrclinksdb"]["url"] % urllib.parse.quote(latin_name)

    c = curl.Curl(url, large = True, silent = False)

    names = [
        re.sub(r'\W|^(?=\d)', '_', n.lower()).replace('class', 'class_')
        for n in next(c.result).split('\t')
    ]
    record = collections.namedtuple("MrclinksdbRaw", field_names=names)
    for line in c.result:
        yield record(*line.split("\t"))


def mrclinksdb_interaction(organism: int | str = 'human') -> Generator[list, None, None]:
    lines = list(mrclinksdb_raw(organism))
    uniprot_id = [row.receptor_uniprot_id for row in lines[1:]]
    all_uniprots = sorted({
        u
        for uu in uniprot_id
        for u in uu.split('_')
    })
    uniprot_locations = {}
    chunk_size = 98

    for i in range(0, len(all_uniprots), chunk_size):

        uniprot_locations.update(
            uniprot.uniprot_locations(
                accession=all_uniprots[i: i + chunk_size],
                organism=None,
                reviewed=None,
            )
        )

    record = collections.namedtuple(
        "MrclinksdbInteractions",
        lines[0]._fields + ("receptor_location",),
    )
    for line in lines:

        uniprot_id = line.receptor_uniprot_id

        if '_' in uniprot_id:

            uniprot_id = intera.Complex(
                components = sorted(uniprot_id.split('_')),
                ncbi_tax_id = 9606,
                sources = 'MRClinksDB',
            )

        location = uniprot_locations.get(uniprot_id)

        yield record(*line, location)


def metabolite_cell() -> Generator[Metabolite_cell, None, None]:
    url = "https://www.cellknowledge.com.cn/mrclinkdb/download/Metabolite-cell%20interaction.txt"
    c = curl.Curl(url)
    metabolite_cell = c.result
    lines = metabolite_cell.strip('\n').split('\n')
    for line in lines:
        hmdb_id, metabolite, interaction, cell_type, experimental_subject, disease, effect, _, pubmed_id = line.split(
            '\t')
        yield Metabolite_cell(hmdb_id, metabolite, interaction, cell_type, experimental_subject, disease, effect,
                              pubmed_id)






