from collections import namedtuple
from collections.abc import Generator
import collections
import urllib
from pypath.share import curl
from pypath.inputs import uniprot
import pypath.resources.urls as urls
import pypath.utils.taxonomy as taxonomy
import re

Metabolite_cell = collections.namedtuple(
    'Metabolite_cell',
    ['hmdb_id', 'metabolite', 'interaction', 'cell_type', 'experimental_subject', 'disease', 'effect', 'pubmed_id'],
)

def mrclinksdb_raw(organism: int | str = 'human') -> Generator[list, None, None]:
    latin_name = taxonomy.ensure_latin_name(organism)
    url = urls.urls["mrclinksdb"]["url"] % urllib.parse.quote(latin_name)

    c = curl.Curl(url, large = True, silent = False)
    for line in c.result:
        yield line.split("\t")

def mrclinksdb_interaction(organism: int | str = 'human') -> Generator[list, None, None]:
    lines = list(mrclinksdb_raw(organism))
    header = lines[0]
    clean_header = [
        "class_" if  h.lower() == "class"
        else re.sub(r'\W|^(?=\d)', '_', h).lower()
        for h in header]
    index = clean_header.index("receptor_uniprot_id")
    uniprot_id = [row[index] for row in lines[1:]]
    all_uniprots = [
        u.split('_', maxsplit=1)[0]
        for u in uniprot_id
        if uniprot.valid_uniprot(u)
    ]
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
    clean_header += ["receptor_location"]
    row = collections.namedtuple("row", field_names=clean_header)
    for line in lines[1:]:
        uniprot_id = line[index].split('.', 1)[0]
        if uniprot_id in uniprot_locations:
            location = uniprot_locations[uniprot_id]
        else:
            location = None

        yield row(*line,location)

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






