from collections.abc import Generator
import collections
from pypath.share import curl
from pypath.resources import urls
from pypath.inputs import uniprot
import pandas as pd

Homo_receptor = collections.namedtuple(
    'Homo_receptor',
    ['hmdb_id', 'metabolite', 'receptor_gene_id', 'receptor_uniprot_id', 'receptor_symbol', 'pubmed_id', 'source_db'],
)
Mouse_receptor = collections.namedtuple(
    'Mouse_receptor',
    ['hmdb_id', 'metabolite', 'receptor_gene_id', 'receptor_uniprot_id', 'receptor_symbol', 'pubmed_id', 'source_db'],
)
Metabolite_cell = collections.namedtuple(
    'Metabolite_cell',
    ['hmdb_id', 'metabolite', 'interaction', 'cell_type', 'experimental_subject', 'disease', 'effect', 'pubmed_id'],
)


def homo_receptor() -> Generator[Homo_receptor, None, None]:
    # homo
    url = "https://www.cellknowledge.com.cn/mrclinkdb/download/Homo%20sapiens%20metabolite%20L-R%20interaction.txt"
    c = curl.Curl(url)
    human_receptor = c.result
    lines = human_receptor.strip('\n').split('\n')
    for line in lines:
        _, hmdb_id, metabolite, _, _, _, _, _, _, receptor_gene_id, receptor_uniprot_id, receptor_symbol, _, pubmed_id, source_db = line.split(
            '\t')
        yield Homo_receptor(hmdb_id, metabolite, receptor_gene_id, receptor_uniprot_id, receptor_symbol, pubmed_id,
                            source_db)


def mouse_receptor() -> Generator[Mouse_receptor, None, None]:
    # mouse
    url = "https://www.cellknowledge.com.cn/mrclinkdb/download/Mus%20musculus%20metabolite%20L-R%20interaction.txt"
    c = curl.Curl(url)
    mouse_receptor = c.result
    lines = mouse_receptor.strip('\n').split('\n')
    for line in lines:
        _, hmdb_id, metabolite, _, _, _, _, _, _, receptor_uniprot_id, receptor_gene_id, receptor_symbol, _, pubmed_id, source_db = line.split(
            '\t')
        yield Mouse_receptor(hmdb_id, metabolite, receptor_gene_id, receptor_uniprot_id, receptor_symbol, pubmed_id,
                             source_db)


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






