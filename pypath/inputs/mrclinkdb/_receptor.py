from collections.abc import Generator
import collections
from pypath.share import curl
from pypath.inputs import uniprot
import pypath.resources.urls as urls

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


def homo_receptor(organism:int | str) -> Generator[Homo_receptor, None, None]:
    # homo
    url = urls.urls["mrclinksdb"]["url"]

    human_receptor = c.result
    lines = human_receptor.strip('\n').split('\n')

    # ID conversion
    # bug: Some uniprot has _xxx_!!
    uniprot = []
    for line in lines:
        _, _, _, _, _, _, _, _, _, _, receptor_uniprot_id, _, _, _, _ = line.split('\t')
        if receptor_uniprot_id:
            uniprot.append(receptor_uniprot_id)

    uniprot = uniprot[1:]
    all_uniprots = set(uniprot)

    all_uniprots = [
        u.split('_', maxsplit=1)[0]
        for u in all_uniprots
        if uniprot.valid_uniprot(u)
    ]

    uniprot_locations = uniprot.uniprot_locations(
        accession=all_uniprots,
        organism=None,
    )

    for line in lines:
        _, hmdb_id, metabolite, _, _, _, _, _, _, receptor_gene_id, receptor_uniprot_id, receptor_symbol, _, pubmed_id, source_db = line.split(
            '\t')
        yield Homo_receptor(hmdb_id, metabolite, receptor_gene_id, receptor_uniprot_id, receptor_symbol, pubmed_id,
                            source_db,
                            uniprot_locations.get(receptor_uniprot_id))


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






