"""
Parse ExoCarta and emit the molecules found in exosomes.

ExoCarta is a manually curated database of the proteins, RNAs and lipids that
were found in exosomes, one entry for each molecule found in a study, with the
study's publication, sample and methods. This module imports the proteins and
RNAs: each is a gene (identified by its Entrez Gene identifier and symbol, which
the table gives) related with ``located_in`` to the vesicle type, here always
exosomes, with the detection method, the study's PubMed identifier, the sample
and the isolation methods of the study. Species is the species of the molecule
as the table writes it. ExoCarta's tables serve Vesiclepedia as well:
the mapping they share is in ``_vesicles.py``, the reader in ``parsers/exocarta.py``.

The tables are read from the Internet Archive's copy of release 5 (September
2015): the current release 6 (March 2025) is behind a Cloudflare challenge and
the Internet Archive holds no capture of its tables (it holds only the top-100
file of release 6; the rescue server has none either). The legacy reader used
release 6, and kept proteins of one species.
"""

from __future__ import annotations

from functools import partial

from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._vesicles import map_detection, map_vesicle_type
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.exocarta import (
    EXOCARTA,
    iter_detections,
    iter_vesicle_types,
)
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV

config = ResourceConfig(
    id=ResourceCv.EXOCARTA,
    name='ExoCarta',
    url='http://exocarta.org/',
    # The site states no license ("ExoCarta Download - 6 (Release date: 20 March 2025)"
    # offers the files, http://exocarta.org/download).
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='26434508',
    topics=(Topic.SUBCELLULAR_LOCALIZATION,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED, KnowledgeOrigin.EXPERIMENTAL),
    summary=(
        'Proteins and RNAs found in exosomes (release 5 of 2015; release 6 of '
        '2025 is not reachable), with species, publication, sample, isolation '
        'and detection methods.'
    ),
    description=(
        'ExoCarta is a manually curated database of exosomal proteins, RNAs and '
        'lipids. Each entry is a molecule identified in exosomes in a published '
        'study, with the species, the sample (cell line, tissue or body fluid), '
        'the isolation method of the exosomes and the method that detected the '
        'molecule. The data are those of release 5 (September 2015).'
    ),
)

# Species of the molecules in the table of ExoCarta (release 5), as written
# there, with the NCBI taxon of the same name. "Drosophila" is the genus: the
# table does not name the species for these 176 proteins, while its experiment
# table names Drosophila melanogaster.
TAXA = {
    'Bos taurus': 'NCBITaxon:9913',
    'Cavia porcellus': 'NCBITaxon:10141',
    'Drosophila': 'NCBITaxon:7215',
    'Homo sapiens': 'NCBITaxon:9606',
    'Human herpesvirus 4': 'NCBITaxon:10376',
    'Mus musculus': 'NCBITaxon:10090',
    'Rattus norvegicus': 'NCBITaxon:10116',
}

VESICLE_TYPE_NAMESPACE = 'exocarta_vesicle_type'


# Release 5 of the tables, from the Internet Archive (captured 2015-09-18).
_ARCHIVE = 'https://web.archive.org/web/20150918122'
experiment_download = Download(
    url=f'{_ARCHIVE}503id_/http://exocarta.org:80/Archive/EXOCARTA_EXPERIMENT_DETAILS_5.txt',
    filename='EXOCARTA_EXPERIMENT_DETAILS_5.txt',
    subfolder='exocarta',
    large=True,
)
detection_download = Download(
    url=f'{_ARCHIVE}556id_/http://exocarta.org:80/Archive/EXOCARTA_PROTEIN_MRNA_DETAILS_5.txt',
    filename='EXOCARTA_PROTEIN_MRNA_DETAILS_5.txt',
    subfolder='exocarta',
    large=True,
)

resource = Resource(
    config,
    terms=Dataset(
        download=experiment_download,
        mapper=partial(
            map_vesicle_type,
            namespace=VESICLE_TYPE_NAMESPACE,
            description='Extracellular vesicles of the type ExoCarta collects.',
        ),
        raw_parser=partial(iter_vesicle_types, database=EXOCARTA),
    ),
    detections=Dataset(
        download=detection_download,
        mapper=partial(
            map_detection, prefix=EXOCARTA, namespace=VESICLE_TYPE_NAMESPACE, taxa=TAXA
        ),
        raw_parser=partial(
            iter_detections, experiment_download=experiment_download, database=EXOCARTA
        ),
    ),
)
