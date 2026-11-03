"""
Parse Vesiclepedia and emit the molecules found in extracellular vesicles.

Vesiclepedia is a manually curated compendium of the proteins, RNAs, lipids and
metabolites found in extracellular vesicles of any kind, from many species. It
uses the tables of ExoCarta, with a vesicle type for each study. Each protein or
RNA found in a study is a gene (identified by its Entrez Gene identifier and
symbol), related with ``located_in`` to each vesicle type of the study with
the detection method, the study's PubMed identifier, the sample and the
isolation methods of the study. The species of the molecule is the species the
table gives, for all of them.

The tables are the release 5.1 files (September 2023) of the Internet Archive's
copy of the site, captured in January 2026. The legacy reader used release 4.1,
and kept the proteins of one species.
"""

from __future__ import annotations

from functools import partial

from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._vesicles import map_detection, map_vesicle_type
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.exocarta import (
    VESICLEPEDIA,
    iter_detections,
    iter_vesicle_types,
)
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV

config = ResourceConfig(
    id=ResourceCv.VESICLEPEDIA,
    name='Vesiclepedia',
    url='http://microvesicles.org/',
    # The download page says "The data is free to all. Please cite the article
    # describing Vesiclepedia if you use the data." (http://microvesicles.org/download),
    # which names no license.
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='37953359',
    topics=(Topic.SUBCELLULAR_LOCALIZATION,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED, KnowledgeOrigin.EXPERIMENTAL),
    summary=(
        'Proteins and RNAs found in extracellular vesicles of many species '
        '(release 5.1), with vesicle type, publication, sample and methods.'
    ),
    description=(
        'Vesiclepedia is a manually curated compendium of the proteins, RNAs, '
        'lipids and metabolites found in extracellular vesicles (exosomes, '
        'microvesicles, apoptotic bodies, outer membrane vesicles and others), '
        'from published studies, with the species, sample, vesicle type, '
        'isolation method and detection method of each entry.'
    ),
)

# Species of the molecules in the protein/mRNA table of release 5.1, as the
# table writes them (without trailing spaces), with the NCBI taxon of the same
# name. A species that is not here raises.
TAXA = {
    'Acinetobacter baumannii': 'NCBITaxon:470',
    'Bos taurus': 'NCBITaxon:9913',
    'Burkholderia pseudomallei': 'NCBITaxon:28450',
    'Candida albicans': 'NCBITaxon:5476',
    'Canis lupus familiaris': 'NCBITaxon:9615',
    'Cavia porcellus': 'NCBITaxon:10141',
    'Citrus clementina': 'NCBITaxon:85681',
    'Citrus sinensis': 'NCBITaxon:2711',
    'Drosophila melanogaster': 'NCBITaxon:7227',
    'Echinostoma caproni': 'NCBITaxon:27848',
    'Edwardsiella tarda': 'NCBITaxon:636',
    'Equus caballus': 'NCBITaxon:9796',
    'Escherichia coli': 'NCBITaxon:562',
    'Fasciola hepatica': 'NCBITaxon:6192',
    # Francisella tularensis subsp. novicida, which NCBI also names Francisella novicida.
    'Francisella novicida': 'NCBITaxon:264',
    'Francisella philomiragia': 'NCBITaxon:28110',
    'Francisella tularensis': 'NCBITaxon:263',
    'Gallus gallus': 'NCBITaxon:9031',
    'Helicobacter pylori': 'NCBITaxon:210',
    'Homo sapiens': 'NCBITaxon:9606',
    'Klebsiella pneumoniae': 'NCBITaxon:573',
    'Listeria monocytogenes': 'NCBITaxon:1639',
    'Mus musculus': 'NCBITaxon:10090',
    'Myxococcus xanthus': 'NCBITaxon:34',
    'Neisseria meningitidis': 'NCBITaxon:487',
    'Ovis aries': 'NCBITaxon:9940',
    'Paracoccidioides brasiliensis': 'NCBITaxon:121759',
    'Pseudomonas aeruginosa': 'NCBITaxon:287',
    'Rattus norvegicus': 'NCBITaxon:10116',
    'Saccharomyces cerevisiae': 'NCBITaxon:4932',
    'Salmo salar': 'NCBITaxon:8030',
    'Staphylococcus aureus': 'NCBITaxon:1280',
    'Sus scrofa': 'NCBITaxon:9823',
    'Trichomonas vaginalis': 'NCBITaxon:5722',
    'Trypanosoma cruzi': 'NCBITaxon:5693',
}

VESICLE_TYPE_NAMESPACE = 'vesiclepedia_vesicle_type'

# Release 5.1 of the tables, from the Internet Archive (captured 2026-01-08 and 2026-01-07).
experiment_download = Download(
    url=(
        'https://web.archive.org/web/20260107233255id_/'
        'http://www.microvesicles.org/Archive/VESICLEPEDIA_EXPERIMENT_DETAILS_5.1.txt'
    ),
    filename='VESICLEPEDIA_EXPERIMENT_DETAILS_5.1.txt',
    subfolder='vesiclepedia',
    large=True,
)
detection_download = Download(
    url=(
        'https://web.archive.org/web/20260108064341id_/'
        'http://www.microvesicles.org/Archive/VESICLEPEDIA_PROTEIN_MRNA_DETAILS_5.1.txt'
    ),
    filename='VESICLEPEDIA_PROTEIN_MRNA_DETAILS_5.1.txt',
    subfolder='vesiclepedia',
    large=True,
)

resource = Resource(
    config,
    terms=Dataset(
        download=experiment_download,
        mapper=partial(
            map_vesicle_type,
            namespace=VESICLE_TYPE_NAMESPACE,
            description='Type of extracellular vesicle of Vesiclepedia.',
        ),
        raw_parser=partial(iter_vesicle_types, database=VESICLEPEDIA),
    ),
    detections=Dataset(
        download=detection_download,
        mapper=partial(
            map_detection, prefix=VESICLEPEDIA, namespace=VESICLE_TYPE_NAMESPACE, taxa=TAXA
        ),
        raw_parser=partial(
            iter_detections, experiment_download=experiment_download, database=VESICLEPEDIA
        ),
    ),
)
