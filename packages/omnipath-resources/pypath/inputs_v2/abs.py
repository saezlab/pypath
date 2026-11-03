"""
Parse the ABS binding sites and emit Relation records of transcription factors and promoters.

ABS (Annotated regulatory Binding Sites from orthologous promoters) collects
the transcription factor binding sites that were verified experimentally in
the promoters of 100 sets of orthologous genes of human, mouse, rat and
chicken, from the papers that describe them. Its GFF file (see
``parsers/abs.py``) has one line for each of the 650 sites: the promoter
sequence in which it lies, the ABS entry of the promoter, the transcription
factor, the positions of the site in the promoter and the species.

The file does not name the gene of a promoter: the promoter is a sequence, a
GenBank record (or, for 58 of the 211, a sequence the authors extracted from a
genome, with a name only ABS knows). The legacy reader treated the sequence
accession as the identifier of a gene. This module does not: the object of
the relation is the promoter, identified by its GenBank accession and by its
name in ABS, with the species the file gives it. The factor is known by the
name ABS writes (``SP1``, ``P53``, ``CEBP``, ``CAAT``), not by a gene
symbol, and the file gives no species for it.

What the source observed is a factor bound to a site in a promoter, with no
mode of regulation, so the relation is ``interacts_with`` without direction or
aspect, as for HTRIdb and ORegAnno (``htri.py``, ``oreganno.py``).
"""

from __future__ import annotations

import re
from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.abs import iter_sites
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.ABS,
    name='ABS',
    full='Annotated regulatory Binding Sites from orthologous promoters',
    url='http://genome.crg.es/datasets/abs2005/index.html',
    # "ABS is under GNU General Public License" (copyright 2005), without a
    # version; version 2 was the current one then, and the legacy metadata says
    # GPLv2.
    license=LicenseCV.GPL_2_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='16381947',
    topics=(Topic.GENE_REGULATION,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Experimentally verified transcription factor binding sites in the '
        'promoters of orthologous genes of human, mouse, rat and chicken.'
    ),
    description=(
        'ABS is a database of experimentally verified transcription factor '
        'binding sites in promoters, collected from the literature by '
        'Blanco et al. (Genome Informatics Research Laboratory, IMIM-UPF-CRG). '
        'It has 650 sites in 211 promoters of 100 sets of orthologous genes '
        'of human, mouse, rat and chicken. The promoter sequences, 500 '
        'nucleotides each, were retrieved with their GenBank accession or '
        'extracted from genomes, and the sites were mapped to them and '
        'checked against the orthologous promoters. The release is of 2005.'
    ),
)

# The file is read from the Internet Archive copy of 2022-08-08. The original
# location, http://genome.crg.es/datasets/abs2005/data/abs.gff, still serves
# it, but the server sends an incomplete certificate chain (the intermediate
# certificate "GEANT TLS RSA 1" is missing), so the build's client cannot
# verify the site. The archive copy is byte for byte the file of the original
# location (compared on 2026-10-12).
download = Download(
    url=(
        'http://web.archive.org/web/20220808195334id_/'
        'https://genome.crg.es/datasets/abs2005/data/abs.gff'
    ),
    filename='abs.gff',
    subfolder='abs',
)

# ``# HS:`` comment of a line -> NCBI Taxonomy ID of the species of the
# promoter. The codes are those of ABS: human, mouse, rat, chicken.
SPECIES = {
    'HS': '9606',
    'MM': '10090',
    'RN': '10116',
    'GG': '9031',
}
# A GenBank accession; the other promoters are called ``extractedNN``.
_GENBANK = re.compile(r'[A-Z]{1,2}\d{5,6}')
_ABS_PROMOTER = 'abs_promoter'


def _taxon(row: dict[str, Any]) -> str:
    species = row.get('species')
    if species not in SPECIES:
        raise ValueError(f'Unknown ABS species: {species!r}')
    return f'NCBITaxon:{SPECIES[species]}'


def _genbank(row: dict[str, Any]) -> str | None:
    return row['sequence'] if _GENBANK.fullmatch(row.get('sequence') or '') else None


f = FieldConfig()

# A row without a factor name or a promoter yields no relation: the
# participant has no identifier, so the builder returns None.
transcription_factor = EntityBuilder(
    entity_type=model.Protein,
    identifiers=IdentifiersBuilder(CV(term=Namespace.NAME, value=f('factor'))),
)
promoter = EntityBuilder(
    entity_type=model.RegulatoryRegion,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.GENBANK, value=_genbank),
        CV(term=_ABS_PROMOTER, value=f('sequence')),
    ),
    annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=_taxon)),
)

sites_schema = RelationBuilder(
    subject=transcription_factor,
    predicate=slots.interacts_with,
    object=promoter,
)

resource = Resource(
    config,
    binding_sites=Dataset(
        download=download,
        mapper=sites_schema,
        raw_parser=iter_sites,
    ),
)
