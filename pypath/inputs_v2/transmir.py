"""
TransmiR: transcription factor - microRNA regulations curated from literature.

Source file: ``db_literaturev3.txt`` of TransmiR v3.0 (August 2024), the
literature curated regulations of all 29 organisms of the database; see
``parsers/transmir.py`` for its columns. The target of a regulation is the
microRNA gene, so the relation is a transcription factor ``affects`` a
microRNA. The action type gives the direction of the effect on the
expression of the microRNA where it is Activation or Repression.

The ChIP-seq derived and the predicted regulations that TransmiR v3.0 also
offers (six million lines) are not part of this module.
"""

from __future__ import annotations

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import PMID_EXTRACT
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.transmir import iter_literature
from pypath.internals.cv_terms import (
    LicenseCV,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

_DIRECTION = model.DirectionQualifierEnum
_ASPECT = model.GeneOrGeneProductOrChemicalEntityAspectEnum

# The action type -> effect on the expression of the microRNA. Activation and
# Repression, with or without a remark in brackets (the feedback remark says
# that the factor is also a target of the microRNA), say the direction. The
# other values say that the factor regulates the microRNA, not how.
_ACTIONS = {
    'Activation': _DIRECTION.increased,
    'Activation(feedback)': _DIRECTION.increased,
    'Activation(a negative feedback loop)': _DIRECTION.increased,
    'Activation(negative regulatory loop)': _DIRECTION.increased,
    'Repression': _DIRECTION.decreased,
    'Repression(feedback)': _DIRECTION.decreased,
    'Regulation': None,
    'Regulation(feedback)': None,
    'Regulation(Double-Negative Feedback Loop)': None,
    'auto-regulatory feedback circuit': None,
    'autoregulatory negative feedback loop(feedback)': None,
    'A negative regulatory loop(feedback)': None,
    'regulatory loop': None,
    'regulatory network': None,
}
# The species of the regulation -> NCBI Taxonomy ID. The viruses are None: a
# regulation of a virus involves a viral and a host participant, and the file
# does not say which is which.
_TAXA = {
    'H.sapiens': '9606',
    'M.musculus': '10090',
    'R.norvegicus': '10116',
    'A.thaliana': '3702',
    'C.elegans': '6239',
    'S.scrofa': '9823',
    'D.melanogaster': '7227',
    'G.gallus': '9031',
    'D.rerio': '7955',
    'G.max': '3847',
    'S.lycopersicum': '4081',
    'O.sativa': '4530',
    'M.notabilis': '981085',
    'C.familiaris': '9615',
    'C.sinensis': '2711',
    'L.chinensis': '151069',
    'H.armigera': '29058',
    'C.idella': '7959',
    'F.vesca': '57918',
    'B.grunniens': '30521',
    'R.chinensis': '74649',
    'M.domestica': '3750',
    'V.vinifera': '29760',
    'C.hircus': '9925',
    'Human herpesvirus 2': None,
    'Epstein-Barr virus (EBV)': None,
    'Human papillomaviruses 16': None,
    'Human herpesvirus 8 type M': None,
    'Human gammaherpesvirus 4': None,
}


def _direction(row: dict[str, object]):
    """The effect of the factor on the expression of the microRNA, if stated."""
    action = row.get('action')
    if action not in _ACTIONS:
        raise ValueError(f'Unknown TransmiR action type: {action!r}')
    return _ACTIONS[action]


def _aspect(row: dict[str, object]):
    return _ASPECT.expression if _direction(row) is not None else None


def _taxon(row: dict[str, object]):
    species = row.get('species')
    if species not in _TAXA:
        raise ValueError(f'Unknown TransmiR species: {species!r}')
    return f'NCBITaxon:{_TAXA[species]}' if _TAXA[species] else None


config = ResourceConfig(
    id=ResourceCv.TRANSMIR,
    name='TransmiR',
    full='Transcription Factor-microRNA Regulation Database',
    url='http://www.cuilab.cn/transmir',
    # The TransmiR site states that the data are free for academic usage
    # (checked when the module was written, October 2026); the legacy registry
    # recorded CC BY-NC 4.0.
    license=LicenseCV.ACADEMIC_FREE,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='39530226',
    topics=(Topic.GENE_REGULATION,),
    content=(Content.CAUSAL_REGULATION,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary=(
        'Literature curated transcription factor - microRNA regulations '
        'with the action type, the publication and the species.'
    ),
    description=(
        'TransmiR is a database of transcription factor - microRNA '
        'regulations of the Cui Lab. Version 3.0 (2024) holds 5,110 '
        'regulations curated from 2,291 publications, in 29 organisms, and '
        'six million regulations derived from ChIP-seq data. This module '
        'imports the literature curated regulations: each gives the '
        'transcription factor, the microRNA, the action type (activation, '
        'repression or regulation), the PubMed ID and the species.'
    ),
)

f = FieldConfig()

# The names are as the curators wrote them: mostly HGNC or species gene
# symbols, but also protein names (HIF-1α) and aliases.
transcription_factor = EntityBuilder(
    entity_type=model.Protein,
    identifiers=IdentifiersBuilder(CV(term=Namespace.NAME, value=f('tf'))),
    annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=_taxon)),
)

# The names are miRBase precursor names (hsa-mir-21) except for the plants
# (ath-MIR156a) and a few viral and other names; none is an accession.
mirna = EntityBuilder(
    entity_type=model.MicroRNA,
    identifiers=IdentifiersBuilder(CV(term=Namespace.NAME, value=f('mirna'))),
    annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=_taxon)),
)

interactions_schema = RelationBuilder(
    subject=transcription_factor,
    predicate=slots.affects,
    object=mirna,
    annotations=AnnotationsBuilder(
        CV(term=slots.object_direction_qualifier, value=_direction),
        CV(term=slots.object_aspect_qualifier, value=_aspect),
        CV(
            term=slots.publications,
            value=f('pmid', extract=PMID_EXTRACT, transform=lambda pmid: f'PMID:{pmid}'),
        ),
    ),
)

resource = Resource(
    config,
    literature=Dataset(
        download=Download(
            url='http://www.cuilab.cn/files/transmir3/download/db_literaturev3.txt',
            filename='db_literaturev3.txt',
            subfolder='transmir/v3.0',
            large=True,
        ),
        mapper=interactions_schema,
        raw_parser=iter_literature,
    ),
)
