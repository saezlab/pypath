"""
Parse the CollecTRI regulons and emit Relation records.

Source file (comma separated, with a header; one transcription factor -
target gene pair per line):

    source,target,weight,TF.category,resources,PMID,sign.decision
    SPI1,BGLAP,1,dbTF,ExTRI,10022617,default activation

It is the table the CollecTRI authors' sign assignment script writes
(``output/CollecTRI/CollecTRI.csv`` of https://github.com/saezlab/CollecTRI),
before any identifier translation. ``weight`` is the mode of regulation
(1 activation, -1 repression) and ``sign.decision`` says how it was chosen:

- ``PMID``: more of the cited publications report one mode than the other;
- ``regulon``: no majority, the mode of most of the TF's other targets is used;
- ``default activation``: no information, activation is filled in.

Only the first is a statement about the pair itself, so only those rows get
an effect direction.
"""

from __future__ import annotations

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import PMID_PATTERN
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_csv
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

# CollecTRI converts every participant to a human gene symbol (PMID 37843125);
# the file has no organism column.
COLLECTRI_TAX_ID = '9606'
# The two TF dimers CollecTRI treats as transcription factors themselves. They
# are names, not gene symbols, and the file does not list their subunits.
COMPLEXES = frozenset({'AP1', 'NFKB'})

_DIRECTION = model.DirectionQualifierEnum
_ASPECT = model.GeneOrGeneProductOrChemicalEntityAspectEnum
# weight -> effect on the target gene's expression.
_MODE_OF_REGULATION = {
    '1': _DIRECTION.increased,
    '-1': _DIRECTION.decreased,
}
# sign.decision -> whether the weight is supported by the row's publications.
_SIGN_DECISIONS = {
    'PMID': True,
    'regulon': False,
    'default activation': False,
}


def _is_complex(row: dict[str, object]) -> bool:
    return row.get('source') in COMPLEXES


def _direction(row: dict[str, object]):
    """The literature supported mode of regulation, None if it was filled in."""
    weight = str(row.get('weight') or '').strip()
    decision = str(row.get('sign.decision') or '').strip()
    if weight not in _MODE_OF_REGULATION:
        raise ValueError(f'Unknown CollecTRI weight: {weight!r}')
    if decision not in _SIGN_DECISIONS:
        raise ValueError(f'Unknown CollecTRI sign decision: {decision!r}')
    return _MODE_OF_REGULATION[weight] if _SIGN_DECISIONS[decision] else None


def _aspect(row: dict[str, object]):
    return _ASPECT.expression if _direction(row) is not None else None


config = ResourceConfig(
    id=ResourceCv.COLLECTRI,
    name='CollecTRI',
    full='Collection of Transcriptional Regulatory Interactions',
    url='https://github.com/saezlab/CollecTRI',
    # Same license as recorded in the legacy registry (resources.json).
    license=LicenseCV.COMPOSITE,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='37843125',
    topics=(Topic.GENE_REGULATION,),
    content=(Content.CAUSAL_REGULATION,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.TEXT_MINED,
    ),
    summary=(
        'Transcription factor–target gene regulation with literature '
        'references, source databases and the mode of regulation where '
        'publications support it.'
    ),
    description=(
        'CollecTRI is a collection of transcriptional regulatory '
        'interactions between human transcription factors and their target '
        'genes, compiled from twelve resources: the text mined ExTRI, HTRI, '
        'TRRUST, TFactS, GOA, IntAct, SIGNOR, CytReg, GEREDB, Pavlidis2021, '
        'the curated level A of DoRothEA and manual curation at NTNU. Almost '
        'every interaction has at least one PubMed reference (185 of 43,536 '
        'have only identifiers of other databases) and a transcription '
        'factor classified as DNA-binding, co-regulatory or general '
        'initiation factor. The regulons give each interaction a mode of '
        'regulation, taken from the references where they agree and '
        'otherwise from the rest of the regulon or set to activation.'
    ),
)

f = FieldConfig()

# A row without a transcription factor or a target yields no relation: the
# participant has no identifier, so the builder returns None.
transcription_factor = EntityBuilder(
    entity_type=lambda row: (
        model.MacromolecularComplex if _is_complex(row) else model.Protein
    ),
    identifiers=IdentifiersBuilder(
        CV(
            term=Namespace.GENESYMBOL,
            value=lambda row: None if _is_complex(row) else row.get('source'),
        ),
        CV(
            term=Namespace.NAME,
            value=lambda row: row.get('source') if _is_complex(row) else None,
        ),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=f'NCBITaxon:{COLLECTRI_TAX_ID}'),
    ),
)

# Protein coding, microRNA and other RNA genes alike are gene symbols here.
target_gene = EntityBuilder(
    entity_type=model.Gene,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.GENESYMBOL, value=f('target')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=f'NCBITaxon:{COLLECTRI_TAX_ID}'),
    ),
)

interactions_schema = RelationBuilder(
    subject=transcription_factor,
    predicate=slots.affects,
    object=target_gene,
    annotations=AnnotationsBuilder(
        # Activation or repression of the target gene's expression, only where
        # the cited publications decide it. An inferred or default sign is not
        # carried over, so those rows have neither qualifier.
        CV(term=slots.object_direction_qualifier, value=_direction),
        CV(term=slots.object_aspect_qualifier, value=_aspect),
        # The reference list also holds a few identifiers of other databases
        # (PAZAR, TRED); only PubMed IDs are publications.
        CV(
            term=slots.publications,
            value=f(
                'PMID',
                delimiter=';',
                extract=rf'^(?:PMID:)?({PMID_PATTERN})$',
                transform=lambda pmid: f'PMID:{pmid}',
            ),
        ),
        # The resources the pair comes from, as the file writes them.
        CV(
            term=slots.supporting_data_source,
            value=f('resources', delimiter=';'),
        ),
    ),
)

resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            url='https://rescued.omnipathdb.org/CollecTRI.csv',
            filename='CollecTRI.csv',
            subfolder='collectri',
            large=True,
        ),
        mapper=interactions_schema,
        raw_parser=iter_csv,
    ),
)
