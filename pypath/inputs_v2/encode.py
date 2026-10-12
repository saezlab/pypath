"""
ENCODE-Nets: transcription factor networks derived from ENCODE ChIP-seq.

The regulatory networks of Gerstein et al., "Architecture of the human
regulatory network derived from ENCODE data" (Nature 2012), as the authors
publish them at encodenets.gersteinlab.org. They are not the ENCODE portal:
the edges are computed from the peaks of 119 transcription-related factors in
458 ChIP-seq experiments of the consortium, pooled over the cell types. Three
files are read; see ``parsers/encode.py`` for their layout.

- ``proximal_filtered`` (enets2): the 26,070 most confident edges between a
  factor and the gene near whose promoter it binds, chosen from 500,542 raw
  edges by a probabilistic model.
- ``distal`` (enets3): 19,258 edges between a factor and a gene through a
  distal regulatory element, found with a statistical model of chromatin
  marks and expression.
- ``tf_mirna`` (enets10): the edges from a factor to a microRNA, from the
  peaks of the factor near the microRNA gene. The 411 edges of the same file
  from a microRNA to a factor are predictions of TargetScan and are left out.

The edges come from the binding of the factor (ChIP-seq peaks) at or near the
gene, at a distal element linked to it, or near the microRNA. The networks do
not say whether a factor activates or represses a target, and binding alone
does not establish a change of expression, so every edge is an
``interacts_with`` between the factor and the gene (or microRNA), not an
``affects``. The network of the edge is on the relation (``encode:edge_type``).
"""

from __future__ import annotations

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.encode import (
    iter_distal,
    iter_proximal_filtered,
    iter_tf_mirna,
)
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

# The networks are of the human cell lines of ENCODE; the file has no species
# column. The microRNA names carry the prefix hsa-.
ENCODE_TAX_ID = '9606'
# The paper the networks are published with.
PUBLICATION = 'PMID:22955619'
# chromatin immunoprecipitation-seq evidence
CHIP_SEQ = 'ECO:0000229'
# Source-prefixed attribute: the kind of edge, as the file writes it
# (``proximal_filtered``, ``distal``, ``(TF-miRNA)``).
EDGE_TYPE = 'encode:edge_type'

config = ResourceConfig(
    id=ResourceCv.ENCODE,
    name='ENCODE',
    full='Encyclopedia of DNA Elements',
    url='http://encodenets.gersteinlab.org/',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='22955619',
    topics=(Topic.GENE_REGULATION,),
    content=(Content.CAUSAL_REGULATION,),
    knowledge_origin=(KnowledgeOrigin.EXPERIMENTAL,),
    summary=(
        'Transcription factor binding to genes and microRNAs, computed from '
        'ENCODE ChIP-seq peaks (ENCODE-Nets).'
    ),
    description=(
        'The ENCODE consortium mapped the binding of 119 transcription-'
        'related factors in 458 ChIP-seq experiments. Gerstein et al. '
        '(Nature 2012) turned the peaks into regulatory networks and '
        'publish them as ENCODE-Nets. This module imports the filtered '
        'promoter-proximal network, the distal network and the edges from '
        'factors to microRNAs. The edges are binding data processed by '
        'models, not literature curation, and have no sign; they state that '
        'a factor binds at or near a gene, not that it regulates it.'
    ),
)

f = FieldConfig()


def _entity(column: str, entity_type, namespace) -> EntityBuilder:
    return EntityBuilder(
        entity_type=entity_type,
        identifiers=IdentifiersBuilder(CV(term=namespace, value=f(column))),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=f'NCBITaxon:{ENCODE_TAX_ID}'),
        ),
    )


def _schema(target: EntityBuilder, evidence: str | None) -> RelationBuilder:
    annotations = [
        CV(term=slots.publications, value=PUBLICATION),
        CV(term=EDGE_TYPE, value=f('edge_type')),
    ]
    if evidence:
        annotations.append(CV(term=slots.has_evidence_of_type, value=evidence))
    return RelationBuilder(
        subject=_entity('source', model.Protein, Namespace.GENESYMBOL),
        predicate=slots.interacts_with,
        object=target,
        annotations=AnnotationsBuilder(*annotations),
    )


# The names are the gene names of the annotation the consortium used
# (GENCODE, 2011); some are older aliases or clone based names, and POL2 stands
# for RNA polymerase II.
gene_target = _entity('target', model.Gene, Namespace.GENESYMBOL)
mirna_target = _entity('target', model.MicroRNA, Namespace.NAME)

# The proximal edges and the edges to microRNAs come from the peaks of the
# factor near the gene. The distal edges come from a model that joins distal
# elements to genes through chromatin marks, so ChIP-seq is not claimed there.
proximal_schema = _schema(gene_target, CHIP_SEQ)
distal_schema = _schema(gene_target, None)
tf_mirna_schema = _schema(mirna_target, CHIP_SEQ)

_BASE_URL = 'http://encodenets.gersteinlab.org/'


def _download(filename: str) -> Download:
    return Download(
        url=f'{_BASE_URL}{filename}',
        filename=filename,
        subfolder='encode',
        large=True,
        encoding='ascii',
    )


resource = Resource(
    config,
    proximal_filtered=Dataset(
        download=_download('enets2.Proximal_filtered.txt'),
        mapper=proximal_schema,
        raw_parser=iter_proximal_filtered,
    ),
    distal=Dataset(
        download=_download('enets3.Distal.txt'),
        mapper=distal_schema,
        raw_parser=iter_distal,
    ),
    tf_mirna=Dataset(
        download=_download('enets10.TF-miRNA.txt'),
        mapper=tf_mirna_schema,
        raw_parser=iter_tf_mirna,
    ),
)
