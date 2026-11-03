"""
Phobius predictions of transmembrane helices and signal peptides for the human proteome.

Phobius (Käll, Krogh and Sonnhammer 2004) predicts the transmembrane topology
and the signal peptide of a protein sequence with one hidden Markov model.
OmniPath ran it on the human Swiss-Prot proteome and published the output
(``phobius_human_proteome.tsv``, last modified 2020-03-30; the Phobius version
and the run date are not recorded in the file). These are predictions made by
OmniPath, not curated or measured annotations.

The module imports:

* ``terms``: two terms of the module, "predicted transmembrane protein" and
  "predicted signal peptide", in the namespace ``phobius_prediction``;
* ``predictions``: a human protein ``associated_with`` the term for each
  prediction it has: the transmembrane term for a protein with at least one
  predicted transmembrane helix (with the number of helices), the signal
  peptide term for a protein with a predicted signal peptide, and on both the
  topology string as Phobius wrote it. A protein with both predictions has two
  statements. A protein with neither (12,337 of the 20,350) has no statement:
  the module does not state that a protein is not transmembrane or has no
  signal peptide.

The predicate is ``associated_with``, not ``located_in``: a predicted helix is
a feature of the sequence, not a location. The legacy input also kept the
number of ``i`` and ``o`` in the topology as the cytoplasmic and
non-cytoplasmic parts; those counts are not carried (the topology string they
are counted from is).
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import OntologyClass, Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.phobius import (
    SIGNAL_PEPTIDE,
    TRANSMEMBRANE,
    iter_predictions,
    iter_terms,
)
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import Annotation, Entity, Identifier, Relation

config = ResourceConfig(
    id=ResourceCv.PHOBIUS,
    name='Phobius',
    full='Phobius transmembrane topology and signal peptide predictions',
    url='https://phobius.sbc.su.se/',
    # The file is OmniPath's own output of the Phobius program, which has terms
    # of its own. The Phobius site (https://phobius.sbc.su.se/download, checked
    # 2026-10-12) says of the program: "Phobius is freely available for local
    # installation for academic use" and "Phobius download via this page is
    # intended for academic users only"; it states nothing about the use of
    # predictions made with it. License as recorded in the legacy registry
    # (resources.json); the source states none today for the file.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    # Käll, Krogh, Sonnhammer: A combined transmembrane topology and signal
    # peptide prediction method. J Mol Biol 338(5):1027-1036, 2004.
    pubmed='15111065',
    topics=(Topic.SUBCELLULAR_LOCALIZATION,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,),
    summary=(
        'Predicted transmembrane proteins and signal peptides of the human '
        'proteome, computed by OmniPath with Phobius.'
    ),
    description=(
        'Phobius is a predictor of transmembrane helices and signal peptides '
        'that combines both in one hidden Markov model. The data are '
        "OmniPath's own run of Phobius on the human Swiss-Prot proteome "
        '(20,350 proteins; the file was last modified on 2020-03-30, and the '
        'Phobius version and the date of the run are not recorded): they are '
        'predictions, not measurements or curated annotations. A protein with '
        'at least one predicted transmembrane helix is associated with the '
        'term "predicted transmembrane protein", with the number of helices '
        'and the predicted topology; a protein with a predicted signal peptide '
        'is associated with "predicted signal peptide", with the topology. A '
        'protein with neither prediction has no statement: the data do not '
        'state that it is not transmembrane or has no signal peptide.'
    ),
)

# ================================== DOWNLOAD ==================================

# OmniPath's rescue copy of the file the legacy input read (the same file over
# http); 795,772 bytes.
download = Download(
    url='https://rescued.omnipathdb.org/phobius_human_proteome.tsv',
    filename='phobius_human_proteome.tsv',
    subfolder='phobius',
    large=False,
)

# =================================== SCHEMA ===================================

# The terms are the module's own (the file lists none), identified in a
# namespace of their own.
TERM_NAMESPACE = 'phobius_prediction'
TERMS = {
    TRANSMEMBRANE: (
        'predicted transmembrane protein',
        'A protein for which Phobius predicts at least one transmembrane '
        'helix. A protein without a predicted helix is not given a term: the '
        'data do not state that it is not transmembrane.',
    ),
    SIGNAL_PEPTIDE: (
        'predicted signal peptide',
        'A protein for which Phobius predicts a signal peptide. A protein '
        'without a predicted signal peptide is not given a term: the data do '
        'not state that it has none.',
    ),
}

# Attributes of a statement: the number of predicted helices (on the
# transmembrane statement) and the Phobius topology string, as written.
_HELICES = 'phobius:transmembrane_helices'
_TOPOLOGY = 'phobius:topology'

# The file is a run on the human proteome and every entry name ends in _HUMAN
# (the parser raises otherwise).
HUMAN = 'NCBITaxon:9606'


def _term(term: str) -> Entity:
    if term not in TERMS:
        raise ValueError(f'Unknown Phobius term: {term!r}')
    return Entity(
        type=OntologyClass,
        identifiers=[Identifier(type=TERM_NAMESPACE, value=term)],
    )


def map_term(row: dict[str, Any]) -> Entity:
    """One of the two terms, with its name and description."""
    name, description = TERMS[row['term']]
    return Entity(
        type=OntologyClass,
        identifiers=[
            Identifier(type=TERM_NAMESPACE, value=row['term']),
            Identifier(type=Namespace.NAME, value=name),
        ],
        annotations=[Annotation(term=slots.description, value=description)],
    )


def map_prediction(row: dict[str, Any]) -> Relation:
    """A protein associated with the term of one of its predictions."""
    annotations = []
    if row['term'] == TRANSMEMBRANE:
        # A count of helices (the TM column), without unit.
        annotations.append(
            Annotation(term=_HELICES, value=measurement(row['helices'], source_field='TM'))
        )
    annotations.append(Annotation(term=_TOPOLOGY, value=row['topology']))
    return Relation(
        subject=Entity(
            type=Protein,
            identifiers=[
                Identifier(type=Namespace.UNIPROT, value=row['accession']),
                Identifier(type=Namespace.UNIPROT_ENTRY, value=row['entry_name']),
            ],
            annotations=[Annotation(term=slots.in_taxon, value=HUMAN)],
        ),
        # A predicted helix or signal peptide is a feature of the sequence,
        # not a location of the protein.
        predicate=slots.associated_with,
        object=_term(row['term']),
        annotations=annotations,
    )


# ================================= RESOURCE ===================================

resource = Resource(
    config,
    terms=Dataset(
        download=None,
        mapper=map_term,
        raw_parser=iter_terms,
    ),
    predictions=Dataset(
        download=download,
        mapper=map_prediction,
        raw_parser=iter_predictions,
    ),
)
