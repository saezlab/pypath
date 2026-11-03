"""
The membranes and topology types of the transmembrane proteins of TOPDB.

TOPDB, the Topology Data Bank of Transmembrane Proteins, collects the
experimentally derived topology of transmembrane proteins: for each protein
the membrane it is in with the compartments on its two sides, its topology
type (alpha-helical polytopic, alpha-helical bitopic or beta-barrel
polytopic) and the segments of its sequence that are inside, outside or in the
membrane, with the experiments behind them.

The module imports:

* ``terms``: the membranes and the topology types of the file, as TOPDB writes
  them (it has no identifiers for them), each in a namespace of its own;
* ``proteins``: for each protein the membranes it is in and its topology type.
  The statement of the topology type carries the number of transmembrane
  segments of the topology and its reliability.

The proteins are the UniProt entry names the file gives as IDs. The file does
not state the species and none is asserted. The segments of the topology
(inside, outside, membrane and others, with their residues) and the experiments
are not imported: see ``parsers/topdb.py``.

The topology of an entry is computed by CCTOP from the experimental constraints
that TOPDB collected for the protein (the site: "for each protein in the
database the most probable topology consistent with the collected experimental
constraints was also calculated using CCTOP"); the number of segments and the
reliability belong to that topology.

The legacy input mapped the entry names to UniProt accessions, filtered one
species, and kept each compartment name of a membrane (``Cytoplasm``) as a
membrane of its own. It exposed the membrane, the topology (the location of
each region), a score (the reliability) and the number of transmembrane regions;
the last two are carried now.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import OntologyClass, Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.topdb import iter_entries, iter_terms
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Association,
    Entity,
    EntityRef,
    Identifier,
    Relation,
)

config = ResourceConfig(
    id=ResourceCv.TOPDB,
    name='TOPDB',
    full='Topology Data Bank of Transmembrane Proteins',
    url='https://topdb.unitmp.org/',
    # The site states no license for the data (home and download pages checked
    # 2026-10-12; the UniTmp pages say "All rights reserved" for the site).
    # License as recorded in the legacy registry (resources.json); the source
    # states none today.
    license=LicenseCV.CC_BY_NC_4_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='25392424',
    topics=(Topic.SUBCELLULAR_LOCALIZATION,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(
        KnowledgeOrigin.EXPERIMENTAL,
        KnowledgeOrigin.LITERATURE_CURATED,
        # The topology of each entry is computed by CCTOP from the
        # experimental constraints collected for the protein.
        KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,
    ),
    summary=(
        'Membranes and topology types of transmembrane proteins of all '
        'species, with the number of transmembrane segments and the '
        'reliability of the topology.'
    ),
    description=(
        'TOPDB (Topology Data Bank of Transmembrane Proteins) is a database of '
        'transmembrane proteins with experimentally derived topology: for '
        'each protein the segments of its sequence that are inside, outside '
        'or in the membrane, collected from the literature and from the '
        'structures of the PDBTM database, with the membrane the protein is '
        'in. The topology of each entry is the most probable one consistent '
        'with the collected experimental constraints, computed with CCTOP. '
        'The module imports the membrane and the topology type of each '
        'protein, with the number of transmembrane segments of the topology '
        'and its reliability (as TOPDB gives it, up to 100). It is part of '
        'the UniTmp resources.'
    ),
)

# ================================== DOWNLOAD ==================================

# Release 3.3, all species (256 MB). The legacy input read the same file.
download = Download(
    url='https://topdb.unitmp.org/data/TOPDB/data/v3.3/sets/topdb_all.xml',
    filename='topdb_all.xml',
    subfolder='topdb',
    large=True,
    encoding=None,
    default_mode='rb',
)

# =================================== SCHEMA ===================================

# TOPDB names its membranes ("Cell membrane") and its topology types
# ("Alpha_Polytopic") and has no identifiers for them: the term is identified by
# the value, as written.
MEMBRANE = 'topdb_membrane'
TYPE = 'topdb_type'
TERM_NAMESPACES = {'membrane': MEMBRANE, 'type': TYPE}

# Attributes of the statement of the topology type. The segment count is the
# entry's ``Numtm`` (a count, a Measurement without unit), equal to its number of Membrane regions (the parser
# raises otherwise). The reliability is the value as TOPDB gives it, observed
# range about 24.6 to 100. It belongs to the topology that TOPDB computed with
# CCTOP from the experimental constraints, and is not a confidence in the
# protein's membrane. The site (https://topdb.unitmp.org/documents/reliability,
# read 2026-10-12): "Reliability is the average of the posterior probabilities
# along the most probably state path (q)", defined as
# R = 100 * (sum of the posterior probabilities along the path) / n, with n the
# sequence length. It has no unit.
SEGMENTS = 'topdb:transmembrane_segments'
RELIABILITY = 'topdb:reliability'

# The two compartments the file names after a membrane
# (``Cell membrane;Cytoplasm;Extracellular``), in the order it writes them.
SIDE_1 = 'topdb:side_1'
SIDE_2 = 'topdb:side_2'


def map_term(row: dict[str, Any]) -> Entity:
    """A membrane, with the compartments on its sides, or a topology type."""
    kind = row['kind']
    if kind not in TERM_NAMESPACES:
        raise ValueError(f'Unknown TOPDB term kind: {kind!r}')
    annotations = None
    if kind == 'membrane':
        annotations = [
            Annotation(SIDE_1, row['side_1']),
            Annotation(SIDE_2, row['side_2']),
        ]
    return Entity(
        type=OntologyClass,
        identifiers=[
            Identifier(TERM_NAMESPACES[kind], row['id']),
            Identifier(Namespace.NAME, row['name']),
        ],
        annotations=annotations,
    )


def exclusion_reason(row: dict[str, Any]) -> str | None:
    """Why a parsed entry gives no statement; None for an entry that is mapped."""
    if row['kind'] != 'uniprot_entry':
        # The entries made from a PDB structure (``2k3c_A``) have no protein
        # identifier, only the PDB code and the chain of the structure.
        return 'no UniProt entry name'
    return None


def map_protein(row: dict[str, Any]) -> Relation | None:
    """A protein associated with its topology type, which carries the number of
    transmembrane segments and the reliability of the topology; the protein is
    located in each of its membranes. None for an entry without a UniProt entry
    name."""
    if exclusion_reason(row):
        return None
    identifiers = [Identifier(Namespace.UNIPROT_ENTRY, row['id'])]
    if row.get('name'):
        identifiers.append(Identifier(Namespace.NAME, row['name']))
    # The file lists the membranes the protein is in: it is located in each
    # (a protein can be in several, as in the endoplasmic reticulum and the
    # plasma membrane). The membranes have nothing attached, so they are
    # associations of the protein, the subject of the type statement.
    protein = Entity(
        type=Protein,
        identifiers=identifiers,
        associations=[
            Association(
                object=EntityRef(OntologyClass, MEMBRANE, membrane['name']),
                predicate=slots.located_in,
            )
            for membrane in row['membranes']
        ],
    )
    # The type classifies the topology; the statement is a relation because the
    # count and the reliability of that topology are attached to it.
    return Relation(
        subject=protein,
        predicate=slots.associated_with,
        object=Entity(
            type=OntologyClass,
            identifiers=[Identifier(TYPE, row['type'])],
        ),
        annotations=[
            # A count of transmembrane segments, without unit.
            Annotation(SEGMENTS, measurement(row['segments'], source_field='Numtm Count')),
            Annotation(
                RELIABILITY,
                measurement(row['reliability'], source_field='Reliability'),
            ),
        ],
    )


# ================================= RESOURCE ===================================

resource = Resource(
    config,
    terms=Dataset(
        download=download,
        mapper=map_term,
        raw_parser=iter_terms,
    ),
    proteins=Dataset(
        download=download,
        mapper=map_protein,
        raw_parser=iter_entries,
    ),
)
