"""
Parse 3did domain-domain and domain-motif interactions and emit Relation records.

3did, the database of three-dimensional interacting domains, lists pairs of
Pfam domains, and a Pfam domain with a linear motif, that are in contact in
structures of the Protein Data Bank. Its flat files describe each structure by
the PDB entry and the chains and residue ranges of the partners in the
numbering of the PDB file. They name no protein: the file does not say which
UniProt protein a chain is, and a mapping of chains to UniProt (SIFTS) is
another resource. The legacy input translated the chains with it and made
protein-level domain interactions; this module does not.

What the files state is an interaction between domain families, so the
partners are the Pfam domains (``ProteinDomain``) and, for the domain-motif
file, 3did's motif class (an ``OntologyClass`` of the source's own vocabulary,
with its pattern). A record is one structural instance: the relation of the
two entries, with the PDB entry as the identifier of the observation and the
chains, ranges and InterPreTS scores as attributes. Instances of one pair
merge into one statement with the structures as its evidence. The residue
contacts listed below each instance of ``3did_flat`` are not read, and 3did's
interface files (``3did_interface_flat``, ``3did_global_interface_flat``) are
not used.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

from biolink_model.datamodel.model import OntologyClass, ProteinDomain, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.threedid import iter_domain_domain, iter_domain_motif
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    IdentifiersBuilder,
    RelationBuilder,
)

config = ResourceConfig(
    id=ResourceCv.THREEDID,
    name='3did',
    full='3did: database of three-dimensional interacting domains',
    url='https://3did.irbbarcelona.org/',
    # The site (download, help and about pages) states no license or terms of use.
    license=LicenseCV.UNSPECIFIED,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='24081580',
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(KnowledgeOrigin.EXPERIMENTAL,),
    summary=(
        'Pfam domain pairs, and Pfam domains with linear motifs, in contact '
        'in PDB structures, with the structures.'
    ),
    description=(
        '3did is a catalogue of high-resolution three-dimensional structural '
        'templates for domain-domain interactions and domain-peptide '
        'interactions, derived from the Protein Data Bank. Its release of '
        'March 2025 (Pfam 37.0, PDB 2024_12) has 20,644 interacting Pfam '
        'domain pairs and 1,657 motifs in interactions of known structure, '
        'which are peptide-binding sites found by the method of Stein and '
        'Aloy (2010). Each interaction comes with its instances in PDB '
        'structures, scored with InterPreTS.'
    ),
)

# The current release; the files are not versioned by name.
download_domain_domain = Download(
    url='https://3did.irbbarcelona.org/download/current/3did_flat.gz',
    filename='3did_flat.gz',
    subfolder='threedid',
    large=True,
    ext='gz',
)
download_domain_motif = Download(
    url='https://3did.irbbarcelona.org/download/current/3did_dmi_flat.gz',
    filename='3did_dmi_flat.gz',
    subfolder='threedid',
    large=True,
    ext='gz',
)

# Pfam's short names of its entries (``SH2``, ``1-cysPrx_C``) are its own
# identifiers; 3did writes the domain-motif domains with them only.
_PFAM_NAME = 'pfam_id'
_PFAM = 'pfam'
# 3did's names of its motif classes (``PCNA_C_LIG_0-0``).
_MOTIF = 'threedid_motif'

# Attributes of the relation; Biolink has no slot for them. The score and the
# Z-score are InterPreTS's (Aloy and Russell 2003) for the interaction in the
# structure. The regions are ``chain:start-end`` in the numbering of the PDB file.
_REGION_1 = 'threedid:domain_1_region'
_REGION_2 = 'threedid:domain_2_region'
_SCORE = 'threedid:interprets_score'
_Z_SCORE = 'threedid:interprets_z_score'
_TOPOLOGY = 'threedid:topology'
_DOMAIN_REGION = 'threedid:domain_region'
_MOTIF_REGION = 'threedid:motif_region'
_PEPTIDE = 'threedid:motif_instance'
_CONTACTS = 'threedid:contextual_contacts'
_PATTERN = 'threedid:motif_pattern'


def _region(region: dict[str, Any]) -> str:
    return f"{region['chain']}:{region['start']}-{region['end']}"


def _domain(
    pfam: Callable[[dict[str, Any]], str | None], name: Callable[[dict[str, Any]], str]
) -> EntityBuilder:
    def identifiers(row: dict[str, Any]) -> list[tuple[str, str]]:
        accession = pfam(row)
        pairs = [(_PFAM, accession)] if accession else []
        return [*pairs, (_PFAM_NAME, name(row)), (Namespace.NAME, name(row))]

    return EntityBuilder(
        entity_type=ProteinDomain,
        identifiers=IdentifiersBuilder(CV.from_pairs(identifiers)),
    )


def _pdb_entry(row: dict[str, Any]) -> str:
    return row['pdb_id']


_structure = IdentifiersBuilder(CV(term=Namespace.PDB, value=_pdb_entry))

# The first domain of the pair of the file is the subject. The domains touch
# each other in the structure: a direct physical interaction, with the contact
# the only thing the file asserts about it.
domain_domain_schema = RelationBuilder(
    subject=_domain(lambda row: row['pfam_1'], lambda row: row['name_1']),
    predicate=slots.directly_physically_interacts_with,
    object=_domain(lambda row: row['pfam_2'], lambda row: row['name_2']),
    identifiers=_structure,
    annotations=AnnotationsBuilder(
        CV(term=_REGION_1, value=lambda row: _region(row['region_1'])),
        CV(term=_REGION_2, value=lambda row: _region(row['region_2'])),
        CV(term=_SCORE, value=lambda row: row['score']),
        CV(term=_Z_SCORE, value=lambda row: row['z_score']),
        CV(term=_TOPOLOGY, value=lambda row: row.get('topology')),
    ),
)

_motif_class = EntityBuilder(
    entity_type=OntologyClass,
    identifiers=IdentifiersBuilder(
        CV(term=_MOTIF, value=lambda row: row['motif']),
        CV(term=Namespace.NAME, value=lambda row: row['motif']),
    ),
    annotations=AnnotationsBuilder(
        CV(term=_PATTERN, value=lambda row: row.get('pattern')),
    ),
)

domain_motif_schema = RelationBuilder(
    subject=_domain(lambda row: None, lambda row: row['domain']),
    predicate=slots.directly_physically_interacts_with,
    object=_motif_class,
    identifiers=_structure,
    annotations=AnnotationsBuilder(
        CV(term=_DOMAIN_REGION, value=lambda row: _region(row['domain_region'])),
        CV(term=_MOTIF_REGION, value=lambda row: _region(row['motif_region'])),
        CV(term=_PEPTIDE, value=lambda row: row['peptide']),
        CV(term=_CONTACTS, value=lambda row: row['contacts']),
        CV(term=_TOPOLOGY, value=lambda row: row['topology']),
    ),
)

resource = Resource(
    config,
    domain_domain=Dataset(
        download=download_domain_domain,
        mapper=domain_domain_schema,
        raw_parser=iter_domain_domain,
    ),
    domain_motif=Dataset(
        download=download_domain_motif,
        mapper=domain_motif_schema,
        raw_parser=iter_domain_motif,
    ),
)
