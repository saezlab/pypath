"""
Parse the human integrin subunits of Takada et al. 2007 and emit the integrin family with its members.

Integrins are the cell adhesion receptors that are alpha-beta heterodimers; the
review "The integrins" lists the 18 alpha and 8 beta subunits of humans in its
Table 1, with their gene symbols and protein accessions. The module emits one
family, the integrins, with each subunit as a member (the ``families`` dataset).
"""

from __future__ import annotations

import re
from typing import Any

from biolink_model.datamodel.model import Protein, ProteinFamily, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.integrins import iter_family
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
)

config = ResourceConfig(
    id=ResourceCv.INTEGRINS,
    name='Integrins',
    full='The integrins (Takada et al. 2007)',
    url='https://doi.org/10.1186/gb-2007-8-5-215',
    # License as recorded in the legacy registry (resources.json); the source states none today
    # (Europe PMC's copy of the article says "Copyright (c) 2007 BioMed Central Ltd" only, and
    # Europe PMC lists the article as open access without naming a license; checked 2026-10-13).
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='17543136',
    topics=(Topic.CELL_CELL_COMMUNICATION,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    knowledge_origin=(KnowledgeOrigin.LITERATURE_CURATED,),
    summary='The family of the human integrin subunits, with each subunit as a member.',
    description=(
        'Integrins are a superfamily of cell adhesion receptors that bind '
        'extracellular matrix ligands, cell-surface ligands and soluble '
        'ligands. They are transmembrane alpha-beta heterodimers; at least 18 '
        'alpha and eight beta subunits are known in humans, generating 24 '
        'heterodimers. Table 1 of the review lists the human subunits with '
        'their gene symbols and protein accessions.'
    ),
)

# Table 1 sits in the article, which Europe PMC serves as XML (the legacy URL, the
# table page of PubMed Central, now asks a script to solve a CAPTCHA).
download = Download(
    url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC1929136/fullTextXML',
    filename='takada2007_integrins_PMC1929136.xml',
    subfolder='integrins',
    default_mode='rb',
    encoding=None,
)

# The family is identified in a namespace of the source, as the other family ports
# do (``kinasedotcom_family``, ``gpcrdb_family``, ``hgnc_group``): the bare name
# "Integrins" would merge with a family of that name from another resource.
FAMILY_NAMESPACE = 'integrins_family'
FAMILY_NAME = 'Integrins'
TAXON = 'NCBITaxon:9606'
# The shape of a UniProt accession. The table's protein accession of ITGA4 is
# a GenBank protein accession (AAB25486).
_UNIPROT = re.compile(
    r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$'
)


def _subunit(row: dict[str, Any]) -> Entity:
    accession = row['accession']
    namespace = Namespace.UNIPROT if _UNIPROT.match(accession) else Namespace.GENBANK
    return Entity(
        type=Protein,
        identifiers=[
            Identifier(type=namespace, value=accession),
            Identifier(type=Namespace.GENESYMBOL, value=row['symbol']),
        ],
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
    )


def map_family(row: dict[str, Any]) -> Entity:
    """The family of the integrin subunits."""
    return Entity(
        type=ProteinFamily,
        identifiers=[
            Identifier(type=FAMILY_NAMESPACE, value=FAMILY_NAME),
            Identifier(type=Namespace.NAME, value=FAMILY_NAME),
        ],
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
        membership=[
            Membership(member=_subunit(member), predicate=slots.has_member)
            for member in row['members']
        ],
    )


resource = Resource(
    config,
    families=Dataset(
        download=download,
        mapper=map_family,
        raw_parser=iter_family,
    ),
)
