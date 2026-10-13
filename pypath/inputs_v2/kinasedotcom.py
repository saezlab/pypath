"""
Parse the human kinome table of kinase.com and emit the kinase families.

kinase.com classifies the human protein kinases of Manning et al. 2002 in
groups (AGC, CAMK, CK1, CMGC, STE, TK, TKL, RGC, Atypical, Other), families
and subfamilies. The module emits each group as a gene family with its
families as members, each family with its subfamilies and the kinases that
the table puts in no subfamily, and each subfamily with its kinases. A kinase
with a second kinase domain is also a member of the family of that domain.

The kinases are genes, identified by Entrez Gene ID and symbol and by their
kinase.com accession. Kinase pseudogenes of the table are not emitted.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import Gene, GeneFamily, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.kinasedotcom import iter_kinase_tree
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
)

config = ResourceConfig(
    id=ResourceCv.KINASEDOTCOM,
    name='kinase.com',
    url='http://kinase.com/',
    # The disclaimer (http://kinase.com/about/Disclaimer.html) says: "User
    # shall use the information obtained from the Website only for his/her own
    # legitimate purposes [...]. User shall not sell any information User
    # obtains from the Website." The legacy registry (resources.json) records
    # the "Kinase.com" license, whose purpose is non-profit.
    license=LicenseCV.NONCOMMERCIAL,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='12471243',
    topics=(Topic.SIGNALING,),
    content=(Content.MOLECULAR_ANNOTATIONS,),
    # The article (Science 2002), abstract: "We have catalogued the protein kinase
    # complement of the human genome (the 'kinome') using public and proprietary
    # genomic, complementary DNA, and expressed sequence tag (EST) sequences", "We
    # identify 518 putative protein kinase genes" and "Classification and
    # comparison with model organism kinomes identified orthologous groups".
    knowledge_origin=(KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,),
    summary='Groups, families and subfamilies of the human protein kinases.',
    description=(
        'kinase.com, hosted by the Salk Institute, is the home of the human '
        'kinome of Manning et al. 2002, the complement of protein kinase '
        'genes of the human genome. Its table classifies each kinase by '
        'sequence of the kinase domain in a group, a family and, where there '
        'is one, a subfamily.'
    ),
)

# kinase.com writes the human kinome.
TAXON = 'NCBITaxon:9606'
# The kinase.com names of the levels have no Namespace member. A name is only
# unique within its level (the family RSK has the subfamily RSK).
NAMESPACES = {
    'group': 'kinasedotcom_group',
    'family': 'kinasedotcom_family',
    'subfamily': 'kinasedotcom_subfamily',
}
# The kinase.com accession of a kinase (``SK018``) has no Namespace member.
KINASE_ACCESSION = 'kinasedotcom'
# Which kinase domain's classification puts a kinase in a family: 1 for the
# first kinase domain, 2 for the second of a kinase with two.
KINASE_DOMAIN = 'kinasedotcom:kinase_domain'


def _kinase(member: dict[str, Any]) -> Entity:
    """A kinase gene. Entrez Gene ID and symbol, where the table has them."""
    identifiers = []
    if member['entrez']:
        identifiers.append(Identifier(type=Namespace.ENTREZ, value=member['entrez']))
    if member['symbol']:
        identifiers.append(Identifier(type=Namespace.GENESYMBOL, value=member['symbol']))
    identifiers.append(Identifier(type=KINASE_ACCESSION, value=member['sk']))
    # kinase.com's own name for the kinase is not always the symbol (CRIK is
    # CIT). It is the label of a kinase without a symbol and a synonym of the
    # others.
    if not member['symbol']:
        identifiers.append(Identifier(type=Namespace.NAME, value=member['name']))
    elif member['name'] != member['symbol']:
        identifiers.append(Identifier(type=Namespace.SYNONYM, value=member['name']))
    return Entity(
        type=Gene,
        identifiers=identifiers,
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
    )


def map_family(node: dict[str, Any]) -> Entity:
    """A group, family or subfamily with the next level and its kinases."""
    return Entity(
        type=GeneFamily,
        identifiers=[
            Identifier(type=NAMESPACES[node['level']], value=node['name']),
            Identifier(type=Namespace.NAME, value=node['name']),
        ],
        annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
        membership=[
            *(Membership(member=map_family(child)) for child in node['children']),
            *(
                Membership(
                    member=_kinase(member),
                    annotations=[
                        Annotation(term=KINASE_DOMAIN, value=member['kinase_domain'])
                    ],
                )
                for member in node['members']
            ),
        ],
    )


resource = Resource(
    config,
    families=Dataset(
        # The same file is on OmniPath's rescue server
        # (https://rescued.omnipathdb.org/Kincat_Hsap.08.02.xls, identical
        # MD5), which the legacy module used.
        download=Download(
            url='http://kinase.com/static/colt/data/human/kinome/tables/Kincat_Hsap.08.02.xls',
            filename='Kincat_Hsap.08.02.xls',
            subfolder='kinasedotcom',
            large=True,
            encoding=None,
            default_mode='rb',
        ),
        mapper=map_family,
        raw_parser=iter_kinase_tree,
    ),
)
