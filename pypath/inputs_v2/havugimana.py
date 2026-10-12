"""
Parse the complexes of Havugimana et al. 2012 and emit Entity records.

"A census of human soluble protein complexes" (Cell 150:1068-1081,
PMID 22939629) identified 622 putative complexes of soluble human proteins.
The authors fractionated cytoplasmic and nuclear extracts of HeLa S3 and
HEK293 cells by several biochemical procedures, measured the co-elution of
proteins by mass spectrometry, scored the physical interactions with a
machine learning procedure that also uses functional association evidence
(correlated evolutionary rates and functional genomics data of four species)
and clustered the high-confidence interactions into complexes.

Table S3A of the supplement lists the complexes with their subunits. The
table calls them "predicted" and "putative". This module reads that table:
each complex with its subunits, as UniProtKB accessions with the entry names.
The comparison of the complexes with the CORUM, Reactome, PINdb, HPRD, GO,
MEMO, fly and yeast complexes, the novelty of the subunits and the
conservation group of each complex are not mapped.

The supplement is read from the OmniPath rescue copy: the publisher's file
sits behind a bot challenge (Cloudflare) and the PubMed Central copy of the
author manuscript behind a proof-of-work page.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel.model import MacromolecularComplex, Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.havugimana import iter_havugimana_complexes
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
)

config = ResourceConfig(
    id=ResourceCv.HAVUGIMANA,
    name='Havugimana2012',
    full='A census of human soluble protein complexes',
    url='https://doi.org/10.1016/j.cell.2012.08.011',
    # The article and its supplement are "(c) 2012 Elsevier Inc. All rights
    # reserved"; they state no license.
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='22939629',
    content=(Content.PROTEIN_COMPLEXES,),
    knowledge_origin=(
        KnowledgeOrigin.EXPERIMENTAL,
        KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,
    ),
    summary='622 putative soluble human protein complexes from biochemical fractionation, with their subunits.',
    description=(
        'Havugimana et al. fractionated cytoplasmic and nuclear extracts of HeLa S3 '
        'and HEK293 cells by complementary biochemical procedures, measured the '
        'co-elution of proteins by mass spectrometry, scored physical interactions '
        'with a machine learning procedure that also uses evolutionary and '
        'functional genomics evidence, and clustered them into 622 putative '
        'complexes of soluble human proteins, many of them new.'
    ),
)

# The complex ids of the table (``C_1``) have no Namespace member yet.
COMPLEX_NAMESPACE = 'havugimana2012'
# The entry names of all subunits are ``_HUMAN``.
TAXON = 'NCBITaxon:9606'


def map_complex(row: dict[str, Any]) -> Entity | None:
    """A complex with its UniProtKB subunits; None without subunits."""
    if not row.get('accessions'):
        return None
    return Entity(
        type=MacromolecularComplex,
        identifiers=[Identifier(type=COMPLEX_NAMESPACE, value=row['complex_id'])],
        membership=[
            Membership(
                member=Entity(
                    type=Protein,
                    identifiers=[
                        Identifier(type=Namespace.UNIPROT, value=accession),
                        Identifier(type=Namespace.UNIPROT_ENTRY, value=entry_name),
                    ],
                    annotations=[Annotation(term=slots.in_taxon, value=TAXON)],
                )
            )
            for accession, entry_name in zip(row['accessions'], row['entry_names'])
        ],
    )


# Supplementary Table S3 of the article (mmc3.xls), kept on the OmniPath rescue server.
resource = Resource(
    config,
    complexes=Dataset(
        download=Download(
            url='https://rescued.omnipathdb.org/havugimana2012-mmc3.xls',
            filename='havugimana2012-mmc3.xls',
            subfolder='havugimana',
            default_mode='rb',
            encoding=None,
        ),
        mapper=map_complex,
        raw_parser=iter_havugimana_complexes,
    ),
)
