#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#  This file is part of the `pypath` python module
#
#  Copyright 2014-2025
#  EMBL, EMBL-EBI, Uniklinik RWTH Aachen, Heidelberg University
#
#  Authors: see the file `README.rst`
#  Contact: Dénes Türei (turei.denes@gmail.com)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      https://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://pypath.omnipathdb.org/
#

"""
Parse Guide to Pharmacology data and emit Entity records.

This module converts Guide to Pharmacology ligand-target interaction data
into Entity records using the schema defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator
import csv

from pypath.share.downloads import download_and_open
from pypath.internals.silver_schema import Entity as SilverEntity, Identifier, Annotation
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    PharmacologicalActionCv,
    LigandTypeCv,
    AffinityUnitCv,
    MoleculeAnnotationsCv,
    ResourceAnnotationCv,
    ResourceCv,
)
from ..internals.tabular_builder import (
    Annotations,
    Column,
    Entity,
    Identifiers,
    Member,
    Members,
)


# Mapping for Action column values to PharmacologicalActionCv
action_cv_mapping = {
    'Activation': PharmacologicalActionCv.ACTIVATION,
    'Agonist': PharmacologicalActionCv.AGONIST,
    'Antagonist': PharmacologicalActionCv.ANTAGONIST,
    'Biased agonist': PharmacologicalActionCv.BIASED_AGONIST,
    'Binding': PharmacologicalActionCv.BINDING,
    'Biphasic': PharmacologicalActionCv.BIPHASIC,
    'Competitive': PharmacologicalActionCv.COMPETITIVE,
    'Feedback inhibition': PharmacologicalActionCv.FEEDBACK_INHIBITION,
    'Full agonist': PharmacologicalActionCv.FULL_AGONIST,
    'Inhibition': PharmacologicalActionCv.INHIBITION,
    'Inverse agonist': PharmacologicalActionCv.INVERSE_AGONIST,
    'Irreversible agonist': PharmacologicalActionCv.IRREVERSIBLE_AGONIST,
    'Irreversible inhibition': PharmacologicalActionCv.IRREVERSIBLE_INHIBITION,
    'Mixed': PharmacologicalActionCv.MIXED,
    'Negative': PharmacologicalActionCv.NEGATIVE,
    'Neutral': PharmacologicalActionCv.NEUTRAL,
    'Non-competitive': PharmacologicalActionCv.NON_COMPETITIVE,
    'None': PharmacologicalActionCv.NONE,
    'Partial agonist': PharmacologicalActionCv.PARTIAL_AGONIST,
    'Pore blocker': PharmacologicalActionCv.PORE_BLOCKER,
    'Positive': PharmacologicalActionCv.POSITIVE,
    'Potentiation': PharmacologicalActionCv.POTENTIATION,
    'Slows inactivation': PharmacologicalActionCv.SLOWS_INACTIVATION,
    'Unknown': PharmacologicalActionCv.UNKNOWN,
    'Voltage-dependent inhibition': PharmacologicalActionCv.VOLTAGE_DEPENDENT_INHIBITION,
}

# Mapping for Type column values to LigandTypeCv
type_cv_mapping = {
    'Activator': LigandTypeCv.ACTIVATOR,
    'Agonist': LigandTypeCv.AGONIST,
    'Allosteric modulator': LigandTypeCv.ALLOSTERIC_MODULATOR,
    'Antagonist': LigandTypeCv.ANTAGONIST,
    'Antibody': LigandTypeCv.ANTIBODY,
    'Channel blocker': LigandTypeCv.CHANNEL_BLOCKER,
    'Fusion protein': LigandTypeCv.FUSION_PROTEIN,
    'Gating inhibitor': LigandTypeCv.GATING_INHIBITOR,
    'Inhibitor': LigandTypeCv.INHIBITOR,
    'None': LigandTypeCv.NONE,
    'Subunit-specific': LigandTypeCv.SUBUNIT_SPECIFIC,
}

# Mapping for Affinity Units column values to AffinityUnitCv
affinity_units_cv_mapping = {
    'pA2': AffinityUnitCv.PA2,
    'pEC50': AffinityUnitCv.PEC50,
    'pIC50': AffinityUnitCv.PIC50,
    'pKB': AffinityUnitCv.PKB,
    'pKd': AffinityUnitCv.PKD,
    'pKi': AffinityUnitCv.PKI,
}

# Mapping for Ligand Type (chemical nature) to EntityTypeCv
ligand_chemical_type_mapping = {
    'Synthetic organic': EntityTypeCv.SYNTHETIC_ORGANIC,
    'Natural product': EntityTypeCv.NATURAL_PRODUCT,
    'Metabolite': EntityTypeCv.METABOLITE,
    'Inorganic': EntityTypeCv.INORGANIC,
    'Peptide': EntityTypeCv.PEPTIDE,
    'Antibody': EntityTypeCv.ANTIBODY,
    'Nucleic acid': EntityTypeCv.NUCLEIC_ACID,
}

# Mapping for Target Type to EntityTypeCv
target_type_mapping = {
    'gpcr': EntityTypeCv.GPCR,
    'lgic': EntityTypeCv.LGIC,
    'vgic': EntityTypeCv.VGIC,
    'other_ic': EntityTypeCv.OTHER_ION_CHANNEL,
    'enzyme': EntityTypeCv.ENZYME,
    'catalytic_receptor': EntityTypeCv.CATALYTIC_RECEPTOR,
    'nhr': EntityTypeCv.NUCLEAR_HORMONE_RECEPTOR,
    'transporter': EntityTypeCv.TRANSPORTER,
    'other_protein': EntityTypeCv.OTHER_PROTEIN,
}

# Mapping for Target Species to NCBI Taxonomy IDs TODO need to be checked again.
species_to_taxid = {
    'Human': '9606',
    'Mouse': '10090',
    'Rat': '10116',
    'Rabbit': '9986',
    'Guinea pig': '10141',
    'Hamster': '10036',
    'Dog': '9615',
    'Pig': '9823',
    'Bovine': '9913',
    'Sheep': '9940',
    'Chicken': '9031',
    'Turkey': '9103',
    'Monkey': '9544',  # Rhesus macaque (common lab monkey)
    'Gorilla': '9593',
    'Ferret': '9669',
    'Honeybee': '7460',
    'spiny starfish': '7609',  # Marthasterias glacialis
    'Escherichia coli': '562',
    'Mycobacterium tuberculosis': '1773',
    'Plasmodium falciparum': '5833',
    'Plasmodium falciparum 3D7': '36329',
    'Plasmodium falciparum 7G8': '57266',
    'Plasmodium falciparum A2': '5833',  # Use parent taxid as strain not in NCBI
    'Plasmodium falciparum D6': '478860',
    'Plasmodium falciparum Dd2': '57267',
    'Plasmodium falciparum FC27': '5833',
    'Plasmodium falciparum FCR3': '5833',  # Use parent taxid as strain not in NCBI
    'Plasmodium falciparum HB3': '137071',
    'Plasmodium falciparum K1': '5839',
    'Plasmodium falciparum NF54': '5833',
    'Plasmodium falciparum TM4': '5833',
    'Plasmodium falciparum TM90C2A': '5833',
    'Plasmodium falciparum TM90C2B': '5833',
    'Plasmodium falciparum TM91C235': '5833',
    'Plasmodium falciparum V1/S': '5833',
    'Plasmodium falciparum W2': '5833',
    'Plasmodium vivax': '5855',
    'Plasmodium berghei': '5821',
    'Plasmodium cynomolgi': '5827',
    'Plasmodium knowlesi': '5850',
    'Plasmodium yoelii': '5861',
    'SARS-CoV': '694009',
    'SARS-CoV-2': '2697049',
    'MERS-CoV': '1335626',
    'HCoV-OC43': '31631',
    'Hepatitis C virus': '3052230',
    'Zika virus': '64320',
}


def guidetopharma() -> Generator[SilverEntity]:
    """
    Yield resource metadata as an Entity record.

    Yields:
        Entity record with type CV_TERM containing Guide to Pharmacology metadata.
    """
    yield SilverEntity(
        type=EntityTypeCv.CV_TERM,
        identifiers=[
            Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=ResourceCv.GUIDETOPHARMA),
            Identifier(type=IdentifierNamespaceCv.NAME, value='Guide to Pharmacology'),
        ],
        annotations=[
            Annotation(term=ResourceAnnotationCv.LICENSE, value=str(LicenseCV.CC_BY_SA_4_0)),
            Annotation(term=ResourceAnnotationCv.UPDATE_CATEGORY, value=str(UpdateCategoryCV.REGULAR)),
            Annotation(term=IdentifierNamespaceCv.PUBMED, value='37953350'),
            Annotation(term=ResourceAnnotationCv.URL, value='https://www.guidetopharmacology.org/'),
            Annotation(term=ResourceAnnotationCv.DESCRIPTION, value=(
                'The IUPHAR/BPS Guide to PHARMACOLOGY is an expert-curated resource '
                'of ligand-activity-target relationships, providing quantitative '
                'information on drug targets and the prescription medicines and '
                'experimental drugs that act on them. It covers G protein-coupled '
                'receptors, voltage-gated ion channels, ligand-gated ion channels, '
                'nuclear hormone receptors, catalytic receptors, enzymes, and '
                'transporters.'
            )),
        ],
    )


def guidetopharma_targets() -> Generator[SilverEntity, None, None]:
    """
    Download and parse Guide to Pharmacology targets as Entity records.

    Yields:
        Entity records representing protein targets with cross-species identifiers.
    """

    # Download the targets table
    url = 'http://www.guidetopharmacology.org/DATA/targets_and_families.csv'
    opener = download_and_open(
        url=url,
        filename='targets_and_families.csv',
        subfolder='guidetopharma',
    )

    # Define the schema mapping
    # TODO: need to think how we deal with this: they all get the same guide to pharma ID but have different species. So guide to pharma ID is not merge-safe. Maybe with our check to not merge when different species this will be solved? Then we need to create separate Entities + Interactions for each species variant.
    schema = Entity(
        entity_type=Column('Type', cv=target_type_mapping),
        identifiers=Identifiers(
            Column('Target id', cv=IdentifierNamespaceCv.GUIDETOPHARMA),
            # Human identifiers
            Column('Human SwissProt', cv=IdentifierNamespaceCv.UNIPROT),
            Column('Human Ensembl Gene', cv=IdentifierNamespaceCv.ENSEMBL),
            Column('Human Entrez Gene', cv=IdentifierNamespaceCv.ENTREZ),
            Column('Human protein RefSeq', cv=IdentifierNamespaceCv.REFSEQ_PROTEIN),
            Column('HGNC id', cv=IdentifierNamespaceCv.HGNC),
            # Rat identifiers
            #Column('Rat SwissProt', cv=IdentifierNamespaceCv.UNIPROT),
            #Column('Rat Ensembl Gene', cv=IdentifierNamespaceCv.ENSEMBL),
            #Column('Rat Entrez Gene', cv=IdentifierNamespaceCv.ENTREZ),
            #Column('Rat protein RefSeq', cv=IdentifierNamespaceCv.REFSEQ_PROTEIN),
            # Mouse identifiers
            #Column('Mouse SwissProt', cv=IdentifierNamespaceCv.UNIPROT),
            #Column('Mouse Ensembl Gene', cv=IdentifierNamespaceCv.ENSEMBL),
            #Column('Mouse Entrez Gene', cv=IdentifierNamespaceCv.ENTREZ),
            #Column('Mouse protein RefSeq', cv=IdentifierNamespaceCv.REFSEQ_PROTEIN),
            # Names
            Column('Target name', cv=IdentifierNamespaceCv.NAME),
            Column('HGNC symbol', cv=IdentifierNamespaceCv.GENE_NAME_PRIMARY),
            Column('synonyms', delimiter='|', cv=IdentifierNamespaceCv.SYNONYM),
            Column('Target systematic name', cv=IdentifierNamespaceCv.SYSTEMATIC_NAME),
            Column('Target abbreviated name', cv=IdentifierNamespaceCv.ABBREVIATED_NAME),

        ),
        annotations=Annotations(
            # Source annotation
            Column('Species', term_cv=IdentifierNamespaceCv.NCBI_TAX_ID, value="9606"),  # Human
            Column('Family name', cv=MoleculeAnnotationsCv.PROTEIN_FAMILY),
            Column('HGNC name', cv=MoleculeAnnotationsCv.DESCRIPTION),
        ),
    )

    # Parse and yield entities
    # Skip the first line (version info)
    next(opener.result)
    reader = csv.DictReader(opener.result)
    for row in reader:
        yield schema(row)


def guidetopharma_ligands() -> Generator[SilverEntity, None, None]:
    """
    Download and parse Guide to Pharmacology ligands as Entity records.

    Yields:
        Entity records representing ligands/small molecules.
    """

    # Download the ligands table
    url = 'http://www.guidetopharmacology.org/DATA/ligands.csv'
    opener = download_and_open(
        url=url,
        filename='ligands.csv',
        subfolder='guidetopharma',
    )

    # Define the schema mapping
    schema = Entity(
        entity_type=Column('Type', cv=ligand_chemical_type_mapping),
        identifiers=Identifiers(
            Column('Ligand ID', cv=IdentifierNamespaceCv.GUIDETOPHARMA),
            Column('PubChem SID', cv=IdentifierNamespaceCv.PUBCHEM_COMPOUND),
            Column('PubChem CID', cv=IdentifierNamespaceCv.PUBCHEM),
            Column('UniProt ID', cv=IdentifierNamespaceCv.UNIPROT),
            Column('Ensembl ID', cv=IdentifierNamespaceCv.ENSEMBL),
            Column('ChEMBL ID', cv=IdentifierNamespaceCv.CHEMBL),
            Column('InChIKey', cv=IdentifierNamespaceCv.STANDARD_INCHI_KEY),
            Column('Synonyms', cv=IdentifierNamespaceCv.SYNONYM),
            Column('Name', cv=IdentifierNamespaceCv.NAME),
            Column('IUPAC name', cv=IdentifierNamespaceCv.IUPAC_NAME),
            Column('SMILES', cv=IdentifierNamespaceCv.SMILES),
            Column('InChI', cv=IdentifierNamespaceCv.STANDARD_INCHI),
            Column('INN', cv=IdentifierNamespaceCv.INN),
            Column('GtoImmuPdb', cv=IdentifierNamespaceCv.GTO_IMMU_PDB),
            Column('GtoMPdb', cv=IdentifierNamespaceCv.GTO_M_PDB),

        ),
        annotations=Annotations(
            # Source annotation
            Column('Species', cv=species_to_taxid, term_cv=IdentifierNamespaceCv.NCBI_TAX_ID),
            Column('Approved', cv={'yes': MoleculeAnnotationsCv.APPROVED}),
            Column('Withdrawn', cv={'yes': MoleculeAnnotationsCv.WITHDRAWN}),
            Column('Labelled', cv={'yes': MoleculeAnnotationsCv.LABELLED}),
            Column('Radioactive', cv={'yes': MoleculeAnnotationsCv.RADIOACTIVE}),
            Column('Antibacterial', cv={'yes': MoleculeAnnotationsCv.ANTIBACTERIAL}),
        ),
    )

    # Parse and yield entities
    # Skip the first line (version info)
    next(opener.result)
    reader = csv.DictReader(opener.result)
    for row in reader:
        yield schema(row)


def guidetopharma_interactions() -> Generator[SilverEntity, None, None]:
    """
    Download and parse Guide to Pharmacology ligand-target interactions as Entity records.

    Yields:
        Entity records with type INTERACTION, representing ligand-target interactions.
    """

    # Download the interactions table
    url = 'http://www.guidetopharmacology.org/DATA/interactions.csv'
    opener = download_and_open(
        url=url,
        filename='interactions.csv',
        subfolder='guidetopharma',
    )

    # Define the schema mapping
    # Note: GuideToPharmacology doesn't provide explicit interaction IDs,
    # so we create a composite identifier from Target ID and Ligand ID
    schema = Entity(
        entity_type=EntityTypeCv.INTERACTION,
        identifiers=Identifiers(
            Column(
                lambda row: f"{row['Target ID']}_{row['Ligand ID']}" if row.get('Target ID') and row.get('Ligand ID') else None,
                cv=IdentifierNamespaceCv.GUIDETOPHARMA
            ),
        ),
        annotations=Annotations(
            # Source annotation
            # Interaction properties
            Column('Action', cv=action_cv_mapping),
            Column('Type', cv=type_cv_mapping),
            Column('Endogenous', cv=MoleculeAnnotationsCv.ENDOGENOUS),
            Column('Affinity High', cv=MoleculeAnnotationsCv.AFFINITY_HIGH),
            Column('Affinity Low', cv=MoleculeAnnotationsCv.AFFINITY_LOW),
            Column('Affinity Median', cv=MoleculeAnnotationsCv.AFFINITY_MEDIAN),
            Column('Affinity Units', cv=affinity_units_cv_mapping),
            Column('Primary Target', cv=MoleculeAnnotationsCv.PRIMARY_TARGET),
            Column('PubMed ID', cv=IdentifierNamespaceCv.PUBMED),
            Column('Target Species', cv=species_to_taxid, term_cv=IdentifierNamespaceCv.NCBI_TAX_ID),
        ),
        members=Members(
            # Ligand
            Member(
                entity=Entity(
                    entity_type=Column('Ligand Type', cv=ligand_chemical_type_mapping),
                    identifiers=Identifiers(
                        Column('Ligand ID', cv=IdentifierNamespaceCv.GUIDETOPHARMA),
                    ),
                ),
            ),
            # Target
            Member(
                entity=Entity(
                    entity_type=EntityTypeCv.PROTEIN,
                    identifiers=Identifiers(
                        Column('Target ID', cv=IdentifierNamespaceCv.GUIDETOPHARMA),
                    ),
                ),
            ),
        ),
    )

    # Parse and yield entities
    # Skip the first line (version info)
    next(opener.result)
    reader = csv.DictReader(opener.result)
    for row in reader:
        yield schema(row)
