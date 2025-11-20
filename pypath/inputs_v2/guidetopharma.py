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
from pypath.internals.silver_schema import Entity, Identifier, Annotation
from pypath.internals.cv_terms import (
    EntityTypeCv,
    ProteinFunctionalClassCv,
    MoleculeSubtypeCv,
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
    AnnotationsBuilder,
    Column,
    CV,
    EntityBuilder,
    IdentifiersBuilder,
    Map,
    Member,
    MembershipBuilder,
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

# Mapping for Ligand Type (chemical nature) to MoleculeSubtypeCv
ligand_chemical_type_mapping = {
    'Synthetic organic': MoleculeSubtypeCv.SYNTHETIC_ORGANIC,
    'Natural product': MoleculeSubtypeCv.NATURAL_PRODUCT,
    'Metabolite': MoleculeSubtypeCv.METABOLITE,
    'Inorganic': MoleculeSubtypeCv.INORGANIC,
    'Peptide': MoleculeSubtypeCv.PEPTIDE,
    'Antibody': MoleculeSubtypeCv.ANTIBODY,
    'Nucleic acid': MoleculeSubtypeCv.NUCLEIC_ACID,
}

# Mapping for Target Type to ProteinFunctionalClassCv
target_type_mapping = {
    'gpcr': ProteinFunctionalClassCv.GPCR,
    'lgic': ProteinFunctionalClassCv.LGIC,
    'vgic': ProteinFunctionalClassCv.VGIC,
    'other_ic': ProteinFunctionalClassCv.OTHER_ION_CHANNEL,
    'enzyme': ProteinFunctionalClassCv.ENZYME,
    'catalytic_receptor': ProteinFunctionalClassCv.CATALYTIC_RECEPTOR,
    'nhr': ProteinFunctionalClassCv.NUCLEAR_HORMONE_RECEPTOR,
    'transporter': ProteinFunctionalClassCv.TRANSPORTER,
    'other_protein': ProteinFunctionalClassCv.OTHER_PROTEIN,
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


def resource() -> Generator[Entity]:
    """
    Yield resource metadata as an Entity record.

    Yields:
        Entity record with type CV_TERM containing Guide to Pharmacology metadata.
    """
    yield Entity(
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


def guidetopharma_targets() -> Generator[Entity, None, None]:
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
    schema = EntityBuilder(
        entity_type=EntityTypeCv.PROTEIN,
        identifiers=IdentifiersBuilder(
            CV(term=IdentifierNamespaceCv.GUIDETOPHARMA, value=Column('Target id')),
            # Human identifiers
            CV(term=IdentifierNamespaceCv.UNIPROT, value=Column('Human SwissProt')),
            CV(term=IdentifierNamespaceCv.ENSEMBL, value=Column('Human Ensembl Gene')),
            CV(term=IdentifierNamespaceCv.ENTREZ, value=Column('Human Entrez Gene')),
            CV(term=IdentifierNamespaceCv.REFSEQ_PROTEIN, value=Column('Human protein RefSeq')),
            CV(term=IdentifierNamespaceCv.HGNC, value=Column('HGNC id')),
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
            CV(term=IdentifierNamespaceCv.NAME, value=Column('Target name')),
            CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=Column('HGNC symbol')),
            CV(term=IdentifierNamespaceCv.SYNONYM, value=Column('synonyms', delimiter='|')),
            CV(term=IdentifierNamespaceCv.SYSTEMATIC_NAME, value=Column('Target systematic name')),
            CV(term=IdentifierNamespaceCv.ABBREVIATED_NAME, value=Column('Target abbreviated name')),

        ),
        annotations=AnnotationsBuilder(
            # Source annotation
            CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value="9606"),  # Human
            CV(term=MoleculeAnnotationsCv.PROTEIN_FAMILY, value=Column('Family name')),
            # Protein functional class as annotation
            CV(term=Map(col=Column('Type'), map=target_type_mapping)),
        ),
    )

    # Parse and yield entities
    # Skip the first line (version info)
    next(opener.result)
    reader = csv.DictReader(opener.result)
    for row in reader:
        yield schema(row)


def guidetopharma_ligands() -> Generator[Entity, None, None]:
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
    def boolean_term(column_name: str, term) -> CV:
        column = Column(column_name)
        return CV(
            term=Map(
                col=column,
                extract=[str.lower],
                map={'yes': term},
            ),
        )

    schema = EntityBuilder(
        entity_type=EntityTypeCv.SMALL_MOLECULE,
        identifiers=IdentifiersBuilder(
            CV(term=IdentifierNamespaceCv.GUIDETOPHARMA, value=Column('Ligand ID')),
            CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=Column('PubChem SID')),
            CV(term=IdentifierNamespaceCv.PUBCHEM, value=Column('PubChem CID')),
            CV(term=IdentifierNamespaceCv.UNIPROT, value=Column('UniProt ID')),
            CV(term=IdentifierNamespaceCv.ENSEMBL, value=Column('Ensembl ID')),
            CV(term=IdentifierNamespaceCv.CHEMBL, value=Column('ChEMBL ID')),
            CV(term=IdentifierNamespaceCv.STANDARD_INCHI_KEY, value=Column('InChIKey')),
            CV(term=IdentifierNamespaceCv.SYNONYM, value=Column('Synonyms')),
            CV(term=IdentifierNamespaceCv.NAME, value=Column('Name')),
            CV(term=IdentifierNamespaceCv.IUPAC_NAME, value=Column('IUPAC name')),
            CV(term=IdentifierNamespaceCv.SMILES, value=Column('SMILES')),
            CV(term=IdentifierNamespaceCv.STANDARD_INCHI, value=Column('InChI')),
            CV(term=IdentifierNamespaceCv.INN, value=Column('INN')),
        ),
        annotations=AnnotationsBuilder(
            # Source annotation
            CV(
                term=IdentifierNamespaceCv.NCBI_TAX_ID,
                value=Map(col=Column('Species'), map=species_to_taxid),
            ),
            boolean_term('Approved', MoleculeAnnotationsCv.APPROVED),
            boolean_term('Withdrawn', MoleculeAnnotationsCv.WITHDRAWN),
            boolean_term('Labelled', MoleculeAnnotationsCv.LABELLED),
            boolean_term('Radioactive', MoleculeAnnotationsCv.RADIOACTIVE),
            boolean_term('Antibacterial', MoleculeAnnotationsCv.ANTIBACTERIAL),
            # Molecule subtype as annotation
            CV(term=Map(col=Column('Type'), map=ligand_chemical_type_mapping)),
        ),
    )

    # Parse and yield entities
    # Skip the first line (version info)
    next(opener.result)
    reader = csv.DictReader(opener.result)
    for row in reader:
        yield schema(row)


def guidetopharma_interactions() -> Generator[Entity, None, None]:
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
    affinity_units_source = Map(
        col=Column('Affinity Units'),
        map=affinity_units_cv_mapping,
    )

    schema = EntityBuilder(
        entity_type=EntityTypeCv.INTERACTION,
        identifiers=IdentifiersBuilder(
            CV(
                term=IdentifierNamespaceCv.GUIDETOPHARMA,
                value=Column(
                    lambda row: f"{row['Target ID']}_{row['Ligand ID']}"
                    if row.get('Target ID') and row.get('Ligand ID')
                    else None
                ),
            ),
        ),
        annotations=AnnotationsBuilder(
            # Source annotation
            # Interaction properties
            CV(term=Map(col=Column('Action'), map=action_cv_mapping)),
            CV(term=Map(col=Column('Type'), map=type_cv_mapping)),
            CV(term=MoleculeAnnotationsCv.ENDOGENOUS, value=Column('Endogenous')),
            CV(
                term=MoleculeAnnotationsCv.AFFINITY_HIGH,
                value=Column('Affinity High'),
                unit=affinity_units_source,
            ),
            CV(
                term=MoleculeAnnotationsCv.AFFINITY_LOW,
                value=Column('Affinity Low'),
                unit=affinity_units_source,
            ),
            CV(
                term=MoleculeAnnotationsCv.AFFINITY_MEDIAN,
                value=Column('Affinity Median'),
                unit=affinity_units_source,
            ),
            CV(term=IdentifierNamespaceCv.PUBMED, value=Column('PubMed ID')),
            CV(
                term=IdentifierNamespaceCv.NCBI_TAX_ID,
                value=Map(col=Column('Target Species'), map=species_to_taxid),
            ),
        ),
        membership=MembershipBuilder(
            # Ligand
            Member(
                entity=EntityBuilder(
                    entity_type=EntityTypeCv.SMALL_MOLECULE,
                    identifiers=IdentifiersBuilder(
                        CV(term=IdentifierNamespaceCv.GUIDETOPHARMA, value=Column('Ligand ID')),
                    ),
                ),
            ),
            # Target
            Member(
                entity=EntityBuilder(
                    entity_type=EntityTypeCv.PROTEIN,
                    identifiers=IdentifiersBuilder(
                        CV(term=IdentifierNamespaceCv.GUIDETOPHARMA, value=Column('Target ID')),
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
