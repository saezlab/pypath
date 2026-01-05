"""
Parse Guide to Pharmacology data and emit Entity records.

This module converts Guide to Pharmacology ligand-target interaction data
into Entity records using the schema defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator
import csv

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
    ResourceCv,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    Member,
    MembershipBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig


# Processing helpers
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

affinity_units_cv_mapping = {
    'pA2': AffinityUnitCv.PA2,
    'pEC50': AffinityUnitCv.PEC50,
    'pIC50': AffinityUnitCv.PIC50,
    'pKB': AffinityUnitCv.PKB,
    'pKd': AffinityUnitCv.PKD,
    'pKi': AffinityUnitCv.PKI,
}

ligand_chemical_type_mapping = {
    'Synthetic organic': MoleculeSubtypeCv.SYNTHETIC_ORGANIC,
    'Natural product': MoleculeSubtypeCv.NATURAL_PRODUCT,
    'Metabolite': MoleculeSubtypeCv.METABOLITE,
    'Inorganic': MoleculeSubtypeCv.INORGANIC,
    'Peptide': MoleculeSubtypeCv.PEPTIDE,
    'Antibody': MoleculeSubtypeCv.ANTIBODY,
    'Nucleic acid': MoleculeSubtypeCv.NUCLEIC_ACID,
}

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
    'Monkey': '9544',
    'Gorilla': '9593',
    'Ferret': '9669',
    'Honeybee': '7460',
    'spiny starfish': '7609',
    'Escherichia coli': '562',
    'Mycobacterium tuberculosis': '1773',
    'Plasmodium falciparum': '5833',
    'Plasmodium falciparum 3D7': '36329',
    'Plasmodium falciparum 7G8': '57266',
    'Plasmodium falciparum A2': '5833',
    'Plasmodium falciparum D6': '478860',
    'Plasmodium falciparum Dd2': '57267',
    'Plasmodium falciparum FC27': '5833',
    'Plasmodium falciparum FCR3': '5833',
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


f = FieldConfig(
    extract={
        'lower': str.lower,
    },
    map={
        'action': action_cv_mapping,
        'ligand_type': type_cv_mapping,
        'affinity_units': affinity_units_cv_mapping,
        'ligand_chemical_type': ligand_chemical_type_mapping,
        'target_type': target_type_mapping,
        'species_taxid': species_to_taxid,
    },
)


def _iter_guidetopharma_csv(opener, **_kwargs: object):
    if not opener or not opener.result:
        return
    next(opener.result, None)
    reader = csv.DictReader(opener.result)
    yield from reader


config = ResourceConfig(
    id=ResourceCv.GUIDETOPHARMA,
    name='Guide to Pharmacology',
    url='https://www.guidetopharmacology.org/',
    license=LicenseCV.CC_BY_SA_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='37953350',
    description=(
        'The IUPHAR/BPS Guide to PHARMACOLOGY is an expert-curated resource '
        'of ligand-activity-target relationships, providing quantitative '
        'information on drug targets and the prescription medicines and '
        'experimental drugs that act on them. It covers G protein-coupled '
        'receptors, voltage-gated ion channels, ligand-gated ion channels, '
        'nuclear hormone receptors, catalytic receptors, enzymes, and '
        'transporters.'
    ),
)

targets_schema = EntityBuilder(
    entity_type=EntityTypeCv.PROTEIN,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.GUIDETOPHARMA, value=f('Target id')),
        CV(term=IdentifierNamespaceCv.UNIPROT, value=f('Human SwissProt')),
        CV(term=IdentifierNamespaceCv.ENSEMBL, value=f('Human Ensembl Gene')),
        CV(term=IdentifierNamespaceCv.ENTREZ, value=f('Human Entrez Gene')),
        CV(term=IdentifierNamespaceCv.REFSEQ_PROTEIN, value=f('Human protein RefSeq')),
        CV(term=IdentifierNamespaceCv.HGNC, value=f('HGNC id')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('Target name')),
        CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=f('HGNC symbol')),
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('synonyms', delimiter='|')),
        CV(term=IdentifierNamespaceCv.SYSTEMATIC_NAME, value=f('Target systematic name')),
        CV(term=IdentifierNamespaceCv.ABBREVIATED_NAME, value=f('Target abbreviated name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value='9606'),
        CV(term=MoleculeAnnotationsCv.PROTEIN_FAMILY, value=f('Family name')),
        CV(term=f('Type', map='target_type')),
    ),
)

ligands_schema = EntityBuilder(
    entity_type=EntityTypeCv.SMALL_MOLECULE,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.GUIDETOPHARMA, value=f('Ligand ID')),
        CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=f('PubChem SID')),
        CV(term=IdentifierNamespaceCv.PUBCHEM, value=f('PubChem CID')),
        CV(term=IdentifierNamespaceCv.UNIPROT, value=f('UniProt ID')),
        CV(term=IdentifierNamespaceCv.ENSEMBL, value=f('Ensembl ID')),
        CV(term=IdentifierNamespaceCv.CHEMBL, value=f('ChEMBL ID')),
        CV(term=IdentifierNamespaceCv.STANDARD_INCHI_KEY, value=f('InChIKey')),
        CV(term=IdentifierNamespaceCv.SYNONYM, value=f('Synonyms')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('Name')),
        CV(term=IdentifierNamespaceCv.IUPAC_NAME, value=f('IUPAC name')),
        CV(term=IdentifierNamespaceCv.SMILES, value=f('SMILES')),
        CV(term=IdentifierNamespaceCv.STANDARD_INCHI, value=f('InChI')),
        CV(term=IdentifierNamespaceCv.INN, value=f('INN')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=IdentifierNamespaceCv.NCBI_TAX_ID,
            value=f('Species', map='species_taxid'),
        ),
        CV(term=f('Approved', extract='lower', map={'yes': MoleculeAnnotationsCv.APPROVED})),
        CV(term=f('Withdrawn', extract='lower', map={'yes': MoleculeAnnotationsCv.WITHDRAWN})),
        CV(term=f('Labelled', extract='lower', map={'yes': MoleculeAnnotationsCv.LABELLED})),
        CV(term=f('Radioactive', extract='lower', map={'yes': MoleculeAnnotationsCv.RADIOACTIVE})),
        CV(term=f('Antibacterial', extract='lower', map={'yes': MoleculeAnnotationsCv.ANTIBACTERIAL})),
        CV(term=f('Type', map='ligand_chemical_type')),
    ),
)

affinity_units_source = f('Affinity Units', map='affinity_units')

interactions_schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    identifiers=IdentifiersBuilder(
        CV(
            term=IdentifierNamespaceCv.GUIDETOPHARMA,
            value=f(
                lambda row: f"{row['Target ID']}_{row['Ligand ID']}"
                if row.get('Target ID') and row.get('Ligand ID')
                else None
            ),
        ),
    ),
    annotations=AnnotationsBuilder(
        CV(term=f('Action', map='action')),
        CV(term=f('Type', map='ligand_type')),
        CV(term=MoleculeAnnotationsCv.ENDOGENOUS, value=f('Endogenous')),
        CV(
            term=MoleculeAnnotationsCv.AFFINITY_HIGH,
            value=f('Affinity High'),
            unit=affinity_units_source,
        ),
        CV(
            term=MoleculeAnnotationsCv.AFFINITY_LOW,
            value=f('Affinity Low'),
            unit=affinity_units_source,
        ),
        CV(
            term=MoleculeAnnotationsCv.AFFINITY_MEDIAN,
            value=f('Affinity Median'),
            unit=affinity_units_source,
        ),
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('PubMed ID')),
        CV(
            term=IdentifierNamespaceCv.NCBI_TAX_ID,
            value=f('Target Species', map='species_taxid'),
        ),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.SMALL_MOLECULE,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.GUIDETOPHARMA, value=f('Ligand ID')),
                ),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.GUIDETOPHARMA, value=f('Target ID')),
                ),
            ),
        ),
    ),
)

resource = Resource(
    config,
    targets=Dataset(
        download=Download(
            url='http://www.guidetopharmacology.org/DATA/targets_and_families.csv',
            filename='targets_and_families.csv',
            subfolder='guidetopharma',
        ),
        mapper=targets_schema,
        raw_parser=_iter_guidetopharma_csv,
    ),
    ligands=Dataset(
        download=Download(
            url='http://www.guidetopharmacology.org/DATA/ligands.csv',
            filename='ligands.csv',
            subfolder='guidetopharma',
        ),
        mapper=ligands_schema,
        raw_parser=_iter_guidetopharma_csv,
    ),
    interactions=Dataset(
        download=Download(
            url='http://www.guidetopharmacology.org/DATA/interactions.csv',
            filename='interactions.csv',
            subfolder='guidetopharma',
        ),
        mapper=interactions_schema,
        raw_parser=_iter_guidetopharma_csv,
    ),
)
