import collections

__all__ = [
    "ChemblTarget",
    "ChemblComponent",
    "ChemblAssay",
    "ChemblParam",
    "ChemblIndication",
]


ChemblTarget = collections.namedtuple(
    "ChemblTarget",
    [
        "chembl_id",
        "target_type",
        "preferred_name",
        "ncbi_taxa_id",
        "organism",
        "components",
        "num_components",
    ]
)

ChemblComponent = collections.namedtuple(
    "ChemblComponent",
    [
        "uniprot_accession",
        "component_type",
        "component_description",
        "component_id",
        "component_relationship",
        "component_number",
    ]
)

ChemblAssay = collections.namedtuple(
    "ChemblAssay",
    [
        "assay_chembl_id",
        "assay_type",
        "assay_type_description",
        "assay_category",
        "target_chembl_id",
        "organism",
        "tax_id",
        "tissue",
        "cell_type",
        "subcellular_fraction",
        "parameters",
        "source_id",
        "confidence_score",
        "confidence_description",
        "document_chembl_id",
        "description",
    ]
)

ChemblParam = collections.namedtuple(
    "ChemblParam",
    [
        "standard_type",
        "standard_value",
        "standard_units",
        "standard_relation",
    ]
)

ChemblActivity = collections.namedtuple(
    "ChemblActivity",
    [
        "activity_id",
        "action_type",
        "standard_relation",
        "standard_value",
        "standard_upper_value",
        "standard_units",
        "standard_type",
        "ligand_efficiency",
        "validity_comment",
        "potential_duplicate",
        "pchembl_value",
        "source_id",
        "molecule_chembl_id",
        "target_chembl_id",
        "target_taxa_id",
        "assay_id",
        "document_chembl_id",
    ]
)

ChemblAction = collections.namedtuple(
    "ChemblAction",
    [
        "action_type",
        "description",
        "parent_type",
        "is_stimulation",
        "is_inhibition",
    ]
)

ChemblDocument = collections.namedtuple(
    "ChemblDocument",
    [
        "document_chembl_id",
        'pubmed_id',
        "patent_id",
        "doc_type",
        "journal",
        "year",
        "doi",
    ]
)

ChemblMolecule = collections.namedtuple(
    "ChemblMolecule",
    [
        "molecule_chembl_id",
        "preferred_name",
        "molecule_type",
        "structure_type",
        "chirality",
        "biotherapeutic",
        "inorganic_flag",
        "natural_flag",
        "polymer_flag",
        "helm_notation",
        "molecule_properties",
        "structure",
    ]
)

ChemblMolProps = collections.namedtuple(
    "ChemblMolProps",
    [
        "mol_formula",
        "full_mwt",
        "monoisotopic_mwt",
        "molecular_species",
        "logp",
        "logd",
        "alogp",
    ]
)

ChemblMolStruct = collections.namedtuple(
    "ChemblMolStruct",
    [
        "canonical_smiles",
        "inchi",
        "inchi_key",
    ]
)

ChemblMechanism = collections.namedtuple(
    "ChemblMechanism",
    [
        "action_type",
        "molecule_chembl_id",
        "target_chembl_id",
        "mechanism_id",
        "drug_rec_id",
        "direct_interaction",
        "variant_sequence",
        "molecular_mechanism",
        "mechanism_refs",
    ]
)

ChemblIndication = collections.namedtuple(
    "ChemblIndication",
    [
        "chembl_id",
        "efo_id",
        "efo_term",
        "mesh_id",
        "mesh_term",
        "max_phase",
        "refs",
    ]
)
