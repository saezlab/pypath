"""
Parse ChEMBL data and emit Entity records.

This module converts ChEMBL ligand-target interaction data
into Entity records using the schema defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from functools import partial
from pathlib import Path

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.chembl import (
    molecules_parser,
    assays_parser,
    documents_parser,
    mechanisms_parser,
)
from pypath.share import cache

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
    MoleculeSubtypeCv,
    ProteinFunctionalClassCv,
    MoleculeAnnotationsCv,
    AssayTypeCv,
    AssayAnnotationsCv,
    CurationCv,
    PharmacologicalActionCv,
    InteractionMetadataCv,
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


VERSION = 36
# Relative path to the SQLite database file within the downloaded archive
DB_REL_PATH = f'chembl_{VERSION}/chembl_{VERSION}_sqlite/chembl_{VERSION}.db'
# Local path to the cached SQLite database file
SQLITE_PATH = Path(cache.get_cachedir()) / f'ChEMBL_SQLite_{VERSION}.sqlite'


def _chembl_url(version: int = VERSION, **_kwargs: object) -> str:
    return ('https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/'
            f'chembl_{version:02d}/chembl_{version:02d}_sqlite.tar.gz')


def _files_needed(version: int = VERSION, **_kwargs: object) -> list[str]:
    return [f'chembl_{version:02d}/chembl_{version:02d}_sqlite/chembl_{version:02d}.db']


config = ResourceConfig(
    id=ResourceCv.CHEMBL,
    name="ChEMBL",
    url="https://www.ebi.ac.uk/chembl/",
    license=LicenseCV.CC_BY_SA_3_0,
    update_category=UpdateCategoryCV.REGULAR,
    description=(
        "ChEMBL is a manually curated chemical database of bioactive molecules "
        "with drug-like properties."
    ),
)

download = Download(
    url=_chembl_url,
    filename=f"chembl_{VERSION}_sqlite.tar.gz",
    subfolder="chembl",
    large=True,
    ext=".tar.gz",
    needed=_files_needed()
)

# Mapping of ChEMBL molecule_type strings to EntityTypeCv
MOLECULE_TYPE_TO_ENTITY_TYPE = {
    'Small molecule': EntityTypeCv.SMALL_MOLECULE,
    'Protein': EntityTypeCv.PROTEIN,
    'Antibody': EntityTypeCv.PROTEIN,
    'Enzyme': EntityTypeCv.PROTEIN,
    'Oligosaccharide': EntityTypeCv.SMALL_MOLECULE,
    'Oligonucleotide': EntityTypeCv.SMALL_MOLECULE,
    'Gene': EntityTypeCv.GENE,
    'Cell': EntityTypeCv.PHYSICAL_ENTITY,
    'Unknown': EntityTypeCv.PHYSICAL_ENTITY,
    'Unclassified': EntityTypeCv.PHYSICAL_ENTITY,
}

# Mapping of ChEMBL molecule_type strings to more specific subclasses
MOLECULE_TYPE_TO_SUBTYPE = {
    'Antibody': MoleculeSubtypeCv.ANTIBODY,
    'Enzyme': ProteinFunctionalClassCv.ENZYME,
    'Small molecule': MoleculeSubtypeCv.SYNTHETIC_ORGANIC,
}

# Mapping of ChEMBL assay_type codes to AssayTypeCv
ASSAY_TYPE_MAP = {
    'B': AssayTypeCv.BINDING,
    'F': AssayTypeCv.FUNCTIONAL,
    'A': AssayTypeCv.ADME,
    'T': AssayTypeCv.TOXICITY,
    'P': AssayTypeCv.PHYSICOCHEMICAL,
}

# Mapping of ChEMBL action_type strings to PharmacologicalActionCv
ACTION_TYPE_MAP = {
    'AGONIST': PharmacologicalActionCv.AGONIST,
    'ANTAGONIST': PharmacologicalActionCv.ANTAGONIST,
    'INHIBITOR': PharmacologicalActionCv.INHIBITION,
    'ACTIVATOR': PharmacologicalActionCv.ACTIVATION,
    'PARTIAL AGONIST': PharmacologicalActionCv.PARTIAL_AGONIST,
    'INVERSE AGONIST': PharmacologicalActionCv.INVERSE_AGONIST,
    'POTENTIATOR': PharmacologicalActionCv.POTENTIATION,
    'POSITIVE ALLosteric MODULATOR': PharmacologicalActionCv.POSITIVE,
    'NEGATIVE ALLosteric MODULATOR': PharmacologicalActionCv.NEGATIVE,
    'BINDING AGENT': PharmacologicalActionCv.BINDING,
}

f = FieldConfig(
    extract={
        'chebi': r'^(?:CHEBI:)?(\d+)$',
    },
    transform={
        'chebi': lambda v: f'CHEBI:{v}' if v else None,
        'bool_to_cv': lambda v, cv: cv if str(v) == '1' else None,
    },
    map={
        'entity_type': MOLECULE_TYPE_TO_ENTITY_TYPE,
        'subtype': MOLECULE_TYPE_TO_SUBTYPE,
        'assay_type': ASSAY_TYPE_MAP,
        'action_type': ACTION_TYPE_MAP,
    },
)

molecules_schema = EntityBuilder(
    entity_type=f('molecule_type', map='entity_type'),
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.CHEMBL, value=f('chembl_id')),
        CV(term=IdentifierNamespaceCv.CHEMBL_INTERNAL_ID, value=f('molregno')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('pref_name')),
        CV(term=IdentifierNamespaceCv.CHEBI, value=f('chebi_id', extract='chebi', transform='chebi')),
        CV(term=IdentifierNamespaceCv.SMILES, value=f('canonical_smiles')),
        CV(term=IdentifierNamespaceCv.STANDARD_INCHI, value=f('standard_inchi')),
        CV(term=IdentifierNamespaceCv.STANDARD_INCHI_KEY, value=f('standard_inchi_key')),
        CV(term=IdentifierNamespaceCv.MOLECULAR_FORMULA, value=f('full_molformula')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.CLINICAL_PHASE, value=f('max_phase')),
        CV(term=f('molecule_type', map='subtype')),
        CV(term=f('therapeutic_flag', transform=lambda v: f.transform['bool_to_cv'](v, MoleculeAnnotationsCv.APPROVED))),
        CV(term=f('withdrawn_flag', transform=lambda v: f.transform['bool_to_cv'](v, MoleculeAnnotationsCv.WITHDRAWN))),
        CV(term=f('polymer_flag', transform=lambda v: f.transform['bool_to_cv'](v, MoleculeAnnotationsCv.POLYMER))),
        CV(term=f('inorganic_flag', transform=lambda v: f.transform['bool_to_cv'](v, MoleculeAnnotationsCv.INORGANIC))),
        CV(term=f('natural_product', transform=lambda v: f.transform['bool_to_cv'](v, MoleculeAnnotationsCv.NATURAL_PRODUCT))),
        CV(term=MoleculeAnnotationsCv.MASS_DALTON, value=f('full_mwt')),
        CV(term=MoleculeAnnotationsCv.MW_MONOISOTOPIC, value=f('mw_monoisotopic')),
        CV(term=MoleculeAnnotationsCv.ALOGP, value=f('alogp')),
        CV(term=MoleculeAnnotationsCv.CX_LOGP, value=f('cx_logp')),
        CV(term=MoleculeAnnotationsCv.CX_LOGD, value=f('cx_logd')),
        CV(term=MoleculeAnnotationsCv.MOLECULAR_SPECIES, value=f('molecular_species')),
    ),
)

assays_schema = EntityBuilder(
    entity_type=EntityTypeCv.ASSAY,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.CHEMBL_ASSAY, value=f('chembl_id')),
        CV(term=IdentifierNamespaceCv.NAME, value=f('description')),
        CV(term=IdentifierNamespaceCv.CHEMBL_TARGET, value=f('target_chembl_id')),
        CV(term=IdentifierNamespaceCv.CHEMBL, value=f('document_chembl_id')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=f('assay_type', map='assay_type')),
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('assay_tax_id')),
        CV(term=AssayAnnotationsCv.CONFIDENCE_SCORE, value=f('confidence_score')),
        CV(term=AssayAnnotationsCv.ASSAY_CATEGORY, value=f('assay_category')),
        CV(term=AssayAnnotationsCv.SUBCELLULAR_FRACTION, value=f('assay_subcellular_fraction')),
        CV(term=AssayAnnotationsCv.TISSUE, value=f('assay_tissue')),
        CV(term=AssayAnnotationsCv.CELL_TYPE, value=f('assay_cell_type')),
        CV(term=AssayAnnotationsCv.DESCRIPTION, value=f('parameters')),
    ),
)

documents_schema = EntityBuilder(
    entity_type=EntityTypeCv.PUBLICATION,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.CHEMBL_DOCUMENT, value=f('chembl_id')),
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('pubmed_id')),
        CV(term=IdentifierNamespaceCv.DOI, value=f('doi')),
        CV(term=IdentifierNamespaceCv.PATENT_NUMBER, value=f('patent_id')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=CurationCv.JOURNAL, value=f('journal')),
        CV(term=CurationCv.YEAR, value=f('year')),
        CV(term=CurationCv.TITLE, value=f('title')),
        CV(term=CurationCv.ABSTRACT, value=f('abstract')),
        CV(term=CurationCv.DOC_TYPE, value=f('doc_type')),
    ),
)

mechanisms_schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.CHEMBL_MECHANISM, value=f('mec_id')),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.SMALL_MOLECULE,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.CHEMBL, value=f('molecule_chembl_id')),
                ),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.CHEMBL_TARGET, value=f('target_chembl_id')),
                ),
            ),
        ),
    ),
    annotations=AnnotationsBuilder(
        CV(term=f('action_type', map='action_type')),
        CV(term=MoleculeAnnotationsCv.DESCRIPTION, value=f('molecular_mechanism')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('mechanism_refs')),
        CV(term=EntityBuilder(
            entity_type=EntityTypeCv.PUBLICATION,
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.CHEMBL, value=f('document_chembl_id')),
                CV(term=IdentifierNamespaceCv.PUBMED, value=f('pubmed_id')),
                CV(term=IdentifierNamespaceCv.DOI, value=f('doi')),
                CV(term=IdentifierNamespaceCv.NAME, value=f('document_title')),
            ),
        )),
    ),
)


resource = Resource(
    config=config,
    molecules=Dataset(
        download=download,
        mapper=molecules_schema,
        raw_parser=partial(
            molecules_parser,
            sqlite_path=SQLITE_PATH,
            db_rel_path=DB_REL_PATH,
        ),
    ),
    assays=Dataset(
        download=download,
        mapper=assays_schema,
        raw_parser=partial(
            assays_parser,
            sqlite_path=SQLITE_PATH,
            db_rel_path=DB_REL_PATH,
        ),
    ),
    documents=Dataset(
        download=download,
        mapper=documents_schema,
        raw_parser=partial(
            documents_parser,
            sqlite_path=SQLITE_PATH,
            db_rel_path=DB_REL_PATH,
        ),
    ),
    mechanisms=Dataset(
        download=download,
        mapper=mechanisms_schema,
        raw_parser=partial(
            mechanisms_parser,
            sqlite_path=SQLITE_PATH,
            db_rel_path=DB_REL_PATH,
        ),
    ),
)
