"""
Parse MRClinksDB metabolite-receptor interaction data and emit Entity records.

MRClinksDB provides curated metabolite–ligand-receptor (L-R) interactions
with compound classification (ClassyFire taxonomy), receptor UniProt IDs,
and supporting literature references.
"""

from __future__ import annotations

from pypath.internals.cv_terms import (
    CurationCv,
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    MoleculeAnnotationsCv,
    ParticipantMetadataCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    Member,
    MembershipBuilder,
    MembersFromList,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.mrclinksdb import iter_mrclinksdb_interactions


# =============================================================================
# Resource Configuration
# =============================================================================

config = ResourceConfig(
    id=ResourceCv.MRCLINKSDB,
    name='MRClinksDB',
    url='https://www.cellknowledge.com.cn/mrclinkdb/',
    license=LicenseCV.ACADEMIC_FREE,
    update_category=UpdateCategoryCV.IRREGULAR,
    primary_category='interactions',
    description=(
        'MRClinksDB is a curated database of metabolite–receptor interactions, '
        'covering ligand-receptor pairs with compound classification based on '
        'the ClassyFire taxonomy, receptor UniProt identifiers, and PubMed '
        'literature support.'
    ),
)


# =============================================================================
# Field Config
# =============================================================================

f = FieldConfig(
    transform={
        # PubChem field is "CID:12345;SID:67890" — extract just the numeric CID
        'pubchem_cid': lambda v: v.split(';')[0].split(':')[-1].strip(),
    },
)


# =============================================================================
# Schema
# =============================================================================

interactions_schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.MRCLINKSDB, value=f('mrid')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('pmid', delimiter=';')),
        CV(term=CurationCv.COMMENT, value=f('other_db', delimiter=';')),
    ),
    membership=MembershipBuilder(
        # Metabolite (small molecule)
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.SMALL_MOLECULE,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.HMDB, value=f('hmdb_id')),
                    CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=f('pubchem_cid_sid', transform='pubchem_cid')),
                    CV(term=IdentifierNamespaceCv.NAME, value=f('metabolite_name')),
                    CV(term=IdentifierNamespaceCv.SMILES, value=f('canonical_smiles')),
                    CV(term=IdentifierNamespaceCv.MOLECULAR_FORMULA, value=f('molecular_formula')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=MoleculeAnnotationsCv.COMPOUND_KINGDOM, value=f('kingdom')),
                    CV(term=MoleculeAnnotationsCv.COMPOUND_SUPERCLASS, value=f('super_class')),
                    CV(term=MoleculeAnnotationsCv.COMPOUND_CLASS, value=f('class')),
                ),
            ),
        ),
        # Receptor — represented as a complex to accommodate heterodimers;
        # for single-protein receptors the complex has exactly one member.
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.COMPLEX,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.NAME, value=f('receptor_symbol')),
                ),
                membership=MembershipBuilder(
                    MembersFromList(
                        entity_type=EntityTypeCv.PROTEIN,
                        identifiers=IdentifiersBuilder(
                            CV(term=IdentifierNamespaceCv.UNIPROT, value=f('receptor_uniprot_id', delimiter='_')),
                            CV(term=IdentifierNamespaceCv.ENTREZ, value=f('receptor_gene_id')),
                            CV(term=IdentifierNamespaceCv.NAME, value=f('protein_name')),
                        ),
                    ),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=ParticipantMetadataCv.TARGET),
            ),
        ),
    ),
)


# =============================================================================
# Resource Definition
# =============================================================================

resource = Resource(
    config,
    human_interactions=Dataset(
        download=Download(
            url='https://www.cellknowledge.com.cn/mrclinkdb/download/Homo%20sapiens%20metabolite%20L-R%20interaction.txt',
            filename='mrclinksdb_human_interactions.txt',
            subfolder='mrclinksdb',
            ext='txt',
        ),
        mapper=interactions_schema,
        raw_parser=iter_mrclinksdb_interactions,
    ),
    mouse_interactions=Dataset(
        download=Download(
            url='https://www.cellknowledge.com.cn/mrclinkdb/download/Mus%20musculus%20metabolite%20L-R%20interaction.txt',
            filename='mrclinksdb_mouse_interactions.txt',
            subfolder='mrclinksdb',
            ext='txt',
        ),
        mapper=interactions_schema,
        raw_parser=iter_mrclinksdb_interactions,
    ),
)
