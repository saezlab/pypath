"""
Parse MRClinksDB metabolite-receptor interaction data and emit Entity records.

MRClinksDB provides curated metabolite–ligand-receptor (L-R) interactions
with compound classification (ClassyFire taxonomy), receptor UniProt IDs,
and supporting literature references.
"""

from __future__ import annotations

import re

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
        # Protein name fields use UniProt format: "Primary name (Alt 1) (Alt 2) ..."
        # Extract the primary name (everything before the first parenthesis)
        'primary_name': lambda v: v.split('(')[0].strip(),
        # Extract all parenthesised alternative names as a list
        'alt_names': lambda v: re.findall(r'\(([^)]+)\)', v),
    },
)


def _if_single(field: str):
    """Return a callable selector that yields *field* only for single-protein receptors.

    Gates on the absence of ``_`` in ``receptor_uniprot_id``. When the row
    represents a heteromeric complex the selector returns ``None``, causing
    the enclosing ``EntityBuilder`` to produce no identifiers and be silently
    skipped.

    The returned callable can be passed directly to ``f()`` alongside any
    transform, e.g. ``f(_if_single('protein_name'), transform='primary_name')``.
    """
    return lambda row, _f=field: row.get(_f) if '_' not in row.get('receptor_uniprot_id', '') else None


def _if_complex(field: str):
    """Return a callable selector that yields *field* only for heteromeric complex receptors.

    Gates on the presence of ``_`` in ``receptor_uniprot_id``. When the row
    represents a single-protein receptor the selector returns ``None``, causing
    the enclosing ``EntityBuilder`` to produce no identifiers and be silently
    skipped.

    The returned callable can be passed directly to ``f()`` alongside any
    transform, e.g. ``f(_if_complex('protein_name'), transform='primary_name')``.
    """
    return lambda row, _f=field: row.get(_f) if '_' in row.get('receptor_uniprot_id', '') else None


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
                    CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND,
                       value=f('pubchem_cid_sid', transform='pubchem_cid')),
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
            annotations=AnnotationsBuilder(
                CV(term=ParticipantMetadataCv.SOURCE),
            ),
        ),
        # Single-protein receptor
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.UNIPROT,  value=f(_if_single('receptor_uniprot_id'))),
                    CV(term=IdentifierNamespaceCv.ENTREZ,   value=f(_if_single('receptor_gene_id'))),
                    CV(term=IdentifierNamespaceCv.NAME,     value=f(_if_single('protein_name'), transform='primary_name')),
                    CV(term=IdentifierNamespaceCv.SYNONYM,  value=f(_if_single('protein_name'), transform='alt_names')),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=ParticipantMetadataCv.TARGET),
            ),
        ),
        # Heteromeric complex receptor
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.COMPLEX,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.NAME,    value=f(_if_complex('protein_name'), transform='primary_name')),
                    CV(term=IdentifierNamespaceCv.SYNONYM, value=f(_if_complex('protein_name'), transform='alt_names')),
                ),
                membership=MembershipBuilder(
                    MembersFromList(
                        entity_type=EntityTypeCv.PROTEIN,
                        identifiers=IdentifiersBuilder(
                            CV(term=IdentifierNamespaceCv.UNIPROT,
                                value=f('receptor_uniprot_id', delimiter='_')),
                            CV(term=IdentifierNamespaceCv.ENTREZ,
                                value=f('receptor_gene_id',    delimiter='_')),
                            CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY,
                               value=f('receptor_symbol',     delimiter='_')),
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
# Transporter Schemas
# =============================================================================

human_transporters_schema = EntityBuilder(
    entity_type=EntityTypeCv.TRANSPORT,
    identifiers=IdentifiersBuilder(
        CV(
            term=IdentifierNamespaceCv.NAME,
            value=f(lambda row: f"{row.get('hmdb_id', '')}_{row.get('uniprot_id', '')}"),
        ),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.SMALL_MOLECULE,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.HMDB, value=f('hmdb_id')),
                    CV(term=IdentifierNamespaceCv.NAME, value=f('metabolite_name')),
                ),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.UNIPROT,           value=f('uniprot_id')),
                    CV(term=IdentifierNamespaceCv.ENTREZ,            value=f('human_geneid')),
                    CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=f('gene_name')),
                    CV(term=IdentifierNamespaceCv.NAME,              value=f('enzyme_name', transform='primary_name')),
                    CV(term=IdentifierNamespaceCv.SYNONYM,           value=f('enzyme_name', transform='alt_names')),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=ParticipantMetadataCv.TARGET),
            ),
        ),
    ),
)

mouse_transporters_schema = EntityBuilder(
    entity_type=EntityTypeCv.TRANSPORT,
    identifiers=IdentifiersBuilder(
        CV(
            term=IdentifierNamespaceCv.NAME,
            value=f(lambda row: f"{row.get('hmdb_id', '')}_{row.get('mouse_uniprot', '')}"),
        ),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.SMALL_MOLECULE,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.HMDB, value=f('hmdb_id')),
                    CV(term=IdentifierNamespaceCv.NAME, value=f('metabolite_name')),
                ),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.UNIPROT,           value=f('mouse_uniprot')),
                    CV(term=IdentifierNamespaceCv.ENTREZ,            value=f('mouse_geneid')),
                    CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=f('mouse_gene_symbol')),
                    CV(term=IdentifierNamespaceCv.NAME,              value=f('enzyme_name', transform='primary_name')),
                    CV(term=IdentifierNamespaceCv.SYNONYM,           value=f('enzyme_name', transform='alt_names')),
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
    human_transporters=Dataset(
        download=Download(
            url='https://www.cellknowledge.com.cn/mrclinkdb/download/Homo%20sapiens%20transporter%20protein.txt',
            filename='mrclinksdb_human_transporters.txt',
            subfolder='mrclinksdb',
            ext='txt',
        ),
        mapper=human_transporters_schema,
        raw_parser=iter_mrclinksdb_interactions,
    ),
    mouse_transporters=Dataset(
        download=Download(
            url='https://www.cellknowledge.com.cn/mrclinkdb/download/Mus%20musculus%20transporter%20protein.txt',
            filename='mrclinksdb_mouse_transporters.txt',
            subfolder='mrclinksdb',
            ext='txt',
        ),
        mapper=mouse_transporters_schema,
        raw_parser=iter_mrclinksdb_interactions,
    ),
)
