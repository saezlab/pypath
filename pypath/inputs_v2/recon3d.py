"""Recon3D genome-scale metabolic model — inputs_v2 module.

Provides a :class:`~pypath.inputs_v2.base.Resource` with five datasets
parsed from the Recon3D BiGG JSON (``Recon3D.json``):

Datasets:
    metabolites: CHEMICAL entities, one per unique base ID (compartment
        suffix stripped and deduplicated).  Cross-references to HMDB, ChEBI,
        KEGG Compound, and MetaNetX are included as identifiers; metabolite
        subtype, molecular formula and charge are stored as annotations.
    reactions: REACTION entities with CHEMICAL sub-members carrying
        stoichiometry and reactant/product role annotations.  Direction is
        derived from the ``lower_bound`` field.
    catalysis: INTERACTION entities linking an enzyme (PROTEIN or COMPLEX)
        as controller to the catalysed REACTION as controlled.  One row is
        emitted per unique (enzyme, reaction) pair derived from the
        ``gene_reaction_rule`` Boolean expression.
    enzyme_complexes: COMPLEX entities for reactions whose
        ``gene_reaction_rule`` contains AND groups (multi-subunit enzymes).
        Each COMPLEX has PROTEIN sub-members identified by Entrez ID.
        Complexes are deduplicated across reactions by their sorted subunit
        composition.
    genes: PROTEIN entities, one per unique Entrez ID.  Recon3D-internal
        isoform suffixes (e.g. ``_AT1``) are stripped and duplicate entries
        collapsed.  The placeholder gene ``'0'`` is discarded.

Note:
    All five datasets share a single :class:`~pypath.inputs_v2.base.Download`
    so the JSON is fetched only once.  The raw parser populates a
    module-level cache on first access and serves subsequent datasets from
    memory.
"""

from __future__ import annotations

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.recon3d import _raw
from pypath.internals.cv_terms import (
    BiologicalRoleCv,
    EntityTypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    LicenseCV,
    MoleculeAnnotationsCv,
    MoleculeSubtypeCv,
    ParticipantMetadataCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    Member,
    MembersFromList,
    MembershipBuilder,
)

config = ResourceConfig(
    id=ResourceCv.RECON3D,
    name='Recon3D',
    url='https://www.vmh.life/',
    license=LicenseCV.BIGG,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='29457794',
    primary_category='metabolic_network',
    description=(
        'Recon3D is a comprehensive human genome-scale metabolic reconstruction '
        'that accounts for three-dimensional molecular structures.'
    ),
)

f = FieldConfig(
    delimiter=';',
    map={
        'entity_type': lambda value: (
            EntityTypeCv.COMPLEX if value == 'complex'
            else EntityTypeCv.PROTEIN
        ),
        'stoich_id': lambda value: value.split(':')[0] if value else None,
        'stoich_comp': lambda value: value.split(':')[1] if value else None,
        'stoich_val': lambda value: value.split(':')[2] if value else None,
        'split_member_values': (
            lambda value: value.split(';;') if value else None
        ),
        'complex_name': lambda value: (
            'complex:' + '-'.join(value.split('||')) if value else ''
        ),
    },
)


download = Download(
    url='http://bigg.ucsd.edu/static/models/Recon3D.json',
    filename='Recon3D.json',
    subfolder='recon3d',
    large=True,
    ext='json',
    default_mode='r',
)

HUMAN_TAXON_ID = '9606'


# ── metabolites ──────────────────────────────────────────────────────────────

metabolites_schema = EntityBuilder(
    entity_type=EntityTypeCv.CHEMICAL,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('name')),
        CV(term=IdentifierNamespaceCv.BIGG_METABOLITE, value=f('bigg_metabolite_id')),
        CV(term=IdentifierNamespaceCv.HMDB, value=f('hmdb')),
        CV(term=IdentifierNamespaceCv.CHEBI, value=f('chebi')),
        CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('kegg_compound')),
        CV(term=IdentifierNamespaceCv.METANETX, value=f('metanetx')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.MOLECULE_SUBTYPE, value=MoleculeSubtypeCv.METABOLITE),
        CV(term=IdentifierNamespaceCv.MOLECULAR_FORMULA, value=f('formula')),
        CV(term=MoleculeAnnotationsCv.MOLECULAR_CHARGE, value=f('charge')),
    ),
)


# ── reactions ────────────────────────────────────────────────────────────────

reactions_schema = EntityBuilder(
    entity_type=EntityTypeCv.REACTION,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('name')),
        CV(term=IdentifierNamespaceCv.BIGG_REACTION, value=f('bigg_reaction_id')),
        CV(term=IdentifierNamespaceCv.EC, value=f('ec')),
        CV(term=IdentifierNamespaceCv.METANETX_REACTION, value=f('metanetx_reaction')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=InteractionMetadataCv.CONVERSION_DIRECTION, value=f('direction')),
        CV(term=MoleculeAnnotationsCv.PATHWAY_PARTICIPATION, value=f('subsystem')),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=EntityTypeCv.CHEMICAL,
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.BIGG_METABOLITE,
                   value=f('reactants', delimiter='||', map='stoich_id', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.NAME,
                   value=f('reactant_name', delimiter='||', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.HMDB,
                   value=f('reactant_hmdb', delimiter='||', map='split_member_values', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.CHEBI,
                   value=f('reactant_chebi', delimiter='||', map='split_member_values', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.KEGG_COMPOUND,
                   value=f('reactant_kegg_compound', delimiter='||', map='split_member_values', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.METANETX,
                   value=f('reactant_metanetx', delimiter='||', map='split_member_values', preserve_indices=True)),
            ),
            annotations=AnnotationsBuilder(
                CV(term=BiologicalRoleCv.REACTANT),
                CV(term=ParticipantMetadataCv.STOICHIOMETRY,
                   value=f('reactants', delimiter='||', map='stoich_val', preserve_indices=True)),
                CV(term=MoleculeAnnotationsCv.SUBCELLULAR_LOCATION,
                   value=f('reactants', delimiter='||', map='stoich_comp', preserve_indices=True)),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(term=MoleculeAnnotationsCv.MOLECULE_SUBTYPE, value=MoleculeSubtypeCv.METABOLITE),
                CV(term=IdentifierNamespaceCv.MOLECULAR_FORMULA,
                   value=f('reactant_formula', delimiter='||', preserve_indices=True)),
                CV(term=MoleculeAnnotationsCv.MOLECULAR_CHARGE,
                   value=f('reactant_charge', delimiter='||', preserve_indices=True)),
            ),
        ),
        MembersFromList(
            entity_type=EntityTypeCv.CHEMICAL,
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.BIGG_METABOLITE,
                   value=f('products', delimiter='||', map='stoich_id', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.NAME,
                   value=f('product_name', delimiter='||', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.HMDB,
                   value=f('product_hmdb', delimiter='||', map='split_member_values', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.CHEBI,
                   value=f('product_chebi', delimiter='||', map='split_member_values', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.KEGG_COMPOUND,
                   value=f('product_kegg_compound', delimiter='||', map='split_member_values', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.METANETX,
                   value=f('product_metanetx', delimiter='||', map='split_member_values', preserve_indices=True)),
            ),
            annotations=AnnotationsBuilder(
                CV(term=BiologicalRoleCv.PRODUCT),
                CV(term=ParticipantMetadataCv.STOICHIOMETRY,
                   value=f('products', delimiter='||', map='stoich_val', preserve_indices=True)),
                CV(term=MoleculeAnnotationsCv.SUBCELLULAR_LOCATION,
                   value=f('products', delimiter='||', map='stoich_comp', preserve_indices=True)),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(term=MoleculeAnnotationsCv.MOLECULE_SUBTYPE, value=MoleculeSubtypeCv.METABOLITE),
                CV(term=IdentifierNamespaceCv.MOLECULAR_FORMULA,
                   value=f('product_formula', delimiter='||', preserve_indices=True)),
                CV(term=MoleculeAnnotationsCv.MOLECULAR_CHARGE,
                   value=f('product_charge', delimiter='||', preserve_indices=True)),
            ),
        ),
    ),
)


# ── transport_reactions ──────────────────────────────────────────────────────

transport_reactions_schema = EntityBuilder(
    entity_type=EntityTypeCv.TRANSPORT,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('name')),
        CV(term=IdentifierNamespaceCv.BIGG_REACTION, value=f('bigg_reaction_id')),
        CV(term=IdentifierNamespaceCv.EC, value=f('ec')),
        CV(term=IdentifierNamespaceCv.METANETX_REACTION, value=f('metanetx_reaction')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=InteractionMetadataCv.CONVERSION_DIRECTION, value=f('direction')),
        CV(term=MoleculeAnnotationsCv.PATHWAY_PARTICIPATION, value=f('subsystem')),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=EntityTypeCv.CHEMICAL,
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.BIGG_METABOLITE,
                   value=f('reactants', delimiter='||', map='stoich_id', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.NAME,
                   value=f('reactant_name', delimiter='||', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.HMDB,
                   value=f('reactant_hmdb', delimiter='||', map='split_member_values', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.CHEBI,
                   value=f('reactant_chebi', delimiter='||', map='split_member_values', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.KEGG_COMPOUND,
                   value=f('reactant_kegg_compound', delimiter='||', map='split_member_values', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.METANETX,
                   value=f('reactant_metanetx', delimiter='||', map='split_member_values', preserve_indices=True)),
            ),
            annotations=AnnotationsBuilder(
                CV(term=BiologicalRoleCv.REACTANT),
                CV(term=ParticipantMetadataCv.STOICHIOMETRY,
                   value=f('reactants', delimiter='||', map='stoich_val', preserve_indices=True)),
                CV(term=MoleculeAnnotationsCv.SUBCELLULAR_LOCATION,
                   value=f('reactants', delimiter='||', map='stoich_comp', preserve_indices=True)),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(term=MoleculeAnnotationsCv.MOLECULE_SUBTYPE, value=MoleculeSubtypeCv.METABOLITE),
                CV(term=IdentifierNamespaceCv.MOLECULAR_FORMULA,
                   value=f('reactant_formula', delimiter='||', preserve_indices=True)),
                CV(term=MoleculeAnnotationsCv.MOLECULAR_CHARGE,
                   value=f('reactant_charge', delimiter='||', preserve_indices=True)),
            ),
        ),
        MembersFromList(
            entity_type=EntityTypeCv.CHEMICAL,
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.BIGG_METABOLITE,
                   value=f('products', delimiter='||', map='stoich_id', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.NAME,
                   value=f('product_name', delimiter='||', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.HMDB,
                   value=f('product_hmdb', delimiter='||', map='split_member_values', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.CHEBI,
                   value=f('product_chebi', delimiter='||', map='split_member_values', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.KEGG_COMPOUND,
                   value=f('product_kegg_compound', delimiter='||', map='split_member_values', preserve_indices=True)),
                CV(term=IdentifierNamespaceCv.METANETX,
                   value=f('product_metanetx', delimiter='||', map='split_member_values', preserve_indices=True)),
            ),
            annotations=AnnotationsBuilder(
                CV(term=BiologicalRoleCv.PRODUCT),
                CV(term=ParticipantMetadataCv.STOICHIOMETRY,
                   value=f('products', delimiter='||', map='stoich_val', preserve_indices=True)),
                CV(term=MoleculeAnnotationsCv.SUBCELLULAR_LOCATION,
                   value=f('products', delimiter='||', map='stoich_comp', preserve_indices=True)),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(term=MoleculeAnnotationsCv.MOLECULE_SUBTYPE, value=MoleculeSubtypeCv.METABOLITE),
                CV(term=IdentifierNamespaceCv.MOLECULAR_FORMULA,
                   value=f('product_formula', delimiter='||', preserve_indices=True)),
                CV(term=MoleculeAnnotationsCv.MOLECULAR_CHARGE,
                   value=f('product_charge', delimiter='||', preserve_indices=True)),
            ),
        ),
    ),
)


# ── catalysis ────────────────────────────────────────────────────────────────

catalysis_schema = EntityBuilder(
    entity_type=EntityTypeCv.CATALYSIS,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('reaction_name')),
        CV(term=IdentifierNamespaceCv.BIGG_REACTION, value=f('reaction_bigg_id')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.PATHWAY_PARTICIPATION, value=f('subsystem')),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=f('enzyme_type', map='entity_type'),
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.ENTREZ, value=f('enzyme_entrez')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=HUMAN_TAXON_ID),
                ),
            ),
            annotations=AnnotationsBuilder(CV(term=BiologicalRoleCv.CONTROLLER)),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.REACTION,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.BIGG_REACTION, value=f('reaction_bigg_id')),
                ),
            ),
            annotations=AnnotationsBuilder(CV(term=BiologicalRoleCv.CONTROLLED)),
        ),
    ),
)


# ── enzyme_complexes ─────────────────────────────────────────────────────────

enzyme_complexes_schema = EntityBuilder(
    entity_type=EntityTypeCv.COMPLEX,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('complex_subunits', map='complex_name')),
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=EntityTypeCv.PROTEIN,
            identifiers=IdentifiersBuilder(
                CV(term=IdentifierNamespaceCv.ENTREZ, value=f('complex_subunits', delimiter='||')),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=HUMAN_TAXON_ID),
            ),
        ),
    ),
)


# ── genes ────────────────────────────────────────────────────────────────────

genes_schema = EntityBuilder(
    entity_type=EntityTypeCv.PROTEIN,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.ENTREZ, value=f('entrez_id')),
        CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=f('name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=HUMAN_TAXON_ID),
    ),
)


# ── resource ─────────────────────────────────────────────────────────────────

resource = Resource(
    config,
    metabolites=Dataset(
        download=download,
        mapper=metabolites_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='metabolites', **kwargs),
    ),
    reactions=Dataset(
        download=download,
        mapper=reactions_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='reactions', **kwargs),
    ),
    transport_reactions=Dataset(
        download=download,
        mapper=transport_reactions_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='transport_reactions', **kwargs),
    ),
    metabolic_reactions=Dataset(
        download=download,
        mapper=reactions_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='metabolic_reactions', **kwargs),
    ),
    catalysis=Dataset(
        download=download,
        mapper=catalysis_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='catalysis', **kwargs),
    ),
    enzyme_complexes=Dataset(
        download=download,
        mapper=enzyme_complexes_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='enzyme_complexes', **kwargs),
    ),
    genes=Dataset(
        download=download,
        mapper=genes_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='genes', **kwargs),
    ),
)
