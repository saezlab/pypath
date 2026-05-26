"""Human-GEM genome-scale metabolic model — inputs_v2 module.

Provides a :class:`~pypath.inputs_v2.base.Resource` with four datasets parsed
from the Human-GEM YAML (``model/Human-GEM.yml``) downloaded directly from
the SysBioChalmers GitHub repository.

Datasets:
    metabolites: SMALL_MOLECULE entities, one per unique base MAM ID
        (compartment suffix stripped and deduplicated). Metabolite subtype,
        molecular formula and charge are stored as annotations.
    reactions: REACTION entities with SMALL_MOLECULE sub-members carrying
        stoichiometry and reactant/product role annotations.  Direction is
        derived from the ``lower_bound`` / ``upper_bound`` bounds.
    catalysis: INTERACTION entities linking an Ensembl gene (PROTEIN) or
        enzyme complex (COMPLEX) as controller to the catalysed REACTION as
        controlled.  One row is emitted per unique (enzyme, reaction) pair
        derived from the ``gene_reaction_rule`` Boolean expression.  Gene
        symbols are included on protein controller entities.
    enzyme_complexes: COMPLEX entities for reactions whose
        ``gene_reaction_rule`` contains AND groups (multi-subunit enzymes).
        Each COMPLEX has PROTEIN sub-members identified by Ensembl ID.
        Complexes are deduplicated across reactions by their sorted subunit
        composition.

Note:
    All four datasets share a single :class:`~pypath.inputs_v2.base.Download`
    so the YAML is fetched only once.

Intentional differences from the legacy ``pypath.inputs.metatlas`` module:

- Only Human-GEM is supported; the multi-model API and git-infrastructure
  helpers are not ported.
- No HMDB / ChEBI / KEGG cross-references: these are absent from the YAML and
  would require additional TSV downloads.
- Transport-network analysis (``metatlas_gem_transport_network``) is not
  ported.
- Gene identifiers are Ensembl ENSG IDs (as in Human-GEM), not Entrez.
- Genes are not a separate dataset; gene symbols are carried on the protein
  controller entities in the catalysis dataset instead.
"""

from __future__ import annotations

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
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    Member,
    MembersFromList,
    MembershipBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.metatlas import _raw


config = ResourceConfig(
    id = ResourceCv.HUMAN_GEM,
    name = 'Human-GEM',
    url = 'https://github.com/SysBioChalmers/Human-GEM',
    license = LicenseCV.CC_BY_4_0,
    update_category = UpdateCategoryCV.IRREGULAR,
    pubmed = '34619098',
    primary_category = 'metabolic_network',
    description = (
        'Human-GEM is a genome-scale metabolic model of human metabolism '
        'maintained by the SysBioChalmers group.'
    ),
)

f = FieldConfig(
    delimiter = ';',
    map = {
        'entity_type': lambda value: (
            EntityTypeCv.COMPLEX if value == 'complex'
            else EntityTypeCv.PROTEIN
        ),
        'stoich_id': lambda value: value.split(':')[0] if value else None,
        'stoich_comp': lambda value: value.split(':')[1] if value else None,
        'stoich_val': lambda value: value.split(':')[2] if value else None,
        'complex_name': lambda value: (
            'complex:' + '-'.join(value.split('||')) if value else ''
        ),
    },
)


download = Download(
    url = 'https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/main/model/Human-GEM.yml',
    filename = 'Human-GEM.yml',
    subfolder = 'metatlas',
    large = True,
    ext = 'yml',
    default_mode = 'r',
)


# ── metabolites ──────────────────────────────────────────────────────────────

metabolites_schema = EntityBuilder(
    entity_type = EntityTypeCv.SMALL_MOLECULE,
    identifiers = IdentifiersBuilder(
        CV(term = IdentifierNamespaceCv.NAME, value = f('name')),
        CV(term = IdentifierNamespaceCv.HUMAN_GEM_METABOLITE, value = f('human_gem_metabolite_id')),
    ),
    annotations = AnnotationsBuilder(
        CV(term = MoleculeAnnotationsCv.MOLECULE_SUBTYPE, value = MoleculeSubtypeCv.METABOLITE),
        CV(term = IdentifierNamespaceCv.MOLECULAR_FORMULA, value = f('formula')),
        CV(term = MoleculeAnnotationsCv.MOLECULAR_CHARGE, value = f('charge')),
    ),
)


# ── reactions ────────────────────────────────────────────────────────────────

reactions_schema = EntityBuilder(
    entity_type = EntityTypeCv.REACTION,
    identifiers = IdentifiersBuilder(
        CV(term = IdentifierNamespaceCv.NAME, value = f('name')),
        CV(term = IdentifierNamespaceCv.HUMAN_GEM_REACTION, value = f('human_gem_reaction_id')),
        CV(term = IdentifierNamespaceCv.EC, value = f('eccodes')),
    ),
    annotations = AnnotationsBuilder(
        CV(term = InteractionMetadataCv.CONVERSION_DIRECTION, value = f('direction')),
        CV(term = MoleculeAnnotationsCv.PATHWAY_PARTICIPATION, value = f('subsystem')),
    ),
    membership = MembershipBuilder(
        MembersFromList(
            entity_type = EntityTypeCv.SMALL_MOLECULE,
            identifiers = IdentifiersBuilder(
                CV(
                    term = IdentifierNamespaceCv.HUMAN_GEM_METABOLITE,
                    value = f('reactants', delimiter = '||', map = 'stoich_id'),
                ),
            ),
            annotations = AnnotationsBuilder(
                CV(term = BiologicalRoleCv.REACTANT),
                CV(
                    term = ParticipantMetadataCv.STOICHIOMETRY,
                    value = f('reactants', delimiter = '||', map = 'stoich_val'),
                ),
                CV(
                    term = MoleculeAnnotationsCv.SUBCELLULAR_LOCATION,
                    value = f('reactants', delimiter = '||', map = 'stoich_comp'),
                ),
            ),
            entity_annotations = AnnotationsBuilder(
                CV(term = MoleculeAnnotationsCv.MOLECULE_SUBTYPE, value = MoleculeSubtypeCv.METABOLITE),
            ),
        ),
        MembersFromList(
            entity_type = EntityTypeCv.SMALL_MOLECULE,
            identifiers = IdentifiersBuilder(
                CV(
                    term = IdentifierNamespaceCv.HUMAN_GEM_METABOLITE,
                    value = f('products', delimiter = '||', map = 'stoich_id'),
                ),
            ),
            annotations = AnnotationsBuilder(
                CV(term = BiologicalRoleCv.PRODUCT),
                CV(
                    term = ParticipantMetadataCv.STOICHIOMETRY,
                    value = f('products', delimiter = '||', map = 'stoich_val'),
                ),
                CV(
                    term = MoleculeAnnotationsCv.SUBCELLULAR_LOCATION,
                    value = f('products', delimiter = '||', map = 'stoich_comp'),
                ),
            ),
            entity_annotations = AnnotationsBuilder(
                CV(term = MoleculeAnnotationsCv.MOLECULE_SUBTYPE, value = MoleculeSubtypeCv.METABOLITE),
            ),
        ),
    ),
)


# ── transport_reactions ──────────────────────────────────────────────────────

transport_reactions_schema = EntityBuilder(
    entity_type = EntityTypeCv.TRANSPORT,
    identifiers = IdentifiersBuilder(
        CV(term = IdentifierNamespaceCv.NAME, value = f('name')),
        CV(term = IdentifierNamespaceCv.HUMAN_GEM_REACTION, value = f('human_gem_reaction_id')),
        CV(term = IdentifierNamespaceCv.EC, value = f('eccodes')),
    ),
    annotations = AnnotationsBuilder(
        CV(term = InteractionMetadataCv.CONVERSION_DIRECTION, value = f('direction')),
        CV(term = MoleculeAnnotationsCv.PATHWAY_PARTICIPATION, value = f('subsystem')),
    ),
    membership = MembershipBuilder(
        MembersFromList(
            entity_type = EntityTypeCv.SMALL_MOLECULE,
            identifiers = IdentifiersBuilder(
                CV(
                    term = IdentifierNamespaceCv.HUMAN_GEM_METABOLITE,
                    value = f('reactants', delimiter = '||', map = 'stoich_id'),
                ),
            ),
            annotations = AnnotationsBuilder(
                CV(term = BiologicalRoleCv.REACTANT),
                CV(
                    term = ParticipantMetadataCv.STOICHIOMETRY,
                    value = f('reactants', delimiter = '||', map = 'stoich_val'),
                ),
                CV(
                    term = MoleculeAnnotationsCv.SUBCELLULAR_LOCATION,
                    value = f('reactants', delimiter = '||', map = 'stoich_comp'),
                ),
            ),
            entity_annotations = AnnotationsBuilder(
                CV(term = MoleculeAnnotationsCv.MOLECULE_SUBTYPE, value = MoleculeSubtypeCv.METABOLITE),
            ),
        ),
        MembersFromList(
            entity_type = EntityTypeCv.SMALL_MOLECULE,
            identifiers = IdentifiersBuilder(
                CV(
                    term = IdentifierNamespaceCv.HUMAN_GEM_METABOLITE,
                    value = f('products', delimiter = '||', map = 'stoich_id'),
                ),
            ),
            annotations = AnnotationsBuilder(
                CV(term = BiologicalRoleCv.PRODUCT),
                CV(
                    term = ParticipantMetadataCv.STOICHIOMETRY,
                    value = f('products', delimiter = '||', map = 'stoich_val'),
                ),
                CV(
                    term = MoleculeAnnotationsCv.SUBCELLULAR_LOCATION,
                    value = f('products', delimiter = '||', map = 'stoich_comp'),
                ),
            ),
            entity_annotations = AnnotationsBuilder(
                CV(term = MoleculeAnnotationsCv.MOLECULE_SUBTYPE, value = MoleculeSubtypeCv.METABOLITE),
            ),
        ),
    ),
)


# ── catalysis ────────────────────────────────────────────────────────────────

catalysis_schema = EntityBuilder(
    entity_type = EntityTypeCv.CATALYSIS,
    identifiers = IdentifiersBuilder(
        CV(term = IdentifierNamespaceCv.NAME, value = f('reaction_name')),
        CV(term = IdentifierNamespaceCv.HUMAN_GEM_REACTION, value = f('reaction_id')),
    ),
    annotations = AnnotationsBuilder(
        CV(term = MoleculeAnnotationsCv.PATHWAY_PARTICIPATION, value = f('subsystem')),
    ),
    membership = MembershipBuilder(
        Member(
            entity = EntityBuilder(
                entity_type = f('enzyme_type', map = 'entity_type'),
                identifiers = IdentifiersBuilder(
                    CV(term = IdentifierNamespaceCv.ENSEMBL, value = f('enzyme_ensembl')),
                    CV(term = IdentifierNamespaceCv.GENE_NAME_PRIMARY, value = f('enzyme_name')),
                    CV(
                        term = IdentifierNamespaceCv.NAME,
                        value = f('complex_subunits', map = 'complex_name'),
                    ),
                ),
            ),
            annotations = AnnotationsBuilder(CV(term = BiologicalRoleCv.CONTROLLER)),
        ),
        Member(
            entity = EntityBuilder(
                entity_type = EntityTypeCv.REACTION,
                identifiers = IdentifiersBuilder(
                    CV(term = IdentifierNamespaceCv.HUMAN_GEM_REACTION, value = f('reaction_id')),
                ),
            ),
            annotations = AnnotationsBuilder(CV(term = BiologicalRoleCv.CONTROLLED)),
        ),
    ),
)


# ── enzyme_complexes ─────────────────────────────────────────────────────────

enzyme_complexes_schema = EntityBuilder(
    entity_type = EntityTypeCv.COMPLEX,
    identifiers = IdentifiersBuilder(
        CV(term = IdentifierNamespaceCv.NAME, value = f('complex_subunits', map = 'complex_name')),
    ),
    membership = MembershipBuilder(
        MembersFromList(
            entity_type = EntityTypeCv.PROTEIN,
            identifiers = IdentifiersBuilder(
                CV(term = IdentifierNamespaceCv.ENSEMBL, value = f('complex_subunits', delimiter = '||')),
            ),
        ),
    ),
)


# ── resource ─────────────────────────────────────────────────────────────────

resource = Resource(
    config,
    metabolites = Dataset(
        download = download,
        mapper = metabolites_schema,
        raw_parser = lambda opener, **kwargs: _raw(opener, data_type = 'metabolites', **kwargs),
    ),
    reactions = Dataset(
        download = download,
        mapper = reactions_schema,
        raw_parser = lambda opener, **kwargs: _raw(opener, data_type = 'reactions', **kwargs),
    ),
    transport_reactions = Dataset(
        download = download,
        mapper = transport_reactions_schema,
        raw_parser = lambda opener, **kwargs: _raw(opener, data_type = 'transport_reactions', **kwargs),
    ),
    metabolic_reactions = Dataset(
        download = download,
        mapper = reactions_schema,
        raw_parser = lambda opener, **kwargs: _raw(opener, data_type = 'metabolic_reactions', **kwargs),
    ),
    catalysis = Dataset(
        download = download,
        mapper = catalysis_schema,
        raw_parser = lambda opener, **kwargs: _raw(opener, data_type = 'catalysis', **kwargs),
    ),
    enzyme_complexes = Dataset(
        download = download,
        mapper = enzyme_complexes_schema,
        raw_parser = lambda opener, **kwargs: _raw(opener, data_type = 'enzyme_complexes', **kwargs),
    ),
)
