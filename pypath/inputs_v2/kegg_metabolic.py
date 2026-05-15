""Parse KEGG metabolic reaction network data into reaction entities."""

from __future__ import annotations

from pypath.internals.cv_terms import (
    BiologicalRoleCv,
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
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
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.kegg import _raw


config = ResourceConfig(
    id=ResourceCv.KEGG_METABOLIC,
    name='KEGG Metabolic Reactions (mouse)',
    url='https://www.genome.jp/kegg/reaction/',
    license=LicenseCV.KEGG_ACADEMIC,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='10592173',
    primary_category='metabolism',
    description=(
        'Mouse (mmu) metabolic reaction network from KEGG. '
        'Reactions are linked to their catalysing enzymes via UniProt IDs '
        'and to substrate / product metabolites via KEGG compound and ChEBI '
        'identifiers. Built from KEGG REST bulk endpoints '
        '(list, conv, link) using the dlmachine download infrastructure.'
    ),
)

f = FieldConfig(
    extract={
        'kegg_cpd': r'^(?:cpd:)?(.+)$',
        'chebi':    r'^(?:chebi:)?(\d+)$',
    },
)

reactions_schema = EntityBuilder(
    entity_type=EntityTypeCv.REACTION,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.KEGG_REACTION, value=f('reaction_id')),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.UNIPROT, value=f('uniprot_id')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('taxon_id')),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=BiologicalRoleCv.CATALYST),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.METABOLITE,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('compound_kegg_id', extract='kegg_cpd')),
                    CV(term=IdentifierNamespaceCv.CHEBI,          value=f('compound_chebi',   extract='chebi')),
                    CV(term=IdentifierNamespaceCv.NAME,            value=f('compound_name')),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=BiologicalRoleCv.SUBSTRATE),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.METABOLITE,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('compound_kegg_id', extract='kegg_cpd')),
                    CV(term=IdentifierNamespaceCv.CHEBI,          value=f('compound_chebi',   extract='chebi')),
                    CV(term=IdentifierNamespaceCv.NAME,            value=f('compound_name')),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=BiologicalRoleCv.PRODUCT),
            ),
        ),
    ),
)

download = Download(
    url={
        'conv_uniprot':  'https://rest.kegg.jp/conv/uniprot/mmu',
        'link_enzyme':   'https://rest.kegg.jp/link/enzyme/mmu',
        'link_reaction': 'https://rest.kegg.jp/link/reaction/enzyme',
        'list_compound': 'https://rest.kegg.jp/list/compound',
        'conv_chebi':    'https://rest.kegg.jp/conv/chebi/compound',
        'list_reaction': 'https://rest.kegg.jp/list/reaction',
    },
    subfolder='kegg_metabolic',
)

resource = Resource(
    config,
    reactions=Dataset(
        download=download,
        mapper=reactions_schema,
        raw_parser=lambda opener, **kwargs: _raw(opener, data_type='reactions', **kwargs),
    ),
)
