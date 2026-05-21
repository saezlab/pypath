"""Parse KEGG metabolic reaction network data into reaction entities."""

from __future__ import annotations

_KEGG_ORGANISM_CODES: dict[int, str] = {
    9606:  'hsa',   # Homo sapiens (human)
    10090: 'mmu',   # Mus musculus (mouse)
    10116: 'rno',   # Rattus norvegicus (rat)
    7955:  'dre',   # Danio rerio (zebrafish)
    6239:  'cel',   # Caenorhabditis elegans
    7227:  'dme',   # Drosophila melanogaster
    4932:  'sce',   # Saccharomyces cerevisiae
}


def kegg_organism_code(ncbi_tax_id: int) -> str | None:
    """Return the KEGG organism code for an NCBI taxonomy ID, or None."""
    return _KEGG_ORGANISM_CODES.get(ncbi_tax_id)


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
                entity_type=EntityTypeCv.SMALL_MOLECULE,
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
                entity_type=EntityTypeCv.SMALL_MOLECULE,
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


def make_kegg_resource(organism_code: str = 'mmu') -> Resource:
    """
    Return a KEGG Resource configured for the given KEGG organism code.

    The two organism-specific endpoints (``/conv/uniprot/{org}`` and
    ``/link/enzyme/{org}``) and the cache subfolder are set per organism.
    The four organism-agnostic endpoints (compound list, ChEBI mapping,
    reaction→enzyme link, reaction list) are shared across all organisms
    but stored in the organism-specific subfolder to keep caches isolated.

    Args:
        organism_code: KEGG organism code (e.g. ``'hsa'``, ``'mmu'``, ``'rno'``).
            Use :func:`kegg_organism_code` to convert an NCBI taxonomy ID.

    Returns:
        :class:`~pypath.inputs_v2.base.Resource` with a ``reactions`` dataset
        whose ``raw()`` method yields flat reaction dicts ready for conversion
        to COSMOS :class:`~omnipath_metabo.datasets.cosmos._record.Interaction`
        records.
    """
    org_download = Download(
        url={
            'conv_uniprot':  f'https://rest.kegg.jp/conv/uniprot/{organism_code}',
            'link_enzyme':   f'https://rest.kegg.jp/link/enzyme/{organism_code}',
            'link_reaction': 'https://rest.kegg.jp/link/reaction/enzyme',
            'list_compound': 'https://rest.kegg.jp/list/compound',
            'conv_chebi':    'https://rest.kegg.jp/conv/chebi/compound',
            'list_reaction': 'https://rest.kegg.jp/list/reaction',
        },
        subfolder=f'kegg_metabolic_{organism_code}',
    )
    return Resource(
        config,
        reactions=Dataset(
            download=org_download,
            mapper=reactions_schema,
            raw_parser=lambda opener, **kwargs: _raw(
                opener,
                data_type='reactions',
                organism=organism_code,
                **kwargs,
            ),
        ),
    )
