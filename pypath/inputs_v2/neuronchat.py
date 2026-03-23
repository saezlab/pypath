"""
Parse NeuronChat data and emit Entity records.

This module converts NeuronChat interaction data into Entity records using the
declarative schema pattern.
"""

from __future__ import annotations

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
    InteractionMetadataCv,
    ParticipantMetadataCv,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembershipBuilder,
    Member,
    MembersFromList,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.neuronchat import iter_neuronchat

config = ResourceConfig(
    id=ResourceCv.NEURONCHAT,
    name='NeuronChat',
    url='https://github.com/Wei-BioMath/NeuronChat',
    license=LicenseCV.GPL_3_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='36800000', # Zhao et al. 2023 DOI: 10.1038/s41467-023-36800-w
    description=(
        'NeuronChat is a manually curated resource of neural-specific '
        'intercellular molecular interactions, designed for inferring '
        'neuron-neuron communication from single-cell and spatial '
        'transcriptomics data.'
    ),
)

f = FieldConfig()

interactions_schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('interaction_name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('ligand_type')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('interaction_type')),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.COMPLEX,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.NAME, value=f('interaction_name', transform=lambda x: f"{x}_source")),
                ),
                membership=MembershipBuilder(
                    MembersFromList(
                        entity_type=EntityTypeCv.PROTEIN,
                        identifiers=IdentifiersBuilder(
                            CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=f('lig_contributor')),
                        ),
                    )
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=ParticipantMetadataCv.SOURCE),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.COMPLEX,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.NAME, value=f('interaction_name', transform=lambda x: f"{x}_target")),
                ),
                membership=MembershipBuilder(
                    MembersFromList(
                        entity_type=EntityTypeCv.PROTEIN,
                        identifiers=IdentifiersBuilder(
                            CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=f('receptor_subunit')),
                        ),
                    )
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=ParticipantMetadataCv.TARGET),
            ),
        ),
    ),
)


def get_neuronchat_dataset(species: str, taxon: int) -> Dataset:
    """
    Helper to create a NeuronChat dataset for a given species.
    """
    return Dataset(
        download=Download(
            url=f'https://github.com/Wei-BioMath/NeuronChat/raw/main/data/interactionDB_{species}.rda',
            filename=f'neuronchat_interactions_{species}.rda',
            subfolder='neuronchat',
            default_mode='rb',
            encoding=None, # Crucial: avoid encoding in binary mode
        ),
        mapper=interactions_schema,
        raw_parser=iter_neuronchat,
    )


resource = Resource(
    config,
    human=get_neuronchat_dataset('human', 9606),
    mouse=get_neuronchat_dataset('mouse', 10090),
)
