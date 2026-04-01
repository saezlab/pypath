"""
Parse STITCH chemical-protein interaction data and emit Entity records.

This module converts STITCH interaction data into Entity records using the
declarative schema pattern. It merges the actions file (directionality, mode)
with the links file (confidence sub-scores) and filters by combined_score.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    LicenseCV,
    PharmacologicalActionCv,
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
from pypath.inputs_v2.parsers.stitch import iter_stitch_interactions


# =============================================================================
# Resource Configuration
# =============================================================================

config = ResourceConfig(
    id=ResourceCv.STITCH,
    name='STITCH',
    url='https://stitch.embl.de/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='26590256',
    description=(
        'STITCH is a database of known and predicted interactions between '
        'chemicals and proteins. Interactions include direct (physical) and '
        'indirect (functional) associations derived from genomic context, '
        'high-throughput experiments, conserved coexpression, and published '
        'literature.'
    ),
)


# =============================================================================
# Custom Download
# =============================================================================

@dataclass(frozen=True)
class _StitchOpener:
    """Holds the two file openers needed to parse a STITCH dataset."""

    links: Any
    actions: Any


@dataclass(frozen=True)
class StitchDownload:
    """
    Combines two Download objects (links and actions) into a single opener.
    """

    links: Download
    actions: Download

    def open(
            self,
            *,
            force_refresh: bool = False,
            **_kwargs: Any,
    ) -> _StitchOpener:
        """
        Open both STITCH files and return a combined opener.

        Args:
            force_refresh: Force re-download even if cached.
            **_kwargs: Ignored; accepted for interface compatibility.

        Returns:
            A _StitchOpener with .links and .actions attributes.
        """
        return _StitchOpener(
            links=self.links.open(force_refresh=force_refresh),
            actions=self.actions.open(force_refresh=force_refresh),
        )


# =============================================================================
# Field Config and Schema
# =============================================================================

f = FieldConfig(
    map={
        'action_cv': {
            'activation': PharmacologicalActionCv.ACTIVATION,
            'inhibition': PharmacologicalActionCv.INHIBITION,
        },
        'mode_cv': {
            'activation': PharmacologicalActionCv.ACTIVATION,
            'inhibition': PharmacologicalActionCv.INHIBITION,
            'binding': PharmacologicalActionCv.BINDING,
        },
    },
)

interactions_schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.NAME, value=f('interaction_name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=InteractionMetadataCv.CONFIDENCE_VALUE, value=f('combined_score')),
        CV(term=InteractionMetadataCv.CONFIDENCE_VALUE, value=f('experimental')),
        CV(term=InteractionMetadataCv.CONFIDENCE_VALUE, value=f('prediction')),
        CV(term=InteractionMetadataCv.CONFIDENCE_VALUE, value=f('database')),
        CV(term=InteractionMetadataCv.CONFIDENCE_VALUE, value=f('textmining')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('mode', map='mode_cv')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('stereospecific')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('chem_is_acting')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('prot_is_acting')),
        CV(term=InteractionMetadataCv.CONTROL_TYPE, value=f('action', map='action_cv')),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.SMALL_MOLECULE,
                identifiers=IdentifiersBuilder(
                    CV(
                        term=IdentifierNamespaceCv.PUBCHEM_COMPOUND,
                        value=f('chemical_id'),
                    ),
                ),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.ENSEMBL, value=f('protein_id')),
                    CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('ncbi_tax_id')),
                ),
            ),
        ),
    ),
)


# =============================================================================
# Resource Definition
# =============================================================================

_SUBFOLDER = 'stitch'
_LINKS_URL = (
    'http://stitch.embl.de/download/'
    'protein_chemical.links.detailed.v5.0/'
    '{taxid}.protein_chemical.links.detailed.v5.0.tsv.gz'
)
_ACTIONS_URL = (
    'http://stitch.embl.de/download/actions.v5.0/'
    '{taxid}.actions.v5.0.tsv.gz'
)


def _stitch_download(ncbi_tax_id: int) -> StitchDownload:
    """
    Build a StitchDownload for the given species.

    Args:
        ncbi_tax_id: NCBI taxonomy ID of the target organism.

    Returns:
        A StitchDownload configured with links and actions downloads.
    """
    return StitchDownload(
        links=Download(
            url=_LINKS_URL.format(taxid=ncbi_tax_id),
            filename=f'{ncbi_tax_id}.protein_chemical.links.detailed.v5.0.tsv.gz',
            subfolder=_SUBFOLDER,
        ),
        actions=Download(
            url=_ACTIONS_URL.format(taxid=ncbi_tax_id),
            filename=f'{ncbi_tax_id}.actions.v5.0.tsv.gz',
            subfolder=_SUBFOLDER,
        ),
    )


resource = Resource(
    config,
    human_interactions=Dataset(
        download=_stitch_download(9606),
        mapper=interactions_schema,
        raw_parser=iter_stitch_interactions,
    ),
    mouse_interactions=Dataset(
        download=_stitch_download(10090),
        mapper=interactions_schema,
        raw_parser=iter_stitch_interactions,
    ),
)
