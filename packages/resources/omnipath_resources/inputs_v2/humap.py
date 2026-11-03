"""
Map the complexes of the three versions of hu.MAP.

hu.MAP (Human Protein Complex Map) is a map of human protein complexes that
the Marcotte and Drew labs predicted by integrating mass spectrometry
experiments: protein interaction scores from a machine learning model are
clustered into complexes. The complexes are computational predictions.

The three versions are different maps: they have other files, identifiers,
members and confidence scales, so each is its own dataset and nothing is mixed.

- ``humap1_complexes``: hu.MAP 1.0 (2017). The members are Entrez Gene IDs
  (a few Ensembl gene IDs) and the file has no identifier or confidence; the
  identifier is the position of the complex in the file (see the parser).
- ``humap2_complexes``: hu.MAP 2.0 (2021). Identifier ``HuMAP2_00000``,
  confidence 1 (extremely high) to 5 (medium), members UniProt accessions
  with their gene names.
- ``humap3_complexes``: hu.MAP 3.0 (2025). Identifier ``huMAP3_00000.1``,
  confidence 1 (extremely high) to 6 (medium), members UniProt accessions.

The confidence level is the one that the file gives; no complex is left out by
it. A complex has its members and no pairs are made. Complexes with the same
members stay separate complexes. The complexes and their members are human, as
the map is.
"""

from __future__ import annotations

import functools
from typing import Any

from biolink_model.datamodel.model import MacromolecularComplex, Protein, slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.humap import iter_humap1, iter_humap_table
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import Annotation, Entity, Identifier, Membership

config = ResourceConfig(
    id=ResourceCv.HUMAP,
    name='hu.MAP',
    url='https://humap3.proteincomplexes.org/',
    # The hu.MAP site states CC0 for its downloads, with a request for
    # attribution (checked when the module was written, October 2026); the
    # legacy registry recorded CC BY 4.0.
    license=LicenseCV.CC0_1_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='40425816',
    content=(Content.PROTEIN_COMPLEXES,),
    knowledge_origin=(KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,),
    summary=(
        'Human protein complexes predicted from mass spectrometry experiments, '
        'hu.MAP 1.0, 2.0 and 3.0 as separate datasets.'
    ),
    description=(
        'hu.MAP is a map of human protein complexes built by integrating '
        'thousands of published mass spectrometry experiments (affinity '
        'purification, co-fractionation, proximity labeling) with a machine '
        'learning model that scores protein pairs, and clustering the scored '
        'network. hu.MAP 1.0 (Drew et al. 2017) has about 4,600 complexes, '
        'hu.MAP 2.0 (Drew et al. 2021) about 7,000 and hu.MAP 3.0 (Fischer et '
        'al. 2025) about 15,000. Each complex of version 2 and 3 has a '
        'confidence level from the clustering that found it.'
    ),
)

_TAXON = 'NCBITaxon:9606'

# The confidence level of a complex -> its name on the download page of the
# version. The scales differ: 2.0 has five levels and 3.0 has six.
_CONFIDENCE = {
    2: {
        '1': 'Extremely High',
        '2': 'Very High',
        '3': 'High',
        '4': 'Medium High',
        '5': 'Medium',
    },
    3: {
        '1': 'Extremely High',
        '2': 'Very High',
        '3': 'High',
        '4': 'Moderate High',
        '5': 'Medium High',
        '6': 'Medium',
    },
}


def _member(member: dict[str, Any]) -> Entity:
    """A protein, identified by the identifiers that the file gives for it."""
    identifiers = [
        Identifier(namespace, member[key])
        for key, namespace in (
            ('uniprot', Namespace.UNIPROT),
            ('entrez', Namespace.ENTREZ),
            ('ensembl', Namespace.ENSG),
            ('entry', Namespace.UNIPROT_ENTRY),
            ('genename', Namespace.GENESYMBOL),
        )
        if member[key]
    ]
    return Entity(
        type=Protein,
        identifiers=identifiers,
        # The map is of human proteins. A UniProt entry name ends with the
        # organism (the one in hu.MAP 1 is of a jellyfish) and states it.
        annotations=None if member['entry'] else [Annotation(slots.in_taxon, _TAXON)],
    )


def map_complex(row: dict[str, Any], version: int) -> Entity | None:
    """A complex with its members and its confidence level; None without members."""
    members = [_member(member) for member in row['members']]
    if not members:
        return None
    annotations = [Annotation(slots.in_taxon, _TAXON)]
    if version > 1:
        level = row['confidence']
        if level not in _CONFIDENCE[version]:
            raise ValueError(f'Unknown hu.MAP {version} confidence level: {level!r}')
        # The level as the file writes it, with its name on the download page.
        annotations.append(
            Annotation(
                slots.has_confidence_level,
                f'{_CONFIDENCE[version][level]} ({level})',
            )
        )
    return Entity(
        type=MacromolecularComplex,
        identifiers=[Identifier(f'humap{version}', row['complex'])],
        annotations=annotations,
        membership=[Membership(member=member) for member in members],
    )


def _dataset(version: int, url: str, filename: str) -> Dataset:
    return Dataset(
        download=Download(
            url=url,
            filename=filename,
            subfolder='humap',
            ext=filename.rsplit('.', 1)[1],
        ),
        mapper=functools.partial(map_complex, version=version),
        raw_parser=(
            iter_humap1
            if version == 1
            else functools.partial(iter_humap_table, version=version)
        ),
    )


resource = Resource(
    config,
    humap1_complexes=_dataset(
        1,
        'https://hu1.proteincomplexes.org/static/downloads/clusters.txt',
        'humap1_clusters.txt',
    ),
    humap2_complexes=_dataset(
        2,
        'https://humap2.proteincomplexes.org/static/downloads/humap2/'
        'humap2_complexes_20200809.txt',
        'humap2_complexes_20200809.csv',
    ),
    humap3_complexes=_dataset(
        3,
        'https://humap3.proteincomplexes.org/static/downloads/humap3/'
        'hu.MAP3.0_complexes_wConfidenceScores_total15326_wGenenames_20240922.csv',
        'humap3_complexes_20240922.csv',
    ),
)
