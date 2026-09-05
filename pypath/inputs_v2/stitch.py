"""
Parse STITCH chemical-protein interaction data and emit Entity records.

This module converts STITCH interaction data into Entity records using the
declarative schema pattern. It merges the actions file (directionality, mode)
with the links file (confidence sub-scores) and filters by combined_score.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Any

from biolink_model.datamodel.model import (
    ChemicalEntity,
    DirectionQualifierEnum,
    Protein,
    QuantityValue,
    slots,
)
from omnipath_core.measurements import Measurement
from omnipath_core.naming import Namespace

from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.stitch import iter_stitch_interactions
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)


def _measurement(value, unit=None, source_field=None, comparator=None):
    """Keep numeric source observations, units and comparison bounds together."""
    if value is None:
        return None
    match = re.fullmatch(
        '\\s*(<=|>=|<|>|=|~)?\\s*([+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?)\\s*',
        str(value),
    )
    if match is None:
        return None
    original_comparator = comparator or match.group(1)
    if original_comparator not in {None, '<', '>', '=', '<=', '>=', '~', '≈'}:
        return None
    relation = {'<': 'less_than', '>': 'greater_than', '=': 'equal_to'}.get(
        original_comparator
    )
    number = float(match.group(2))
    if not math.isfinite(number):
        return None
    quantity = QuantityValue(
        has_numeric_value=number,
        has_unit=unit or None,
        has_binary_relation=relation,
    )
    return Measurement(
        quantity=quantity,
        source_field=source_field,
        comparator=original_comparator,
    )


config = ResourceConfig(
    id=ResourceCv.STITCH,
    name='STITCH',
    url='https://stitch.embl.de/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='26590256',
    primary_category='interactions',
    description='STITCH is a database of known and predicted interactions between chemicals and proteins. Interactions include direct (physical) and indirect (functional) associations derived from genomic context, high-throughput experiments, conserved coexpression, and published literature.',
)


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
        self, *, force_refresh: bool = False, **_kwargs: Any
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


f = FieldConfig()
ACTION_DIRECTION = {
    'activation': DirectionQualifierEnum.increased,
    'inhibition': DirectionQualifierEnum.decreased,
}


def stitch_direction(row):
    if row.get('chem_role') != 'source' and row.get('prot_role') != 'source':
        return None
    return ACTION_DIRECTION.get(row.get('action')) or ACTION_DIRECTION.get(
        row.get('mode')
    )


def stitch_predicate(row):
    if stitch_direction(row) is not None:
        return slots.affects
    if row.get('mode') == 'binding':
        return slots.physically_interacts_with
    return slots.associated_with


chemical_builder = EntityBuilder(
    entity_type=ChemicalEntity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.PUBCHEM, value=f('chemical_id'))
    ),
    annotations=AnnotationsBuilder(),
)
protein_builder = EntityBuilder(
    entity_type=Protein,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.ENSP, value=f('protein_id'))
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.in_taxon,
            value=f(
                'ncbi_tax_id',
                transform=lambda v: 'NCBITaxon:'
                + str(v).removeprefix('NCBITaxon:'),
            ),
        )
    ),
)
interactions_schema = RelationBuilder(
    subject=lambda row: (
        protein_builder
        if row.get('prot_role') == 'source'
        else chemical_builder
    ).build(row),
    predicate=stitch_predicate,
    object=lambda row: (
        chemical_builder
        if row.get('prot_role') == 'source'
        else protein_builder
    ).build(row),
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.NAME, value=f('interaction_name'))
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.object_direction_qualifier, value=stitch_direction),
        CV(
            term=slots.has_confidence_score,
            value=f(
                'combined_score',
                transform=lambda v: _measurement(
                    v, source_field='combined_score'
                ),
            ),
        ),
        CV(
            term=slots.has_confidence_score,
            value=f(
                'action_score',
                transform=lambda v: _measurement(
                    v, source_field='action_score'
                ),
            ),
        ),
    ),
)
_LINKS_URL = 'http://stitch.embl.de/download/protein_chemical.links.detailed.v5.0/{taxid}.protein_chemical.links.detailed.v5.0.tsv.gz'
_ACTIONS_URL = (
    'http://stitch.embl.de/download/actions.v5.0/{taxid}.actions.v5.0.tsv.gz'
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
            subfolder='stitch',
            ext='gz',
        ),
        actions=Download(
            url=_ACTIONS_URL.format(taxid=ncbi_tax_id),
            filename=f'{ncbi_tax_id}.actions.v5.0.tsv.gz',
            subfolder='stitch',
            ext='gz',
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
