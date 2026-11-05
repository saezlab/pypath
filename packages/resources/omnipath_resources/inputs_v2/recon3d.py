"""Recon3D metabolic activities, chemical participants and model genes.

A small module over the shared BiGG factory (:mod:`omnipath_resources.inputs_v2._bigg`):
states Recon3D's taxon and GPR identifier prefix, and nothing the factory
doesn't already know. See :mod:`omnipath_resources.inputs_v2._bigg` for what the four
datasets carry and :mod:`omnipath_resources.inputs_v2.imm1415` for the other model built
on the same factory.
"""

from __future__ import annotations

from omnipath_core.resource_metadata import Content, Topic

from omnipath_resources.inputs_v2._bigg import bigg_resource
from omnipath_resources.inputs_v2.base import Download, Resource, ResourceConfig
from omnipath_resources.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from omnipath_resources.resources import urls

HUMAN_TAXON_ID = 'NCBITaxon:9606'
GPR_PREFIX = 'recon3d.gpr_and'

config = ResourceConfig(
    id=ResourceCv.RECON3D,
    name='Recon3D',
    url='https://www.vmh.life/',
    license=LicenseCV.BIGG,
    update_category=UpdateCategoryCV.IRREGULAR,
    pubmed='29457794',
    topics=(Topic.METABOLISM,),
    content=(Content.BIOCHEMICAL_REACTIONS,),
    summary='Human metabolic activities, chemical participants and model gene associations.',
    description=(
        'Recon3D is a comprehensive human genome-scale metabolic reconstruction '
        'that accounts for three-dimensional molecular structures.'
    ),
)

download = Download(
    url=urls.urls['recon3d']['bigg_json'],
    filename='Recon3D.json',
    subfolder='recon3d',
    large=True,
    ext='json',
    default_mode='r',
)

_bigg_datasets = bigg_resource(
    download=download, taxon=HUMAN_TAXON_ID, gpr_prefix=GPR_PREFIX, orphan_filter=True,
)
metabolites_schema = _bigg_datasets['metabolites'].mapper
reactions_schema = _bigg_datasets['reactions'].mapper
enzyme_complexes_schema = _bigg_datasets['enzyme_complexes'].mapper
genes_schema = _bigg_datasets['genes'].mapper

resource = Resource(config, **_bigg_datasets)
