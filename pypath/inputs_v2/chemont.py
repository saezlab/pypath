"""ChemOnt chemical taxonomy ontology resource."""

from __future__ import annotations

from pypath.internals.cv_terms import (
    LicenseCV,
    OntologyCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.inputs_v2.base import (
    Dataset,
    Download,
    Resource,
    ResourceConfig,
    ontology_entity_mapper,
)
from pypath.inputs_v2.parsers.obo import iter_obo, obo_record_to_term


config = ResourceConfig(
    id=ResourceCv.CHEMONT,
    name='ChemOnt',
    url='http://classyfire.wishartlab.com/',
    license=LicenseCV.ACADEMIC_FREE,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='27867422',
    primary_category='ontologies',
    annotation_ontologies=(OntologyCv.CHEMONT,),
    resource_kind='ontology',
    description=(
        'ChemOnt is the chemical taxonomy used by ClassyFire for hierarchical '
        'chemical classification of compounds.'
    ),
)


terms_schema = ontology_entity_mapper(obo_record_to_term, ontology_id='chemont')


resource = Resource(
    config,
    terms=Dataset(
        download=Download(
            url=(
                'http://classyfire.wishartlab.com/system/downloads/1_0/'
                'chemont/ChemOnt_2_1.obo.zip'
            ),
            filename='ChemOnt_2_1.obo.zip',
            subfolder='chemont',
            large=True,
            ext='zip',
            needed=['ChemOnt_2_1.obo'],
        ),
        mapper=terms_schema,
        raw_parser=iter_obo,
    ),
)
