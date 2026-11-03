"""Evidence and Conclusion Ontology (ECO): the types of evidence that support an assertion."""

from __future__ import annotations

from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content

from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.inputs_v2.base import (
    Dataset,
    Download,
    Resource,
    ResourceConfig,
    ontology_entity_mapper,
)
from pypath.inputs_v2.parsers.obo import iter_obo, obo_record_to_term


config = ResourceConfig(
    id=ResourceCv.ECO,
    name='ECO',
    full='Evidence and Conclusion Ontology',
    url='https://www.evidenceontology.org/',
    # eco.obo: "remark: ECO (https://github.com/evidenceontology/evidenceontology) is released
    # into the public domain under CC0 1.0 Universal (CC0 1.0)." and "property_value:
    # dcterms:license https://creativecommons.org/publicdomain/zero/1.0/"; the OBO Foundry
    # registry entry (obofoundry.org/ontology/eco.html) gives "CC0 1.0",
    # https://creativecommons.org/publicdomain/zero/1.0/.
    license=LicenseCV.CC0_1_0,
    update_category=UpdateCategoryCV.REGULAR,
    # ECO: the Evidence and Conclusion Ontology, an update for 2022 (Nucleic Acids Res 2022),
    # the publication the OBO Foundry registry marks as preferred.
    pubmed='34986598',
    content=(Content.ONTOLOGIES,),
    summary='Evidence and Conclusion Ontology: the types of scientific evidence (assay, computational, curated) that support an assertion.',
    description=(
        'The Evidence and Conclusion Ontology (ECO) describes the types of scientific '
        'evidence that arise from laboratory experiments, computational methods, literature '
        'curation and other means, and the assertion methods (manual or automatic) that '
        'use them. Resources cite its terms (ECO:0000178, ECO:0000068) as the type of '
        'evidence behind an interaction, a modification or an annotation; this module '
        'provides their names, synonyms, definitions and hierarchy.'
    ),
)

terms_schema = ontology_entity_mapper(
    obo_record_to_term,
    ontology_id='eco',
    identifier_type=Namespace.ECO,
    id_prefix='ECO:',
)


resource = Resource(
    config,
    terms=Dataset(
        download=Download(
            # eco.obo (data-version releases/2026-07-10) holds ECO's own terms alone; the PURL
            # of eco-base.obo and eco-basic.obo does not resolve, and eco-basic.obo on GitHub
            # carries the same terms.
            url='https://purl.obolibrary.org/obo/eco.obo',
            filename='eco.obo',
            subfolder='eco',
            large=True,
        ),
        mapper=terms_schema,
        raw_parser=iter_obo,
    ),
)
