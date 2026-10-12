"""
Parse the HTRIdb transcription factor binding observations and emit Relation records.

Source file (semicolon separated, with a header and a ``;`` ending every line;
one experimental observation per line):

    OID;GENEID_TF;SYMBOL_TF;GENEID_TG;SYMBOL_TG;TECHNIQUE;PUBMED_ID;
    32399;142;PARP1;675;BRCA2;Electrophoretic Mobility Shift Assay;18990703;
    28907;196;AHR;1543;CYP1A1;Chromatin Immunoprecipitation;16972788;

HTRIdb (Human Transcriptional Regulation Interactions database) lists human
transcription factors (TF) and the genes (target genes, TG) they were found to
act on, each with the technique that detected it and the PubMed ID of the
article. The techniques are binding assays (ChIP, EMSA, footprinting, ...):
what the source observed is a transcription factor bound to the target gene,
not a change of its expression, and the file has no mode of regulation. The
relation is therefore ``interacts_with`` without direction or aspect, in
contrast to the regulons that state an effect (``affects``, see
``trrust.py``).

The file is the copy OmniPath rescued from the HTRIdb site, which has been
offline for years. It has 18,630 lines, not the 51,871 interactions the
publication reports.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._publications import PMID_EXTRACT
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.htri import iter_observations
from pypath.internals.cv_terms import (
    LicenseCV,
    OntologyCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

# HTRIdb is a database of human interactions; the file has no organism column.
HTRI_TAX_ID = '9606'
# Technique -> PSI-MI term of the same method, where PSI-MI has one. A
# technique mapped to None is known but has no exact PSI-MI term; its name
# stays in the parsed row only. Any other value is an error.
_TECHNIQUES = {
    'Chromatin Immunoprecipitation': 'MI:0402',
    # PSI-MI names this method "chromatin immunoprecipitation array"
    # (synonym "chip-chip").
    'Chromatin Immunoprecipitation coupled with microarray': 'MI:0225',
    'Electrophoretic Mobility Shift Assay': 'MI:0413',
    'Dnase I Footprinting': 'MI:0606',  # DNase I footprinting
    'Surface Plasmon Resonance': 'MI:0107',
    'Yeast One-Hybrid Assay': 'MI:0432',  # one hybrid
    'Southwestern Blotting': 'MI:2283',
    'CpG Chromatin Immunoprecipitation': None,
    'Concatenate Chromatin Immunoprecipitation': None,
    'Streptavidin Chromatin Immunoprecipitation': None,
    'Avidin/Biotin-conjugated DNA Binding Assay': None,
    'DNA Affinity Chromatography': None,
    'DNA Affinity Precipitation Assay': None,
}


def _evidence_types(row: dict[str, Any]) -> list[str]:
    technique = row.get('TECHNIQUE')
    if technique not in _TECHNIQUES:
        raise ValueError(f'Unknown HTRIdb technique: {technique!r}')
    term = _TECHNIQUES[technique]
    return [term] if term else []


config = ResourceConfig(
    id=ResourceCv.HTRI,
    name='HTRIdb',
    full='Human Transcriptional Regulation Interactions database',
    url='http://www.lbbc.ibb.unesp.br/htri',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.LGPL_3_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='22900683',
    topics=(Topic.GENE_REGULATION,),
    content=(Content.PHYSICAL_INTERACTIONS,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.EXPERIMENTAL,
    ),
    annotation_ontologies=(OntologyCv.PSI_MI,),
    summary=(
        'Human transcription factor binding to target genes, with the '
        'detection technique and the PubMed ID of each observation.'
    ),
    description=(
        'HTRIdb (Human Transcriptional Regulation Interactions database) is '
        'a database of experimentally verified interactions between human '
        'transcription factors and their target genes, from small and '
        'mid-scale assays such as chromatin immunoprecipitation, '
        'electrophoretic mobility shift assay and DNase I footprinting, and '
        'from high-throughput ChIP-chip. The articles and their abstracts '
        'were searched by text mining for gene and technique names and the '
        'interactions were then curated; each interaction has the Entrez '
        'Gene IDs and symbols of the factor and the gene, the technique and '
        'the PubMed ID. The website, at the Universidade Estadual Paulista, '
        'is no longer online.'
    ),
)

f = FieldConfig(
    extract={'pmid': PMID_EXTRACT},
    transform={'pmid_curie': lambda v: f'PMID:{v}'},
)


def _participant(entity_type, role: str) -> EntityBuilder:
    return EntityBuilder(
        entity_type=entity_type,
        identifiers=IdentifiersBuilder(
            CV(term=Namespace.GENESYMBOL, value=f(f'SYMBOL_{role}')),
            CV(term=Namespace.ENTREZ, value=f(f'GENEID_{role}')),
        ),
        annotations=AnnotationsBuilder(
            CV(term=slots.in_taxon, value=f'NCBITaxon:{HTRI_TAX_ID}'),
        ),
    )


# A row without a transcription factor or a target yields no relation: the
# participant has no identifier, so the builder returns None.
interactions_schema = RelationBuilder(
    subject=_participant(model.Protein, 'TF'),
    predicate=slots.interacts_with,
    object=_participant(model.Gene, 'TG'),
    annotations=AnnotationsBuilder(
        CV(term=slots.has_evidence_of_type, value=_evidence_types),
        CV(
            term=slots.publications,
            value=f('PUBMED_ID', extract='pmid', transform='pmid_curie'),
        ),
    ),
)

resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            url='https://rescued.omnipathdb.org/HTRIdb.csv',
            filename='HTRIdb.csv',
            subfolder='htri',
            large=True,
        ),
        mapper=interactions_schema,
        raw_parser=iter_observations,
    ),
)
