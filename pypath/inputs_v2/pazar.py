"""
PAZAR: transcription factors and the genes they were assayed on.

Source file: ``tftargets.zip``, the transcription factor - target gene lists the
PAZAR site offered (release of 2012-05-22). One line is one assay of a
transcription factor on a regulatory sequence of a gene; see
``parsers/pazar.py`` for the columns. The factor is named by an Ensembl
transcript (a few by gene or by a FlyBase or worm identifier) and has no
species of its own; the species is that of the target gene.

The assays are mostly binding assays (ChIP, EMSA, footprinting) and records of
binding sites (the imports of ABS, ORegAnno and other collections of sites,
which give no assay or a curation label, and a few predictions and conserved
sites); a few hundred are reporter gene and expression assays. What a line
states follows its assay:

- a reporter gene assay or an RNA or protein expression assay shows the effect
  of the factor on expression, without a sign in the file: ``affects``,
  without direction or aspect;
- every other line records that the factor binds a site of the gene, with a
  binding assay, a prediction, a conservation of the site or no assay named:
  ``interacts_with``. Binding does not show an effect on expression.

The assay of the line is on the relation (``pazar:method``), so that this
rule can be applied again.
"""

from __future__ import annotations

from typing import Any

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace
from omnipath_core.resource_metadata import Content, KnowledgeOrigin, Topic

from pypath.inputs_v2._mitab import ensembl_namespace
from pypath.inputs_v2._publications import PMID_EXTRACT
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.pazar import iter_tf_targets
from pypath.internals.cv_terms import (
    LicenseCV,
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

# The species of the target gene, as the file writes it.
_TAXA = {
    'Homo sapiens': '9606',
    'Mus musculus': '10090',
    'Rattus norvegicus': '10116',
    'Canis familiaris': '9615',
    'Drosophila melanogaster': '7227',
    'Gallus gallus': '9031',
    'Bos taurus': '9913',
    'Caenorhabditis elegans': '6239',
    'Danio rerio': '7955',
}

# The assay of a line -> an Evidence & Conclusion Ontology term, where one
# names the assay. The assays of the same kind are one term; those with no term
# are None. An assay that is not listed is an error.
_ECO_CHIP = 'ECO:0000226'  # chromatin immunoprecipitation evidence
_ECO_EMSA = 'ECO:0000096'  # electrophoretic mobility shift assay evidence
_ECO_LUCIFERASE = 'ECO:0001805'  # luciferase reporter gene assay evidence
_ECO_FOOTPRINT = 'ECO:0001810'  # DNAse footprinting evidence
_ECO_YEAST_ONE_HYBRID = 'ECO:0000066'  # yeast one-hybrid evidence
_ECO_WESTERN = 'ECO:0000112'  # qualitative western immunoblotting evidence
_METHODS = {
    'PROTEIN BINDING ASSAY::CHROMATIN IMMUNOPRECIPITATION (CHIP)': _ECO_CHIP,
    'PROTEIN BINDING ASSAY::CHROMATIN IMMUNOPRECIPITATION WITH TAG SEQUENCING (CHIP-TS)': _ECO_CHIP,
    'UNKNOWN::CHROMATIN IMMUNOPRECIPITATION (CHIP)': _ECO_CHIP,
    'ELECTROPHORETIC MOBILITY SHIFT ASSAY (EMSA)::SUPERSHIFT': _ECO_EMSA,
    'ELECTROPHORETIC MOBILITY SHIFT ASSAY (EMSA)::DIRECT GEL SHIFT': _ECO_EMSA,
    'ELECTROPHORETIC MOBILITY SHIFT ASSAY (EMSA)::GEL SHIFT COMPETITION': _ECO_EMSA,
    'ELECTROPHORETIC MOBILITY SHIFT ASSAY (EMSA)::UNKNOWN': _ECO_EMSA,
    'ELECTROPHORETIC MOBILITY SHIFT ASSAY (EMSA)::DIRECT GEL SHIFT WITH WESTERN BLOT DETECTION': _ECO_EMSA,
    'EMSA': _ECO_EMSA,
    'REPORTER GENE ASSAY::CHLORAMPHENICOL ACETYLTRANSFERASE (CAT) ASSAY': 'ECO:0001803',
    'REPORTER GENE ASSAY::TRANSIENT TRANSFECTION LUCIFERASE ASSAY': _ECO_LUCIFERASE,
    'REPORTER GENE ASSAY::DUAL LUCIFERASE REPORTER GENE ASSAY': _ECO_LUCIFERASE,
    'REPORTER GENE ASSAY::IN-VIVO LUCIFERASE EXPRESSION ASSAY': _ECO_LUCIFERASE,
    'REPORTER GENE ASSAY::STABLE TRANSFECTION LUCIFERASE ASSAY': _ECO_LUCIFERASE,
    'REPORTER GENE ASSAY::IN-VIVO GFP EXPRESSION ASSAY': 'ECO:0001816',
    'REPORTER GENE ASSAY::IN-VIVO LACZ EXPRESSION ASSAY': 'ECO:0001802',
    # The reporter gene assay in general (the reporter is the globin gene).
    'REPORTER GENE ASSAY::RABBIT BETA-GLOBIN EXPRESSION ASSAY': 'ECO:0000049',
    'PROTEIN BINDING ASSAY::DNASE FOOTPRINTING ASSAY': _ECO_FOOTPRINT,
    'PROTEIN BINDING ASSAY::DNASE FOOTPRINTING COMPETITION ASSAY': _ECO_FOOTPRINT,
    'DNASE I FOOTPRINTING': _ECO_FOOTPRINT,
    'YEAST ONE-HYBRID': _ECO_YEAST_ONE_HYBRID,
    'PROTEIN BINDING ASSAY::YEAST 1-HYBRID ASSAY': _ECO_YEAST_ONE_HYBRID,
    'RNA EXPRESSION ASSAY::REVERSE TRANSCRIPTASE POLYMERASE CHAIN REACTION (RT-PCR)': 'ECO:0000109',
    'RNA EXPRESSION ASSAY::NORTHERN BLOT': 'ECO:0000106',
    'RNA EXPRESSION ASSAY::RNASE PROTECTION ASSAY (RPA)': 'ECO:0000110',
    'PROTEIN BINDING ASSAY::WESTERN BLOT ASSAY': _ECO_WESTERN,
    'PROTEIN EXPRESSION ASSAY::WESTERN BLOT ASSAY': _ECO_WESTERN,
    'PROTEIN EXPRESSION ASSAY::ENZYME-LINKED IMMUNOSORBENT ASSAY (ELISA)': 'ECO:0000267',
    'COMPUTATIONAL PREDICTION': 'ECO:0007672',
    'SEQUENCE CONSERVATION::CONSERVATION FOUND BY ALIGNMENT': 'ECO:0000200',
    'SEQUENCE CONSERVATION::ORTHOLOGOUS GENE CONSERVATION': 'ECO:0000201',
    # No term for the assay, or no assay given. The last is a PubMed ID that
    # sits in the column of one line.
    'SEQUENCE CONSERVATION::CO-EXPRESSED GENE CONSERVATION': None,
    'PROTEIN BINDING ASSAY::DNA-PROTEIN PRECIPITATION ASSAY': None,
    'PROTEIN BINDING ASSAY::DNA-PROTEIN PRECIPITATION ASSAY WITH ELISA DETECTION': None,
    'PROTEIN BINDING ASSAY::METHYLATION INTERFERENCE ASSAY': None,
    'PROTEIN BINDING ASSAY::METHYLATION INTERFERENCE': None,
    'PROTEIN BINDING ASSAY::IN-VIVO FOOTPRINTING ASSAY': None,
    'PROTEIN BINDING ASSAY::UNKNOWN': None,
    'MUTAGENESIS::SITE-DIRECTED': None,
    'MUTAGENESIS::OLIGONUCLEOTIDE SYNTHESIS': None,
    'MUTAGENESIS::UNKNOWN': None,
    'EXPERT CURATED::LITERATURE DERIVED': None,
    'UNKNOWN::UNKNOWN': None,
    'NA': None,
    '': None,
    '15688014': None,
}
# The classes ("class::method") of the assays that measure expression: the
# reporter gene assays and the RNA and protein expression assays. A line with
# one of these shows an effect of the factor on expression, without a sign.
# Every other class of the file (protein binding, EMSA, footprinting, ChIP,
# mutagenesis, sequence conservation, computational prediction, expert
# curation, unknown) is a binding site record of the factor.
_EXPRESSION_ASSAY_CLASSES = frozenset(
    {
        'REPORTER GENE ASSAY',
        'RNA EXPRESSION ASSAY',
        'PROTEIN EXPRESSION ASSAY',
    }
)
# The assay column holds no assay in these lines: "NA", empty, and a PubMed ID
# that sits in the column of one line.
_NO_ASSAY = frozenset({'NA', '', '15688014'})
_METHOD = 'pazar:method'
# The gene identifier types a target gene can have; a transcript or protein
# identifier in the column of the gene does not name a gene.
_GENE_NAMESPACES = frozenset({Namespace.ENSG, Namespace.FLYBASE, Namespace.NAME})


def _identifier(value: object) -> tuple[Any, str] | None:
    """The Ensembl or FlyBase identifier, else a name; None for digits only."""
    text = str(value or '').strip()
    if not text:
        return None
    namespace = ensembl_namespace(text)
    if namespace is not None:
        return namespace, text
    if text.startswith(('FBgn', 'FBtr')) and text[4:].isdigit():
        return Namespace.FLYBASE, text
    # A line has the digits of a transcript identifier without its prefix.
    if text.isdigit():
        return None
    return Namespace.NAME, text


def _tf_identifiers(row: dict[str, Any]) -> list[tuple[Any, str]]:
    pairs = []
    pair = _identifier(row.get('tf_transcript'))
    if pair:
        pairs.append(pair)
    name = str(row.get('tf_name') or '').strip()
    if name:
        pairs.append((Namespace.NAME, name))
    return pairs


def _gene_identifiers(row: dict[str, Any]) -> list[tuple[Any, str]]:
    pair = _identifier(row.get('gene'))
    return [pair] if pair and pair[0] in _GENE_NAMESPACES else []


def _taxon(row: dict[str, Any]) -> str:
    species = row.get('species')
    if species not in _TAXA:
        raise ValueError(f'Unknown PAZAR species: {species!r}')
    return f'NCBITaxon:{_TAXA[species]}'


def _evidence_types(row: dict[str, Any]) -> list[str]:
    method = row.get('method')
    if method not in _METHODS:
        raise ValueError(f'Unknown PAZAR assay: {method!r}')
    return [_METHODS[method]] if _METHODS[method] else []


def _predicate(row: dict[str, Any]):
    """``affects`` for an expression assay, ``interacts_with`` otherwise."""
    _evidence_types(row)  # an assay that is not listed is an error
    if row['method'].partition('::')[0] in _EXPRESSION_ASSAY_CLASSES:
        return slots.affects
    return slots.interacts_with


def _method(row: dict[str, Any]) -> str | None:
    """The assay as the file writes it; None where the line names none."""
    _evidence_types(row)
    return None if row['method'] in _NO_ASSAY else row['method']


config = ResourceConfig(
    id=ResourceCv.PAZAR,
    name='PAZAR',
    full='Public Database of Transcription Factor and Regulatory Sequence Annotation',
    url='http://www.pazar.info/',
    # License as recorded in the legacy registry (resources.json); the source states none today.
    license=LicenseCV.LGPL_3_0,
    update_category=UpdateCategoryCV.STATIC,
    pubmed='18971253',
    topics=(Topic.GENE_REGULATION,),
    content=(Content.CAUSAL_REGULATION,),
    knowledge_origin=(
        KnowledgeOrigin.LITERATURE_CURATED,
        KnowledgeOrigin.EXPERIMENTAL,
        KnowledgeOrigin.COMPUTATIONALLY_PREDICTED,
    ),
    summary=(
        'Transcription factors and the genes whose regulatory sequences they '
        'were assayed on, from PAZAR projects, with the assay, the '
        'publication and the project of each line; binding records are '
        'interactions, reporter and expression assays are regulation.'
    ),
    description=(
        'PAZAR is a public database of transcription factor and regulatory '
        'sequence annotation, run by the Wasserman Lab. Its transcription '
        'factor - target gene lists, as the 2012 export of the site holds '
        'them, give for each assay of a factor on the regulatory sequence of '
        'a gene the Ensembl identifiers of the factor and of the gene, the '
        'species of the gene, the PAZAR project and the PubMed ID. The '
        'projects are literature annotations (EMSA, footprinting, reporter '
        'assays; the imports of ORegAnno, TFe and ABS) and ChIP-seq data '
        'sets, each of which the export cuts to 200 lines. The site has been '
        'offline for years; this module reads the copy of the export on the '
        'OmniPath rescue server. No line states the direction of the effect '
        'of the factor on the expression of the gene; a line with a reporter '
        'gene or expression assay shows that there is an effect, every other '
        'line records the binding of the factor to a site of the gene.'
    ),
)

f = FieldConfig()

# The factor has no species of its own in the file, so it has no taxon.
transcription_factor = EntityBuilder(
    entity_type=model.Protein,
    identifiers=IdentifiersBuilder(CV.from_pairs(_tf_identifiers)),
)

# A line whose target column holds no gene identifier yields no relation: the
# target has no identifier, so the builder returns None.
target_gene = EntityBuilder(
    entity_type=model.Gene,
    identifiers=IdentifiersBuilder(CV.from_pairs(_gene_identifiers)),
    annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=_taxon)),
)

interactions_schema = RelationBuilder(
    subject=transcription_factor,
    predicate=_predicate,
    object=target_gene,
    annotations=AnnotationsBuilder(
        # Only the form of a PubMed ID: one line has nine digits.
        CV(
            term=slots.publications,
            value=f('pmid', extract=PMID_EXTRACT, transform=lambda pmid: f'PMID:{pmid}'),
        ),
        CV(term=slots.has_evidence_of_type, value=_evidence_types),
        CV(term=_METHOD, value=_method),
        # The PAZAR project of the line, as the file writes it.
        CV(term=slots.supporting_data_source, value=f('project')),
    ),
)

resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            url='https://rescued.omnipathdb.org/pazar_tftargets.zip',
            filename='pazar_tftargets.zip',
            subfolder='pazar',
            large=True,
        ),
        mapper=interactions_schema,
        raw_parser=iter_tf_targets,
    ),
)
