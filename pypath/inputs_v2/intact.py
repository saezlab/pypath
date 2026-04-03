"""
Parse IntAct data and emit Entity records.

This module converts IntAct MITAB data into Entity records using the schema
defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator
import csv

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    ParticipantMetadataCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    Column,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    Member,
    MembershipBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig


_IDENTIFIER_CV_MAPPING = {
    'bind smid': IdentifierNamespaceCv.BIND,
    'cas registry number': IdentifierNamespaceCv.CAS,
    'chebi': IdentifierNamespaceCv.CHEBI,
    'chembl': IdentifierNamespaceCv.CHEMBL,
    'chembl compound': IdentifierNamespaceCv.CHEMBL_COMPOUND,
    'ddbj/embl/genbank': IdentifierNamespaceCv.REFSEQ,
    'dip': IdentifierNamespaceCv.DIP,
    'ensembl': IdentifierNamespaceCv.ENSEMBL,
    'ensemblgenomes': IdentifierNamespaceCv.ENSEMBL_GENOMES,
    'entrezgene/locuslink': IdentifierNamespaceCv.ENTREZ,
    'flybase': IdentifierNamespaceCv.FLYBASE,
    'genbank identifier': IdentifierNamespaceCv.GENBANK_IDENTIFIER,
    'genbank_nucl_gi': IdentifierNamespaceCv.GENBANK_NUCL_GI,
    'genbank_protein_gi': IdentifierNamespaceCv.GENBANK_PROTEIN_GI,
    'hgnc': IdentifierNamespaceCv.HGNC,
    'genesymbol': IdentifierNamespaceCv.GENE_NAME_PRIMARY,
    'imex': IdentifierNamespaceCv.IMEX,
    'intact': IdentifierNamespaceCv.INTACT,
    'ipi': IdentifierNamespaceCv.IPI,
    'mint': IdentifierNamespaceCv.MINT,
    'mirbase': IdentifierNamespaceCv.MIRBASE,
    'psi-mi': IdentifierNamespaceCv.CV_TERM_ACCESSION,
    'refseq': IdentifierNamespaceCv.REFSEQ,
    'rfam': IdentifierNamespaceCv.RFAM,
    'rnacentral': IdentifierNamespaceCv.RNACENTRAL,
    'uniparc': IdentifierNamespaceCv.UNIPARC,
    'uniprotkb': IdentifierNamespaceCv.UNIPROT,
}

f = FieldConfig(
    extract={
        'prefix_lower': [r'^([^:]+):', str.lower],
        'value': r'^[^:]+:([^|"]+)',
        'mi': r'(MI:\d+)',
        # Tax fields can look like:
        # taxid:9606(human)|taxid:9606(Homo sapiens)
        # Extract the first signed integer appearing in the field.
        'tax': r'(-?\d+)',
        'pubmed': r'(?i)pubmed:(\d+)',
        'intact': r'intact:([^|"]+)',
    },
    map={
        'identifier_cv': _IDENTIFIER_CV_MAPPING,
    },
    delimiter='|',
)


def _normalize_identifier_prefix(prefix: str, value: str) -> str | None:
    prefix = prefix.lower()
    value = value.strip().strip('"')
    if prefix == 'intact' and value.startswith('MINT-'):
        return 'mint'
    if prefix == 'intact' and not (value.startswith('EBI-') or value.startswith('IM-')):
        return None
    if prefix == 'hgnc' and not value.isdigit() and not value.upper().startswith('HGNC:'):
        return 'genesymbol'
    return prefix


def _parse_identifier_pairs(raw: object) -> list[tuple[object, str]]:
    if raw is None:
        return []
    text = str(raw).strip()
    if not text or text == '-':
        return []
    pairs: list[tuple[object, str]] = []
    for item in text.split('|'):
        item = item.strip().strip('"')
        if not item or ':' not in item:
            continue
        prefix, value = item.split(':', 1)
        value = value.strip().strip('"')
        prefix = _normalize_identifier_prefix(prefix, value)
        if prefix is None:
            continue
        if prefix == 'uniprotkb' and '-PRO_' in value:
            value = value.split('-PRO_', 1)[0]
        mapped = _IDENTIFIER_CV_MAPPING.get(prefix)
        if mapped is not None and value and value != '-':
            pairs.append((mapped, value))
    return pairs


def parsed_identifier_terms(column_name: str):
    return lambda row: [term for term, _ in _parse_identifier_pairs(row.get(column_name))]


def parsed_identifier_values(column_name: str):
    return lambda row: [value for _, value in _parse_identifier_pairs(row.get(column_name))]


def _intact_raw(opener, organism: int = 9606, **_kwargs: object):
    if organism != 9606:
        raise ValueError('Currently only human (9606) is supported for IntAct')
    if not opener or not opener.result:
        return
    for file_handle in opener.result.values():
        reader = csv.DictReader(file_handle, delimiter='\t')
        yield from reader
        break


config = ResourceConfig(
    id=ResourceCv.INTACT,
    name='IntAct',
    url='https://www.ebi.ac.uk/intact/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='37953288',
    primary_category='interactions',
    description=(
        'IntAct provides a freely available, open source database system '
        'and analysis tools for molecular interaction data. All interactions '
        'are derived from literature curation or direct user submissions and '
        'are freely available in PSI-MITAB format. The database includes '
        'protein-protein, protein-small molecule and protein-nucleic acid '
        'interactions with detailed experimental evidence.'
    ),
)

interactions_schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    identifiers=IdentifiersBuilder(
        CV(
            term=IdentifierNamespaceCv.INTACT,
            value=f('Interaction identifier(s)', extract='intact'),
        ),
    ),
    annotations=AnnotationsBuilder(
        CV(term=f('Interaction type(s)', extract='mi')),
        CV(term=f('Interaction detection method(s)', extract='mi')),
        CV(term=f('Source database(s)', extract='mi')),
        CV(term=InteractionMetadataCv.CONFIDENCE_VALUE, value=f('Confidence value(s)')),
        CV(term=InteractionMetadataCv.EXPANSION_METHOD, value=f('Expansion method(s)')),
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('Publication Identifier(s)', extract='pubmed')),
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('Host organism(s)', extract='tax')),
        CV(term=InteractionMetadataCv.INTERACTION_PARAMETER, value=f('Interaction parameter(s)')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('Interaction annotation(s)')),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=parsed_identifier_terms('#ID(s) interactor A'), value=parsed_identifier_values('#ID(s) interactor A')),
                    CV(term=parsed_identifier_terms('Alt. ID(s) interactor A'), value=parsed_identifier_values('Alt. ID(s) interactor A')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('Taxid interactor A', extract='tax')),
                    CV(term=ParticipantMetadataCv.ALIAS, value=f('Alias(es) interactor A')),
                    CV(term=ParticipantMetadataCv.PARTICIPANT_ANNOTATION, value=f('Annotation(s) interactor A')),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=f('Biological role(s) interactor A', extract='mi')),
                CV(term=f('Experimental role(s) interactor A', extract='mi')),
                CV(term=f('Type(s) interactor A', extract='mi')),
                CV(term=ParticipantMetadataCv.PARTICIPANT_FEATURE, value=f('Feature(s) interactor A')),
                CV(term=ParticipantMetadataCv.STOICHIOMETRY, value=f('Stoichiometry(s) interactor A')),
                CV(term=f('Identification method participant A', extract='mi')),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=parsed_identifier_terms('ID(s) interactor B'), value=parsed_identifier_values('ID(s) interactor B')),
                    CV(term=parsed_identifier_terms('Alt. ID(s) interactor B'), value=parsed_identifier_values('Alt. ID(s) interactor B')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('Taxid interactor B', extract='tax')),
                    CV(term=ParticipantMetadataCv.ALIAS, value=f('Alias(es) interactor B')),
                    CV(term=ParticipantMetadataCv.PARTICIPANT_ANNOTATION, value=f('Annotation(s) interactor B')),
                ),
            ),
            annotations=AnnotationsBuilder(
                CV(term=f('Biological role(s) interactor B', extract='mi')),
                CV(term=f('Experimental role(s) interactor B', extract='mi')),
                CV(term=f('Type(s) interactor B', extract='mi')),
                CV(term=ParticipantMetadataCv.PARTICIPANT_FEATURE, value=f('Feature(s) interactor B')),
                CV(term=ParticipantMetadataCv.STOICHIOMETRY, value=f('Stoichiometry(s) interactor B')),
                CV(term=f('Identification method participant B', extract='mi')),
            ),
        ),
    ),
)

resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            url='https://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/species/human.zip',
            filename='human.zip',
            subfolder='intact',
            large=True,
            ext='zip',
        ),
        mapper=interactions_schema,
        raw_parser=_intact_raw,
    ),
)
