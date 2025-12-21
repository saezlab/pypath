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
    'imex': IdentifierNamespaceCv.IMEX,
    'intact': IdentifierNamespaceCv.INTACT,
    'ipi': IdentifierNamespaceCv.IPI,
    'mint': IdentifierNamespaceCv.MINT,
    'mirbase': IdentifierNamespaceCv.MIRBASE,
    'pdbe': IdentifierNamespaceCv.PDB,
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
        'tax': r'taxid:([-\d]+)',
        'pubmed': r'(?i)pubmed:(\d+)',
        'intact': r'intact:([^|"]+)',
    },
    map={
        'identifier_cv': _IDENTIFIER_CV_MAPPING,
    },
    delimiter='|',
)


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
        CV(term=InteractionMetadataCv.INTERACTION_XREF, value=f('Interaction Xref(s)')),
        CV(term=InteractionMetadataCv.INTERACTION_ANNOTATION, value=f('Interaction annotation(s)')),
        CV(term=InteractionMetadataCv.INTERACTION_CHECKSUM, value=f('Interaction Checksum(s)')),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=f('#ID(s) interactor A', extract='prefix_lower', map='identifier_cv'), value=f('#ID(s) interactor A', extract='value')),
                    CV(term=f('Alt. ID(s) interactor A', extract='prefix_lower', map='identifier_cv'), value=f('Alt. ID(s) interactor A', extract='value')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('Taxid interactor A', extract='tax')),
                    CV(term=ParticipantMetadataCv.ALIAS, value=f('Alias(es) interactor A')),
                    CV(term=ParticipantMetadataCv.PARTICIPANT_XREF, value=f('Xref(s) interactor A')),
                    CV(term=ParticipantMetadataCv.PARTICIPANT_ANNOTATION, value=f('Annotation(s) interactor A')),
                    CV(term=ParticipantMetadataCv.PARTICIPANT_CHECKSUM, value=f('Checksum(s) interactor A')),
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
                    CV(term=f('ID(s) interactor B', extract='prefix_lower', map='identifier_cv'), value=f('ID(s) interactor B', extract='value')),
                    CV(term=f('Alt. ID(s) interactor B', extract='prefix_lower', map='identifier_cv'), value=f('Alt. ID(s) interactor B', extract='value')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('Taxid interactor B', extract='tax')),
                    CV(term=ParticipantMetadataCv.ALIAS, value=f('Alias(es) interactor B')),
                    CV(term=ParticipantMetadataCv.PARTICIPANT_XREF, value=f('Xref(s) interactor B')),
                    CV(term=ParticipantMetadataCv.PARTICIPANT_ANNOTATION, value=f('Annotation(s) interactor B')),
                    CV(term=ParticipantMetadataCv.PARTICIPANT_CHECKSUM, value=f('Checksum(s) interactor B')),
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

interactions = resource.interactions
