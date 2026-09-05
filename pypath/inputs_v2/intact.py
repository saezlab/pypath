"""
Parse IntAct data and emit Entity records.

This module converts IntAct MITAB data into Entity records using the schema
defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from biolink_model.datamodel import model
from biolink_model.datamodel.model import slots
from omnipath_core.naming import Namespace

import csv
import re

from pypath.internals.cv_terms import (
    LicenseCV,
    OntologyCv,
    UpdateCategoryCV,
    ResourceCv,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig


_IDENTIFIER_CV_MAPPING = {
    'bind smid': Namespace.BIND,
    'cas registry number': Namespace.CAS,
    'chebi': Namespace.CHEBI,
    'chembl': Namespace.CHEMBL,
    'chembl compound': Namespace.CHEMBL,
    'chembl internal id': Namespace.CHEMBL_INTERNAL_ID,
    'ddbj/embl/genbank': Namespace.GENBANK,
    'dip': Namespace.DIP,
    'ensembl': Namespace.ENSEMBL,
    'ensemblgenomes': Namespace.ENSEMBL_GENOMES,
    'entrezgene/locuslink': Namespace.ENTREZ,
    'flybase': Namespace.FLYBASE,
    'genbank identifier': Namespace.GENBANK,
    'genbank_nucl_gi': Namespace.GENBANK_GI,
    'genbank_protein_gi': Namespace.GENBANK_GI,
    'hgnc': Namespace.HGNC,
    'genesymbol': Namespace.GENESYMBOL,
    'imex': Namespace.IMEX,
    'intact': Namespace.INTACT,
    'ipi': Namespace.IPI,
    'mint': Namespace.MINT,
    'mirbase': Namespace.MIRBASE,
    'psi-mi': Namespace.MI,
    'refseq': Namespace.REFSEQ,
    'rfam': Namespace.RFAM,
    'rnacentral': Namespace.RNACENTRAL,
    'uniparc': Namespace.UNIPARC,
    'uniprotkb': Namespace.UNIPROT,
}

_INTERACTOR_TYPE_MAPPING = {
    'MI:0326': model.Protein,
    'MI:0314': model.MacromolecularComplex,
    'MI:0328': model.ChemicalEntity,
    'MI:0320': model.RNAProduct,
    'MI:0250': model.Gene,
    'MI:0681': model.NucleicAcidEntity,
    'MI:0318': model.NucleicAcidEntity,
    'MI:0317': model.MolecularEntity,
}

_ENSEMBL_RE = re.compile(r'^ENS[A-Z0-9]*\d+(?:\.\d+)?$')
_REFSEQ_RE = re.compile(r'^[A-Z]{2}_[0-9]+(?:\.\d+)?$')


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
    if value.startswith('MINT-') and prefix in {'intact', 'psi-mi'}:
        return 'mint'
    if prefix == 'intact' and not (
        value.startswith('EBI-') or value.startswith('IM-')
    ):
        return None
    if (
        prefix == 'hgnc'
        and not value.isdigit()
        and not value.upper().startswith('HGNC:')
    ):
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
        if prefix == 'chebi':
            match = re.fullmatch(
                r'(?:CHEBI:)?(\d+)', value, flags=re.IGNORECASE
            )
            value = match.group(1) if match else None
        if prefix == 'chembl' and value.upper().startswith('CHEMBL'):
            value = f'CHEMBL{value[len("CHEMBL") :]}'
        if prefix == 'chembl compound' and value.isdigit():
            prefix = 'chembl internal id'
        if prefix == 'refseq' and not _REFSEQ_RE.fullmatch(value):
            continue
        if prefix == 'entrezgene/locuslink' and not value.isdigit():
            continue
        if prefix == 'ensembl' and not _ENSEMBL_RE.fullmatch(value):
            continue
        if (
            prefix in {'genbank_nucl_gi', 'genbank_protein_gi'}
            and not value.isdigit()
        ):
            continue
        mapped = _IDENTIFIER_CV_MAPPING.get(prefix)
        if mapped is not None and value and value != '-':
            pairs.append((mapped, value))
    return pairs


def parsed_identifier_terms(column_name: str):
    return lambda row: [
        term for term, _ in _parse_identifier_pairs(row.get(column_name))
    ]


def parsed_identifier_values(column_name: str):
    return lambda row: [
        value for _, value in _parse_identifier_pairs(row.get(column_name))
    ]


def _interactor_entity_type(suffix: str):
    column = f'Type(s) interactor {suffix}'

    def _value(row: dict[str, object]) -> type:
        match = re.search('(MI:\\d+)', str(row.get(column) or ''))
        if match:
            return _INTERACTOR_TYPE_MAPPING.get(
                match.group(1), model.PhysicalEntity
            )
        return model.PhysicalEntity

    return _value


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
    pubmed='34761267',
    annotation_ontologies=(OntologyCv.PSI_MI,),
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


def intact_predicate(row):
    # MITAB describes an interaction; arbitrary PSI-MI types are preserved below.
    return slots.interacts_with


interactor_a_builder = EntityBuilder(
    entity_type=_interactor_entity_type('A'),
    identifiers=IdentifiersBuilder(
        CV(
            term=parsed_identifier_terms('#ID(s) interactor A'),
            value=parsed_identifier_values('#ID(s) interactor A'),
        ),
        CV(
            term=parsed_identifier_terms('Alt. ID(s) interactor A'),
            value=parsed_identifier_values('Alt. ID(s) interactor A'),
        ),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.in_taxon,
            value=f(
                'Taxid interactor A',
                extract='tax',
                transform=lambda v: 'NCBITaxon:'
                + str(v).removeprefix('NCBITaxon:')
                if str(v).removeprefix('NCBITaxon:').isdigit()
                and int(str(v).removeprefix('NCBITaxon:')) > 0
                else None,
            ),
        ),
        CV(
            term=slots.has_topic,
            value=f('Biological role(s) interactor A', extract='mi'),
        ),
        CV(
            term=slots.has_topic,
            value=f('Experimental role(s) interactor A', extract='mi'),
        ),
        CV(term=slots.stoichiometry, value=f('Stoichiometry(s) interactor A')),
        CV(
            term=slots.has_topic,
            value=f('Identification method participant A', extract='mi'),
        ),
    ),
)

interactor_b_builder = EntityBuilder(
    entity_type=_interactor_entity_type('B'),
    identifiers=IdentifiersBuilder(
        CV(
            term=parsed_identifier_terms('ID(s) interactor B'),
            value=parsed_identifier_values('ID(s) interactor B'),
        ),
        CV(
            term=parsed_identifier_terms('Alt. ID(s) interactor B'),
            value=parsed_identifier_values('Alt. ID(s) interactor B'),
        ),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.in_taxon,
            value=f(
                'Taxid interactor B',
                extract='tax',
                transform=lambda v: 'NCBITaxon:'
                + str(v).removeprefix('NCBITaxon:')
                if str(v).removeprefix('NCBITaxon:').isdigit()
                and int(str(v).removeprefix('NCBITaxon:')) > 0
                else None,
            ),
        ),
        CV(
            term=slots.has_topic,
            value=f('Biological role(s) interactor B', extract='mi'),
        ),
        CV(
            term=slots.has_topic,
            value=f('Experimental role(s) interactor B', extract='mi'),
        ),
        CV(term=slots.stoichiometry, value=f('Stoichiometry(s) interactor B')),
        CV(
            term=slots.has_topic,
            value=f('Identification method participant B', extract='mi'),
        ),
    ),
)

interactions_schema = RelationBuilder(
    subject=interactor_a_builder,
    predicate=intact_predicate,
    object=interactor_b_builder,
    identifiers=IdentifiersBuilder(
        CV(
            term=Namespace.INTACT,
            value=f('Interaction identifier(s)', extract='intact'),
        )
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.intact_confidence_value,
            value=f(
                'Confidence value(s)',
                extract=r'(?:^|\|)intact-miscore:([0-9.]+)',
            ),
        ),
        CV(term=slots.has_topic, value=f('Interaction type(s)', extract='mi')),
        CV(
            term=slots.has_topic,
            value=f('Interaction detection method(s)', extract='mi'),
        ),
        CV(term=slots.has_topic, value=f('Source database(s)', extract='mi')),
        CV(
            term=slots.publications,
            value=f(
                'Publication Identifier(s)',
                extract='pubmed',
                transform=lambda v: 'PMID:' + str(v).removeprefix('PMID:')
                if str(v).removeprefix('PMID:').isdigit()
                and int(str(v).removeprefix('PMID:')) > 0
                else None,
            ),
        ),
        CV(
            term=slots.in_taxon,
            value=f(
                'Host organism(s)',
                extract='tax',
                transform=lambda v: 'NCBITaxon:'
                + str(v).removeprefix('NCBITaxon:')
                if str(v).removeprefix('NCBITaxon:').isdigit()
                and int(str(v).removeprefix('NCBITaxon:')) > 0
                else None,
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
