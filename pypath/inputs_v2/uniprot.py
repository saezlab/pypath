"""Parse UniProt data and emit Entity records.

This module converts UniProt protein data into Entity records using the schema
defined in pypath.internals.silver_schema. It also exposes narrow row-oriented
translation datasets for reference identifiers and secondary accessions.
"""

from __future__ import annotations

from collections.abc import Iterable
import re

from pypath.inputs_v2.base import (
    Dataset,
    Download,
    Resource,
    ResourceConfig,
    ontology_entity_mapper,
)
from pypath.inputs_v2.parsers.base import _first_handle, iter_tsv
from pypath.inputs_v2.parsers.obo import iter_obo, obo_record_to_term
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    MoleculeAnnotationsCv,
    OntologyCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.tabular_builder import (
    AssociationBuilder,
    AssociationsBuilder,
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.resources import urls
from pypath.share.downloads import dm

# UniProt REST API URL for protein data.
# Currently hardcoded for human (9606), mouse (10090), and rat (10116).
# Cross-reference identifiers are loaded through reference_id_translation and
# restored on canonical entities from the resolver, so protein evidence rows
# only need the primary UniProt accession.
UNIPROT_DATA_URL = (
    'https://rest.uniprot.org/uniprotkb/stream'
    '?compressed=true'
    '&format=tsv'
    '&query=(taxonomy_id:9606 OR taxonomy_id:10090 OR taxonomy_id:10116) AND reviewed:true'
    '&fields=accession,length,mass,sequence,organism_id,cc_disease,ft_mutagen,'
    'cc_subcellular_location,cc_ptm,lit_pubmed_id,cc_function,cc_pathway,'
    'cc_activity_regulation,keywordid,ec,go_id,ft_transmem,protein_families'
)

UNIPROT_KEYWORDS_OBO_URL = (
    'https://rest.uniprot.org/keywords/stream'
    '?compressed=true&format=obo&query=%28*%29'
)
UNIPROT_SPROT_FLATFILE_URL = (
    'https://ftp.uniprot.org/pub/databases/uniprot/current_release/'
    'knowledgebase/complete/uniprot_sprot.dat.gz'
)

config = ResourceConfig(
    id=ResourceCv.UNIPROT,
    name='UniProt',
    url='https://www.uniprot.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='33237286',
    primary_category='proteins',
    annotation_ontologies=(OntologyCv.GENE_ONTOLOGY, OntologyCv.UNIPROT_KEYWORDS),
    description=(
        'UniProt is a comprehensive resource for protein sequence and '
        'functional information. It provides high-quality, manually annotated '
        'protein data including function, structure, localization, interactions, '
        'disease associations, and cross-references to other databases.'
    ),
)

f = FieldConfig()
_UNIPROT_KEYWORDS_ONTOLOGY_ID = 'uniprot_keywords'

PROTEIN_REFERENCE_KEY_TYPES: tuple[str, ...] = (
    'MI:1097:Uniprot',
    'OM:0221:Uniprot Entry Name',
    'OM:0200:Gene Name Primary',
    'OM:0201:Gene Name Synonym',
    'MI:0477:Entrez',
    'MI:0476:Ensembl',
    'MI:1095:HGNC',
)

_GN_NAME_RE = re.compile(r'(?:^|; )Name=([^;{]+)')
_GN_SYNONYMS_RE = re.compile(r'(?:^|; )Synonyms=([^;]+)')
_OX_TAXONOMY_RE = re.compile(r'NCBI_TaxID=(\d+)')


def _split_semicolon_field(value: str | None) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in value.split(';') if item.strip()]


def _split_gene_synonym_field(value: str | None) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in value.split() if item.strip()]


def _versioned_and_unversioned(identifier: str | None) -> tuple[str, ...]:
    if not identifier:
        return ()

    identifier = identifier.strip()
    if not identifier:
        return ()
    if '.' not in identifier:
        return (identifier,)

    return (identifier, identifier.split('.', maxsplit=1)[0])


def _extract_hgnc_ids(value: str | None) -> list[str]:
    if not value:
        return []
    hgnc_ids = []
    for item in _split_semicolon_field(value):
        match = re.match(r'^(?:HGNC:)?(\d+)$', item.strip())
        if match:
            hgnc_ids.append(match.group(1))
    return sorted(set(hgnc_ids))


def _iter_uniprot_flatfile_records(opener):
    source = _first_handle(opener)
    if not source:
        return

    record: list[str] = []
    for line in source:
        if line.startswith('//'):
            if record:
                yield record
                record = []
            continue

        record.append(line.rstrip('\n'))

    if record:
        yield record


def _parse_uniprot_flatfile_record(lines: list[str]) -> dict[str, object]:
    record: dict[str, object] = {
        'primary_uniprot': None,
        'entry_name': None,
        'taxonomy_id': None,
        'gene_primary': None,
        'gene_synonyms': [],
        'ensembl': [],
        'entrez': [],
        'hgnc': [],
    }
    gene_text = []

    for line in lines:
        if line.startswith('ID   '):
            parts = line[5:].split()
            record['entry_name'] = parts[0] if parts else None
        elif line.startswith('AC   '):
            accessions = _split_semicolon_field(line[5:])
            if accessions and record['primary_uniprot'] is None:
                record['primary_uniprot'] = accessions[0]
        elif line.startswith('OX   '):
            match = _OX_TAXONOMY_RE.search(line)
            if match:
                record['taxonomy_id'] = match.group(1)
        elif line.startswith('GN   '):
            gene_text.append(line[5:].strip())
        elif line.startswith('DR   Ensembl;'):
            record['ensembl'].extend(_parse_ensembl_dr_line(line))
        elif line.startswith('DR   GeneID;'):
            parts = _parse_dr_parts(line)
            if parts:
                record['entrez'].append(parts[0])
        elif line.startswith('DR   HGNC;'):
            parts = _parse_dr_parts(line)
            if parts:
                record['hgnc'].extend(_extract_hgnc_ids(parts[0]))

    _parse_gene_text(' '.join(gene_text), record)

    return record


def _parse_dr_parts(line: str) -> list[str]:
    value = line[5:].rstrip('.')
    parts = [part.strip() for part in value.split(';')]
    return [part for part in parts[1:] if part and part != '-']


def _parse_ensembl_dr_line(line: str) -> list[str]:
    identifiers: set[str] = set()

    for part in _parse_dr_parts(line):
        identifiers.update(_versioned_and_unversioned(part))

    return sorted(identifiers)


def _parse_gene_text(text: str, record: dict[str, object]) -> None:
    if not text:
        return

    name_match = _GN_NAME_RE.search(text)
    if name_match:
        record['gene_primary'] = _strip_flatfile_evidence(name_match.group(1))

    synonyms_match = _GN_SYNONYMS_RE.search(text)
    if synonyms_match:
        record['gene_synonyms'] = [
            _strip_flatfile_evidence(value)
            for value in synonyms_match.group(1).split(',')
            if _strip_flatfile_evidence(value)
        ]


def _strip_flatfile_evidence(value: str) -> str:
    return re.sub(r'\s*\{[^}]*\}', '', value).strip()


def _reference_id_translation_raw(opener, max_records: int | None = None, **kwargs):
    emitted = 0
    taxonomy_ids = kwargs.get('taxonomy_ids')
    taxonomy_ids = {str(int(t)) for t in taxonomy_ids} if taxonomy_ids else None

    for lines in _iter_uniprot_flatfile_records(opener):
        rec = _parse_uniprot_flatfile_record(lines)
        primary_uniprot = rec.get('primary_uniprot')
        if not primary_uniprot:
            continue
        primary_uniprot = str(primary_uniprot)

        entry_name = rec.get('entry_name')
        taxonomy_id = rec.get('taxonomy_id')
        taxonomy_id = str(taxonomy_id) if taxonomy_id is not None else None
        if taxonomy_ids is not None and taxonomy_id not in taxonomy_ids:
            continue

        rows = [
            {
                'key_type': 'MI:1097:Uniprot',
                'key_value': primary_uniprot,
                'taxonomy_id': taxonomy_id,
                'primary_uniprot': primary_uniprot,
            },
        ]

        if entry_name:
            rows.append({
                'key_type': 'OM:0221:Uniprot Entry Name',
                'key_value': str(entry_name),
                'taxonomy_id': taxonomy_id,
                'primary_uniprot': primary_uniprot,
            })

        if rec.get('gene_primary'):
            rows.append({
                'key_type': 'OM:0200:Gene Name Primary',
                'key_value': rec.get('gene_primary'),
                'taxonomy_id': taxonomy_id,
                'primary_uniprot': primary_uniprot,
            })

        rows.extend(
            {
                'key_type': 'OM:0201:Gene Name Synonym',
                'key_value': synonym,
                'taxonomy_id': taxonomy_id,
                'primary_uniprot': primary_uniprot,
            }
            for synonym in rec.get('gene_synonyms', [])
        )
        rows.extend(
            {
                'key_type': 'MI:0477:Entrez',
                'key_value': entrez,
                'taxonomy_id': taxonomy_id,
                'primary_uniprot': primary_uniprot,
            }
            for entrez in rec.get('entrez', [])
        )
        rows.extend(
            {
                'key_type': 'MI:0476:Ensembl',
                'key_value': ensembl,
                'taxonomy_id': taxonomy_id,
                'primary_uniprot': primary_uniprot,
            }
            for ensembl in rec.get('ensembl', [])
        )
        rows.extend(
            {
                'key_type': 'MI:1095:HGNC',
                'key_value': hgnc,
                'taxonomy_id': taxonomy_id,
                'primary_uniprot': primary_uniprot,
            }
            for hgnc in rec.get('hgnc', [])
        )

        for row in rows:
            yield row
            emitted += 1
            if max_records is not None and emitted >= max_records:
                return


def _secondary_to_primary_raw(_opener, max_records: int | None = None, **_kwargs):
    emitted = 0
    path = dm.download(urls.urls['uniprot_sec']['url'])

    with open(path) as handle:
        for i, line in enumerate(handle):
            if i < 30:
                continue

            parts = line.split()
            if len(parts) != 2:
                continue

            secondary, primary = parts
            if set(secondary) == {'_'} or set(primary) == {'_'}:
                continue

            yield {
                'secondary_uniprot': secondary,
                'primary_uniprot': primary,
            }
            emitted += 1
            if max_records is not None and emitted >= max_records:
                break


proteins_schema = EntityBuilder(
    entity_type=EntityTypeCv.PROTEIN,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.UNIPROT, value=f('Entry')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=MoleculeAnnotationsCv.SEQUENCE_LENGTH, value=f('Length')),
        CV(term=MoleculeAnnotationsCv.MASS_DALTON, value=f('Mass')),
        CV(term=MoleculeAnnotationsCv.FUNCTION, value=f('Function [CC]')),
        CV(term=MoleculeAnnotationsCv.SUBCELLULAR_LOCATION, value=f('Subcellular location [CC]')),
        CV(term=MoleculeAnnotationsCv.POST_TRANSLATIONAL_MODIFICATION, value=f('Post-translational modification')),
        CV(term=MoleculeAnnotationsCv.DISEASE_INVOLVEMENT, value=f('Involvement in disease')),
        CV(term=MoleculeAnnotationsCv.PATHWAY_PARTICIPATION, value=f('Pathway')),
        CV(term=MoleculeAnnotationsCv.ACTIVITY_REGULATION, value=f('Activity regulation')),
        CV(term=MoleculeAnnotationsCv.MUTAGENESIS, value=f('Mutagenesis')),
        CV(term=MoleculeAnnotationsCv.TRANSMEMBRANE_REGION, value=f('Transmembrane')),
        CV(term=MoleculeAnnotationsCv.PROTEIN_FAMILY, value=f('Protein families', delimiter=',')),
        CV(term=MoleculeAnnotationsCv.EC_NUMBER, value=f('EC number', delimiter=';')),
        CV(term=MoleculeAnnotationsCv.AMINO_ACID_SEQUENCE, value=f('Sequence')),
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('Organism (ID)')),
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('PubMed ID', delimiter=';')),
    ),
    associations=AssociationsBuilder(
        AssociationBuilder(
            object_entity_type=EntityTypeCv.CV_TERM,
            object_identifier_type=IdentifierNamespaceCv.CV_TERM_ACCESSION,
            object_identifier=f('Gene Ontology IDs', delimiter=';'),
        ),
        AssociationBuilder(
            object_entity_type=EntityTypeCv.CV_TERM,
            object_identifier_type=IdentifierNamespaceCv.CV_TERM_ACCESSION,
            object_identifier=f('Keyword ID', delimiter=';'),
        ),
    ),
)


keyword_terms_schema = ontology_entity_mapper(
    obo_record_to_term,
    ontology_id=_UNIPROT_KEYWORDS_ONTOLOGY_ID,
)

resource = Resource(
    config,
    keywords=Dataset(
        download=Download(
            url=UNIPROT_KEYWORDS_OBO_URL,
            filename='uniprot_keywords.obo.gz',
            subfolder='uniprot',
            large=True,
            encoding='utf-8',
            default_mode='r',
            ext='gz',
        ),
        mapper=keyword_terms_schema,
        raw_parser=iter_obo,
    ),
    proteins=Dataset(
        download=Download(
            url=UNIPROT_DATA_URL,
            filename='uniprot_proteins_slim_9606_10090_10116.tsv.gz',
            subfolder='uniprot',
            large=True,
            encoding='utf-8',
            default_mode='r',
            ext='gz',
        ),
        mapper=proteins_schema,
        raw_parser=iter_tsv,
    ),
    reference_id_translation=Dataset(
        download=Download(
            url=UNIPROT_SPROT_FLATFILE_URL,
            filename='uniprot_sprot.dat.gz',
            subfolder='uniprot',
            large=True,
            encoding='utf-8',
            default_mode='r',
            ext='gz',
        ),
        mapper=lambda row: row,
        raw_parser=_reference_id_translation_raw,
        kind='id_translation',
    ),
    secondary_to_primary=Dataset(
        download=None,
        mapper=lambda row: row,
        raw_parser=_secondary_to_primary_raw,
        kind='id_translation',
    ),
)
