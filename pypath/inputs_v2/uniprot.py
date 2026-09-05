"""Parse UniProt data and emit Entity records.

This module converts UniProt protein data into Entity records using the schema
defined in pypath.internals.silver_schema. It also exposes narrow row-oriented
translation datasets for reference identifiers and secondary accessions.
"""

from __future__ import annotations

from collections.abc import Iterable
import re

from biolink_model.datamodel.model import OntologyClass, Protein, slots
from omnipath_core.naming import Namespace

from pypath.internals.silver_schema import EntityRef, OntologyRelation
from pypath.inputs_v2.base import (
    Dataset,
    Download,
    Resource,
    ResourceConfig,
)
from pypath.inputs_v2.parsers.base import _first_handle, iter_tsv
from pypath.inputs_v2.parsers.obo import iter_obo
from pypath.internals.cv_terms import (
    LicenseCV,
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
    annotation_ontologies=(
        OntologyCv.GENE_ONTOLOGY,
        OntologyCv.UNIPROT_KEYWORDS,
    ),
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
    Namespace.UNIPROT,
    Namespace.UNIPROT_ENTRY,
    Namespace.GENESYMBOL,
    Namespace.GENESYMBOL_SYN,
    Namespace.ENTREZ,
    Namespace.ENST,
    Namespace.ENSP,
    Namespace.ENSG,
    Namespace.HGNC,
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


def _parse_ensembl_dr_line(line: str) -> list[tuple[Namespace, str]]:
    identifiers: set[tuple[Namespace, str]] = set()

    parts = [part.strip() for part in line[5:].rstrip('.').split(';')[1:4]]
    for namespace, part in zip(
        (Namespace.ENST, Namespace.ENSP, Namespace.ENSG), parts
    ):
        if not part or part == '-':
            continue
        identifier = part.split(' [', 1)[0]
        identifiers.update(
            (namespace, value)
            for value in _versioned_and_unversioned(identifier)
        )

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


def _reference_id_translation_raw(
    opener, max_records: int | None = None, **kwargs
):
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
                'key_type': Namespace.UNIPROT,
                'key_value': primary_uniprot,
                'taxonomy_id': taxonomy_id,
                'primary_uniprot': primary_uniprot,
            },
        ]

        if entry_name:
            rows.append(
                {
                    'key_type': Namespace.UNIPROT_ENTRY,
                    'key_value': str(entry_name),
                    'taxonomy_id': taxonomy_id,
                    'primary_uniprot': primary_uniprot,
                }
            )

        if rec.get('gene_primary'):
            rows.append(
                {
                    'key_type': Namespace.GENESYMBOL,
                    'key_value': rec.get('gene_primary'),
                    'taxonomy_id': taxonomy_id,
                    'primary_uniprot': primary_uniprot,
                }
            )

        rows.extend(
            {
                'key_type': Namespace.GENESYMBOL_SYN,
                'key_value': synonym,
                'taxonomy_id': taxonomy_id,
                'primary_uniprot': primary_uniprot,
            }
            for synonym in rec.get('gene_synonyms', [])
        )
        rows.extend(
            {
                'key_type': Namespace.ENTREZ,
                'key_value': entrez,
                'taxonomy_id': taxonomy_id,
                'primary_uniprot': primary_uniprot,
            }
            for entrez in rec.get('entrez', [])
        )
        rows.extend(
            {
                'key_type': namespace,
                'key_value': ensembl,
                'taxonomy_id': taxonomy_id,
                'primary_uniprot': primary_uniprot,
            }
            for namespace, ensembl in rec.get('ensembl', [])
        )
        rows.extend(
            {
                'key_type': Namespace.HGNC,
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


def _secondary_to_primary_raw(
    _opener, max_records: int | None = None, **_kwargs
):
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


# Source narrative headings are retained verbatim; none imply graph predicates.
# Numeric sequence length/mass and family classifications stay in raw payloads.
_PROTEIN_DESCRIPTION_FIELDS = (
    'Function [CC]',
    'Subcellular location [CC]',
    'Post-translational modification',
    'Involvement in disease',
    'Pathway',
    'Activity regulation',
    'Mutagenesis',
    'Transmembrane',
)


def _protein_descriptions(row):
    return [
        f'{field}: {row[field]}'
        for field in _PROTEIN_DESCRIPTION_FIELDS
        if row.get(field)
    ]


proteins_schema = EntityBuilder(
    entity_type=Protein,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.UNIPROT, value=f('Entry')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.has_biological_sequence, value=f('Sequence')),
        CV(term=slots.description, value=_protein_descriptions),
        CV(
            term=slots.in_taxon,
            value=lambda row: (
                f'NCBITaxon:{row["Organism (ID)"]}'
                if row.get('Organism (ID)')
                else None
            ),
        ),
        CV(
            term=slots.publications,
            value=lambda row: [
                f'PMID:{value}'
                for value in _split_semicolon_field(row.get('PubMed ID'))
            ],
        ),
        CV(
            term=slots.has_topic,
            value=lambda row: [
                f'EC:{value}'
                for value in _split_semicolon_field(row.get('EC number'))
            ],
        ),
    ),
    associations=AssociationsBuilder(
        AssociationBuilder(
            predicate=slots.associated_with,
            object_entity_type=OntologyClass,
            object_identifier_type=Namespace.GO,
            object_identifier=f('Gene Ontology IDs', delimiter=';'),
        ),
        AssociationBuilder(
            predicate=slots.associated_with,
            object_entity_type=OntologyClass,
            object_identifier_type=Namespace.UNIPROT_KEYWORD,
            object_identifier=f('Keyword ID', delimiter=';'),
        ),
    ),
)


_keyword_builder = EntityBuilder(
    entity_type=OntologyClass,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.UNIPROT_KEYWORD, value=f('id')),
        CV(term=Namespace.NAME, value=f('name')),
        CV(term=Namespace.SYNONYM, value=lambda row: row.get('synonyms', [])),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.description, value=f('definition')),
        CV(
            term=slots.xref,
            value=lambda row: [x.split()[0] for x in row.get('xrefs', [])],
        ),
        CV(
            term=slots.has_topic,
            value=lambda row: [
                f'UniProtKB-KW:{rel["target"]}'
                for rel in row.get('relationships', [])
                if rel['type'] == 'category'
            ],
        ),
    ),
    ontology_relations=lambda row: [
        OntologyRelation(
            predicate=slots.subclass_of,
            object=EntityRef(OntologyClass, Namespace.UNIPROT_KEYWORD, parent),
            ontology_id=_UNIPROT_KEYWORDS_ONTOLOGY_ID,
        )
        for parent in row.get('is_a', [])
    ],
)


def keyword_terms_schema(row):
    """Represent source vocabulary concepts uniformly as ontology classes."""
    if row.get('is_obsolete') or not row.get('id'):
        return None
    return _keyword_builder(row)


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
