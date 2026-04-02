"""
Parse UniProt data and emit Entity records.

This module converts UniProt protein data into Entity records using the schema
defined in pypath.internals.silver_schema. It also exposes narrow row-oriented
translation datasets for reference identifiers and secondary accessions.
"""

from __future__ import annotations

from collections.abc import Generator, Iterable
import re
from urllib.parse import quote_plus

from pypath.inputs import uniprot as uniprot_input

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    MoleculeAnnotationsCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceCv,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.base import iter_tsv

# UniProt REST API URL for comprehensive protein data
# Currently hardcoded for human (9606), mouse (10090), and rat (10116)
UNIPROT_DATA_URL = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?compressed=true"
    "&format=tsv"
    "&query=(taxonomy_id:9606 OR taxonomy_id:10090 OR taxonomy_id:10116) AND reviewed:true"
    "&fields=accession,id,protein_name,length,mass,sequence,gene_primary,gene_synonym,"
    "organism_id,cc_disease,ft_mutagen,cc_subcellular_location,cc_ptm,lit_pubmed_id,"
    "cc_function,xref_ensembl,xref_geneid,xref_kegg,cc_pathway,cc_activity_regulation,keywordid,"
    "ec,go_id,ft_transmem,protein_families,xref_refseq,xref_alphafolddb,"
    "xref_chembl,xref_phosphositeplus,xref_signor,xref_pathwaycommons,"
    "xref_biogrid,xref_complexportal"
)

config = ResourceConfig(
    id=ResourceCv.UNIPROT,
    name='UniProt',
    url='https://www.uniprot.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='33237286',
    description=(
        'UniProt is a comprehensive resource for protein sequence and '
        'functional information. It provides high-quality, manually annotated '
        'protein data including function, structure, localization, interactions, '
        'disease associations, and cross-references to other databases.'
    ),
)

f = FieldConfig(
    extract={
        'protein_name': r'^([^(]+)',
        'protein_synonym': r'^([^)]+)\)',
        'ensembl_id': r'^(ENS[A-Z0-9]*\d+(?:\.\d+)?)',
        'complexportal_id': r'^(CPX-\d+)',
    },
)

PROTEIN_REFERENCE_KEY_TYPES: tuple[str, ...] = (
    'MI:1097:Uniprot',
    'OM:0221:Uniprot Entry Name',
    'OM:0200:Gene Name Primary',
    'OM:0201:Gene Name Synonym',
    'MI:0477:Entrez',
    'MI:0476:Ensembl',
)

_REFERENCE_FIELDS = (
    'accession',
    'id',
    'gene_primary',
    'gene_synonym',
    'organism_id',
    'xref_ensembl',
    'xref_geneid',
)

_UNIPROT_ENSEMBL_RE = re.compile(r'(ENS[A-Z0-9]*\d+(?:\.\d+)?)')

protein_name_column = f('Protein names', extract='protein_name')
protein_synonym_column = f('Protein names', delimiter='(', extract='protein_synonym')

def _build_uniprot_reference_query(taxonomy_ids: Iterable[int | str] | None = None) -> str:
    if taxonomy_ids is None:
        return 'reviewed:true'

    tax_terms = [f'taxonomy_id:{int(t)}' for t in taxonomy_ids]
    if not tax_terms:
        return 'reviewed:true'
    if len(tax_terms) == 1:
        return f'{tax_terms[0]} AND reviewed:true'
    return f"({' OR '.join(tax_terms)}) AND reviewed:true"


def _build_uniprot_reference_url(taxonomy_ids: Iterable[int | str] | None = None, **_kwargs: object) -> str:
    query = quote_plus(_build_uniprot_reference_query(taxonomy_ids))
    fields = ','.join(_REFERENCE_FIELDS)
    return (
        'https://rest.uniprot.org/uniprotkb/stream'
        f'?compressed=true&format=tsv&query={query}&fields={fields}'
    )


def _reference_filename(taxonomy_ids: Iterable[int | str] | None = None, **_kwargs: object) -> str:
    if taxonomy_ids is None:
        suffix = 'all_reviewed_swissprot'
    else:
        values = [str(int(t)) for t in taxonomy_ids]
        suffix = '_'.join(values) if values else 'all_reviewed_swissprot'
    return f'uniprot_reference_{suffix}.tsv.gz'


def _split_semicolon_field(value: str | None) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in value.split(';') if item.strip()]


def _split_gene_synonym_field(value: str | None) -> list[str]:
    if not value:
        return []
    return [item.strip() for item in value.split() if item.strip()]


def _extract_ensembl_ids(value: str | None) -> list[str]:
    if not value:
        return []
    return sorted({match.group(1) for match in _UNIPROT_ENSEMBL_RE.finditer(value)})


def _reference_id_translation_raw(opener, max_records: int | None = None, **kwargs):
    emitted = 0
    for rec in iter_tsv(opener, **kwargs):
        primary_uniprot = (rec.get('Entry') or '').strip()
        if not primary_uniprot:
            continue

        entry_name = (rec.get('Entry Name') or '').strip()
        taxonomy_id = (rec.get('Organism (ID)') or '').strip() or None
        gene_primary = (rec.get('Gene Names (primary)') or '').strip()

        rows = [
            {
                'key_type': 'MI:1097:Uniprot',
                'key_value': primary_uniprot,
                'taxonomy_id': None,
                'primary_uniprot': primary_uniprot,
            },
        ]

        if entry_name:
            rows.append({
                'key_type': 'OM:0221:Uniprot Entry Name',
                'key_value': entry_name,
                'taxonomy_id': None,
                'primary_uniprot': primary_uniprot,
            })
        if gene_primary:
            rows.append({
                'key_type': 'OM:0200:Gene Name Primary',
                'key_value': gene_primary,
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
            for synonym in _split_gene_synonym_field(rec.get('Gene Names (synonym)'))
        )
        rows.extend(
            {
                'key_type': 'MI:0477:Entrez',
                'key_value': entrez,
                'taxonomy_id': taxonomy_id,
                'primary_uniprot': primary_uniprot,
            }
            for entrez in _split_semicolon_field(rec.get('GeneID'))
        )
        rows.extend(
            {
                'key_type': 'MI:0476:Ensembl',
                'key_value': ensembl,
                'taxonomy_id': taxonomy_id,
                'primary_uniprot': primary_uniprot,
            }
            for ensembl in _extract_ensembl_ids(rec.get('Ensembl'))
        )

        for row in rows:
            yield row
            emitted += 1
            if max_records is not None and emitted >= max_records:
                return


def _secondary_to_primary_raw(_opener, max_records: int | None = None, **_kwargs):
    emitted = 0
    for secondary, primary in uniprot_input.get_uniprot_sec(organism=None):
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
        CV(term=IdentifierNamespaceCv.UNIPROT_ENTRY_NAME, value=f('Entry Name')),
        CV(term=IdentifierNamespaceCv.GENE_NAME_PRIMARY, value=f('Gene Names (primary)')),
        CV(
            term=IdentifierNamespaceCv.GENE_NAME_SYNONYM,
            value=f('Gene Names (synonym)', delimiter=' '),
        ),
        CV(
            term=IdentifierNamespaceCv.NAME,
            value=protein_name_column,
        ),
        CV(
            term=IdentifierNamespaceCv.SYNONYM,
            value=protein_synonym_column,
        ),
        CV(term=IdentifierNamespaceCv.ENSEMBL, value=f('Ensembl', delimiter=';', extract='ensembl_id')),
        CV(term=IdentifierNamespaceCv.ENTREZ, value=f('GeneID', delimiter=';')),
        #CV(term=IdentifierNamespaceCv.REFSEQ, value=f('RefSeq', delimiter=';')),
        #CV(term=IdentifierNamespaceCv.ALPHAFOLDDB, value=f('AlphaFoldDB', delimiter=';')),
        CV(term=IdentifierNamespaceCv.KEGG, value=f('KEGG', delimiter=';')),
        CV(term=IdentifierNamespaceCv.CHEMBL, value=f('ChEMBL', delimiter=';')),
        CV(term=IdentifierNamespaceCv.SIGNOR, value=f('SIGNOR', delimiter=';')),
        CV(term=IdentifierNamespaceCv.BIOGRID, value=f('BioGRID', delimiter=';')),
        CV(term=IdentifierNamespaceCv.COMPLEXPORTAL, value=f('ComplexPortal', delimiter=';', extract='complexportal_id')),
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
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('Gene Ontology IDs', delimiter=';')),
        CV(term=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=f('Keywords IDs', delimiter=';')),
        CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=f('Organism (ID)')),
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('PubMed ID', delimiter=';')),
    ),
)

resource = Resource(
    config,
    proteins=Dataset(
        download=Download(
            url=UNIPROT_DATA_URL,
            filename='uniprot_proteins_9606_10090_10116.tsv.gz',
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
            url=_build_uniprot_reference_url,
            filename=_reference_filename,
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
