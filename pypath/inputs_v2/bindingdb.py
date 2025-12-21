"""
Parse BindingDB data and emit Entity records.

This module converts BindingDB binding affinity data into Entity records using
the schema defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Generator
import csv

from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    InteractionParameterCv,
    AffinityUnitCv,
    CurationCv,
    ResourceCv,
)
from pypath.internals.tabular_builder import (
    AnnotationsBuilder,
    CV,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    Member,
    MembershipBuilder,
)
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig




def _bindingdb_url(dataset: str = 'All', **_kwargs: object) -> str:
    return f'https://bindingdb.org/rwd/bind/downloads/BindingDB_{dataset}_202506_tsv.zip'


def _bindingdb_filename(dataset: str = 'All', **_kwargs: object) -> str:
    return f'{dataset}_table.zip'


def _bindingdb_raw(opener, max_lines: int | None = None, **_kwargs: object):
    if not opener or not opener.result:
        return
    for file_handle in opener.result.values():
        header_line = file_handle.readline().strip()
        header = header_line.split('\t')

        columns_to_keep = min(49, len(header))
        filtered_header = header[:columns_to_keep]

        def filtered_rows():
            for line in file_handle:
                columns = line.strip().split('\t')
                yield '\t'.join(columns[:columns_to_keep])

        reader = csv.DictReader(filtered_rows(), fieldnames=filtered_header, delimiter='\t')

        if max_lines:
            for i, row in enumerate(reader):
                if i >= max_lines:
                    break
                yield row
        else:
            yield from reader
        break


config = ResourceConfig(
    id=ResourceCv.BINDINGDB,
    name='BindingDB',
    url='https://www.bindingdb.org/',
    license=LicenseCV.CC_BY_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='26481362',
    description=(
        'BindingDB is a public, web-accessible database of measured binding '
        'affinities, focusing chiefly on the interactions of proteins considered '
        'to be drug-targets with small, drug-like molecules. It contains binding '
        'data for over 2 million protein-ligand complexes with experimental '
        'measurements including Ki, Kd, IC50, and EC50 values.'
    ),
)

f = FieldConfig(
    extract={
        'chembl': r'^(CHEMBL\d+)$',
        'zinc': r'^(ZINC\d+)$',
        'cas': r'^(\d{2,7}-\d{2}-\d)$',
        'chebi': r'^CHEBI[:\s]?(\d+)$',
        'pubchem_cid': r'^CID[:\s]?(\d+)$',
        'kegg': r'^(C\d{5})$',
        'tax': r'(\d+)',
    },
)

tax_value = f('Target Source Organism According to Curator or DataSource', extract='tax')

interactions_schema = EntityBuilder(
    entity_type=EntityTypeCv.INTERACTION,
    identifiers=IdentifiersBuilder(
        CV(term=IdentifierNamespaceCv.BINDINGDB, value=f('BindingDB Reactant_set_id')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=InteractionParameterCv.KI, value=f('Ki (nM)'), unit=AffinityUnitCv.NANOMOLAR),
        CV(term=InteractionParameterCv.KD, value=f('Kd (nM)'), unit=AffinityUnitCv.NANOMOLAR),
        CV(term=InteractionParameterCv.IC50, value=f('IC50 (nM)'), unit=AffinityUnitCv.NANOMOLAR),
        CV(term=InteractionParameterCv.EC50, value=f('EC50 (nM)'), unit=AffinityUnitCv.NANOMOLAR),
        CV(term=InteractionParameterCv.KON, value=f('kon (M-1-s-1)'), unit=AffinityUnitCv.PER_MOLAR_PER_SECOND),
        CV(term=InteractionParameterCv.KOFF, value=f('koff (s-1)'), unit=AffinityUnitCv.PER_SECOND),
        CV(term=InteractionParameterCv.PH, value=f('pH')),
        CV(
            term=InteractionParameterCv.TEMPERATURE_CELSIUS,
            value=f('Temp (C)'),
            unit=AffinityUnitCv.DEGREE_CELSIUS,
        ),
        CV(term=IdentifierNamespaceCv.PUBMED, value=f('PMID')),
        CV(term=IdentifierNamespaceCv.DOI, value=f('Article DOI')),
        CV(term=IdentifierNamespaceCv.PATENT_NUMBER, value=f('Patent Number')),
        CV(term=CurationCv.COMMENT, value=f('Curation/DataSource')),
    ),
    membership=MembershipBuilder(
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.SMALL_MOLECULE,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.BINDINGDB, value=f('BindingDB MonomerID')),
                    CV(term=IdentifierNamespaceCv.CHEMBL_COMPOUND, value=f('BindingDB Ligand Name', delimiter='::', extract='chembl')),
                    CV(term=IdentifierNamespaceCv.ZINC, value=f('BindingDB Ligand Name', delimiter='::', extract='zinc')),
                    CV(term=IdentifierNamespaceCv.CAS, value=f('BindingDB Ligand Name', delimiter='::', extract='cas')),
                    CV(term=IdentifierNamespaceCv.CHEBI, value=f('BindingDB Ligand Name', delimiter='::', extract='chebi')),
                    CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=f('BindingDB Ligand Name', delimiter='::', extract='pubchem_cid')),
                    CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('BindingDB Ligand Name', delimiter='::', extract='kegg')),
                    CV(term=IdentifierNamespaceCv.STANDARD_INCHI_KEY, value=f('Ligand InChI Key')),
                    CV(term=IdentifierNamespaceCv.STANDARD_INCHI, value=f('Ligand InChI')),
                    CV(term=IdentifierNamespaceCv.SMILES, value=f('Ligand SMILES')),
                    CV(term=IdentifierNamespaceCv.PUBCHEM_COMPOUND, value=f('PubChem CID')),
                    CV(term=IdentifierNamespaceCv.PUBCHEM, value=f('PubChem SID')),
                    CV(term=IdentifierNamespaceCv.CHEBI, value=f('ChEBI ID of Ligand')),
                    CV(term=IdentifierNamespaceCv.CHEMBL_COMPOUND, value=f('ChEMBL ID of Ligand')),
                    CV(term=IdentifierNamespaceCv.DRUGBANK, value=f('DrugBank ID of Ligand')),
                    CV(term=IdentifierNamespaceCv.KEGG_COMPOUND, value=f('KEGG ID of Ligand')),
                    CV(term=IdentifierNamespaceCv.ZINC, value=f('ZINC ID of Ligand')),
                    CV(term=IdentifierNamespaceCv.PDB, value=f('Ligand HET ID in PDB')),
                ),
            ),
        ),
        Member(
            entity=EntityBuilder(
                entity_type=EntityTypeCv.PROTEIN,
                identifiers=IdentifiersBuilder(
                    CV(term=IdentifierNamespaceCv.NAME, value=f('Target Name')),
                    CV(term=IdentifierNamespaceCv.UNIPROT, value=f('UniProt (SwissProt) Primary ID of Target Chain')),
                    CV(term=IdentifierNamespaceCv.NAME, value=f('UniProt (SwissProt) Recommended Name of Target Chain')),
                    CV(term=IdentifierNamespaceCv.UNIPROT_TREMBL, value=f('UniProt (TrEMBL) Primary ID of Target Chain')),
                    CV(term=IdentifierNamespaceCv.NAME, value=f('UniProt (TrEMBL) Submitted Name of Target Chain')),
                ),
                annotations=AnnotationsBuilder(
                    CV(term=IdentifierNamespaceCv.NCBI_TAX_ID, value=tax_value),
                ),
            ),
        ),
    ),
)

resource = Resource(
    config,
    interactions=Dataset(
        download=Download(
            url=_bindingdb_url,
            filename=_bindingdb_filename,
            subfolder='bindingdb',
            large=True,
            encoding='utf-8',
            ext='zip',
        ),
        mapper=interactions_schema,
        raw_parser=_bindingdb_raw,
    ),
)

interactions = resource.interactions
