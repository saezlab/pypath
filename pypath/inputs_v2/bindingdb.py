"""
Parse BindingDB data and emit Entity records.

This module converts BindingDB binding affinity data into Entity records using
the schema defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from biolink_model.datamodel.model import (
    ChemicalEntity,
    MacromolecularComplex,
    Protein,
    slots,
)
from omnipath_core.naming import Namespace

from pypath.inputs_v2._measurements import measurement as _measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.bindingdb import _raw
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.silver_schema import (
    Entity,
    Identifier,
    Membership,
)
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)


def _bindingdb_url(bindingdb_subset: str = 'All', **_kwargs: object) -> str:
    return f'https://bindingdb.org/rwd/bind/downloads/BindingDB_{bindingdb_subset}_202605_tsv.zip'


def _bindingdb_filename(bindingdb_subset: str = 'All', **_kwargs: object) -> str:
    return f'BindingDB_{bindingdb_subset}_202605_tsv.zip'


config = ResourceConfig(
    id=ResourceCv.BINDINGDB,
    name='BindingDB',
    url='https://www.bindingdb.org/',
    license=LicenseCV.CC_BY_3_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='26481362',
    primary_category='interactions',
    description='BindingDB is a public, web-accessible database of measured binding affinities, focusing chiefly on the interactions of proteins considered to be drug-targets with small, drug-like molecules. It contains binding data for over 2 million protein-ligand complexes with experimental measurements including Ki, Kd, IC50, and EC50 values.',
)
f = FieldConfig(
    extract={
        'chembl': '^(CHEMBL\\d+)$',
        'zinc': '^(ZINC\\d+)$',
        'cas': '^(\\d{2,7}-\\d{2}-\\d)$',
        'chebi': '^(?:CHEBI[:\\s]?)?(\\d+)$',
        'pubchem_cid': '^CID[:\\s]?(\\d+)$',
        'kegg': '^(C\\d{5})$',
        'tax': '(\\d+)',
        'uniprot': '^([A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}|[OPQ][0-9][A-Z0-9]{3}[0-9])$',
    }
)
tax_value = lambda row: {
    'Homo sapiens': 'NCBITaxon:9606',
    'Mus musculus': 'NCBITaxon:10090',
    'Rattus norvegicus': 'NCBITaxon:10116',
}.get(row.get('Target Source Organism According to Curator or DataSource'))
chemical_builder = EntityBuilder(
    entity_type=ChemicalEntity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.BINDINGDB, value=f('BindingDB MonomerID')),
        CV(
            term=Namespace.CHEMBL,
            value=f('BindingDB Ligand Name', delimiter='::', extract='chembl'),
        ),
        CV(
            term=Namespace.ZINC,
            value=f('BindingDB Ligand Name', delimiter='::', extract='zinc'),
        ),
        CV(
            term=Namespace.CAS,
            value=f('BindingDB Ligand Name', delimiter='::', extract='cas'),
        ),
        CV(
            term=Namespace.PUBCHEM,
            value=f(
                'BindingDB Ligand Name', delimiter='::', extract='pubchem_cid'
            ),
        ),
        CV(
            term=Namespace.KEGG,
            value=f('BindingDB Ligand Name', delimiter='::', extract='kegg'),
        ),
        CV(term=Namespace.INCHIKEY, value=f('Ligand InChI Key')),
        CV(term=Namespace.INCHI, value=f('Ligand InChI')),
        CV(term=Namespace.SMILES, value=f('Ligand SMILES')),
        CV(term=Namespace.PUBCHEM, value=f('PubChem CID')),
        CV(term=Namespace.PUBCHEM_SUBSTANCE, value=f('PubChem SID')),
        CV(
            term=Namespace.CHEBI, value=f('ChEBI ID of Ligand', extract='chebi')
        ),
        CV(
            term=Namespace.CHEMBL,
            value=f('ChEMBL ID of Ligand', delimiter='::', extract='chembl'),
        ),
        CV(term=Namespace.DRUGBANK, value=f('DrugBank ID of Ligand')),
        CV(term=Namespace.KEGG, value=f('KEGG ID of Ligand')),
        CV(term=Namespace.ZINC, value=f('ZINC ID of Ligand')),
        CV(
            term=Namespace.SYNONYM,
            value=f('BindingDB Ligand Name', delimiter='::'),
        ),
    ),
    annotations=AnnotationsBuilder(),
)
target_builder = EntityBuilder(
    entity_type=Protein,
    identifiers=IdentifiersBuilder(
        CV(
            term=Namespace.UNIPROT,
            value=f(
                'UniProt (SwissProt) Primary ID of Target Chain 1',
                delimiter=' ',
                extract='uniprot',
            ),
        ),
        CV(
            term=Namespace.UNIPROT_TREMBL,
            value=f(
                'UniProt (TrEMBL) Primary ID of Target Chain 1',
                delimiter=' ',
                extract='uniprot',
            ),
        ),
        CV(term=Namespace.NAME, value=lambda row: row.get('Target Name') if not (row.get('UniProt (SwissProt) Primary ID of Target Chain 1') or row.get('UniProt (TrEMBL) Primary ID of Target Chain 1')) else None),
        CV(
            term=Namespace.NAME,
            value=f('UniProt (SwissProt) Recommended Name of Target Chain 1'),
        ),
        CV(
            term=Namespace.NAME,
            value=f('UniProt (TrEMBL) Submitted Name of Target Chain 1'),
        ),
    ),
    annotations=AnnotationsBuilder(CV(term=slots.in_taxon, value=tax_value)),
)


def _target(row):
    count = int(
        row.get(
            'Number of Protein Chains in Target (>1 implies a multichain complex)'
        )
        or 1
    )
    if count <= 1:
        return target_builder.build(row)
    members = []
    for index in range(1, count + 1):
        chain = dict(row)
        # The whole target name identifies the complex, not an unnamed chain.
        chain['Target Name'] = None
        for stem in (
            'UniProt (SwissProt) Primary ID of Target Chain ',
            'UniProt (TrEMBL) Primary ID of Target Chain ',
            'UniProt (SwissProt) Recommended Name of Target Chain ',
            'UniProt (TrEMBL) Submitted Name of Target Chain ',
        ):
            chain[stem + '1'] = row.get(stem + str(index))
        if (protein := target_builder.build(chain)) is not None:
            members.append(Membership(member=protein, predicate=slots.has_part))
    return Entity(
        type=MacromolecularComplex,
        identifiers=[
            Identifier(type=Namespace.NAME, value=row.get('Target Name'))
        ],
        membership=members,
        annotations=target_builder.annotations.build(row),
    )


interactions_schema = RelationBuilder(
    subject=chemical_builder,
    predicate=slots.associated_with,
    object=_target,
    annotations=AnnotationsBuilder(
        CV(term=slots.original_object, value=f('Target Name')),
        CV(
            term=slots.has_quantitative_value,
            value=f(
                'pchembl_value',
                transform=lambda v: _measurement(
                    v, source_field='pchembl_value'
                ),
            ),
        ),
        CV(
            term=slots.has_quantitative_value,
            value=f(
                'pH', transform=lambda v: _measurement(v, source_field='pH')
            ),
        ),
        CV(
            term=slots.has_quantitative_value,
            value=f(
                'Temp (C)',
                transform=lambda v: _measurement(v, 'Cel', 'Temp (C)'),
            ),
        ),
        CV(
            term=slots.publications,
            value=f(
                'PMID',
                transform=lambda v: 'PMID:' + str(v).removeprefix('PMID:'),
            ),
        ),
        CV(
            term=slots.publications,
            value=f(
                'Article DOI',
                transform=lambda v: 'doi:' + str(v).removeprefix('doi:'),
            ),
        ),
        CV(
            term='BAO:0000192',
            value=f(
                'Ki (nM)', transform=lambda v: _measurement(v, 'nM', 'Ki (nM)')
            ),
        ),
        CV(
            term='BAO:0000034',
            value=f(
                'Kd (nM)', transform=lambda v: _measurement(v, 'nM', 'Kd (nM)')
            ),
        ),
        CV(
            term='BAO:0000190',
            value=f(
                'IC50 (nM)',
                transform=lambda v: _measurement(v, 'nM', 'IC50 (nM)'),
            ),
        ),
        CV(
            term='BAO:0000188',
            value=f(
                'EC50 (nM)',
                transform=lambda v: _measurement(v, 'nM', 'EC50 (nM)'),
            ),
        ),
        CV(
            term='BAO:0000480',
            value=f(
                'kon (M-1-s-1)',
                transform=lambda v: _measurement(v, '/M/s', 'kon (M-1-s-1)'),
            ),
        ),
        CV(
            term='BAO:0000479',
            value=f(
                'koff (s-1)',
                transform=lambda v: _measurement(v, '/s', 'koff (s-1)'),
            ),
        ),
    ),
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.BINDINGDB, value=f('BindingDB Reactant_set_id')),
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
        raw_parser=_raw,
    ),
)
