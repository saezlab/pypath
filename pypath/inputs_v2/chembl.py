"""
Parse ChEMBL data and emit Entity records.

This module converts ChEMBL ligand-target interaction data
into Entity records using the schema defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from functools import partial
from pathlib import Path
import re

from biolink_model.datamodel.model import (
    AnatomicalEntity,
    Cell,
    CellLine,
    ChemicalEntity,
    DirectionQualifierEnum,
    Gene,
    ProteinFamily,
    MacromolecularComplex,
    MolecularActivity,
    MolecularEntity,
    NamedThing,
    NucleicAcidEntity,
    OrganismTaxon,
    Protein,
    RNAProduct,
    slots,
)
from omnipath_core.naming import Namespace

from pypath.inputs_v2._measurements import measurement as _measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.inputs_v2.parsers.chembl import (
    activities_parser,
    mechanisms_parser,
    molecules_parser,
)
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    MembersFromList,
    MembershipBuilder,
    RelationBuilder,
)
from pypath.share import cache


VERSION = 36
DB_REL_PATH = f'chembl_{VERSION}/chembl_{VERSION}_sqlite/chembl_{VERSION}.db'
SQLITE_PATH = Path(cache.get_cachedir()) / f'ChEMBL_SQLite_{VERSION}.sqlite'


def _chembl_url(version: int = VERSION, **_kwargs: object) -> str:
    return f'https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_{version:02d}/chembl_{version:02d}_sqlite.tar.gz'


def _files_needed(version: int = VERSION, **_kwargs: object) -> list[str]:
    return [
        f'chembl_{version:02d}/chembl_{version:02d}_sqlite/chembl_{version:02d}.db'
    ]


config = ResourceConfig(
    id=ResourceCv.CHEMBL,
    name='ChEMBL',
    url='https://www.ebi.ac.uk/chembl/',
    license=LicenseCV.CC_BY_SA_3_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='21948594',
    primary_category='interactions',
    description='ChEMBL is a manually curated chemical database of bioactive molecules with drug-like properties.',
)
download = Download(
    url=_chembl_url,
    filename=f'chembl_{VERSION}_sqlite.tar.gz',
    subfolder='chembl',
    large=True,
    ext='.tar.gz',
    needed=_files_needed(),
)
MOLECULE_TYPE_TO_ENTITY_TYPE = {
    'Small molecule': ChemicalEntity,
    'Protein': Protein,
    'Antibody': Protein,
    'Enzyme': Protein,
    'Oligosaccharide': ChemicalEntity,
    'Oligonucleotide': ChemicalEntity,
    'Gene': Gene,
    'Cell': Cell,
    'Unknown': ChemicalEntity,
    'Unclassified': ChemicalEntity,
}
TARGET_TYPE_MAP = {
    'SINGLE PROTEIN': Protein,
    'PROTEIN COMPLEX': MacromolecularComplex,
    'PROTEIN FAMILY': ProteinFamily,
    'PROTEIN-PROTEIN INTERACTION': MolecularActivity,
    'SELECTIVITY GROUP': NamedThing,
    'NUCLEIC-ACID': NucleicAcidEntity,
    'ORGANISM': OrganismTaxon,
    'CELL LINE': CellLine,
    'SUBCELLULAR': NamedThing,
    'MACROMOLECULE': MolecularEntity,
    'TISSUE': AnatomicalEntity,
}
COMPONENT_TYPE_MAP = {
    'PROTEIN': Protein,
    'RNA': RNAProduct,
    'DNA': NucleicAcidEntity,
}
AFFINITY_TERMS = {
    'IC50': 'BAO:0000190',
    'EC50': 'BAO:0000188',
    'KI': 'BAO:0000192',
    'KD': 'BAO:0000034',
    'KON': 'BAO:0000480',
    'KOFF': 'BAO:0000479',
}


def _ensembl_namespace(value):
    match = re.fullmatch('ENS[A-Z]*([GTP])\\d+(?:\\.\\d+)?', str(value))
    return (
        {'G': Namespace.ENSG, 'T': Namespace.ENST, 'P': Namespace.ENSP}.get(
            match.group(1)
        )
        if match
        else None
    )


f = FieldConfig(
    map={
        'entity_type': MOLECULE_TYPE_TO_ENTITY_TYPE,
        'target_type': TARGET_TYPE_MAP,
        'component_type': COMPONENT_TYPE_MAP,
    }
)


def _split_chembl_list(value: object) -> list[str]:
    if value is None:
        return []
    return [
        item.strip() for item in str(value).split(',') if item and item.strip()
    ]


def _target_component_values(row: dict[str, object], key: str) -> list[str]:
    if row.get('target_type') != 'SINGLE PROTEIN':
        return []
    return _split_chembl_list(row.get(key))


molecules_schema = EntityBuilder(
    entity_type=f('molecule_type', map='entity_type'),
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.CHEMBL, value=f('chembl_id')),
        CV(term=Namespace.SMILES, value=f('canonical_smiles')),
        CV(term=Namespace.INCHI, value=f('standard_inchi')),
        CV(term=Namespace.INCHIKEY, value=f('standard_inchi_key')),
        CV(term=Namespace.NAME, value=f('pref_name')),
    ),
    annotations=AnnotationsBuilder(
        CV(term='chemrof:mass', value=f('full_mwt'))
    ),
)
targets_schema = EntityBuilder(
    entity_type=f('target_type', map='target_type'),
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.CHEMBL_TARGET, value=f('chembl_id')),
        CV(term=Namespace.NAME, value=f('pref_name')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.in_taxon,
            value=f(
                'tax_id',
                transform=lambda v: 'NCBITaxon:'
                + str(v).removeprefix('NCBITaxon:'),
            ),
        )
    ),
    membership=MembershipBuilder(
        MembersFromList(
            entity_type=f(
                'component_types', delimiter=',', map='component_type'
            ),
            identifiers=IdentifiersBuilder(
                CV(
                    term=Namespace.UNIPROT,
                    value=f(
                        'component_uniprot_accessions',
                        delimiter=',',
                        preserve_indices=True,
                    ),
                ),
                CV(
                    term=f(
                        'component_ensembl_accessions',
                        delimiter=',',
                        preserve_indices=True,
                        transform=_ensembl_namespace,
                    ),
                    value=f(
                        'component_ensembl_accessions',
                        delimiter=',',
                        preserve_indices=True,
                    ),
                ),
            ),
            entity_annotations=AnnotationsBuilder(
                CV(
                    term=slots.in_taxon,
                    value=f(
                        'tax_id',
                        transform=lambda v: 'NCBITaxon:'
                        + str(v).removeprefix('NCBITaxon:'),
                    ),
                ),
                CV(
                    term=slots.description,
                    value=f('component_descriptions', delimiter=','),
                ),
            ),
        )
    ),
)
ACTION_DIRECTION = {
    'AGONIST': DirectionQualifierEnum.increased,
    'PARTIAL AGONIST': DirectionQualifierEnum.increased,
    'ACTIVATOR': DirectionQualifierEnum.increased,
    'POTENTIATOR': DirectionQualifierEnum.increased,
    'POSITIVE ALLOSTERIC MODULATOR': DirectionQualifierEnum.increased,
    'ANTAGONIST': DirectionQualifierEnum.decreased,
    'INVERSE AGONIST': DirectionQualifierEnum.decreased,
    'INHIBITOR': DirectionQualifierEnum.decreased,
    'NEGATIVE ALLOSTERIC MODULATOR': DirectionQualifierEnum.decreased,
}


def chembl_predicate(row):
    action = str(row.get('action_type') or '').strip().upper()
    if action in ACTION_DIRECTION:
        return slots.affects
    if action == 'BINDING AGENT':
        return slots.physically_interacts_with
    return slots.associated_with


molecule_builder = EntityBuilder(
    entity_type=ChemicalEntity,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.CHEMBL, value=f('molecule_chembl_id'))
    ),
)
target_builder = EntityBuilder(
    entity_type=f('target_type', map='target_type', default=NamedThing),
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.CHEMBL_TARGET, value=f('target_chembl_id')),
        CV(
            term=Namespace.UNIPROT,
            value=lambda row: _target_component_values(
                row, 'target_component_uniprot_accessions'
            ),
        ),
        CV(
            term=lambda row: [
                _ensembl_namespace(v)
                for v in _target_component_values(
                    row, 'target_component_ensembl_accessions'
                )
            ],
            value=lambda row: _target_component_values(
                row, 'target_component_ensembl_accessions'
            ),
        ),
        CV(term=Namespace.NAME, value=f('target_pref_name')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.in_taxon,
            value=f(
                'target_tax_id',
                transform=lambda v: 'NCBITaxon:'
                + str(v).removeprefix('NCBITaxon:'),
            ),
        )
    ),
)
activities_schema = RelationBuilder(
    subject=molecule_builder,
    predicate=chembl_predicate,
    object=target_builder,
    annotations=AnnotationsBuilder(
        CV(term=slots.source_record_urls, value=f('assay_chembl_id', transform=lambda v: f'https://www.ebi.ac.uk/chembl/explore/assay/{v}')),
        CV(term=slots.source_record_urls, value=f('document_chembl_id', transform=lambda v: f'https://www.ebi.ac.uk/chembl/explore/document/{v}')),
        CV(
            term=slots.has_quantitative_value,
            value=f(
                'pchembl_value',
                transform=lambda v: _measurement(
                    v, source_field='pchembl_value'
                ),
            ),
        ),
        CV(term=slots.description, value=f('data_validity_comment')),
        CV(term=slots.description, value=f('action_description')),
        CV(
            term=slots.in_taxon,
            value=f(
                'assay_tax_id',
                transform=lambda v: 'NCBITaxon:'
                + str(v).removeprefix('NCBITaxon:'),
            ),
        ),
        CV(
            term=slots.publications,
            value=f(
                'pubmed_id',
                transform=lambda v: 'PMID:' + str(v).removeprefix('PMID:'),
            ),
        ),
        CV(
            term=slots.publications,
            value=f(
                'doi', transform=lambda v: 'doi:' + str(v).removeprefix('doi:')
            ),
        ),
        CV(term=slots.description, value=f('mechanism_of_action')),
        CV(term=slots.description, value=f('mechanism_comment')),
        CV(term=slots.description, value=f('selectivity_comment')),
        CV(term=slots.description, value=f('binding_site_comment')),
        CV(
            term=slots.object_direction_qualifier,
            value=lambda row: ACTION_DIRECTION.get(
                str(row.get('action_type') or '').strip().upper()
            ),
        ),
        CV(
            term=lambda row: AFFINITY_TERMS.get(
                str(row.get('standard_type') or '').upper()
            ),
            value=lambda row: _measurement(
                row.get('standard_value'),
                row.get('standard_units'),
                row.get('standard_type'),
                row.get('standard_relation'),
            ),
        ),
        CV(term=slots.chembl_confidence_score, value=f('confidence_score')),
    ),
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.CHEMBL_ACTIVITY, value=f('activity_id')),
        CV(term=Namespace.CHEMBL_MECHANISM, value=f('mec_id')),
    ),
)
resource = Resource(
    config=config,
    molecules=Dataset(
        download=download,
        mapper=molecules_schema,
        raw_parser=partial(
            molecules_parser, sqlite_path=SQLITE_PATH, db_rel_path=DB_REL_PATH
        ),
    ),
    activities=Dataset(
        download=download,
        mapper=activities_schema,
        raw_parser=partial(
            activities_parser, sqlite_path=SQLITE_PATH, db_rel_path=DB_REL_PATH
        ),
    ),
    mechanisms=Dataset(
        download=download,
        mapper=activities_schema,
        raw_parser=partial(
            mechanisms_parser, sqlite_path=SQLITE_PATH, db_rel_path=DB_REL_PATH
        ),
    ),
)
