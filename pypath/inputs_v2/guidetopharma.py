"""
Parse Guide to Pharmacology data and emit Entity records.

This module converts Guide to Pharmacology ligand-target interaction data
into Entity records using the schema defined in pypath.internals.silver_schema.
"""

from __future__ import annotations

from collections.abc import Mapping
import csv
import re

from biolink_model.datamodel.model import (
    CausalMechanismQualifierEnum,
    ChemicalEntity,
    DirectionQualifierEnum,
    Gene,
    MacromolecularComplex,
    NamedThing,
    NucleicAcidEntity,
    Polypeptide,
    Protein,
    RNAProduct,
    slots,
)
from omnipath_core.naming import Namespace

from pypath.inputs_v2._measurements import measurement
from pypath.inputs_v2.base import Dataset, Download, Resource, ResourceConfig
from pypath.internals.cv_terms import LicenseCV, ResourceCv, UpdateCategoryCV
from pypath.internals.tabular_builder import (
    CV,
    AnnotationsBuilder,
    EntityBuilder,
    FieldConfig,
    IdentifiersBuilder,
    RelationBuilder,
)

ligand_chemical_type_mapping = {
    'Synthetic organic': 'synthetic organic',
    'Natural product': 'natural product',
    'Metabolite': 'metabolite',
    'Inorganic': 'inorganic',
    'Peptide': 'peptide',
    'Antibody': 'antibody',
    'Nucleic acid': 'nucleic acid',
}
species_to_taxid = {
    'Human': '9606',
    'Mouse': '10090',
    'Rat': '10116',
    'Rabbit': '9986',
    'Guinea pig': '10141',
    'Hamster': '10036',
    'Dog': '9615',
    'Pig': '9823',
    'Bovine': '9913',
    'Sheep': '9940',
    'Chicken': '9031',
    'Turkey': '9103',
    'Monkey': '9544',
    'Gorilla': '9593',
    'Ferret': '9669',
    'Honeybee': '7460',
    'spiny starfish': '7609',
    'Escherichia coli': '562',
    'Mycobacterium tuberculosis': '1773',
    'Plasmodium falciparum': '5833',
    'Plasmodium falciparum 3D7': '36329',
    'Plasmodium falciparum 7G8': '57266',
    'Plasmodium falciparum A2': '5833',
    'Plasmodium falciparum D6': '478860',
    'Plasmodium falciparum Dd2': '57267',
    'Plasmodium falciparum FC27': '5833',
    'Plasmodium falciparum FCR3': '5833',
    'Plasmodium falciparum HB3': '137071',
    'Plasmodium falciparum K1': '5839',
    'Plasmodium falciparum NF54': '5833',
    'Plasmodium falciparum TM4': '5833',
    'Plasmodium falciparum TM90C2A': '5833',
    'Plasmodium falciparum TM90C2B': '5833',
    'Plasmodium falciparum TM91C235': '5833',
    'Plasmodium falciparum V1/S': '5833',
    'Plasmodium falciparum W2': '5833',
    'Plasmodium vivax': '5855',
    'Plasmodium berghei': '5821',
    'Plasmodium cynomolgi': '5827',
    'Plasmodium knowlesi': '5850',
    'Plasmodium yoelii': '5861',
    'SARS-CoV': '694009',
    'SARS-CoV-2': '2697049',
    'MERS-CoV': '1335626',
    'HCoV-OC43': '31631',
    'Hepatitis C virus': '3052230',
    'Zika virus': '64320',
}
f = FieldConfig(
    extract={'lower': str.lower},
    map={
        'species_taxid': {
            name: f'NCBITaxon:{taxon}'
            for name, taxon in species_to_taxid.items()
        }
    },
)
_NUCLEIC_ACID_LIGAND_SUBTYPES = {'nucleic acid'}
_MACROMOLECULE_LIGAND_SUBTYPES = {'peptide', 'antibody'}
_TAXON_SCOPED_ENTITY_TYPES = {
    Protein,
    Gene,
    RNAProduct,
    NucleicAcidEntity,
    Polypeptide,
    NucleicAcidEntity,
}


def _split_pipe_values(raw: object) -> list[str]:
    if raw is None:
        return []
    return [
        part.strip() for part in str(raw).split('|') if part and part.strip()
    ]


def _filter_refseq(raw: object, prefixes: tuple[str, ...]) -> list[str]:
    values = []
    for value in _split_pipe_values(raw):
        if value.upper().startswith(prefixes):
            values.append(value)
    return values


_UNIPROT_ACCESSION_RE = re.compile(
    '^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\\d+)?$'
)
_UNIPROT_ENTRY_NAME_RE = re.compile('^[A-Z0-9]+_[A-Z0-9]+$')


def _normalize_chembl(raw: object) -> list[str]:
    values = []
    for value in _split_pipe_values(raw):
        upper = value.upper()
        if upper.startswith('CHEMBL'):
            values.append(f'CHEMBL{value[len("CHEMBL") :]}')
    return values


def _filter_uniprot_accessions(raw: object) -> list[str]:
    return [
        value
        for value in _split_pipe_values(raw)
        if _UNIPROT_ACCESSION_RE.fullmatch(value)
    ]


def _filter_uniprot_entry_names(raw: object) -> list[str]:
    return [
        value
        for value in _split_pipe_values(raw)
        if _UNIPROT_ENTRY_NAME_RE.fullmatch(value)
    ]


def _iter_guidetopharma_csv(opener, **_kwargs: object):
    if not opener or not opener.result:
        return
    next(opener.result, None)
    reader = csv.DictReader(opener.result)
    yield from reader


def _guidetopharma_rows_by_id(
    download: Download, key: str, *, force_refresh: bool = False
) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    opener = download.open(force_refresh=force_refresh)
    for row in _iter_guidetopharma_csv(opener):
        value = (row.get(key) or '').strip()
        if value:
            rows[value] = row
    return rows


def _iter_guidetopharma_interactions(opener, **kwargs: object):
    ligands = _guidetopharma_rows_by_id(
        _ligands_download,
        'Ligand ID',
        force_refresh=bool(kwargs.get('force_refresh', False)),
    )
    targets = _guidetopharma_rows_by_id(
        _targets_download,
        'Target id',
        force_refresh=bool(kwargs.get('force_refresh', False)),
    )
    for row in _iter_guidetopharma_csv(opener):
        row = dict(row)
        row['_ligand'] = ligands.get((row.get('Ligand ID') or '').strip(), {})
        row['_target'] = targets.get((row.get('Target ID') or '').strip(), {})
        row['_target_ligand'] = ligands.get(
            (row.get('Target Ligand ID') or '').strip(), {}
        )
        yield row


def _nested_row(row: Mapping[str, object], key: str) -> Mapping[str, object]:
    value = row.get(key)
    return value if isinstance(value, Mapping) else {}


def _ligand_row(row: Mapping[str, object]) -> Mapping[str, object]:
    return _nested_row(row, '_ligand') or row


def _interaction_target_row(row: Mapping[str, object]) -> Mapping[str, object]:
    return _nested_row(row, '_target')


def _target_ligand_row(row: Mapping[str, object]) -> Mapping[str, object]:
    return _nested_row(row, '_target_ligand')


def _guidetopharma_ligand_entity_type(row):
    row = _ligand_row(row)
    ligand_subtype = ligand_chemical_type_mapping.get(
        (row.get('Type') or '').strip()
    )
    if _filter_uniprot_accessions(row.get('UniProt ID')):
        return Protein
    if ligand_subtype in _NUCLEIC_ACID_LIGAND_SUBTYPES:
        return NucleicAcidEntity
    if ligand_subtype in _MACROMOLECULE_LIGAND_SUBTYPES:
        return Polypeptide
    return ChemicalEntity


def _guidetopharma_ligand_taxon_id(row):
    row = _ligand_row(row)
    entity_type = _guidetopharma_ligand_entity_type(row)
    taxon_id = species_to_taxid.get((row.get('Species') or '').strip())
    return (
        f'NCBITaxon:{taxon_id}'
        if taxon_id and entity_type in _TAXON_SCOPED_ENTITY_TYPES
        else None
    )


def _guidetopharma_target_entity_type(row):
    if row.get('_entity_type') is not None:
        return row['_entity_type']
    if _split_pipe_values(row.get('Subunit id')):
        return MacromolecularComplex
    if any(
        len(_split_pipe_values(row.get(field))) > 1
        for field in (
            'Human SwissProt',
            'Human Ensembl Gene',
            'Human Entrez Gene',
            'HGNC id',
        )
    ):
        return NamedThing
    if any(
        (
            _split_pipe_values(row.get(field))
            for field in (
                'Human SwissProt',
                'Human Ensembl Gene',
                'Human Entrez Gene',
                'HGNC id',
                'HGNC symbol',
            )
        )
    ):
        return Protein
    return NamedThing


config = ResourceConfig(
    id=ResourceCv.GUIDETOPHARMA,
    name='Guide to Pharmacology',
    url='https://www.guidetopharmacology.org/',
    license=LicenseCV.CC_BY_SA_4_0,
    update_category=UpdateCategoryCV.REGULAR,
    pubmed='41160876',
    primary_category='interactions',
    description='The IUPHAR/BPS Guide to PHARMACOLOGY is an expert-curated resource of ligand-activity-target relationships, providing quantitative information on drug targets and the prescription medicines and experimental drugs that act on them. It covers G protein-coupled receptors, voltage-gated ion channels, ligand-gated ion channels, nuclear hormone receptors, catalytic receptors, enzymes, and transporters.',
)
_targets_download = Download(
    url='http://www.guidetopharmacology.org/DATA/targets_and_families.csv',
    filename='targets_and_families.csv',
    subfolder='guidetopharma',
)
_ligands_download = Download(
    url='http://www.guidetopharmacology.org/DATA/ligands.csv',
    filename='ligands.csv',
    subfolder='guidetopharma',
)
_interactions_download = Download(
    url='http://www.guidetopharmacology.org/DATA/interactions.csv',
    filename='interactions.csv',
    subfolder='guidetopharma',
)


class _TargetIdentifiersBuilder(IdentifiersBuilder):
    """Component accessions are not aliases of a source target group."""

    def build(self, row, cache=None):
        identifiers = super().build(row, cache)
        if _guidetopharma_target_entity_type(row) != Protein:
            # Component membership needs a join on source Subunit id. Until that
            # is available, retain the group without asserting false identity.
            identifiers = [
                identifier
                for identifier in identifiers
                if identifier.type
                in {
                    Namespace.GUIDETOPHARMA_TARGET,
                    Namespace.NAME,
                    Namespace.SYNONYM,
                }
            ]
        return identifiers


targets_schema = EntityBuilder(
    entity_type=_guidetopharma_target_entity_type,
    identifiers=_TargetIdentifiersBuilder(
        CV(term=Namespace.GUIDETOPHARMA_TARGET, value=f('Target id')),
        CV(term=Namespace.UNIPROT, value=f('Human SwissProt', delimiter='|')),
        CV(term=Namespace.ENSG, value=f('Human Ensembl Gene', delimiter='|')),
        CV(term=Namespace.ENTREZ, value=f('Human Entrez Gene', delimiter='|')),
        CV(
            term=Namespace.REFSEQ,
            value=lambda row: _filter_refseq(
                row.get('Human nucleotide RefSeq'),
                (
                    'NM_',
                    'XM_',
                    'NR_',
                    'XR_',
                    'NG_',
                    'NT_',
                    'NW_',
                    'NC_',
                    'AC_',
                    'NP_',
                ),
            ),
        ),
        CV(
            term=Namespace.REFSEQ_PROTEIN,
            value=lambda row: _filter_refseq(
                row.get('Human protein RefSeq'),
                ('NP_', 'XP_', 'YP_', 'WP_', 'AP_', 'ZP_'),
            ),
        ),
        CV(
            term=Namespace.REFSEQ,
            value=lambda row: _filter_refseq(
                row.get('Human protein RefSeq'), ('NM_', 'XM_', 'NR_', 'XR_')
            ),
        ),
        CV(term=Namespace.HGNC, value=f('HGNC id', delimiter='|')),
        CV(term=Namespace.GENESYMBOL, value=f('HGNC symbol', delimiter='|')),
        CV(term=Namespace.NAME, value=f('Target name')),
        CV(term=Namespace.SYNONYM, value=f('synonyms', delimiter='|')),
        CV(term=Namespace.SYNONYM, value=f('Target systematic name')),
        CV(term=Namespace.SYNONYM, value=f('Target abbreviated name')),
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.in_taxon,
            value=lambda row: f'NCBITaxon:{row["_taxon_id"]}'
            if row.get('_taxon_id')
            else 'NCBITaxon:9606'
            if any(
                (
                    row.get(k)
                    for k in (
                        'Human SwissProt',
                        'Human Ensembl Gene',
                        'Human Entrez Gene',
                        'HGNC id',
                    )
                )
            )
            else None,
        )
    ),
)
ligands_schema = EntityBuilder(
    entity_type=_guidetopharma_ligand_entity_type,
    identifiers=IdentifiersBuilder(
        CV(term=Namespace.GUIDETOPHARMA, value=f('Ligand ID')),
        CV(term=Namespace.PUBCHEM_SUBSTANCE, value=f('PubChem SID')),
        CV(term=Namespace.PUBCHEM, value=f('PubChem CID')),
        CV(
            term=Namespace.UNIPROT,
            value=lambda row: _filter_uniprot_accessions(row.get('UniProt ID')),
        ),
        CV(
            term=Namespace.UNIPROT_ENTRY,
            value=lambda row: _filter_uniprot_entry_names(
                row.get('UniProt ID')
            ),
        ),
        CV(term=Namespace.ENSG, value=f('Ensembl ID', delimiter='|')),
        CV(
            term=Namespace.CHEMBL,
            value=lambda row: _normalize_chembl(row.get('ChEMBL ID')),
        ),
        CV(term=Namespace.INCHIKEY, value=f('InChIKey')),
        CV(term=Namespace.SMILES, value=f('SMILES')),
        CV(term=Namespace.INCHI, value=f('InChI')),
        CV(term=Namespace.SYNONYM, value=f('Synonyms')),
        CV(term=Namespace.NAME, value=f('Name')),
        CV(term=Namespace.SYNONYM, value=f('IUPAC name')),
        CV(term=Namespace.SYNONYM, value=f('INN')),
    ),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=_guidetopharma_ligand_taxon_id)
    ),
)


class _GuidetopharmaLigandMemberBuilder:
    def __call__(self, row: Mapping[str, object]) -> object | None:
        return self.build(row)

    def build(self, row: Mapping[str, object]) -> object | None:
        member_row = _ligand_row(row)
        return ligands_schema.build(member_row) if member_row else None


class _GuidetopharmaTargetMemberBuilder:
    def __call__(self, row: Mapping[str, object]) -> object | None:
        return self.build(row)

    def build(self, row: Mapping[str, object]) -> object | None:
        if row.get('Target ID'):
            target = _interaction_target_row(row)
            if target:
                target = dict(target)
                taxon = species_to_taxid.get(
                    str(row.get('Target Species') or '').strip()
                )
                if taxon:
                    target['_taxon_id'] = taxon
                if taxon and taxon != '9606':
                    # Human reference aliases cannot identify a nonhuman assay
                    # target. Preserve the taxon-scoped source target accession.
                    target['_entity_type'] = _guidetopharma_target_entity_type(
                        target
                    )
                    for field in list(target):
                        if field.startswith('Human ') or field in (
                            'HGNC id',
                            'HGNC symbol',
                        ):
                            target.pop(field)
                return targets_schema.build(target)
            # Retain an explicit source target accession even when a cached
            # metadata table does not contain its newly added record.
            return targets_schema.build(
                {
                    'Target id': row['Target ID'],
                    'Target name': row.get('Target'),
                }
            )
        if row.get('Target Ligand ID'):
            return ligands_schema.build(_target_ligand_row(row))
        return None


# Exact source actions map to schema-defined mechanisms. Broader source
# qualifiers such as "Competitive" do not alone determine inhibition.
_ACTION_MECHANISM = {
    'Activation': CausalMechanismQualifierEnum.activation,
    'Agonist': CausalMechanismQualifierEnum.agonism,
    'Full agonist': CausalMechanismQualifierEnum.agonism,
    'Partial agonist': CausalMechanismQualifierEnum.partial_agonism,
    'Biased agonist': CausalMechanismQualifierEnum.biased_agonism,
    'Antagonist': CausalMechanismQualifierEnum.antagonism,
    'Inhibition': CausalMechanismQualifierEnum.inhibition,
    'Feedback inhibition': CausalMechanismQualifierEnum.feedback_inhibition,
    'Irreversible inhibition': CausalMechanismQualifierEnum.irreversible_inhibition,
    'Inverse agonist': CausalMechanismQualifierEnum.inverse_agonism,
    'Positive': CausalMechanismQualifierEnum.positive_modulation,
    'Negative': CausalMechanismQualifierEnum.negative_modulation,
    'Potentiation': CausalMechanismQualifierEnum.potentiation,
    'Pore blocker': CausalMechanismQualifierEnum.molecular_channel_blockage,
    'Binding': CausalMechanismQualifierEnum.binding,
    'Mixed': CausalMechanismQualifierEnum.modulation,
    'Biphasic': CausalMechanismQualifierEnum.modulation,
}

_AFFINITY_SLOT = {
    'pIC50': slots.pIC50,
    'pEC50': slots.pEC50,
    'pKi': slots.pKi,
    'pKd': slots.pKd,
}


def _affinity_measurements(row):
    measure = str(row.get('Affinity Units') or '').strip()
    return [
        value
        for field in ('Affinity High', 'Affinity Low', 'Affinity Median')
        if (
            value := measurement(
                row.get(field),
                source_field=f'{measure or "unspecified measure"} {field}',
            )
        )
        is not None
    ]


def _original_affinity_measurements(row):
    measure = str(row.get('Original Affinity Units') or '').strip()
    return [
        value
        for field in (
            'Original Affinity High nm',
            'Original Affinity Low nm',
            'Original Affinity Median nm',
        )
        if (
            value := measurement(
                row.get(field),
                # These source columns explicitly express concentrations in nM.
                unit='nM',
                source_field=f'{measure or "unspecified measure"} {field}',
                comparator=row.get('Original Affinity Relation'),
            )
        )
        is not None
    ]


_ACTION_DIRECTION = {
    'Activation': DirectionQualifierEnum.increased,
    'Agonist': DirectionQualifierEnum.increased,
    'Full agonist': DirectionQualifierEnum.increased,
    'Partial agonist': DirectionQualifierEnum.increased,
    'Inhibition': DirectionQualifierEnum.decreased,
    'Antagonist': DirectionQualifierEnum.decreased,
    'Feedback inhibition': DirectionQualifierEnum.decreased,
    'Irreversible inhibition': DirectionQualifierEnum.decreased,
    'Inverse agonist': DirectionQualifierEnum.decreased,
    'Voltage-dependent inhibition': DirectionQualifierEnum.decreased,
}


def guidetopharma_predicate(row: Mapping[str, object]):
    action = str(row.get('Action') or '').strip()
    if action in _ACTION_DIRECTION or (
        action in _ACTION_MECHANISM and action != 'Binding'
    ):
        return slots.affects
    if action == 'Binding':
        return slots.directly_physically_interacts_with
    return slots.interacts_with


interactions_schema = RelationBuilder(
    subject=_GuidetopharmaLigandMemberBuilder(),
    predicate=guidetopharma_predicate,
    object=_GuidetopharmaTargetMemberBuilder(),
    identifiers=IdentifiersBuilder(
        CV(
            term=Namespace.NAME,
            value=f(
                lambda row: f'{row["Target ID"]}_{row["Ligand ID"]}'
                if row.get('Target ID') and row.get('Ligand ID')
                else None
            ),
        )
    ),
    annotations=AnnotationsBuilder(
        CV(
            term=slots.has_quantitative_value,
            value=_original_affinity_measurements,
        ),
        CV(term=slots.description, value=f('Action comment')),
        CV(term=slots.description, value=f('Assay Description')),
        CV(
            term=slots.causal_mechanism_qualifier,
            value=lambda row: _ACTION_MECHANISM.get(
                str(row.get('Action') or '').strip()
            ),
        ),
        CV(
            term=lambda row: _AFFINITY_SLOT.get(
                str(row.get('Affinity Units') or '').strip(),
                slots.has_quantitative_value,
            ),
            value=_affinity_measurements,
        ),
        CV(
            term=slots.object_direction_qualifier,
            value=lambda row: _ACTION_DIRECTION.get(
                str(row.get('Action') or '').strip()
            ),
        ),
        CV(
            term=slots.publications,
            value=f(
                'PubMed ID',
                delimiter='|',
                transform=lambda v: f'PMID:{v.removeprefix("PMID:")}',
            ),
        ),
        CV(term=slots.in_taxon, value=f('Target Species', map='species_taxid')),
    ),
)
resource = Resource(
    config,
    interactions=Dataset(
        download=_interactions_download,
        mapper=interactions_schema,
        raw_parser=_iter_guidetopharma_interactions,
    ),
)
