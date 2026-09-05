"""Regression coverage for numeric and provenance review findings."""

import importlib

from omnipath_build.silver import SilverExtractor
from omnipath_core.measurements import quantity_dict
from omnipath_core.naming import Namespace
import pytest

from pypath.inputs_v2._measurements import measurement


@pytest.mark.parametrize('comparison', ['<', '>', '=', '<=', '>=', '~', '≈'])
def test_numeric_comparisons_preserved(comparison):
    result = measurement(f'{comparison} -1.25e2', 'nM', 'IC50')
    assert result.quantity.has_numeric_value == -125
    assert result.quantity.has_unit == 'nM'
    assert result.source_field == 'IC50'
    assert result.comparator == comparison
    assert quantity_dict(result)['has_binary_relation'] == {
        '<': 'less_than',
        '>': 'greater_than',
        '=': 'equal_to',
    }.get(comparison)


@pytest.mark.parametrize(
    'raw', [None, '', '-', 'NaN', 'inf', '1e999', 'unknown', '1-2']
)
def test_invalid_numbers_are_omitted(raw):
    assert measurement(raw) is None


def test_conflicting_comparisons_are_not_silently_overwritten():
    assert measurement('<3', comparator='>') is None
    assert measurement('3', comparator=' ≤ ') is None
    assert measurement('≈3', comparator=' ≈ ').comparator == '≈'
    assert measurement(0).quantity.has_numeric_value == 0


@pytest.mark.parametrize(
    'module',
    [
        'bindingdb',
        'chembl',
        'drugcentral',
        'foodb',
        'phenol_explorer',
        'ptfi',
        'stitch',
    ],
)
def test_modules_use_shared_numeric_parser(module):
    assert (
        importlib.import_module(f'pypath.inputs_v2.{module}')._measurement
        is measurement
    )


def _extract(module, row):
    module = importlib.import_module(f'pypath.inputs_v2.{module}')
    mapped = module.interactions_schema(row)
    assert mapped is not None
    result = SilverExtractor(module.config.id.value, 'fixture')
    result.process_record(mapped, row, 'fixture:1', 1)
    return result


def test_bindingdb_retains_occurrence_id_and_approximate_measurement():
    out = _extract(
        'bindingdb',
        {
            'BindingDB MonomerID': '12',
            'BindingDB Reactant_set_id': '34',
            'UniProt (SwissProt) Primary ID of Target Chain 1': 'P00533',
            'IC50 (nM)': '≈ 3',
        },
    )
    assert out.relations[0].upstream_id == '34'
    quantity = next(a['quantity'] for a in out.relations[0].annotations)
    assert quantity['comparator'] == '≈'


def test_drugcentral_preserves_parallel_identifier_positions():
    from pypath.inputs_v2.drugcentral import _target_entity

    target = _target_entity(
        {
            'ACCESSION': 'P00533|P04626|P21860',
            'GENE': 'EGFR||ERBB3',
            'SWISSPROT': '|ERBB2_HUMAN|ERBB3_HUMAN',
            'TARGET_NAME': 'ERBB group',
        }
    )
    identifiers = [
        {i.type: i.value for i in membership.member.identifiers}
        for membership in target.membership
    ]
    assert identifiers[0][Namespace.GENESYMBOL] == 'EGFR'
    assert Namespace.UNIPROT_ENTRY not in identifiers[0]
    assert Namespace.GENESYMBOL not in identifiers[1]
    assert identifiers[1][Namespace.UNIPROT_ENTRY] == 'ERBB2_HUMAN'
    assert identifiers[2][Namespace.GENESYMBOL] == 'ERBB3'


def test_drugcentral_sources_are_serving_annotations():
    out = _extract(
        'drugcentral',
        {
            'STRUCT_ID': '1',
            'ACCESSION': 'P00533',
            'ACT_SOURCE': 'ChEMBL',
            'MOA_SOURCE': 'DrugBank',
            'ACT_SOURCE_URL': 'https://example.org/assay',
            'MOA_SOURCE_URL': 'https://example.org/mechanism',
        },
    )
    annotations = {
        (a['term'], a['value']) for a in out.relations[0].annotations
    }
    assert ('supporting_data_source', 'ChEMBL') in annotations
    assert ('supporting_data_source', 'DrugBank') in annotations
    assert ('source_record_urls', 'https://example.org/assay') in annotations
    assert (
        'source_record_urls',
        'https://example.org/mechanism',
    ) in annotations


@pytest.mark.parametrize('raw', ['unknown', 'NaN', 'inf', '1e999', '-'])
def test_guidetopharma_invalid_quantity_does_not_drop_relation(raw):
    out = _extract(
        'guidetopharma',
        {
            'Target ID': '1',
            'Ligand ID': '1',
            'Action': 'Inhibition',
            '_ligand': {'Ligand ID': '1', 'Type': 'Synthetic organic'},
            '_target': {'Target id': '1', 'Human SwissProt': 'P00533'},
            'Affinity Units': 'pIC50',
            'Affinity High': raw,
            'Affinity Median': '0',
            'Original Affinity Units': 'IC50',
            'Original Affinity High nm': raw,
            'Original Affinity Median nm': '≈3',
        },
    )
    quantities = [
        a['quantity'] for a in out.relations[0].annotations if a.get('quantity')
    ]
    assert len(quantities) == 2
    assert {q['has_numeric_value'] for q in quantities} == {0, 3}
    assert (
        next(q for q in quantities if q['has_numeric_value'] == 3)['comparator']
        == '≈'
    )


@pytest.mark.parametrize(
    'module,row',
    [
        (
            'bindingdb',
            {
                'Number of Protein Chains in Target (>1 implies a multichain complex)': '2',
                'Target Name': 'complex',
                'UniProt (SwissProt) Primary ID of Target Chain 1': 'P00533',
                'UniProt (SwissProt) Primary ID of Target Chain 2': 'P04626',
                'Target Source Organism According to Curator or DataSource': 'Homo sapiens',
            },
        ),
        (
            'drugcentral',
            {
                'ACCESSION': 'P00533|P04626',
                'TARGET_NAME': 'group',
                'ORGANISM': 'Homo sapiens',
            },
        ),
    ],
)
def test_multicomponent_target_retains_taxon(module, row):
    module = importlib.import_module(f'pypath.inputs_v2.{module}')
    target = getattr(
        module, '_target', getattr(module, '_target_entity', None)
    )(row)
    assert any(a.value == 'NCBITaxon:9606' for a in target.annotations)


@pytest.mark.parametrize(
    'subunits,expected', [('2|3', 'macromolecular_complex'), ('', 'named_thing')]
)
def test_guidetopharma_group_accessions_do_not_alias_components(
    subunits, expected
):
    from pypath.inputs_v2.guidetopharma import targets_schema

    entity = targets_schema.build(
        {
            'Target id': '1',
            'Target name': 'receptor group',
            'Subunit id': subunits,
            'Human SwissProt': 'P00533|P04626',
            'Human Ensembl Gene': 'ENSG00000146648|ENSG00000141736',
        }
    )
    assert entity.type == expected
    assert {i.type for i in entity.identifiers} == {
        Namespace.GUIDETOPHARMA_TARGET,
        Namespace.NAME,
    }
    assert any(a.value == 'NCBITaxon:9606' for a in entity.annotations)
