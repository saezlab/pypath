"""Native chemical input contracts; fixtures are source-shaped, no network."""

import importlib
import pytest

from omnipath_build.silver import SilverExtractor


def project(module, mapper, row):
    module = importlib.import_module('pypath.inputs_v2.' + module)
    record = getattr(module, mapper)(row)
    assert record is not None
    ex = SilverExtractor(module.config.id.value, 'fixture')
    ex.process_record(record, row, 'fixture:1', 1)
    return ex


@pytest.mark.parametrize(
    'module,mapper,row',
    [
        (
            'bindingdb',
            'interactions_schema',
            {
                'BindingDB MonomerID': '1',
                'UniProt (SwissProt) Primary ID of Target Chain 1': 'P00533',
                'Target Source Organism According to Curator or DataSource': 'Homo sapiens',
                'PMID': '123',
                'Ki (nM)': '10',
            },
        ),
        (
            'chembl',
            'molecules_schema',
            {
                'chembl_id': 'CHEMBL1',
                'molecule_type': 'Small molecule',
                'pref_name': 'compound',
                'therapeutic_flag': 1,
                'max_phase': 4,
            },
        ),
        (
            'chembl',
            'activities_schema',
            {
                'molecule_chembl_id': 'CHEMBL1',
                'target_chembl_id': 'CHEMBL2',
                'target_type': 'SINGLE PROTEIN',
                'action_type': 'INHIBITOR',
                'standard_value': '10',
                'standard_type': 'IC50',
                'standard_units': 'nM',
                'pubmed_id': '123',
            },
        ),
        (
            'drugcentral',
            'interactions_schema',
            {
                'STRUCT_ID': '1',
                'ACCESSION': 'P00533',
                'ACTION_TYPE': 'INHIBITOR',
                'ORGANISM': 'Homo sapiens',
                'ACT_VALUE': '10',
                'ACT_TYPE': 'IC50',
                'ACT_UNIT': 'nM',
            },
        ),
        (
            'foodb',
            'foods_schema',
            {
                'public_id': 'FOOD00001',
                'name': 'Apple',
                'member_compound_public_id': 'FDB000001',
                'member_content': '1',
                'member_unit': 'mg/100g',
                'ncbi_taxonomy_id': '3750',
            },
        ),
        (
            'hmdb',
            'metabolites_schema',
            {
                'accession': 'HMDB0000001',
                'name': 'metabolite',
                'chemont_ids': 'CHEMONTID:0000001',
                'pubmed_ids': '123',
                'cellular_locations': 'Cytoplasm',
            },
        ),
        (
            'lipidmaps',
            'lipids_schema',
            {
                'LM_ID': 'LMFA0000001',
                'COMMON_NAME': 'lipid',
                'FORMULA': 'C2H6O',
                'EXACT_MASS': '46.1',
            },
        ),
        (
            'metatlas',
            'reactions_schema',
            {
                'human_gem_reaction_id': 'MAR00001',
                'name': 'reaction',
                'reactants': 'MAM00001:c:1',
                'products': 'MAM00002:c:2',
                'enzyme_ensembl': 'ENSG000001',
                'eccodes': '1.1.1.1',
            },
        ),
        (
            'metatlas',
            'enzyme_complexes_schema',
            {'complex_subunits': 'ENSG000001||ENSG000002'},
        ),
        (
            'phenol_explorer',
            'foods_schema',
            {
                'id': '1',
                'name': 'Apple',
                'member_compound_id': '1',
                'member_compound_name': 'compound',
                'member_formula': 'C2H6O',
                'member_mean': '10',
                'member_units': 'mg/100g',
            },
        ),
        (
            'ptfi',
            'foods_schema',
            {
                'specimen_id': 'specimen_FOODON_000001',
                'name': 'Food',
                'member_ptfi_id': 'compound1',
                'member_name': 'molecule',
                'member_formula': 'C2H6O',
                'member_value': '10',
                'member_units': 'mg',
            },
        ),
        (
            'rampdb',
            'analytehaspathway_schema',
            {'rampId': 'RAMP_C_000000001', 'pathwayRampId': 'RAMP_P_000000001'},
        ),
        (
            'refmet',
            'schema',
            {
                'refmet_id': 'RM0108606',
                'refmet_name': 'Acutumidine',
                'pubchem_cid': '442840',
                'super_class': 'Alkaloids',
                'formula': 'C18H22ClNO6',
                'exactmass': '383.113567',
            },
        ),
        (
            'stitch',
            'interactions_schema',
            {
                'chemical_id': '123',
                'protein_id': 'ENSP000001',
                'ncbi_tax_id': '9606',
                'mode': 'activation',
                'action': 'activation',
                'chem_role': 'source',
                'prot_role': 'target',
                'combined_score': '900',
            },
        ),
        (
            'swisslipids',
            'lipids_schema',
            {
                'Lipid ID': 'SLM:000001',
                'Name': 'lipid',
                'Parent': 'SLM:000002',
                'PMID': '123',
                'Charge (pH7.3)': '0',
            },
        ),
    ],
)
def test_source_rows_reach_silver(module, mapper, row):
    ex = project(module, mapper, row)
    assert ex.entities
    for entity in ex.entities.values():
        assert not entity.entity_type.startswith('OM:')
        assert all(
            not str(i['ns']).startswith('OM:') for i in entity.identifiers
        )
    assert all(not r.predicate.startswith('OM:') for r in ex.relations)


def test_unknown_actions_do_not_become_effects_and_measurements_are_not_descriptions():
    row = {
        'molecule_chembl_id': 'CHEMBL1',
        'target_chembl_id': 'CHEMBL2',
        'target_type': 'SINGLE PROTEIN',
        'action_type': 'UNREVIEWED',
        'standard_value': '10',
        'standard_type': 'IC50',
        'standard_units': 'nM',
    }
    ex = project('chembl', 'activities_schema', row)
    assert ex.relations[-1].predicate == 'associated_with'
    assert not any(
        a['term'] in {'object_direction_qualifier', 'description'}
        for a in ex.relations[-1].annotations
    )


def test_stitch_direction_uses_source_endpoint_and_requires_orientation():
    row = {
        'chemical_id': '123',
        'protein_id': 'ENSP000001',
        'mode': 'activation',
        'action': 'activation',
        'chem_role': 'target',
        'prot_role': 'source',
    }
    ex = project('stitch', 'interactions_schema', row)
    relation = ex.relations[-1]
    assert ex.entities[relation.subject_entity_key].entity_type == 'protein'
    assert relation.predicate == 'affects'
    ex = project(
        'stitch',
        'interactions_schema',
        {**row, 'chem_role': 'undirected', 'prot_role': 'undirected'},
    )
    assert ex.relations[-1].predicate == 'associated_with'


def test_ramp_pathway_id_is_a_distinct_endpoint():
    ex = project(
        'rampdb',
        'analytehaspathway_schema',
        {'rampId': 'RAMP_C_000000001', 'pathwayRampId': 'RAMP_P_000000001'},
    )
    assert len(ex.entities) == 2
    assert ex.relations[0].predicate == 'participates_in'


def test_hmdb_locations_are_not_inferred_graph_links():
    ex = project(
        'hmdb',
        'metabolites_schema',
        {
            'accession': 'HMDB0000001',
            'cellular_locations': 'Cytoplasm',
            'diseases': 'Cancer',
        },
    )
    assert not ex.relations


def test_bindingdb_multichain_target_is_not_chain_one():
    row = {
        'BindingDB MonomerID': '1',
        'Number of Protein Chains in Target (>1 implies a multichain complex)': '2',
        'UniProt (SwissProt) Primary ID of Target Chain 1': 'P00533',
        'UniProt (SwissProt) Primary ID of Target Chain 2': 'P04637',
    }
    ex = project('bindingdb', 'interactions_schema', row)
    target = ex.entities[ex.relations[-1].object_entity_key]
    assert target.entity_type == 'macromolecular_complex'
    assert sum(r.predicate == 'has_part' for r in ex.relations) == 2


def test_bindingdb_ic50_bound_is_typed_and_preserved():
    ex = project(
        'bindingdb',
        'interactions_schema',
        {
            'BindingDB MonomerID': '1',
            'UniProt (SwissProt) Primary ID of Target Chain 1': 'P00533',
            'IC50 (nM)': '<=10',
        },
    )
    ann = next(
        a for a in ex.relations[-1].annotations if a['term'] == 'BAO:0000190'
    )
    assert ann['quantity']['has_numeric_value'] == 10
    assert ann['quantity']['has_unit'] == 'nM'
    assert ann['quantity']['comparator'] == '<='
    assert ann['quantity']['source_field'] == 'IC50 (nM)'


def test_food_measurements_keep_member_alignment_units_and_statistic():
    ex = project(
        'foodb',
        'foods_schema',
        {
            'public_id': 'FOOD1',
            'member_compound_public_id': 'FDB1||FDB2',
            'member_content': '1||2',
            'member_min': '0.5||1.5',
            'member_unit': 'mg/100g||ug/g',
        },
    )
    quantities = [[a['quantity'] for a in r.annotations] for r in ex.relations]
    assert [[q['has_numeric_value'] for q in qs] for qs in quantities] == [
        [1, 0.5],
        [2, 1.5],
    ]
    assert [[q['has_unit'] for q in qs] for qs in quantities] == [
        ['mg/100g', 'mg/100g'],
        ['ug/g', 'ug/g'],
    ]
    assert [[q['source_field'] for q in qs] for qs in quantities] == [
        ['member_content', 'member_min']
    ] * 2


def test_ptfi_does_not_infer_formula_from_ordinary_name():
    from pypath.inputs_v2.parsers.ptfi import parse_compound_name

    assert 'formula' not in parse_compound_name('Yohimbine')
    assert parse_compound_name('C28H34O4-470 (PTF28874)') == {
        'formula': 'C28H34O4',
        'ptfi_id': 'PTF28874',
    }


def test_chembl_ensembl_identifiers_retain_molecule_domain():
    from pypath.inputs_v2.chembl import _ensembl_namespace

    assert _ensembl_namespace('ENSP000001') == 'ensp'
    assert _ensembl_namespace('ENSG000001') == 'ensg'
    assert _ensembl_namespace('ENST000001') == 'enst'
    assert _ensembl_namespace('unreviewed') is None


def test_invalid_quantity_does_not_drop_relation():
    ex = project(
        'bindingdb',
        'interactions_schema',
        {
            'BindingDB MonomerID': '1',
            'UniProt (SwissProt) Primary ID of Target Chain 1': 'P00533',
            'IC50 (nM)': '1e999',
        },
    )
    assert len(ex.relations) == 1
    assert not any(
        a['term'] == 'BAO:0000190' for a in ex.relations[0].annotations
    )


def test_drugcentral_accession_group_does_not_guess_a_complex():
    ex = project(
        'drugcentral',
        'interactions_schema',
        {
            'STRUCT_ID': '1',
            'ACCESSION': 'P00533|P04637',
            'TARGET_NAME': 'a complex-sounding name',
            'ACTION_TYPE': 'INHIBITOR',
        },
    )
    assert sum(e.entity_type == 'protein' for e in ex.entities.values()) == 2
    assert (
        ex.entities[ex.relations[-1].object_entity_key].entity_type
        == 'named_thing'
    )
    assert sum(r.predicate == 'has_member' for r in ex.relations) == 2


@pytest.mark.parametrize(
    'module, row, expected',
    [
        (
            'bindingdb',
            {
                'BindingDB MonomerID': '1',
                'UniProt (SwissProt) Primary ID of Target Chain 1': 'P00533',
                'pchembl_value': '7.25',
                'pH': '7.4',
                'Temp (C)': '25',
            },
            {
                'pchembl_value': (7.25, None),
                'pH': (7.4, None),
                'Temp (C)': (25, 'Cel'),
            },
        ),
        (
            'chembl',
            {
                'molecule_chembl_id': 'CHEMBL1',
                'target_chembl_id': 'CHEMBL2',
                'target_type': 'SINGLE PROTEIN',
                'pchembl_value': '6.5',
            },
            {'pchembl_value': (6.5, None)},
        ),
    ],
)
def test_source_defined_assay_quantities_remain_numeric(module, row, expected):
    mapper = (
        'activities_schema' if module == 'chembl' else 'interactions_schema'
    )
    ex = project(module, mapper, row)
    actual = {
        a['quantity']['source_field']: (
            a['quantity']['has_numeric_value'],
            a['quantity']['has_unit'],
        )
        for a in ex.relations[-1].annotations
        if a['term'] == 'has_quantitative_value'
    }
    assert actual == expected
    assert ex.relations[-1].predicate == 'associated_with'


def test_human_gem_parser_and_mapper_preserve_flux_bounds(monkeypatch):
    from pypath.inputs_v2.parsers import metatlas

    monkeypatch.setattr(metatlas, '_metabolite_xrefs', lambda: {})
    parsed = next(
        metatlas._reactions(
            {
                'reactions': [
                    {
                        'id': 'MAR1',
                        'name': 'reaction',
                        'metabolites': {},
                        'lower_bound': -1000.0,
                        'upper_bound': 500.0,
                    }
                ]
            }
        )
    )
    assert parsed['lower_bound'] == -1000.0
    assert parsed['upper_bound'] == 500.0
    for mapper in ('reactions_schema', 'transport_reactions_schema'):
        ex = project('metatlas', mapper, parsed)
        entity = next(iter(ex.entities.values()))
        quantities = {
            a['quantity']['source_field']: a['quantity']
            for a in entity.annotations
            if a['term'] == 'has_quantitative_value'
        }
        assert quantities['lower_bound']['has_numeric_value'] == -1000.0
        assert quantities['upper_bound']['has_numeric_value'] == 500.0
        assert all(q['has_unit'] is None for q in quantities.values())


def test_human_gem_invalid_bound_omits_only_quantity():
    ex = project(
        'metatlas',
        'reactions_schema',
        {
            'human_gem_reaction_id': 'MAR1',
            'lower_bound': float('inf'),
            'upper_bound': 0,
        },
    )
    entity = next(iter(ex.entities.values()))
    quantities = [
        a['quantity']
        for a in entity.annotations
        if a['term'] == 'has_quantitative_value'
    ]
    assert len(quantities) == 1
    assert quantities[0]['source_field'] == 'upper_bound'
    assert quantities[0]['has_numeric_value'] == 0


def test_swisslipids_names_are_supplemental_identifiers():
    from pypath.inputs_v2.swisslipids import lipids_schema

    record = lipids_schema(
        {
            'Lipid ID': 'SLM:000001',
            'Name': 'A lipid',
            'Synonyms*': 'First synonym;Second synonym',
            'Abbreviation*': 'AL',
        }
    )
    ids = [(i.type, i.value) for i in record.identifiers]
    assert ids[0] == ('swisslipids', 'SLM:000001')
    assert ('name', 'A lipid') in ids
    assert ('synonym', 'First synonym') in ids
    assert ('synonym', 'Second synonym') in ids
    assert ('synonym', 'AL') in ids
    assert all(
        a.term not in {'name', 'synonym'} for a in record.annotations or []
    )


def test_nested_and_direct_names_use_identifier_namespaces():
    from pypath.inputs_v2.phenol_explorer import foods_schema
    from pypath.inputs_v2.drugcentral import _target_entity

    food = foods_schema(
        {
            'id': '1',
            'name': 'Food',
            'member_compound_id': '2',
            'member_compound_name': 'Compound',
            'member_synonyms': 'Alias',
        }
    )
    assert food.identifiers[0].type == 'phenol_explorer_food'
    member = food.membership[0].member
    assert member.identifiers[0].type == 'phenol_explorer_compound'
    assert ('name', 'Compound') in [
        (i.type, i.value) for i in member.identifiers
    ]
    assert ('synonym', 'Alias') in [
        (i.type, i.value) for i in member.identifiers
    ]
    group = _target_entity(
        {'ACCESSION': 'P00533|P04637', 'TARGET_NAME': 'Target group'}
    )
    assert group.identifiers[0].type == 'drugcentral.target_group'
    assert (group.identifiers[1].type, group.identifiers[1].value) == ('name', 'Target group')
    assert len(group.membership) == 2


def test_complex_name_does_not_become_an_unidentified_chain_name():
    from pypath.inputs_v2.bindingdb import _target

    target = _target(
        {
            'Target Name': 'Whole complex',
            'Number of Protein Chains in Target (>1 implies a multichain complex)': '2',
            'UniProt (SwissProt) Primary ID of Target Chain 1': 'P00533',
        }
    )
    assert [(i.type, i.value) for i in target.identifiers] == [
        ('name', 'Whole complex')
    ]
    assert len(target.membership) == 1


def test_chembl_parser_accepts_pipeline_dataset_context(monkeypatch):
    from pypath.inputs_v2.parsers import chembl
    queries = []
    def rows(opener, *, query, **kwargs):
        queries.append(query)
        yield {'chembl_id': 'CHEMBL1'}
    monkeypatch.setattr(chembl, 'iter_sqlite', rows)
    assert list(chembl.molecules_parser(None, dataset='molecules', use_parquet_cache=False, max_records=10)) == [{'chembl_id': 'CHEMBL1'}]
    assert queries[0].endswith('LIMIT 10')


def test_chembl_staging_cache_requires_measurement_units(tmp_path):
    import duckdb
    from pypath.inputs_v2.parsers.chembl import _duckdb_cache_compatible
    path = tmp_path / 'chembl.duckdb'
    assert not _duckdb_cache_compatible(path)
    with duckdb.connect(str(path)) as con:
        con.execute('CREATE TABLE activities (standard_value DOUBLE)')
    assert not _duckdb_cache_compatible(path)
    with duckdb.connect(str(path)) as con:
        con.execute('ALTER TABLE activities ADD COLUMN standard_units VARCHAR')
    assert _duckdb_cache_compatible(path)


def test_chembl_sample_does_not_materialize_full_prepared_cache(monkeypatch):
    from pypath.inputs_v2.parsers import chembl
    calls = []
    def full_cache(*args, **kwargs):
        raise AssertionError('Sample attempted complete cache preparation')
    def sqlite_rows(opener, **kwargs):
        calls.append(kwargs['query'])
        yield {'standard_value': 0.5, 'standard_units': 'nM'}
    monkeypatch.setattr(chembl, '_ensure_duckdb_cache', full_cache)
    monkeypatch.setattr(chembl, 'iter_sqlite', sqlite_rows)
    rows = list(chembl._iter_chembl_prepared('activities', None, sqlite_path='source.db', max_records=10))
    assert rows == [{'standard_value': 0.5, 'standard_units': 'nM'}]
    assert calls[0].endswith('LIMIT 10')
