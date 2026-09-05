"""Source-specific Biolink contracts for communication resource inputs."""

from importlib import import_module
from types import SimpleNamespace
import io

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import entity_type, slot_name, predicate
from omnipath_core.naming import Namespace


MODULES = 'cellchat cellinker cellphonedb connectomedb guidetopharma icellnet imm1415 mebocost mrclinksdb neuronchat nichenet'.split()


def extract(module, dataset, row):
    resource = import_module('pypath.inputs_v2.' + module).resource
    record = resource.datasets()[dataset].mapper(row)
    assert record is not None
    out = SilverExtractor(module, dataset)
    out.process_record(record, row, 'fixture:1', 1)
    return record, out


@pytest.mark.parametrize('module', MODULES)
def test_module_imports_with_native_schema(module):
    assert import_module('pypath.inputs_v2.' + module).resource.datasets()


@pytest.mark.parametrize(
    'module,dataset,row',
    [
        (
            'cellchat',
            'human_interactions',
            {
                'ligand_name': 'TGFB1',
                'ligand_genes': ['TGFB1'],
                'receptor_name': 'TGFBR1+TGFBR2',
                'receptor_genes': ['TGFBR1', 'TGFBR2'],
                'taxon_id': '9606',
                'interaction_name': 'TGFB1_TGFBR',
                'evidence': 'PMID:123',
            },
        ),
        (
            'cellinker',
            'human_protein_interactions',
            {
                'LRID': 'LRH1',
                'Ligand_Uniprot': 'P04637',
                'ligand_symbol': 'TP53',
                'Receptor_Uniprot': 'P00533',
                'Receptor_symbol': 'EGFR',
                'PMID': '123;456',
            },
        ),
        (
            'cellphonedb',
            'interactions',
            {
                'partner_a': 'P04637',
                'partner_b': 'P00533',
                'source': 'PMID:123;PMC456',
                'interactors': 'TP53-EGFR',
                'directionality': 'ligand-receptor',
            },
        ),
        (
            'connectomedb',
            'interactions',
            {
                'Species': 'human',
                'Ligand Symbols': 'TGFB1',
                'Receptor Symbols': 'TGFBR1',
                'Interaction ID': 'CDB25:1',
                'Ligand Location': 'Not secreted; near membrane',
                'Evidence': 'not direct',
            },
        ),
        (
            'guidetopharma',
            'interactions',
            {
                'Target ID': '1',
                'Ligand ID': '1',
                'Target Species': 'Human',
                'Action': 'Activation',
                '_ligand': {
                    'Ligand ID': '1',
                    'Type': 'Synthetic organic',
                    'PubChem SID': '42',
                    'PubChem CID': '43',
                    'Name': 'Example',
                },
                '_target': {
                    'Target id': '1',
                    'Human SwissProt': 'P00533',
                    'Target name': 'EGFR',
                },
                'PubMed ID': '123',
            },
        ),
        (
            'icellnet',
            'interactions',
            {
                'Ligand 1': 'TIMP1',
                'Ligand 2': 'MMP9',
                'Receptor 1': 'LRP1',
                'Reference': '123',
            },
        ),
        (
            'imm1415',
            'metabolites',
            {
                'name': 'chemical',
                'chebi': ['CHEBI:123'],
                'kegg.drug': ['D00001'],
                'sbo': 'SBO:0000247',
                'compartment': 'cytosol',
            },
        ),
        (
            'mebocost',
            'human',
            {
                'HMDB_ID': 'HMDB0000012',
                'standard_metName': 'Deoxyuridine',
                'metName': 'Deoxyuridine',
                'Gene_name': 'CDA',
                'Protein_name': 'Cytidine deaminase',
                'Evidence': '123; HMDB',
            },
        ),
        (
            'mrclinksdb',
            'human_interactions',
            {
                'hmdb_id': 'HMDB0000012',
                'pubchem_cid_sid': 'SID:999;CID:123',
                'metabolite_name': 'Example',
                'receptor_uniprot_id': 'P00533',
                'receptor_gene_id': '1956',
                'taxon_id': '9606',
                'pmid': '123',
            },
        ),
        (
            'neuronchat',
            'human_interactions',
            {
                'interaction_name': 'GABA_GABRA1',
                'receptor_subunit': ['GABRA1'],
                'ligand_type': 'Small molecule',
                'interaction_type': 'ligand-receptor',
            },
        ),
        (
            'nichenet',
            'human_interactions',
            {'ligand': 'TGFB1', 'receptor': 'TGFBR1', 'taxon_id': '9606'},
        ),
    ],
)
def test_representative_rows_reach_silver_extractor(module, dataset, row):
    record, out = extract(module, dataset, row)
    assert out.entities
    assert all(
        a['term'] not in {'name', 'synonym'}
        for entity in out.entities.values()
        for a in entity.annotations
    )
    assert all(
        a['term'] not in {'name', 'synonym'}
        for relation in out.relations
        for a in relation.annotations
    )
    assert all(
        not e.entity_type.startswith('OM:') for e in out.entities.values()
    )
    assert all(predicate(r.predicate) == r.predicate for r in out.relations)
    for e in out.entities.values():
        assert e.taxon is None or e.taxon.isdigit()


def test_cellchat_unknown_label_does_not_invent_complex_members_or_gene_ids():
    mod = import_module('pypath.inputs_v2.cellchat')
    entity = mod._protein_group_entity(name='A_B_2', genes=[], taxon_id='9606')
    assert entity_type(entity.type) == 'named_thing'
    assert not entity.membership
    assert [(i.type, i.value) for i in entity.identifiers] == [
        (Namespace.NAME, 'A_B_2')
    ]


def test_cellchat_cofactor_effect_is_qualified_and_unknown_effect_rejected():
    mod = import_module('pypath.inputs_v2.cellchat')
    row = {
        'taxon_id': '9606',
        'cofactor_gene': 'FST',
        'target_name': 'TGFB1',
        'target_genes': ['TGFB1'],
        'effect': 'inhibition',
    }
    relation = mod.map_cellchat_cofactor_interaction(row)
    assert slot_name(relation.predicate) == 'affects'
    assert any(
        slot_name(a.term) == 'object_direction_qualifier'
        for a in relation.annotations
    )
    with pytest.raises(ValueError, match='Unknown CellChat'):
        mod.map_cellchat_cofactor_interaction(dict(row, effect='ambiguous'))


def test_guidetopharma_id_domains_and_unknown_action():
    mod = import_module('pypath.inputs_v2.guidetopharma')
    ligand = mod.ligands_schema(
        {
            'Ligand ID': '1',
            'PubChem SID': '2',
            'PubChem CID': '3',
            'Type': 'Synthetic organic',
        }
    )
    target = mod.targets_schema({'Target id': '1', 'Human SwissProt': 'P00533'})
    assert (Namespace.GUIDETOPHARMA, '1') in ligand.identifiers
    assert (Namespace.GUIDETOPHARMA_TARGET, '1') in target.identifiers
    assert (Namespace.PUBCHEM_SUBSTANCE, '2') in ligand.identifiers
    assert (Namespace.PUBCHEM, '3') in ligand.identifiers
    assert (
        slot_name(mod.guidetopharma_predicate({'Action': 'Competitive'}))
        == 'interacts_with'
    )
    assert (
        slot_name(mod.guidetopharma_predicate({'Action': 'Binding'}))
        == 'directly_physically_interacts_with'
    )


def test_mrclinksdb_never_promotes_pubchem_substance_to_compound():
    mod = import_module('pypath.inputs_v2.mrclinksdb')
    for value in ['SID:999', 'SID:999;CID:123', 'garbage']:
        entity = mod.metabolite_builder(
            {'hmdb_id': 'HMDB0000012', 'pubchem_cid_sid': value}
        )
        ids = [
            x.value for x in entity.identifiers if x.type == Namespace.PUBCHEM
        ]
        assert ids == (['123'] if 'CID:123' in value else [])


def test_transporter_is_association_not_fabricated_transport_event():
    _, out = extract(
        'mrclinksdb',
        'human_transporters',
        {'hmdb_id': 'HMDB0000012', 'uniprot_id': 'P00533', 'taxon_id': '9606'},
    )
    assert [r.predicate for r in out.relations] == ['associated_with']
    assert {e.entity_type for e in out.entities.values()} == {
        'protein',
        'chemical_entity',
    }


def test_neuronchat_chemical_has_no_organism_taxonomy_or_inferred_effect():
    relation, out = extract(
        'neuronchat',
        'human_interactions',
        {
            'interaction_name': 'GABA_GABRA1',
            'receptor_subunit': ['GABRA1'],
            'interaction_type': 'inhibition',
        },
    )
    assert slot_name(relation.predicate) == 'interacts_with'
    chemical = next(
        e for e in out.entities.values() if e.entity_type == 'chemical_entity'
    )
    assert chemical.taxon is None
    assert not any(
        slot_name(a.term) == 'object_direction_qualifier'
        for a in relation.annotations
    )


def test_missing_participants_do_not_emit_nichenet_or_icellnet_statements():
    assert (
        import_module('pypath.inputs_v2.nichenet').map_nichenet_interaction(
            {'ligand': 'A'}
        )
        is None
    )
    assert (
        import_module('pypath.inputs_v2.icellnet').map_icellnet_interaction(
            {'Ligand 1': 'A'}
        )
        is None
    )


def test_imm1415_missing_optional_annotations_and_compound_namespace():
    mod = import_module('pypath.inputs_v2.imm1415')
    rows = list(
        mod.parser(
            SimpleNamespace(
                result=io.StringIO(
                    '{"compartments": {}, "metabolites": [{"id":"x","name":"X"}]}'
                )
            )
        )
    )
    assert rows == [{'id': 'x', 'name': 'X', 'compartment': ''}]
    entity = mod.schema({'chebi': ['CHEBI:123'], 'kegg.drug': ['D00001']})
    assert (Namespace.KEGG_DRUG, 'D00001') in entity.identifiers
    assert (Namespace.KEGG, 'D00001') not in entity.identifiers


@pytest.mark.parametrize(
    'measure,term',
    [
        ('pIC50', 'pic50'),
        ('pKi', 'pki'),
        ('pKd', 'pkd'),
        ('pEC50', 'pec50'),
        ('pA2', 'has_quantitative_value'),
        ('pKB', 'has_quantitative_value'),
    ],
)
def test_guidetopharma_quantities_preserve_reported_high_low_median(
    measure, term
):
    _, out = extract(
        'guidetopharma',
        'interactions',
        {
            'Target ID': '1',
            'Ligand ID': '1',
            'Action': 'Inhibition',
            '_ligand': {'Ligand ID': '1', 'Type': 'Synthetic organic'},
            '_target': {'Target id': '1', 'Human SwissProt': 'P00533'},
            'Affinity Units': measure,
            'Affinity High': '8.2',
            'Affinity Low': '6.0',
            'Affinity Median': '7.1',
        },
    )
    annotations = out.relations[0].annotations
    quantities = [a['quantity'] for a in annotations if a['term'] == term]
    assert len(quantities) == 3
    assert {
        (q['source_field'], q['has_numeric_value']) for q in quantities
    } == {
        (f'{measure} Affinity High', 8.2),
        (f'{measure} Affinity Low', 6.0),
        (f'{measure} Affinity Median', 7.1),
    }
    assert all(q['has_unit'] is None for q in quantities)
    assert any(
        a['term'] == 'causal_mechanism_qualifier' and a['value'] == 'inhibition'
        for a in annotations
    )


def test_guidetopharma_original_ic50_and_comparator_remain_numeric():
    _, out = extract(
        'guidetopharma',
        'interactions',
        {
            'Target ID': '1',
            'Ligand ID': '1',
            'Action': 'Inhibition',
            '_ligand': {'Ligand ID': '1', 'Type': 'Synthetic organic'},
            '_target': {'Target id': '1', 'Human SwissProt': 'P00533'},
            'Original Affinity Units': 'IC50',
            'Original Affinity Median nm': '30.0',
            'Original Affinity Relation': '<',
        },
    )
    q = next(
        a['quantity']
        for a in out.relations[0].annotations
        if a['term'] == 'has_quantitative_value'
    )
    assert q['has_numeric_value'] == 30.0
    assert q['has_unit'] == 'nM'
    assert q['comparator'] == '<'
    assert q['source_field'] == 'IC50 Original Affinity Median nm'


def test_guidetopharma_nonhuman_assay_does_not_reuse_human_gene_aliases():
    mod = import_module('pypath.inputs_v2.guidetopharma')
    target = mod._GuidetopharmaTargetMemberBuilder().build(
        {
            'Target ID': '1',
            'Target Species': 'Mouse',
            '_target': {
                'Target id': '1',
                'Human SwissProt': 'P00533',
                'HGNC symbol': 'EGFR',
                'Target name': 'EGFR',
            },
        }
    )
    assert all(
        i.type not in {Namespace.UNIPROT, Namespace.HGNC, Namespace.GENESYMBOL}
        for i in target.identifiers
    )
    out = SilverExtractor('guidetopharma', 'target')
    out.process_record(target, {}, '1', 1)
    assert next(iter(out.entities.values())).taxon == '10090'


def test_cellinker_missing_participants_are_not_anonymous_protein_records():
    mod = import_module('pypath.inputs_v2.cellinker')
    assert (
        mod.resource.human_protein_interactions.mapper(
            {'Ligand_Uniprot': 'P00533'}
        )
        is None
    )


def test_cellphonedb_keeps_fifth_source_complex_subunit():
    mod = import_module('pypath.inputs_v2.cellphonedb')
    row = {
        'complex_name': 'five_member_complex',
        **{f'uniprot_{i}': f'P0000{i}' for i in range(1, 6)},
    }
    entity = mod.complexes_schema(row)
    assert len(entity.membership) == 5


def test_names_and_synonyms_are_supplemental_identifiers():
    mod = import_module('pypath.inputs_v2.guidetopharma')
    record = mod.ligands_schema(
        {
            'Ligand ID': '42',
            'Name': 'Example',
            'Synonyms': 'Alias',
            'IUPAC name': 'Systematic example',
            'INN': 'Short example',
            'Type': 'Synthetic organic',
        }
    )
    assert record.identifiers[0] == (Namespace.GUIDETOPHARMA, '42')
    assert (Namespace.NAME, 'Example') in record.identifiers
    assert (Namespace.SYNONYM, 'Alias') in record.identifiers
    assert (Namespace.SYNONYM, 'Systematic example') in record.identifiers
    assert (Namespace.SYNONYM, 'Short example') in record.identifiers
    assert not record.annotations


def test_direct_cellchat_single_protein_uses_gene_name_not_group_role():
    mod = import_module('pypath.inputs_v2.cellchat')
    record = mod._protein_group_entity(
        name='Source label', genes=['TGFB1'], taxon_id='9606'
    )
    assert record.identifiers == [
        (Namespace.GENESYMBOL, 'TGFB1'),
        (Namespace.NAME, 'TGFB1'),
    ]
    assert all(
        slot_name(a.term) not in {'name', 'synonym'} for a in record.annotations
    )


def test_tabular_download_error_page_is_not_a_successful_empty_resource():
    from io import StringIO
    from types import SimpleNamespace
    import pytest
    from pypath.inputs_v2.parsers.base import iter_tsv
    opener = SimpleNamespace(result=StringIO('<html>\n<head><title>404 Not Found</title></head>\n</html>'))
    with pytest.raises(ValueError, match='HTML document'):
        list(iter_tsv(opener))
    valid = SimpleNamespace(result=StringIO('LRID\tLigand_geneid\nLRH1\t1\n'))
    assert list(iter_tsv(valid)) == [{'LRID': 'LRH1', 'Ligand_geneid': '1'}]
