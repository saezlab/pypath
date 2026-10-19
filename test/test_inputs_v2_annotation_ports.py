"""Source fields the move to Biolink had dropped; fixtures are source-shaped, no network."""

import importlib

from omnipath_build.silver import SilverExtractor


def project(module, mapper, row):
    module = importlib.import_module('pypath.inputs_v2.' + module)
    record = getattr(module, mapper)(row)
    assert record is not None
    ex = SilverExtractor(module.config.id.value, 'fixture')
    ex.process_record(record, row, 'fixture:1', 1)
    return ex


def statements(ex):
    """(subject id, predicate, object namespace, object id) of every relation."""
    out = set()
    for r in ex.relations:
        subject = ex.entities[r.subject_entity_key].identifiers[0]['id']
        obj = ex.entities[r.object_entity_key].identifiers[0]
        out.add((subject, r.predicate, str(obj['ns']), obj['id']))
    return out


def annotations(ex):
    return {
        (a['term'], a['value'])
        for holder in [*ex.entities.values(), *ex.relations]
        for a in (holder.annotations or [])
    }


def test_lipidmaps_classes_are_coded_terms_of_every_level():
    from pypath.inputs_v2 import _class_terms, lipidmaps

    row = {
        'LM_ID': 'LMFA01010001',
        'CATEGORY': 'Fatty Acyls [FA]',
        'MAIN_CLASS': 'Fatty Acids and Conjugates [FA01]',
        'SUB_CLASS': 'Straight chain fatty acids [FA0101]',
    }
    ex = project('lipidmaps', 'lipids_schema', row)
    assert {s[2:] for s in statements(ex)} == {
        ('lipidmaps_class', 'FA'),
        ('lipidmaps_class', 'FA01'),
        ('lipidmaps_class', 'FA0101'),
    }
    terms = list(_class_terms.iter_class_terms([row, row], lipidmaps.CLASS_LEVELS))
    assert [(t['id'], t['name'], t['parent']) for t in terms] == [
        ('FA', 'Fatty Acyls', None),
        ('FA01', 'Fatty Acids and Conjugates', 'FA'),
        ('FA0101', 'Straight chain fatty acids', 'FA01'),
    ]


def test_class_term_is_below_its_parent():
    from pypath.inputs_v2._class_terms import map_class_term

    record = map_class_term(
        {
            'namespace': 'refmet_main_class',
            'id': 'Fatty acids',
            'name': 'Fatty acids',
            'parent_namespace': 'refmet_super_class',
            'parent': 'Fatty Acyls',
        }
    )
    ex = SilverExtractor('refmet', 'terms')
    ex.process_record(record, {}, 'terms:1', 1)
    assert statements(ex) == {
        ('Fatty acids', 'subclass_of', 'refmet_super_class', 'Fatty Acyls')
    }


def test_refmet_levels_have_namespaces_of_their_own_and_empty_levels_are_skipped():
    ex = project(
        'refmet',
        'schema',
        {'refmet_id': 'RM1', 'super_class': 'Fatty Acyls', 'main_class': 'Fatty acids', 'sub_class': ''},
    )
    assert {s[2:] for s in statements(ex)} == {
        ('refmet_super_class', 'Fatty Acyls'),
        ('refmet_main_class', 'Fatty acids'),
    }


def test_swisslipids_class_components_and_level():
    ex = project(
        'swisslipids',
        'lipids_schema',
        {
            'Lipid ID': 'SLM:000000178',
            'Level': 'Isomeric subspecies',
            'Lipid class*': 'SLM:000000002',
            'Parent': 'SLM:000392021',
            'Components*': 'SLM:000000418 (sn2) / SLM:000000510 (sn1)',
        },
    )
    assert {s[1:] for s in statements(ex)} == {
        ('subclass_of', 'swisslipids', 'SLM:000392021'),
        ('subclass_of', 'swisslipids', 'SLM:000000002'),
        ('has_part', 'swisslipids', 'SLM:000000418'),
        ('has_part', 'swisslipids', 'SLM:000000510'),
    }
    assert {('swisslipids:position', 'sn1'), ('swisslipids:position', 'sn2')} <= annotations(ex)
    assert ('swisslipids:level', 'Isomeric subspecies') in annotations(ex)


def test_classyfire_names_are_chemont_terms_only_when_unambiguous(monkeypatch):
    from pypath.inputs_v2 import _class_terms

    monkeypatch.setattr(
        _class_terms,
        'chemont_labels',
        lambda opener: {'fatty acyls': 'CHEMONTID:0000262', 'lipids and lipid-like molecules': 'CHEMONTID:0000012'},
    )
    monkeypatch.setattr(_class_terms, 'chemont_download', type('D', (), {'open': lambda self, **kw: None})())
    parser = _class_terms.with_chemont(
        lambda opener, **kw: iter([{'superklass': 'Lipids and lipid-like molecules', 'klass': 'Fatty Acyls', 'subklass': 'Unknown class'}]),
        ('superklass', 'klass', 'subklass'),
    )
    row = next(parser(None))
    assert row['chemont_ids'] == ['CHEMONTID:0000012', 'CHEMONTID:0000262']
    ex = project('foodb', 'compounds_schema', {'public_id': 'FDB000001', **row})
    assert {s[2:] for s in statements(ex)} == {
        ('chemont', 'CHEMONTID:0000012'),
        ('chemont', 'CHEMONTID:0000262'),
    }


def test_foodb_food_groups():
    ex = project(
        'foodb',
        'foods_schema',
        {'public_id': 'FOOD00001', 'name': 'Apple', 'food_group': 'Fruits', 'food_subgroup': 'Pomes'},
    )
    assert {s[2:] for s in statements(ex)} == {
        ('foodb_food_group', 'Fruits'),
        ('foodb_food_subgroup', 'Pomes'),
    }


def test_phenol_explorer_member_classes_and_measurement_context():
    ex = project(
        'phenol_explorer',
        'foods_schema',
        {
            'id': '1',
            'name': 'Apple',
            'food_group': 'Fruits and fruit products',
            'food_subgroup': 'Fruits - Pomes',
            'botanical_family': 'Rosaceae',
            'member_compound_id': '10||20',
            'member_compound_class': 'Flavonoids||Phenolic acids',
            'member_compound_subclass': 'Flavonols||Hydroxybenzoic acids',
            'member_aglycones': 'Quercetin // Kaempferol||-',
            'member_n': '3||1',
            'member_N': '5||1',
            'member_experimental_method': 'Chromatography||Folin assay',
            'member_mean': '1.0||2.0',
            'member_units': 'mg/100 g||mg/100 g',
        },
    )
    s = statements(ex)
    assert {
        ('10', 'associated_with', 'phenol_explorer_compound_class', 'Flavonoids'),
        ('10', 'associated_with', 'phenol_explorer_aglycone', 'Quercetin'),
        ('10', 'associated_with', 'phenol_explorer_aglycone', 'Kaempferol'),
        ('20', 'associated_with', 'phenol_explorer_compound_subclass', 'Hydroxybenzoic acids'),
        ('1', 'associated_with', 'phenol_explorer_botanical_family', 'Rosaceae'),
    } <= s
    assert not any(x[3] == '-' for x in s)
    by_member = {
        ex.entities[r.object_entity_key].identifiers[0]['id']: {
            a['term']: a['value'] for a in r.annotations if a['term'] != 'PATO:0000033'
        }
        for r in ex.relations
        if r.predicate == 'has_part'
    }
    assert by_member == {
        '10': {
            'supporting_study_size': '3',
            'phenol_explorer:data_point_count': '5',
            'supporting_study_method_description': 'Chromatography',
        },
        '20': {
            'supporting_study_size': '1',
            'phenol_explorer:data_point_count': '1',
            'supporting_study_method_description': 'Folin assay',
        },
    }


def test_ptfi_categories_skip_uncategorized():
    ex = project(
        'ptfi',
        'foods_schema',
        {
            'specimen_id': 'S1',
            'num_samples': '3',
            'member_ptfi_id': 'PTF1||PTF2',
            'member_category': 'Triacylglycerol (TG)||Uncategorized',
            'member_n_samples': '2||1',
        },
    )
    assert ('PTF1', 'associated_with', 'ptfi_compound_category', 'Triacylglycerol (TG)') in statements(ex)
    assert not any(x[3] == 'Uncategorized' for x in statements(ex))
    assert ('ptfi:num_samples', '3') in annotations(ex)
    assert ('supporting_study_size', '2') in annotations(ex)


def test_ramp_class_levels_and_hmdb_status():
    from pypath.inputs_v2 import rampdb

    rows = list(
        rampdb._class_rows(
            [
                {
                    'ramp_id': 'RAMP_C_1',
                    'class_source_id': 'hmdb:HMDB0000001',
                    'class_level_name': 'LipidMaps_main_class',
                    'class_name': 'Fatty Acids and Conjugates [FA01]',
                    'source': 'hmdb',
                },
                {
                    'ramp_id': 'RAMP_C_1',
                    'class_source_id': 'hmdb:HMDB0000001',
                    'class_level_name': 'ClassyFire_class',
                    'class_name': 'Fatty Acyls',
                    'source': 'hmdb',
                    'chemont_ids': ['CHEMONTID:0000262'],
                },
            ]
        )
    )
    assert {s[2:] for s in statements(project('rampdb', 'metabolite_class_schema', rows[0]))} == {
        ('lipidmaps_class', 'FA01')
    }
    assert {s[2:] for s in statements(project('rampdb', 'metabolite_class_schema', rows[1]))} == {
        ('chemont', 'CHEMONTID:0000262')
    }
    source = {
        'rampId': 'RAMP_C_1',
        'IDtype': 'hmdb',
        'sourceId': 'hmdb:HMDB0000001',
        'geneOrCompound': 'compound',
    }
    assert ('rampdb:priority_hmdb_status', 'quantified') in annotations(
        project('rampdb', 'source_schema', {**source, 'priorityHMDBStatus': 'quantified'})
    )
    assert not annotations(
        project('rampdb', 'source_schema', {**source, 'priorityHMDBStatus': 'no_HMDB_status'})
    )


def test_mrclinksdb_metabolite_classes_are_chemont_terms():
    ex = project(
        'mrclinksdb',
        'interactions_schema',
        {
            'hmdb_id': 'HMDB0000012',
            'metabolite_name': 'Deoxyuridine',
            'receptor_uniprot_id': 'Q14542',
            'taxon_id': '9606',
            'chemont_ids': ['CHEMONTID:0000495'],
        },
    )
    assert ('HMDB0000012', 'associated_with', 'chemont', 'CHEMONTID:0000495') in statements(ex)


def test_chembl_flags_phase_and_withdrawal():
    ex = project(
        'chembl',
        'molecules_schema',
        {
            'chembl_id': 'CHEMBL1',
            'molecule_type': 'Small molecule',
            'max_phase': 0.5,
            'therapeutic_flag': 1,
            'withdrawn_flag': 1,
            'polymer_flag': 0,
            'inorganic_flag': -1,
            'natural_product': 0,
            'alogp': 3.25,
        },
    )
    assert annotations(ex) >= {
        ('chembl:max_phase', '0.5'),
        ('chembl:therapeutic_flag', 'true'),
        ('chembl:polymer_flag', 'false'),
        ('chembl:natural_product', 'false'),
        ('chembl:alogp', '3.25'),
        ('drug_regulatory_status_world_wide', 'post_approval_withdrawal'),
    }
    assert not any(term == 'chembl:inorganic_flag' for term, _ in annotations(ex))


def test_chembl_molecule_never_withdrawn_has_no_regulatory_status():
    ex = project('chembl', 'molecules_schema', {'chembl_id': 'CHEMBL1', 'withdrawn_flag': 0})
    assert not any(term == 'drug_regulatory_status_world_wide' for term, _ in annotations(ex))


def test_guidetopharma_ligand_properties_and_antibacterial_role():
    ex = project(
        'guidetopharma',
        'ligand_entities_schema',
        {
            'Ligand ID': '1',
            'Name': 'Penicillin',
            'Type': 'Synthetic organic',
            'Approved': 'yes',
            'Withdrawn': '',
            'Labelled': '',
            'Radioactive': '',
            'Antibacterial': 'yes',
        },
    )
    assert annotations(ex) >= {
        ('guidetopharma:approved', 'true'),
        ('guidetopharma:labelled', 'false'),
        ('guidetopharma:radioactive', 'false'),
    }
    assert not any(term == 'drug_regulatory_status_world_wide' for term, _ in annotations(ex))
    assert {s[1:] for s in statements(ex)} == {('has_chemical_role', 'chebi', 'CHEBI:33282')}


def test_guidetopharma_family_has_its_targets():
    from pypath.inputs_v2 import guidetopharma

    family = {
        'id': '1',
        'name': '5-Hydroxytryptamine receptors',
        'targets': [{'Type': 'gpcr', 'Target id': '1', 'Human SwissProt': 'P08908'}],
    }
    ex = SilverExtractor('guidetopharma', 'families')
    ex.process_record(guidetopharma.map_family(family), family, 'families:1', 1)
    assert statements(ex) == {('1', 'has_member', 'uniprot', 'P08908')}


def test_stitch_stereospecific_flag():
    from pypath.inputs_v2 import stitch

    cv = next(
        cv for cv in stitch.interactions_schema.annotations.cvs
        if str(getattr(cv.term_source, 'value', '')) == 'stitch:stereospecific'
    )
    assert cv.value_source.extract({'stereospecific': 'True'}) == ['true']


def _record_statements(module, mapper, row, dataset='fixture'):
    module = importlib.import_module('pypath.inputs_v2.' + module)
    ex = SilverExtractor(module.config.id.value, dataset)
    ex.process_record(getattr(module, mapper)(row), row, f'{dataset}:1', 1)
    return statements(ex)


def test_tcdb_family_has_its_proteins():
    assert _record_statements(
        'tcdb', 'map_family', {'id': '2.A.1', 'name': 'MFS', 'members': ['P11166']}
    ) == {('2.A.1', 'has_member', 'uniprot', 'P11166')}


def test_slc_family_has_its_transporters_and_prose_stays_text():
    row = {
        'SLC_family': 'SLC1',
        'SLC name': 'SLC1A1',
        'Protein name': 'EAAC1, EAAT3',
        'Transport type*': 'C / Na+, H+, K+',
        'Substrates': 'L-Glu, D/L-Asp',
        'Tissue and cellular expression': 'brain (neurons), intestine',
    }
    assert _record_statements('slctables', 'map_family', {'family': 'SLC1', 'rows': [row]}) == {
        ('SLC1', 'has_member', 'genesymbol', 'SLC1A1')
    }
    ex = project('slctables', 'schema', row)
    assert not ex.relations
    assert annotations(ex) >= {
        ('slctables:transport_type', 'C / Na+, H+, K+'),
        ('slctables:substrates', 'L-Glu, D/L-Asp'),
        ('slctables:tissue_expression', 'brain (neurons), intestine'),
    }


def test_corum_funcat_categories():
    row = {
        'ComplexID': '1',
        'subunits(UniProt IDs)': 'P41182',
        'FunCat ID': '10.01.09.05;42.10.03',
        'FunCat description': 'DNA conformation modification (e.g. chromatin);organization of chromosome structure',
    }
    assert {s[2:] for s in statements(project('corum', 'complexes_schema', row)) if s[1] == 'associated_with'} == {
        ('funcat', '10.01.09.05'),
        ('funcat', '42.10.03'),
    }
    from pypath.inputs_v2 import corum

    assert [(str(i.type), i.value) for i in corum.map_funcat({'id': '42.10.03', 'name': 'organization of chromosome structure'}).identifiers] == [
        ('funcat', '42.10.03'),
        ('name', 'organization of chromosome structure'),
    ]
