"""OPM: readers of the tables, mappers and extractor on rows of the real exports."""

import inspect
import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import opm
from pypath.inputs_v2.parsers import opm as opm_parser

# Columns of the tables that the reader uses, with rows copied from the exports.
MEMBRANES = """\
id,name,short_name,topology_in,topology_out,description
1,Archaebacterial membrane,Archaebac.,cytoplasmic side ,extracellular side,
3,Bacterial Gram-negative outer membrane,Gram-neg. outer,periplasm,extracellular side,
20,Secreted,Secreted,inner side,outer side,
"""
TYPES = """\
id,name,description
1,Transmembrane,
2,Monotopic/peripheral,
"""
CLASSES = """\
id,name,description,type_id
1,Alpha-helical polytopic,,1
2,Beta-barrel transmembrane,,1
"""
SUPERFAMILIES = """\
id,name,description,comment,tcdb,pfam,type_id,classtype_id
2,Photosystems,,,,,1,1
3,Electron transport chain complex II superfamily,,,5.A.3,CL0328,1,1
109,Viral membrane fusion glycoprotein,,,1.G,,1,11
"""
FAMILIES = """\
id,name,description,tcdb,pfam,type_id,classtype_id,superfamily_id,interpro
1,Light-harvesting complexes from bacteria,,,PF00556,1,1,2,IPR000066
3,Photosynthetic reaction centers from bacteria,,3.E.2,PF00124,1,1,2,IPR000484
960,EF-hand,,8.A.82,PF0036,2,4,3,IPR028846
"""
# The class table of the real export has the 11 classes; the fixture keeps the
# two that the superfamilies above refer to, and class 11 for the 1.G row.
CLASSES += "11,Bitopic proteins,,1\n"

STRUCTURES = '''\
id,pdbid,name,uniprotcode,membrane_id,family_id,thickness,thicknesserror,tilt,tilterror,gibbs,tau
1,"=""1qjp""","Outer membrane protein A (OMPA), disordered loops",OMPA_ECOLI,3,34,25.4,1.5,11,1,-29.5,
30,"=""1ots""",ClC chloride transporter,CLCA_ECOLI HVM38_MOUSE IGKC_MOUSE,2,18,29.7,0.8,0,0,-115.3,
99,"=""1grm""","Gramicidin A, head-to-head dimer, right-handed",,20,80,23.2,4.0,7,9,-28.1,
28,"=""7ahl""",Alpha-hemolysin,HLA_STAAU,20,46,23.5,0.9,0,1,-49.3,
34,"=""2oau""",Mechanosensitive channel MscS open,MSCS_ECOLI,2,19,31.1,0.8,0,1,-181.5,80
'''


def opener(text):
    return SimpleNamespace(result=io.StringIO(text))


class Table:
    """A download that gives a fixed table."""

    def __init__(self, text):
        self.text = text

    def open(self, **_kwargs):
        return opener(self.text)


def terms(membranes=MEMBRANES, superfamilies=SUPERFAMILIES, families=FAMILIES):
    return list(
        opm_parser.iter_terms(
            opener(membranes),
            types=Table(TYPES),
            classes=Table(CLASSES),
            superfamilies=Table(superfamilies),
            families=Table(families),
        )
    )


def structures(text=STRUCTURES):
    return list(opm_parser.iter_structures(opener(text)))


def extract(mapper, rows, dataset):
    out = SilverExtractor('opm', dataset)
    for index, row in enumerate(rows):
        out.process_record(mapper(row), row, f'fixture:{index}', index)
    return out


def relations(out):
    return [r for r in out.relations if r.predicate != 'has_member']


def test_term_reader_gives_membranes_and_the_four_levels_in_order():
    rows = terms()
    assert [r['kind'] for r in rows] == (
        ['membrane'] * 3 + ['type'] * 2 + ['class'] * 3 + ['superfamily'] * 3
        + ['family'] * 3
    )
    assert rows[0] == {
        'kind': 'membrane', 'id': '1', 'name': 'Archaebacterial membrane',
        'short_name': 'Archaebac.', 'inside': 'cytoplasmic side',
        'outside': 'extracellular side',
    }
    superfamily = next(r for r in rows if r['kind'] == 'superfamily' and r['id'] == '3')
    assert superfamily['parent_kind'] == 'class' and superfamily['parent_id'] == '1'
    assert superfamily['xrefs'] == {'pfam': 'CL0328', 'tcdb': '5.A.3'}


def test_term_reader_repairs_trailing_characters_and_leaves_out_a_bad_accession():
    csv = FAMILIES.replace('PF00556', 'PF00556﻿').replace(',8.A.82,', ',8.A.82.,')
    families = {r['id']: r for r in terms(families=csv) if r['kind'] == 'family'}
    assert families['1']['xrefs'] == {'pfam': 'PF00556', 'interpro': 'IPR000066'}
    # "8.A.82." ends with a dot in two rows of the export; "PF0036" has a digit
    # too few in one and is not repaired.
    assert families['960']['xrefs'] == {'interpro': 'IPR028846', 'tcdb': '8.A.82'}
    # A TCDB class number (1.G) is a valid accession.
    superfamily = next(r for r in terms() if r['kind'] == 'superfamily' and r['id'] == '109')
    assert superfamily['xrefs'] == {'tcdb': '1.G'}


def test_term_reader_rejects_a_family_without_its_superfamily_and_a_missing_column():
    families = FAMILIES + '5,Orphan,,,,1,1,9999,\n'
    with pytest.raises(ValueError, match='not in the superfamily table'):
        terms(families=families)
    with pytest.raises(ValueError, match='short_name'):
        terms(membranes='id,name,topology_in,topology_out\n1,A,x,y\n')


def test_structure_reader_reads_the_pdb_formula_and_the_uniprot_codes():
    rows = structures()
    assert rows[0] == {
        'id': '1', 'pdb_id': '1qjp',
        'name': 'Outer membrane protein A (OMPA), disordered loops',
        'entries': ['OMPA_ECOLI'], 'membrane_id': '3', 'family_id': '34',
        'thickness': '25.4', 'thickness_error': '1.5', 'tilt': '11',
        'tilt_error': '1', 'transfer_energy': '-29.5', 'curvature_radius': None,
    }
    assert rows[1]['entries'] == ['CLCA_ECOLI', 'HVM38_MOUSE', 'IGKC_MOUSE']
    assert rows[2]['entries'] == []
    assert rows[4]['curvature_radius'] == '80'


def test_structure_reader_rejects_an_unexpected_pdb_code():
    with pytest.raises(ValueError, match='PDB code'):
        structures(STRUCTURES.replace('"=""1qjp"""', '1qjp'))
    with pytest.raises(ValueError, match='gibbs'):
        structures('id,pdbid,name\n1,x,y\n')


def test_a_term_has_its_namespace_name_and_parent():
    out = extract(opm.map_term, terms(), 'terms')
    by_id = {
        (e.namespace, e.identifier): e for e in out.entities.values() if e.identifiers
    }
    family = by_id[(opm.FAMILY, '3')]
    assert family.entity_type == 'ontology_class'
    assert {i['id'] for i in family.identifiers if i['ns'] == Namespace.NAME} == {
        'Photosynthetic reaction centers from bacteria'
    }
    assert {a['value'] for a in family.annotations if a['term'] == 'xref'} == {
        'Pfam:PF00124', 'InterPro:IPR000484', 'TCDB:3.E.2',
    }
    parents = {
        (
            out.entities[r.subject_entity_key].namespace,
            out.entities[r.subject_entity_key].identifier,
            out.entities[r.object_entity_key].namespace,
            out.entities[r.object_entity_key].identifier,
        )
        for r in out.relations
        if r.predicate == 'subclass_of'
    }
    assert (opm.FAMILY, '3', opm.SUPERFAMILY, '2') in parents
    assert (opm.SUPERFAMILY, '2', opm.CLASS, '1') in parents
    assert (opm.CLASS, '1', opm.TYPE, '1') in parents
    # A membrane has no parent; its sides and short name are kept.
    membrane = by_id[(opm.MEMBRANE, '1')]
    assert not any(
        r.subject_entity_key == membrane.entity_key and r.predicate == 'subclass_of'
        for r in out.relations
    )
    assert {a['term']: a['value'] for a in membrane.annotations} == {
        'opm:topology_in': 'cytoplasmic side',
        'opm:topology_out': 'extracellular side',
    }
    assert any(i['ns'] == Namespace.SYNONYM and i['id'] == 'Archaebac.' for i in membrane.identifiers)


def test_unknown_term_kind_raises():
    with pytest.raises(ValueError, match='Unknown OPM term kind'):
        opm.map_term({'kind': 'order', 'id': '1', 'name': 'x'})


def test_membrane_statement_of_a_structure_of_one_protein():
    out = extract(opm.map_membrane, structures()[:1], 'membrane_localizations')
    (relation,) = relations(out)
    subject = out.entities[relation.subject_entity_key]
    assert (subject.entity_type, subject.namespace, subject.identifier) == (
        'protein', Namespace.UNIPROT_ENTRY, 'OMPA_ECOLI',
    )
    obj = out.entities[relation.object_entity_key]
    assert (obj.entity_type, obj.namespace, obj.identifier) == (
        'ontology_class', opm.MEMBRANE, '3',
    )
    assert relation.predicate == 'located_in'
    annotations = {a['term']: a for a in relation.annotations}
    assert annotations['opm:pdb_id']['value'] == '1qjp'
    thickness = annotations['opm:hydrophobic_thickness']['quantity']
    assert (thickness['has_numeric_value'], thickness['has_unit']) == (25.4, 'Å')
    assert thickness['source_field'] == 'thickness'
    assert annotations['opm:tilt_angle']['quantity']['has_unit'] == '°'
    assert annotations['opm:transfer_energy']['quantity']['has_unit'] == 'kcal/mol'
    # The structure has no radius of curvature.
    assert 'opm:curvature_radius' not in annotations
    # No species is asserted: OPM gives a name, not a taxon.
    assert not any(a['term'] == 'in_taxon' for a in subject.annotations)


def test_structure_of_several_proteins_is_the_pdb_entry_with_its_members():
    out = extract(opm.map_membrane, structures()[1:2], 'membrane_localizations')
    (relation,) = relations(out)
    subject = out.entities[relation.subject_entity_key]
    assert (subject.entity_type, subject.namespace, subject.identifier) == (
        'macromolecular_complex', Namespace.PDB, '1ots',
    )
    members = {
        out.entities[r.object_entity_key].identifier
        for r in out.relations
        if r.predicate == 'has_member'
    }
    assert members == {'CLCA_ECOLI', 'HVM38_MOUSE', 'IGKC_MOUSE'}
    # The membrane is not stated of the members.
    assert {r.subject_entity_key for r in relations(out)} == {subject.entity_key}


def test_structure_without_a_uniprot_code_is_identified_by_its_pdb_code():
    out = extract(opm.map_membrane, structures()[2:3], 'membrane_localizations')
    (relation,) = relations(out)
    subject = out.entities[relation.subject_entity_key]
    assert (subject.entity_type, subject.namespace, subject.identifier) == (
        'protein', Namespace.PDB, '1grm',
    )


def test_secreted_and_undefined_are_associations_not_locations():
    out = extract(opm.map_membrane, structures()[3:4], 'membrane_localizations')
    (relation,) = relations(out)
    assert relation.predicate == 'associated_with'
    assert out.entities[relation.object_entity_key].identifier == '20'
    assert not any(r.predicate == 'located_in' for r in out.relations)


def test_family_statement_goes_to_the_family_only():
    out = extract(opm.map_family, structures()[:1], 'family_annotations')
    (relation,) = relations(out)
    assert relation.predicate == 'associated_with'
    obj = out.entities[relation.object_entity_key]
    assert (obj.entity_type, obj.namespace, obj.identifier) == (
        'ontology_class', opm.FAMILY, '34',
    )
    assert {a['term']: a['value'] for a in relation.annotations} == {'opm:pdb_id': '1qjp'}


def test_statement_object_is_the_entity_of_the_terms_dataset():
    term_out = extract(opm.map_term, terms(), 'terms')
    keys = {
        (e.entity_type, e.namespace, e.identifier)
        for e in term_out.entities.values()
    }
    rows = structures()
    for mapper, dataset, namespace, field in (
        (opm.map_membrane, 'membrane_localizations', opm.MEMBRANE, 'membrane_id'),
        (opm.map_family, 'family_annotations', opm.FAMILY, 'family_id'),
    ):
        out = extract(mapper, rows[:1], dataset)
        (relation,) = relations(out)
        obj = out.entities[relation.object_entity_key]
        reference = (obj.entity_type, obj.namespace, obj.identifier)
        assert reference == ('ontology_class', namespace, rows[0][field])
    # The membrane of the first structure is in the fixture's terms.
    assert ('ontology_class', opm.MEMBRANE, '3') in keys


def test_unknown_uniprot_code_and_unparsable_measurement_raise():
    row = structures()[0]
    with pytest.raises(ValueError, match='UniProt code'):
        opm.map_membrane({**row, 'entries': ['not a code']})
    with pytest.raises(ValueError, match='thickness'):
        opm.map_membrane({**row, 'thickness': '25.4 +- 1'})


def test_uniprot_accession_is_not_an_entry_name():
    out = extract(opm.map_family, [{**structures()[0], 'entries': ['Q9RP17']}], 'family_annotations')
    (relation,) = relations(out)
    subject = out.entities[relation.subject_entity_key]
    assert (subject.namespace, subject.identifier) == (Namespace.UNIPROT, 'Q9RP17')


def test_resource_datasets_and_import_does_not_download():
    # The vocabulary dataset is ``terms``; the statements are named by what they state.
    datasets = opm.resource.datasets()
    assert set(datasets) == {'terms', 'membrane_localizations', 'family_annotations'}
    assert datasets['membrane_localizations'].mapper is opm.map_membrane
    assert datasets['family_annotations'].mapper is opm.map_family
    assert opm.config.name == 'OPM'


def test_module_has_no_reading_code_and_the_terms_dataset_binds_the_parser():
    # The vocabulary tables are read by the parser; the module only names them.
    assert 'opener' not in inspect.getsource(opm)
    parser = opm.resource.terms._raw_parser
    assert parser.func is opm_parser.iter_terms
    assert parser.keywords == {
        'types': opm.download_types,
        'classes': opm.download_classes,
        'superfamilies': opm.download_superfamilies,
        'families': opm.download_families,
    }
