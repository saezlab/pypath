"""Interactome3D: reader, mapper and extractor on lines of the human representative table."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import i3d
from pypath.inputs_v2.parsers.i3d import COLUMNS

HEADER = '\t'.join(COLUMNS) + '\n'
# Lines copied from interactions.dat (release 2024_12): a homology model, an
# experimental structure, a model on a domain-domain template, and a structure
# of a protein with itself.
MODEL = (
    'A0A024RBG1\tQ8NFP7\t1\t0\tModel\t5ltu\t1\tB\t0\t92.2\t77.9\t8\t148\t-\tA\t0\t87.1\t84.8\t8\t146\t-\t'
    'A0A024RBG1-Q8NFP7-MDL-5ltu.pdb1-B-0-A-0.pdb\n'
)
STRUCTURE = (
    'A0A075B5G3\tQ8NBP7\t1\t0\tStructure\t4ov6\t1\tF\t0\t100.0\t100.0\t1\t99\t-\tB\t0\t95.6\t42.5\t153\t446\t-\t'
    'A0A075B5G3-Q8NBP7-EXP-4ov6.pdb1-F-0-B-0.pdb\n'
)
DOMAIN_MODEL = (
    'A0AV96\tQ8N6W0\t1\t0\tDom_dom_model\t6nty\t-\tB\t-\t22.1\t10.8\t74\t137\tRRM_1\tA\t-\t30.4\t14.0\t403\t471\tRRM_1\t'
    'A0AV96-Q8N6W0-MDD-RRM_1-RRM_1-6nty-B-23-92-A-23-92.pdb\n'
)
HOMODIMER = (
    'A0A075B6T6\tA0A075B6T6\t1\t0\tStructure\t1ao7\t2\tD\t2\t98.9\t79.6\t23\t112\t-\tD\t3\t98.9\t79.6\t23\t112\t-\t'
    'A0A075B6T6-A0A075B6T6-EXP-1ao7.pdb2-D-2-D-3.pdb\n'
)
TEXT = HEADER + MODEL + STRUCTURE + DOMAIN_MODEL + HOMODIMER


def read(text=TEXT):
    opener = SimpleNamespace(result=io.StringIO(text))
    return list(i3d.resource.interactions._raw_parser(opener))


def extract(rows):
    out = SilverExtractor('i3d', 'interactions')
    for index, row in enumerate(rows):
        record = i3d.resource.interactions.mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def values(relation, term, scope='relation'):
    return [
        a['value'] for a in relation.annotations
        if a['term'] == term and a['scope'] == scope
    ]


def quantities(relation, term, scope):
    return [
        a['quantity'] for a in relation.annotations
        if a['term'] == term and a['scope'] == scope
    ]


def test_module_import_defines_the_dataset_without_downloading():
    assert list(i3d.resource.datasets()) == ['interactions']
    assert i3d.resource.interactions.download.url.endswith(
        '/human/download/representative/interactions.dat'
    )


def test_reader_yields_numbers_and_none_for_the_dash():
    rows = read()
    assert len(rows) == 4
    assert rows[2] == {
        'rank_major': '1',
        'rank_minor': '0',
        'type': 'Dom_dom_model',
        'pdb_id': '6nty',
        'bio_unit': None,
        'filename': 'A0AV96-Q8N6W0-MDD-RRM_1-RRM_1-6nty-B-23-92-A-23-92.pdb',
        'uniprot_a': 'A0AV96',
        'chain_a': 'B',
        'model_a': None,
        'seq_identity_a': '22.1',
        'coverage_a': '10.8',
        'seq_begin_a': 74,
        'seq_end_a': 137,
        'domain_a': 'RRM_1',
        'uniprot_b': 'Q8N6W0',
        'chain_b': 'A',
        'model_b': None,
        'seq_identity_b': '30.4',
        'coverage_b': '14.0',
        'seq_begin_b': 403,
        'seq_end_b': 471,
        'domain_b': 'RRM_1',
    }


def test_reader_rejects_a_table_that_is_not_the_expected_one():
    with pytest.raises(ValueError, match='lacks expected columns'):
        read('PROT1\tPROT2\nA0A024RBG1\tQ8NFP7\n')
    with pytest.raises(ValueError, match='no rows'):
        read(HEADER)
    with pytest.raises(ValueError, match='ends before it begins'):
        read(HEADER + MODEL.replace('\t8\t148\t', '\t148\t8\t'))
    with pytest.raises(ValueError, match='not a UniProt accession'):
        read(HEADER + MODEL.replace('A0A024RBG1', 'NOT_AN_ACC', 1))
    with pytest.raises(ValueError, match='percentage'):
        read(HEADER + MODEL.replace('92.2', '192.2'))


def test_experimental_structure_is_a_direct_interaction():
    (relation,) = extract(read()[1:2]).relations
    assert relation.predicate == 'directly_physically_interacts_with'
    assert values(relation, 'has_confidence_level') == ['experimental structure']
    assert values(relation, 'i3d:pdb_id') == ['4ov6']
    assert values(relation, 'i3d:chain', 'subject') == ['F']
    assert values(relation, 'i3d:chain', 'object') == ['B']


def test_models_do_not_get_the_standing_of_a_structure():
    model, domain_model = extract([read()[0], read()[2]]).relations
    assert model.predicate == 'interacts_with'
    assert values(model, 'has_confidence_level') == ['homology model from a global template']
    assert domain_model.predicate == 'interacts_with'
    assert values(domain_model, 'has_confidence_level') == [
        'model from a domain-domain template'
    ]
    # The PDB entry of a model is the template, which is not the proteins' own.
    assert values(model, 'i3d:pdb_id') == ['5ltu']


def test_sequence_identity_and_coverage_belong_to_each_protein():
    (relation,) = extract(read()[1:2]).relations
    (identity_a,) = quantities(relation, 'i3d:sequence_identity', 'subject')
    (identity_b,) = quantities(relation, 'i3d:sequence_identity', 'object')
    assert (identity_a['has_numeric_value'], identity_b['has_numeric_value']) == (100.0, 95.6)
    assert identity_b['has_unit'] == '%'
    assert identity_b['source_field'] == 'SEQ_IDENT2'
    (coverage_b,) = quantities(relation, 'i3d:coverage', 'object')
    assert coverage_b['has_numeric_value'] == 42.5
    assert coverage_b['source_field'] == 'COVERAGE2'


def test_covered_residues_are_a_region_of_each_protein():
    mapper = i3d.resource.interactions.mapper
    structure = mapper(read()[1])
    assert structure.object.molecular_form['regions'] == [
        {
            'type': 'region covered by the structure',
            'identifier': None,
            'position': 153,
            'end_position': 446,
            'coordinate_reference': {
                'identifier': None,
                'coordinate_system': 'protein',
                'position_base': 1,
            },
            'description': None,
        }
    ]
    model = mapper(read()[2])
    (region,) = model.subject.molecular_form['regions']
    assert region['type'] == 'region covered by the model'
    assert (region['position'], region['end_position']) == (74, 137)
    # A domain-domain model names the domain, as the source does.
    assert region['description'] == 'RRM_1'


def test_participants_are_human_proteins_with_a_uniprot_accession():
    out = extract(read())
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert {e.taxon for e in out.entities.values()} == {'9606'}
    assert {e.namespace for e in out.entities.values()} == {Namespace.UNIPROT}


def test_structure_of_a_protein_with_itself_is_kept_as_a_self_relation():
    (relation,) = extract(read()[3:4]).relations
    # Both chains are the same protein over the same residues: one participant.
    assert relation.subject_entity_key == relation.object_entity_key
    assert relation.predicate == 'directly_physically_interacts_with'


def test_unknown_structure_type_raises_and_no_publication_or_method_is_invented():
    row = dict(read()[0], type='Docking')
    with pytest.raises(ValueError, match='Unknown Interactome3D structure type'):
        i3d.resource.interactions.mapper(row)
    (relation,) = extract(read()[:1]).relations
    terms = {a['term'] for a in relation.annotations}
    assert 'publications' not in terms
    assert 'has_evidence_of_type' not in terms
