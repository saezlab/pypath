"""PhosphoNetworks high-resolution network: parser, mapper and extractor on source lines."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import phosphonetworks
from pypath.inputs_v2.parsers.phosphonetworks import iter_network

# Lines copied from the source file: substrate lines start with ">", a site
# near the start of the sequence has a window padded with "-".
SOURCE = (
    '>ABCA1\n'
    'SSYGISETTLEEIFL\tT1242\tCSNK2A1\t1.0162\n'
    '>ABI2\n'
    'PPTQKPPSPPMSGKG\tS177\tCDK5\t2.0987\n'
    'PPTQKPPSPPMSGKG\tS177\tDYRK4\t3.5732\n'
    'PPVVPNDYVPSPTRN\tY207\tMET\t1.4211\n'
    '>ACTB\n'
    '-----MDDDIAALVV\tD1\tFAKE1\t1.0\n'
)


def parse(text=SOURCE):
    return list(iter_network(SimpleNamespace(result=io.StringIO(text))))


def extract(rows):
    out = SilverExtractor('phosphonetworks', 'enzyme_substrate')
    for index, row in enumerate(rows):
        record = phosphonetworks.resource.enzyme_substrate.mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def annotations(relation, term):
    return [a for a in relation.annotations if a['term'] == term]


def test_parser_adds_the_substrate_line_to_the_kinase_site_lines():
    rows = parse()
    assert len(rows) == 5
    assert rows[0] == {
        'substrate': 'ABCA1',
        'window': 'SSYGISETTLEEIFL',
        'site': 'T1242',
        'residue': 'T',
        'position': 1242,
        'enzyme': 'CSNK2A1',
        'score': '1.0162',
    }
    assert [(r['substrate'], r['enzyme']) for r in rows[1:4]] == [
        ('ABI2', 'CDK5'),
        ('ABI2', 'DYRK4'),
        ('ABI2', 'MET'),
    ]


def test_parser_rejects_changed_layout():
    with pytest.raises(ValueError, match='comes before any substrate'):
        parse('SSYGISETTLEEIFL\tT1242\tCSNK2A1\t1.0162\n')
    with pytest.raises(ValueError, match='not a kinase-site line'):
        parse('>ABCA1\nSSYGISETTLEEIFL\tT1242\tCSNK2A1\n')
    with pytest.raises(ValueError, match='window does not centre'):
        parse('>ABCA1\nSSYGISETTLEEIFL\tS1242\tCSNK2A1\t1.0\n')
    with pytest.raises(ValueError, match='unreadable site'):
        parse('>ABCA1\nSSYGISETTLEEIFL\tT\tCSNK2A1\t1.0\n')
    with pytest.raises(ValueError, match='no kinase-site lines'):
        parse('>ABCA1\n')
    with pytest.raises(ValueError, match='not a kinase-site line'):
        parse('<html>\n<body>404</body>\n')


def test_kinase_increases_substrate_phosphorylation_at_the_site():
    record = phosphonetworks.resource.enzyme_substrate.mapper(parse()[0])
    assert slot_name(record.predicate) == 'affects'
    assert record.subject.identifiers == [(Namespace.GENESYMBOL, 'CSNK2A1')]
    assert record.object.identifiers == [(Namespace.GENESYMBOL, 'ABCA1')]
    assert record.subject.molecular_form is None
    assert record.object.molecular_form['modifications'] == [
        {
            'term': 'phosphorylated residue',
            'residue': 'T',
            'position': 1242,
            'end_position': 1242,
            'coordinate_reference': {
                'identifier': None,
                'coordinate_system': 'protein',
                'position_base': 1,
            },
            'description': 'T1242 (SSYGISETTLEEIFL)',
        }
    ]

    out = extract(parse()[:1])
    (relation,) = out.relations
    assert relation.predicate == 'affects'
    (aspect,) = annotations(relation, 'object_aspect_qualifier')
    (direction,) = annotations(relation, 'object_direction_qualifier')
    assert (aspect['value'], direction['value']) == ('phosphorylation', 'increased')
    (score,) = annotations(relation, 'has_confidence_score')
    assert score['quantity']['has_numeric_value'] == 1.0162
    assert score['quantity']['source_field'] == 'score'
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert {e.taxon for e in out.entities.values()} == {'9606'}
    # The file cites no publication per line.
    assert annotations(relation, 'publications') == []


def test_each_kinase_site_line_is_its_own_relation():
    out = extract(parse()[1:3])
    # Same substrate site, two kinases: two relations and one substrate form.
    assert len(out.relations) == 2
    assert out.relations[0].object_entity_key == out.relations[1].object_entity_key
    assert out.relations[0].relation_key != out.relations[1].relation_key
    assert [annotations(r, 'has_confidence_score')[0]['quantity']['has_numeric_value'] for r in out.relations] == [2.0987, 3.5732]


def test_padded_window_is_kept_as_written_and_nothing_is_invented():
    rows = parse()
    out = extract(rows)
    assert len(out.relations) == 5
    form = next(
        e.molecular_form
        for e in out.entities.values()
        if e.molecular_form and e.molecular_form['modifications'][0]['position'] == 1
    )
    assert form['modifications'][0]['description'] == 'D1 (-----MDDDIAALVV)'
    for relation in out.relations:
        # Phosphorylation is not an effect on activity, and no reference, mechanism or level.
        assert [a['value'] for a in annotations(relation, 'object_aspect_qualifier')] == ['phosphorylation']
        assert annotations(relation, 'publications') == []
        assert annotations(relation, 'has_confidence_level') == []
    # A row without a kinase gives no relation.
    assert phosphonetworks.resource.enzyme_substrate.mapper(dict(rows[0], enzyme=None)) is None
