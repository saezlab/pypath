"""Source contract of the hippocampal CA1 neuron signaling network input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.ca1 import resource
from pypath.inputs_v2.parsers.ca1 import COLUMNS, HEADER, iter_links

# Lines of S1.txt.
INHIBITION = 'CAMKII Q9UQM7 P11798 Kinase Cytosol SYNAPSIN P17600 O88935 Vesicle Vesicles _ Phosphorylation 7890745'
ACTIVATION = 'CAMKII Q9UQM7 P11798 Kinase Cytosol MAP2 P11137 P20357 Cytoskeleton Cytosol + Phosphorylation 7890745'
NEUTRAL = 'GABAAR P14867 P18504 Receptor Membrane GABARAP Q9BY60 P60518 Adapter Cytosol 0 Binding 9892355'
# A small molecule has no accession.
NO_ACCESSION = 'GLYCINE NA NA Ligand Extracellular NMDAR Q12879 P35436 Receptor Membrane + Binding 7870047'
# A protein with a mouse accession only.
MOUSE_ONLY = 'MGLUR1 Q13255 P97772 Receptor Membrane TAMALIN NA Q8R4T5 Adapter Cytosol + Binding 15173175'
# The source's location is missing from this line.
MALFORMED = 'CALMODULIN P62158 P62204 Activator ADDUCIN P35611 Q9QYB5 Adapter Cytosol + Binding 8810272'

mapper = resource.datasets()['interactions'].mapper


def table(*lines, header=HEADER):
    text = '\r\n'.join((header, *lines)) + '\r\n'
    return SimpleNamespace(result=io.BytesIO(text.encode('ascii')))


def read(line):
    (row,) = iter_links(table(line))
    return row


def extract(source_row):
    record = mapper(source_row)
    assert record is not None
    out = SilverExtractor('ca1', 'interactions')
    out.process_record(record, source_row, 'fixture:1', 1)
    (relation,) = out.relations
    return (
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
    )


def identifiers(entity):
    return {(i['ns'], i['id']) for i in entity.identifiers}


def annotations(relation, scope='relation'):
    return sorted(
        (a['term'], a['value'])
        for a in relation.annotations
        if a['scope'] == scope
    )


def test_reader_names_the_repeated_columns_by_their_side():
    assert read(INHIBITION) == dict(
        zip(
            COLUMNS,
            (
                'CAMKII', 'Q9UQM7', 'P11798', 'Kinase', 'Cytosol',
                'SYNAPSIN', 'P17600', 'O88935', 'Vesicle', 'Vesicles',
                '_', 'Phosphorylation', '7890745',
            ),
            strict=True,
        )
    )


def test_inhibitory_link_decreases_the_target():
    relation, source, target = extract(read(INHIBITION))
    assert relation.predicate == 'affects'
    assert identifiers(source) == {
        (Namespace.UNIPROT, 'Q9UQM7'),
        (Namespace.NAME, 'CAMKII'),
    }
    assert identifiers(target) == {
        (Namespace.UNIPROT, 'P17600'),
        (Namespace.NAME, 'SYNAPSIN'),
    }
    assert (source.namespace, source.identifier) == ('uniprot', 'Q9UQM7')
    assert {source.entity_type, target.entity_type} == {'protein'}
    assert {source.taxon, target.taxon} == {'9606'}
    assert annotations(relation) == [
        ('object_direction_qualifier', 'decreased'),
        ('publications', 'PMID:7890745'),
    ]


def test_activating_link_increases_the_target():
    relation, _, target = extract(read(ACTIVATION))
    assert relation.predicate == 'affects'
    assert (Namespace.UNIPROT, 'P11137') in identifiers(target)
    assert annotations(relation) == [
        ('object_direction_qualifier', 'increased'),
        ('publications', 'PMID:7890745'),
    ]


def test_neutral_link_is_an_interaction_without_direction():
    relation, _, _ = extract(read(NEUTRAL))
    assert relation.predicate == 'interacts_with'
    assert annotations(relation) == [('publications', 'PMID:9892355')]


def test_mechanism_is_not_turned_into_an_aspect():
    # The sign is the functional effect on the target: an inhibitory
    # phosphorylation must not read as decreased phosphorylation.
    relation, _, _ = extract(read(INHIBITION))
    assert {a['term'] for a in relation.annotations} == {
        'in_taxon',
        'object_direction_qualifier',
        'publications',
    }


def test_component_class_location_and_mouse_accession_are_not_emitted():
    relation, source, target = extract(read(INHIBITION))
    values = {value for _, value in identifiers(source) | identifiers(target)}
    values |= {a['value'] for a in relation.annotations}
    assert not values & {'P11798', 'O88935', 'Kinase', 'Vesicle', 'Cytosol', 'Vesicles'}


@pytest.mark.parametrize('line', [NO_ACCESSION, MOUSE_ONLY])
def test_no_relation_for_a_component_without_a_human_accession(line):
    row = read(line)
    assert mapper(row) is None


def test_missing_accession_is_read_as_none():
    row = read(NO_ACCESSION)
    assert (row['Source'], row['Source Human'], row['Source Mouse']) == (
        'GLYCINE', None, None,
    )


def test_line_lacking_a_field_is_kept_apart_and_not_mapped():
    row = read(MALFORMED)
    assert row['Malformed'] == MALFORMED
    assert all(row[column] is None for column in COLUMNS)
    assert mapper(row) is None


def test_unknown_effect_gets_no_default():
    with pytest.raises(ValueError, match="'-'"):
        mapper(dict(read(ACTIVATION), Effects='-'))


def test_reader_rejects_a_table_with_another_header():
    with pytest.raises(ValueError, match='expected header'):
        list(iter_links(table(ACTIVATION, header=HEADER.replace('Effects', 'Effect'))))


def test_reader_rejects_an_empty_download():
    with pytest.raises(ValueError, match='no rows'):
        list(iter_links(table()))


def test_reader_rejects_an_unknown_effect_on_any_line():
    # The line has no accession and is never mapped; the value is still checked.
    unmapped = NO_ACCESSION.replace(' + Binding', ' - Binding')
    with pytest.raises(ValueError, match="Unknown CA1 effect: '-'"):
        list(iter_links(table(ACTIVATION, unmapped)))
    # The three known values, and the line that lacks a field, are accepted.
    assert len(list(iter_links(table(INHIBITION, ACTIVATION, NEUTRAL, MALFORMED)))) == 4


@pytest.mark.parametrize('number', ['0123456', '123456789'])
def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication(number):
    relation, _, _ = extract(read(ACTIVATION.replace('7890745', number)))
    assert relation.predicate == 'affects'
    assert annotations(relation) == [('object_direction_qualifier', 'increased')]
