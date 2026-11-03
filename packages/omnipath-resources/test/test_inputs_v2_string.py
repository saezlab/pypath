"""STRING: the links reader and the three network mappers on lines of the v12.5 files."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.string import (
    CONFIDENCE,
    FUNCTIONAL_CHANNELS,
    PHYSICAL_CHANNELS,
    REGULATORY_CHANNELS,
    REGULATORY_SIGNS,
    functional_schema,
    physical_schema,
    regulatory_direction,
    regulatory_schema,
    resource,
)
from pypath.inputs_v2.parsers.string import iter_links

ARF5, M6PR, FKBP4 = (
    '9606.ENSP00000000233',
    '9606.ENSP00000000412',
    '9606.ENSP00000001008',
)

# 9606.protein.physical.links.detailed.v12.5.txt: both orientations of two pairs.
PHYSICAL = [
    'protein1 protein2 experimental database textmining combined_score',
    f'{ARF5} 9606.ENSP00000040663 191 0 0 191',
    f'{ARF5} 9606.ENSP00000239690 902 0 0 902',
    f'{FKBP4} 9606.ENSP00000231509 976 800 340 996',
    f'9606.ENSP00000040663 {ARF5} 191 0 0 191',
    f'9606.ENSP00000231509 {FKBP4} 976 800 340 996',
    f'9606.ENSP00000239690 {ARF5} 902 0 0 902',
]
# 9606.protein.links.detailed.v12.5.txt.
FUNCTIONAL = [
    'protein1 protein2 neighborhood fusion cooccurence coexpression experimental '
    'database textmining combined_score',
    f'{ARF5} 9606.ENSP00000040663 0 0 0 83 191 0 0 226',
    f'{ARF5} 9606.ENSP00000158762 0 0 0 63 46 0 835 839',
    '9606.ENSP00000205402 9606.ENSP00000227868 61 200 290 836 999 979 777 999',
    f'9606.ENSP00000040663 {ARF5} 0 0 0 83 191 0 0 226',
    f'9606.ENSP00000158762 {ARF5} 0 0 0 63 46 0 835 839',
    '9606.ENSP00000227868 9606.ENSP00000205402 61 200 290 836 999 979 777 999',
]
# 9606.protein.regulatory.links.detailed.v12.5.txt: directed, one row per direction.
REGULATORY = [
    'protein1 protein2 experimental database textmining positive negative combined_score',
    f'{M6PR} 9606.ENSP00000216124 0 750 0 0 0 750',
    f'{FKBP4} 9606.ENSP00000219548 0 0 182 0 182 182',
    f'{FKBP4} 9606.ENSP00000231509 0 0 569 517 260 569',
    f'{FKBP4} 9606.ENSP00000262719 0 0 330 329 0 330',
    f'9606.ENSP00000231509 {FKBP4} 0 0 203 203 0 203',
]
REGULATORY_COLUMNS = (*REGULATORY_CHANNELS, *REGULATORY_SIGNS)


def table(lines):
    return SimpleNamespace(result=io.StringIO('\n'.join(lines) + '\n'))


def extract(schema, row):
    record = schema(row)
    assert record is not None
    out = SilverExtractor('string', 'fixture')
    out.process_record(record, row, 'fixture:1', 1)
    (relation,) = out.relations
    return (
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
    )


def identifiers(entity):
    return {(i['ns'], i['id']) for i in entity.identifiers}


def scores(relation):
    return {
        a['quantity']['source_field']: a['quantity']['has_numeric_value']
        for a in relation.annotations
        if a['term'] == 'has_confidence_score'
    }


def terms(relation, scope='relation'):
    return {a['term'] for a in relation.annotations if a['scope'] == scope}


def test_reader_yields_each_symmetric_pair_once():
    rows = list(iter_links(table(PHYSICAL), columns=PHYSICAL_CHANNELS))
    assert [(row['protein1'], row['protein2']) for row in rows] == [
        (ARF5, '9606.ENSP00000040663'),
        (ARF5, '9606.ENSP00000239690'),
        (FKBP4, '9606.ENSP00000231509'),
    ]
    assert rows[0] == {
        'protein1': ARF5,
        'protein2': '9606.ENSP00000040663',
        'experimental': '191',
        'database': '0',
        'textmining': '0',
        'combined_score': '191',
    }


def test_reader_applies_only_the_named_minimum_score():
    def pairs(minimum_score):
        return [
            row['combined_score']
            for row in iter_links(
                table(PHYSICAL), columns=PHYSICAL_CHANNELS, minimum_score=minimum_score
            )
        ]

    # No threshold unless one is asked for.
    assert pairs(0) == ['191', '902', '996']
    assert pairs(CONFIDENCE['high']) == ['902', '996']
    assert pairs(CONFIDENCE['highest']) == ['902', '996']
    assert pairs(903) == ['996']
    assert CONFIDENCE == {'low': 150, 'medium': 400, 'high': 700, 'highest': 900}


def test_reader_keeps_both_directions_of_a_directed_file():
    rows = list(iter_links(table(REGULATORY), columns=REGULATORY_COLUMNS, directed=True))
    assert len(rows) == 5
    assert (FKBP4, '9606.ENSP00000231509') in {(r['protein1'], r['protein2']) for r in rows}
    assert ('9606.ENSP00000231509', FKBP4) in {(r['protein1'], r['protein2']) for r in rows}


def test_reader_rejects_a_symmetric_file_with_one_orientation_missing():
    with pytest.raises(ValueError, match='not symmetric: 3 pairs.*2 in the other'):
        list(iter_links(table(PHYSICAL[:-1]), columns=PHYSICAL_CHANNELS))


def test_reader_rejects_changed_columns_and_an_empty_file():
    with pytest.raises(ValueError, match='experimental, database, textmining'):
        list(
            iter_links(
                table(['protein1 protein2 combined_score', f'{ARF5} {M6PR} 200']),
                columns=PHYSICAL_CHANNELS,
            )
        )
    with pytest.raises(ValueError, match='no rows'):
        list(iter_links(table(PHYSICAL[:1]), columns=PHYSICAL_CHANNELS))


def test_datasets_name_their_network_and_confidence_level():
    datasets = resource.datasets()
    assert set(datasets) == {
        'human_physical_high_confidence',
        'human_functional_high_confidence',
        'human_regulatory_high_confidence',
    }
    for name, dataset in datasets.items():
        assert dataset._raw_parser.keywords['minimum_score'] == 700
        assert dataset._raw_parser.keywords['directed'] is ('regulatory' in name)
        assert '/9606.protein.' in dataset.download.url
        assert dataset.download.url.endswith('.detailed.v12.5.txt.gz')
    assert 'protein.physical.links.' in datasets['human_physical_high_confidence'].download.url
    assert 'protein.links.' in datasets['human_functional_high_confidence'].download.url
    assert 'protein.regulatory.links.' in datasets['human_regulatory_high_confidence'].download.url
    assert datasets['human_physical_high_confidence'].mapper is physical_schema
    assert datasets['human_functional_high_confidence'].mapper is functional_schema
    assert datasets['human_regulatory_high_confidence'].mapper is regulatory_schema


def test_physical_link_is_an_interaction_with_its_channel_scores():
    rows = list(iter_links(table(PHYSICAL), columns=PHYSICAL_CHANNELS))
    relation, first, second = extract(physical_schema, rows[2])
    assert relation.predicate == 'interacts_with'
    assert identifiers(first) == {
        (Namespace.STRING, '9606.ENSP00000001008'),
        (Namespace.ENSP, 'ENSP00000001008'),
    }
    assert (first.namespace, first.identifier) == ('string', '9606.ENSP00000001008')
    assert (second.namespace, second.identifier) == ('string', '9606.ENSP00000231509')
    assert {first.entity_type, second.entity_type} == {'protein'}
    assert {first.taxon, second.taxon} == {'9606'}
    assert scores(relation) == {
        'combined_score': 996.0,
        'experimental': 976.0,
        'database': 800.0,
        'textmining': 340.0,
    }
    # A symmetric link has no direction, sign or mechanism.
    assert terms(relation) == {'has_confidence_score'}


@pytest.mark.parametrize('value', ['', 'x', '19.1', '-1', '1001', '²', '1e2'])
@pytest.mark.parametrize('column', ['experimental', 'combined_score'])
def test_reader_rejects_a_malformed_score_instead_of_leaving_it_out(column, value):
    # A channel score the mapper's transform cannot read would otherwise be
    # dropped without a trace: the builder swallows the error.
    lines = PHYSICAL[:3] + PHYSICAL[4:]
    index = PHYSICAL[0].split().index(column)
    cells = lines[2].split()
    cells[index] = value
    lines[2] = ' '.join(cells)
    lines = [line for line in lines if line]
    with pytest.raises(ValueError, match=f'{column} score'):
        list(iter_links(table(lines), columns=PHYSICAL_CHANNELS))


def test_reader_checks_the_scores_of_rows_it_does_not_yield():
    # The mirrored orientation and a row below the minimum are checked too.
    mirrored = PHYSICAL[:-1] + [f'9606.ENSP00000239690 {ARF5} 902 0 x 902']
    with pytest.raises(ValueError, match='textmining score'):
        list(iter_links(table(mirrored), columns=PHYSICAL_CHANNELS))
    low = PHYSICAL[:-2] + [f'9606.ENSP00000231509 {FKBP4} 976 800 340 996', f'9606.ENSP00000239690 {ARF5} 902 0 0 902']
    low[1] = f'{ARF5} 9606.ENSP00000040663 191 y 0 191'
    with pytest.raises(ValueError, match='database score'):
        list(iter_links(table(low), columns=PHYSICAL_CHANNELS, minimum_score=700))


def test_channel_without_evidence_has_no_score():
    rows = list(iter_links(table(PHYSICAL), columns=PHYSICAL_CHANNELS))
    relation, _, _ = extract(physical_schema, rows[1])
    # The combined score equals the only channel's score; both are kept.
    assert scores(relation) == {'combined_score': 902.0, 'experimental': 902.0}


def test_functional_link_is_an_association_not_an_interaction():
    rows = list(iter_links(table(FUNCTIONAL), columns=FUNCTIONAL_CHANNELS))
    relation, first, _ = extract(functional_schema, rows[1])
    assert relation.predicate == 'associated_with'
    assert (Namespace.ENSP, 'ENSP00000000233') in identifiers(first)
    assert scores(relation) == {
        'combined_score': 839.0,
        'coexpression': 63.0,
        'experimental': 46.0,
        'textmining': 835.0,
    }
    relation, _, _ = extract(functional_schema, rows[2])
    assert scores(relation) == {
        'combined_score': 999.0,
        'neighborhood': 61.0,
        'fusion': 200.0,
        'cooccurence': 290.0,
        'coexpression': 836.0,
        'experimental': 999.0,
        'database': 979.0,
        'textmining': 777.0,
    }
    assert terms(relation) == {'has_confidence_score'}


def regulatory_rows():
    return list(iter_links(table(REGULATORY), columns=REGULATORY_COLUMNS, directed=True))


def direction(relation):
    return [
        a['value'] for a in relation.annotations if a['term'] == 'object_direction_qualifier'
    ]


def test_regulatory_link_points_from_the_first_protein_to_the_second():
    relation, regulator, target = extract(regulatory_schema, regulatory_rows()[3])
    assert relation.predicate == 'affects'
    assert regulator.identifier == FKBP4
    assert target.identifier == '9606.ENSP00000262719'
    assert direction(relation) == ['increased']
    assert scores(relation) == {
        'combined_score': 330.0,
        'textmining': 330.0,
        'positive': 329.0,
    }
    # What is regulated is not in the file: no aspect or mechanism.
    assert terms(relation) == {'has_confidence_score', 'object_direction_qualifier'}
    # The opposite direction is a separate row and a separate relation.
    reverse, regulator, target = extract(regulatory_schema, regulatory_rows()[4])
    assert (regulator.identifier, target.identifier) == ('9606.ENSP00000231509', FKBP4)
    assert reverse.relation_key != relation.relation_key


def test_negative_regulation_is_a_decrease():
    relation, _, _ = extract(regulatory_schema, regulatory_rows()[1])
    assert direction(relation) == ['decreased']
    assert scores(relation) == {
        'combined_score': 182.0,
        'textmining': 182.0,
        'negative': 182.0,
    }


def test_no_sign_without_a_score_for_exactly_one_sign():
    unsigned, both = regulatory_rows()[0], regulatory_rows()[2]
    relation, _, _ = extract(regulatory_schema, unsigned)
    assert direction(relation) == []
    assert scores(relation) == {'combined_score': 750.0, 'database': 750.0}
    # Scores for both signs: the two scores are kept and no sign is chosen.
    relation, _, _ = extract(regulatory_schema, both)
    assert direction(relation) == []
    assert scores(relation) == {
        'combined_score': 569.0,
        'textmining': 569.0,
        'positive': 517.0,
        'negative': 260.0,
    }
    assert regulatory_direction({'positive': '517', 'negative': '260'}) is None
    assert regulatory_direction({'positive': '0', 'negative': '0'}) is None
    assert regulatory_direction({}) is None


def test_identifier_of_another_genome_is_kept_without_an_ensembl_identifier():
    # Escherichia coli K-12: the local identifier is a locus tag.
    row = {
        'protein1': '511145.b0002',
        'protein2': '511145.b0003',
        'experimental': '0',
        'database': '900',
        'textmining': '0',
        'combined_score': '900',
    }
    _, first, second = extract(physical_schema, row)
    assert identifiers(first) == {(Namespace.STRING, '511145.b0002')}
    assert {first.taxon, second.taxon} == {'511145'}


@pytest.mark.parametrize('column', ['protein1', 'protein2'])
@pytest.mark.parametrize('value', ['', 'ENSP00000000233', 'ARF5'])
def test_no_relation_without_a_string_identifier(column, value):
    row = dict(
        zip(
            PHYSICAL[0].split(),
            PHYSICAL[2].split(),
            strict=True,
        )
    )
    assert physical_schema(dict(row, **{column: value})) is None
