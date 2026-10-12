"""Source-specific Biolink contract for the DoRothEA input."""

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import dorothea as mod
from pypath.inputs_v2.parsers import dorothea as parser

_COLUMNS = (
    'tf',
    'target',
    'mor',
    'confidence',
    'is_evidence_chip_seq',
    'is_evidence_curated',
    'is_evidence_inferred',
    'is_evidence_tfbs',
    'which_chip_seq',
    'which_curated',
    'which_inferred',
    'which_tfbs',
    'pubmed_id',
)
# Rows of the entire_database data frame, as the shared R data reader yields
# them.
_ROWS = [
    ('AHR', 'CYP1A1', 1.0, 'A', False, True, False, True, 'none', 'htri_db,trrust', 'none', 'hocomoco', '9698073,17785579,10359656,16972788,12213390,9865727'),
    ('AR', 'ACP3', 0.0, 'A', False, True, False, False, 'none', 'int_act,tred_via_reg_network', 'none', 'none', 'imex:IM-18734|pubmed:23518348'),
    ('CTCF', 'EPX', 1.0, 'A', False, True, False, False, 'none', 'oreganno,pazar', 'none', 'none', 'N/A'),
    ('AHR', 'ASAP1', 0.0, 'B', True, False, True, True, 'remap', 'none', 'gtex', 'hocomoco', '-'),
    ('AHR', 'CTNNA2', 0.0, 'C', True, False, False, True, 'remap', 'none', 'none', 'hocomoco', '-'),
    ('AHR', 'CA9', -1.0, 'D', False, True, False, False, 'none', 'trrust', 'none', 'none', '19154183'),
    ('ADNP', 'ABCC1', 0.0, 'D', True, False, False, False, 'remap', 'none', 'none', 'none', '-'),
    ('ADNP', 'AASDH', 0.0, 'E', False, False, True, False, 'none', 'none', 'gtex', 'none', '-'),
    ('AHR', 'ASAP2', 0.0, 'E', False, False, True, True, 'none', 'none', 'gtex', 'hocomoco', '-'),
]
ROWS = {(row[0], row[1]): dict(zip(_COLUMNS, row)) for row in _ROWS}


def extract(row):
    record = mod.interactions_schema(row)
    assert record is not None
    out = SilverExtractor('dorothea', 'level_a')
    out.process_record(record, row, 'fixture:1', 1)
    return record, out


def annotations(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def test_module_has_one_dataset_per_confidence_level():
    assert list(mod.resource.datasets()) == [
        'level_a',
        'level_b',
        'level_c',
        'level_d',
        'level_e',
    ]


def test_activation_with_curated_and_motif_evidence():
    _, out = extract(ROWS['AHR', 'CYP1A1'])
    relation = out.relations[0]
    assert relation.predicate == 'affects'
    assert annotations(relation, 'object_direction_qualifier') == ['increased']
    assert annotations(relation, 'object_aspect_qualifier') == ['expression']
    assert annotations(relation, 'has_confidence_level') == ['A']
    assert annotations(relation, 'has_evidence_of_type') == ['ECO:0005536']
    assert annotations(relation, 'supporting_data_source') == [
        'htri_db',
        'trrust',
        'hocomoco',
    ]
    assert annotations(relation, 'publications') == [
        'PMID:9698073',
        'PMID:17785579',
        'PMID:10359656',
        'PMID:16972788',
        'PMID:12213390',
        'PMID:9865727',
    ]
    subject = out.entities[relation.subject_entity_key]
    target = out.entities[relation.object_entity_key]
    assert (subject.entity_type, subject.namespace, subject.identifier) == (
        'protein',
        Namespace.GENESYMBOL,
        'AHR',
    )
    assert (target.entity_type, target.namespace, target.identifier) == (
        'gene',
        Namespace.GENESYMBOL,
        'CYP1A1',
    )
    assert subject.taxon == target.taxon == '9606'


def test_repression_is_decreased_expression():
    _, out = extract(ROWS['AHR', 'CA9'])
    relation = out.relations[0]
    assert relation.predicate == 'affects'
    assert annotations(relation, 'object_direction_qualifier') == ['decreased']
    assert annotations(relation, 'object_aspect_qualifier') == ['expression']
    assert annotations(relation, 'has_confidence_level') == ['D']
    assert not annotations(relation, 'has_evidence_of_type')


@pytest.mark.parametrize(
    'pair,predicate',
    [
        # Curated: a regulatory statement without a usable sign.
        (('AR', 'ACP3'), 'affects'),
        # ChIP-seq or a motif match, with or without inference: binding.
        (('AHR', 'ASAP1'), 'interacts_with'),
        (('AHR', 'CTNNA2'), 'interacts_with'),
        (('ADNP', 'ABCC1'), 'interacts_with'),
        (('AHR', 'ASAP2'), 'interacts_with'),
        # Inference from expression alone: an association.
        (('ADNP', 'AASDH'), 'associated_with'),
    ],
)
def test_unknown_mode_of_regulation_gets_no_direction(pair, predicate):
    record, out = extract(ROWS[pair])
    relation = out.relations[0]
    assert relation.predicate == predicate
    assert not annotations(relation, 'object_direction_qualifier')
    assert not annotations(relation, 'object_aspect_qualifier')
    # A source-prefixed term (dorothea:evidence) is a string, not a slot.
    assert all(
        'qualifier' not in (a.term if isinstance(a.term, str) else slot_name(a.term))
        for a in record.annotations
    )


@pytest.mark.parametrize(
    'pair,predicate,kinds',
    [
        (('AHR', 'CYP1A1'), 'affects', ['curated', 'tfbs']),
        (('AHR', 'CA9'), 'affects', ['curated']),
        (('AR', 'ACP3'), 'affects', ['curated']),
        # Curated evidence wins over binding and inference.
        (dict(ROWS['AHR', 'ASAP1'], is_evidence_curated=True), 'affects',
         ['curated', 'chip_seq', 'tfbs', 'inferred']),
        # Binding wins over inference.
        (('AHR', 'ASAP1'), 'interacts_with', ['chip_seq', 'tfbs', 'inferred']),
        (('ADNP', 'ABCC1'), 'interacts_with', ['chip_seq']),
        (('AHR', 'ASAP2'), 'interacts_with', ['tfbs', 'inferred']),
        (('ADNP', 'AASDH'), 'associated_with', ['inferred']),
    ],
)
def test_predicate_follows_the_strongest_line_of_evidence(pair, predicate, kinds):
    row = pair if isinstance(pair, dict) else ROWS[pair]
    _, out = extract(row)
    relation = out.relations[0]
    assert relation.predicate == predicate
    # The flags travel with the relation, so the rule can be applied again.
    assert annotations(relation, 'dorothea:evidence') == kinds
    assert annotations(relation, 'has_confidence_level') == [row['confidence']]


def test_binding_and_inferred_pairs_keep_confidence_sources_and_eco_terms():
    _, out = extract(ROWS['AHR', 'ASAP1'])
    relation = out.relations[0]
    assert annotations(relation, 'has_evidence_of_type') == ['ECO:0000229', 'ECO:0005536']
    assert annotations(relation, 'supporting_data_source') == ['remap', 'hocomoco', 'gtex']
    _, out = extract(ROWS['ADNP', 'AASDH'])
    assert annotations(out.relations[0], 'supporting_data_source') == ['gtex']


def test_a_sign_without_curated_evidence_is_an_error_not_dropped():
    # The table gives mor only to curated pairs; a signed binding-only pair
    # would be a statement the rule cannot place.
    with pytest.raises(ValueError, match='without curated evidence'):
        mod.interactions_schema(dict(ROWS['ADNP', 'ABCC1'], mor=1.0))


def test_row_without_any_line_of_evidence_is_an_error():
    row = dict(
        ROWS['ADNP', 'AASDH'],
        is_evidence_inferred=False,
    )
    with pytest.raises(ValueError, match='without a line of evidence'):
        mod.interactions_schema(row)


@pytest.mark.parametrize(
    'pair,level,evidence,sources',
    [
        (('AHR', 'ASAP1'), 'B', ['ECO:0000229', 'ECO:0005536'], ['remap', 'hocomoco', 'gtex']),
        (('AHR', 'CTNNA2'), 'C', ['ECO:0000229', 'ECO:0005536'], ['remap', 'hocomoco']),
        (('ADNP', 'ABCC1'), 'D', ['ECO:0000229'], ['remap']),
        (('ADNP', 'AASDH'), 'E', [], ['gtex']),
        (('AHR', 'ASAP2'), 'E', ['ECO:0005536'], ['hocomoco', 'gtex']),
    ],
)
def test_level_evidence_kinds_and_sources_are_kept(pair, level, evidence, sources):
    _, out = extract(ROWS[pair])
    relation = out.relations[0]
    assert annotations(relation, 'has_confidence_level') == [level]
    assert annotations(relation, 'has_evidence_of_type') == evidence
    assert annotations(relation, 'supporting_data_source') == sources
    assert not annotations(relation, 'publications')


def test_only_pubmed_ids_become_publications():
    _, out = extract(ROWS['AR', 'ACP3'])
    assert annotations(out.relations[0], 'publications') == ['PMID:23518348']
    _, out = extract(ROWS['CTCF', 'EPX'])
    assert not annotations(out.relations[0], 'publications')
    row = dict(ROWS['AHR', 'CA9'], pubmed_id='S7501,mint:MINT-8429205, ,10551783')
    _, out = extract(row)
    assert annotations(out.relations[0], 'publications') == ['PMID:10551783']


def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication():
    # The file has ``0`` and fused numbers of 16 digits.
    row = dict(
        ROWS['AHR', 'CA9'],
        pubmed_id='0,0123456,123456789,1451920411085504,pubmed:10551783',
    )
    _, out = extract(row)
    assert annotations(out.relations[0], 'publications') == ['PMID:10551783']


@pytest.mark.parametrize('missing', ['tf', 'target'])
def test_row_without_a_participant_emits_no_relation(missing):
    assert mod.interactions_schema(dict(ROWS['AHR', 'CA9'], **{missing: None})) is None


@pytest.mark.parametrize(
    'field,value,message',
    [
        ('mor', 0.5, 'Unknown DoRothEA mode of regulation'),
        ('mor', None, 'Unknown DoRothEA mode of regulation'),
        ('mor', True, 'Unknown DoRothEA mode of regulation'),
        ('confidence', 'F', 'Unknown DoRothEA confidence level'),
        ('confidence', None, 'Unknown DoRothEA confidence level'),
        ('is_evidence_tfbs', 'TRUE', 'is not a boolean'),
        ('is_evidence_chip_seq', None, 'is not a boolean'),
    ],
)
def test_unknown_categorical_values_are_rejected(field, value, message):
    with pytest.raises(ValueError, match=message):
        mod.interactions_schema(dict(ROWS['AHR', 'CA9'], **{field: value}))


def test_level_datasets_partition_the_table(monkeypatch):
    seen = {}

    def frame(opener, name=None, **_kwargs):
        seen['name'] = name
        yield from ROWS.values()

    monkeypatch.setattr(parser, 'iter_rdata', frame)
    by_level = {
        level: list(parser.iter_level(None, level=level))
        for level in parser.CONFIDENCE_LEVELS
    }
    assert seen['name'] == 'entire_database'
    assert {level: len(rows) for level, rows in by_level.items()} == {
        'A': 3,
        'B': 1,
        'C': 1,
        'D': 2,
        'E': 2,
    }
    assert all(
        row['confidence'] == level
        for level, rows in by_level.items()
        for row in rows
    )
    assert sum(map(len, by_level.values())) == len(ROWS)


def test_row_with_an_unknown_level_is_not_silently_left_out(monkeypatch):
    rows = [ROWS['AHR', 'CA9'], dict(ROWS['AHR', 'CA9'], confidence='F')]
    monkeypatch.setattr(parser, 'iter_rdata', lambda opener, **_kwargs: iter(rows))
    with pytest.raises(ValueError, match='Unknown DoRothEA confidence level'):
        list(parser.iter_level(None, level='A'))
