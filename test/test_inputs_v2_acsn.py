"""Source contract of the ACSN binary relations input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.acsn import INTERACTION_TYPES, STATEMENTS, resource
from pypath.inputs_v2.parsers.acsn import COLUMNS, hugo_names, join_names
from pypath.inputs_v2.parsers.base import iter_gmt

# Lines of ACSN2_HUGO_Correspondence.gmt.
GMT = [
    'VEGFA\tna\tVEGFA\t',
    'FLT1\tna\tFLT1\t',
    'SH2D2A\tna\tSH2D2A\t',
    'LCK\tna\tLCK\t',
    'PTP4A3\tna\tPTP4A3\t',
    'ARF1\tna\tARF1\t',
    'RAC1\tna\tRAC1\t',
    'NCKAP1\tna\tNCKAP1\t',
    'TNF\tna\tTNF\t',
    'NF-_kappa_B_p65*\tna\tRELA\t',
    'IGF1R\tna\tIGF1R\t',
    'IGF1/2*\tna\tIGF1\tIGF2\t',
    'CDKN2A\tna\tCDKN2A\t',
    'p16INK4A*\tna\tCDKN2A\t',
    'MAPK14\tna\tMAPK14\t',
    'EGFR\tna\tEGFR\tEGFRHUGO:EGFR\t',
    'HB-EGF*\tna\tHBEGF\t',
    'HLA-C\tna\tHLA-C\t',
    'KIR2DL2\tna\t:KIR2DL2\t',
    'SNAI2\tna\tSNAI2\t',
    'MIR200*\tna\tMIR200A\tMIR200B\tMIR200C\tMIR141\tMIR429\t',
    'MIR149\tna\tMIR149\t',
    'PTGER2\tna\tPTGER2\t',
    'TGFBR2\tna\tTGFBR2\t',
    'IQGAP1 \tna\tIQGAP1\t',
    'IQGAP1\tna\tIQGAP1\t',
    'CD36\tna\tCD36\t',
    'TSP*\tna\tTHBS1\t',
    'CYFIP2\tna\tCYFIP2\t',
    'BNIP2\tna\tBNIP2\t',
]


def table(lines):
    return SimpleNamespace(result=io.StringIO('\n'.join(lines) + '\n'))


NAMES = hugo_names(iter_gmt(table(GMT)))


def statements(*values):
    """The parsed rows of one line of the relations table."""
    return list(
        join_names(
            [dict(zip(COLUMNS, values, strict=True))],
            NAMES,
            INTERACTION_TYPES,
            STATEMENTS,
        )
    )


# Rows of ACSN2_binary_relations_between_proteins_with_PMID.txt. Long PMID
# cells are shortened.
ASSOCIATION_ROW = (
    'VEGFA', 'HETERODIMER_ASSOCIATION;PROTEIN_INTERACTION', 'FLT1',
    '10022831;11166270;8837607',
)
INHIBITION_ROW = ('SH2D2A', 'INHIBITION', 'LCK', '10975797;27896837')
MIXED_ROW = (
    'PTP4A3', 'CATALYSIS;COMPLEX_EXPANSION;HETERODIMER_ASSOCIATION', 'ARF1',
    '22595524;23691193;24376839',
)
NEGATIVE_INFLUENCE_ROW = ('RAC1', 'NEGATIVE_INFLUENCE', 'NCKAP1', '12181570')
POSITIVE_INFLUENCE_ROW = (
    'TNF', 'POSITIVE_INFLUENCE', 'NF-_kappa_B_p65*', '19411070',
)
UNKNOWN_INFLUENCE_ROW = ('CYFIP2', 'UNKNOWN_NEGATIVE_INFLUENCE', 'BNIP2', 'N/A')
FAMILY_ROW = (
    'IGF1R', 'HETERODIMER_ASSOCIATION;PROTEIN_INTERACTION', 'IGF1/2*',
    '10086337;10826997',
)
TRANSCRIPTION_ROW = ('CDKN2A', 'TRANSCRIPTION', 'p16INK4A*', 'N/A')
ALIAS_ROW = ('MAPK14', 'CATALYSIS', 'p16INK4A*', '17254968;23850759')
EXTRA_NAME_ROW = (
    'EGFR', 'COMPLEX_EXPANSION;HETERODIMER_ASSOCIATION', 'HB-EGF*',
    '16525676;18372913',
)
NO_HUGO_ROW = (
    'HLA-C', 'HETERODIMER_ASSOCIATION;PROTEIN_INTERACTION', 'KIR2DL2',
    '11513144;21339333;9605119',
)
MIRNA_FAMILY_ROW = (
    'SNAI2', 'HETERODIMER_ASSOCIATION;PROTEIN_INTERACTION', 'MIR200*',
    '22370643',
)
MIRNA_ROW = (
    'MIR149', 'HETERODIMER_ASSOCIATION;INHIBITION;PROTEIN_INTERACTION',
    'PTGER2', '21209948;22025564;25916550',
)
TRAILING_SPACE_ROW = (
    'TGFBR2', 'HETERODIMER_ASSOCIATION;PROTEIN_INTERACTION', 'IQGAP1 ',
    '14557817;16474430',
)
LONG_NUMBER_ROW = (
    'CD36', 'COMPLEX_EXPANSION;HETERODIMER_ASSOCIATION', 'TSP*',
    '19909365;20150200015;21149606',
)


def extract(parsed):
    record = resource.datasets()['interactions'].mapper(parsed)
    assert record is not None
    out = SilverExtractor('acsn', 'interactions')
    out.process_record(record, parsed, 'fixture:1', 1)
    (relation,) = [r for r in out.relations if r.predicate != 'has_member']
    return (
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
        out,
    )


def identifiers(entity):
    return {(i['ns'], i['id']) for i in entity.identifiers}


def annotations(relation, scope='relation'):
    return sorted(
        (a['term'], a['value'])
        for a in relation.annotations
        if a['scope'] == scope
    )


def members(out, entity_key):
    return {
        identifier['id']
        for relation in out.relations
        if relation.predicate == 'has_member'
        and relation.subject_entity_key == entity_key
        for identifier in out.entities[relation.object_entity_key].identifiers
    }


def test_association_is_an_undirected_interaction_with_publications():
    (parsed,) = statements(*ASSOCIATION_ROW)
    assert parsed['STATEMENT_TYPES'] == [
        'HETERODIMER_ASSOCIATION',
        'PROTEIN_INTERACTION',
    ]
    relation, first, second, _ = extract(parsed)
    assert relation.predicate == 'interacts_with'
    assert identifiers(first) == {(Namespace.GENESYMBOL, 'VEGFA')}
    assert identifiers(second) == {(Namespace.GENESYMBOL, 'FLT1')}
    assert {first.entity_type, second.entity_type} == {'protein'}
    assert {first.taxon, second.taxon} == {'9606'}
    assert annotations(relation) == [
        ('publications', 'PMID:10022831'),
        ('publications', 'PMID:11166270'),
        ('publications', 'PMID:8837607'),
    ]


def test_inhibition_of_a_reaction_gives_no_direction_on_the_second_entity():
    (parsed,) = statements(*INHIBITION_ROW)
    relation, first, second, _ = extract(parsed)
    assert relation.predicate == 'affects'
    assert identifiers(first) == {(Namespace.GENESYMBOL, 'SH2D2A')}
    assert identifiers(second) == {(Namespace.GENESYMBOL, 'LCK')}
    # The map says SH2D2A inhibits a reaction of LCK, not that LCK decreases.
    assert {a['term'] for a in relation.annotations} == {
        'in_taxon',
        'publications',
    }


def test_row_with_two_kinds_of_type_gives_two_statements():
    association, modulation = statements(*MIXED_ROW)
    assert association['STATEMENT_TYPES'] == [
        'COMPLEX_EXPANSION',
        'HETERODIMER_ASSOCIATION',
    ]
    assert modulation['STATEMENT_TYPES'] == ['CATALYSIS']
    assert extract(association)[0].predicate == 'interacts_with'
    relation = extract(modulation)[0]
    assert relation.predicate == 'affects'
    # Both statements carry the publications of the row.
    assert annotations(relation) == annotations(extract(association)[0])
    assert len(annotations(relation)) == 3


@pytest.mark.parametrize(
    'row,direction',
    [
        (NEGATIVE_INFLUENCE_ROW, 'decreased'),
        (POSITIVE_INFLUENCE_ROW, 'increased'),
        (UNKNOWN_INFLUENCE_ROW, 'decreased'),
    ],
)
def test_influence_arrow_gives_the_direction(row, direction):
    (parsed,) = statements(*row)
    relation, _, _, _ = extract(parsed)
    assert relation.predicate == 'affects'
    assert ('object_direction_qualifier', direction) in annotations(relation)
    assert not any(
        term == 'object_aspect_qualifier' for term, _ in annotations(relation)
    )


def test_pair_without_publication_has_none():
    (parsed,) = statements(*UNKNOWN_INFLUENCE_ROW)
    relation, _, _, _ = extract(parsed)
    assert not any(term == 'publications' for term, _ in annotations(relation))


def test_alias_node_is_the_gene_of_the_correspondence_file():
    (parsed,) = statements(*POSITIVE_INFLUENCE_ROW)
    _, _, second, _ = extract(parsed)
    assert second.entity_type == 'protein'
    assert identifiers(second) == {
        (Namespace.GENESYMBOL, 'RELA'),
        (Namespace.NAME, 'NF-_kappa_B_p65*'),
    }
    (parsed,) = statements(*ALIAS_ROW)
    _, _, second, _ = extract(parsed)
    assert identifiers(second) == {
        (Namespace.GENESYMBOL, 'CDKN2A'),
        (Namespace.NAME, 'p16INK4A*'),
    }


def test_generic_node_is_a_family_with_its_listed_members():
    (parsed,) = statements(*FAMILY_ROW)
    relation, first, second, out = extract(parsed)
    assert second.entity_type == 'protein_family'
    assert identifiers(second) == {(Namespace.NAME, 'IGF1/2*')}
    assert members(out, relation.object_entity_key) == {'IGF1', 'IGF2'}
    # One relation with the family; nothing is expanded to its members.
    assert [r.predicate for r in out.relations].count('interacts_with') == 1
    assert identifiers(first) == {(Namespace.GENESYMBOL, 'IGF1R')}


def test_reactant_to_product_pair_gives_no_relation():
    (parsed,) = statements(*TRANSCRIPTION_ROW)
    assert parsed['STATEMENT'] is None
    assert parsed['STATEMENT_TYPES'] == ['TRANSCRIPTION']
    assert resource.datasets()['interactions'].mapper(parsed) is None


def test_malformed_hugo_name_is_not_an_identifier():
    (parsed,) = statements(*EXTRA_NAME_ROW)
    _, first, _, _ = extract(parsed)
    assert identifiers(first) == {(Namespace.GENESYMBOL, 'EGFR')}
    # The only HUGO name of KIR2DL2 in the file is ``:KIR2DL2``.
    (parsed,) = statements(*NO_HUGO_ROW)
    _, _, second, _ = extract(parsed)
    assert identifiers(second) == {(Namespace.NAME, 'KIR2DL2')}


def test_microrna_nodes_are_not_proteins():
    (parsed,) = statements(*MIRNA_FAMILY_ROW)
    relation, _, second, out = extract(parsed)
    assert second.entity_type == 'gene_family'
    family_members = {
        out.entities[r.object_entity_key].entity_type
        for r in out.relations
        if r.predicate == 'has_member'
    }
    assert family_members == {'microrna'}
    assert members(out, relation.object_entity_key) == {
        'MIR200A', 'MIR200B', 'MIR200C', 'MIR141', 'MIR429',
    }
    association, modulation = statements(*MIRNA_ROW)
    _, first, second, _ = extract(modulation)
    assert first.entity_type == 'microrna'
    assert second.entity_type == 'protein'


def test_space_after_a_name_is_dropped():
    (parsed,) = statements(*TRAILING_SPACE_ROW)
    assert parsed['INTERACTOR2'] == 'IQGAP1'
    _, _, second, _ = extract(parsed)
    assert identifiers(second) == {(Namespace.GENESYMBOL, 'IQGAP1')}


def test_number_with_a_leading_zero_is_not_a_publication():
    # The table has ``0643644`` in six rows; a PubMed ID has no leading zero.
    row = (*LONG_NUMBER_ROW[:3], '0643644;0123456;123456789;21149606')
    (parsed,) = statements(*row)
    relation, _, _, _ = extract(parsed)
    assert annotations(relation) == [('publications', 'PMID:21149606')]


def test_number_longer_than_a_pmid_is_not_a_publication():
    (parsed,) = statements(*LONG_NUMBER_ROW)
    relation, _, second, _ = extract(parsed)
    assert annotations(relation) == [
        ('publications', 'PMID:19909365'),
        ('publications', 'PMID:21149606'),
    ]
    assert (Namespace.GENESYMBOL, 'THBS1') in identifiers(second)


def test_node_missing_from_the_correspondence_file_keeps_only_its_name():
    (parsed,) = statements('VEGFA', 'CATALYSIS', 'XYZ*', 'N/A')
    assert parsed['HUGO2'] == []
    _, _, second, _ = extract(parsed)
    assert identifiers(second) == {(Namespace.NAME, 'XYZ*')}


def test_unknown_interaction_type_is_an_error():
    with pytest.raises(ValueError, match='ACTIVATION'):
        statements('VEGFA', 'CATALYSIS;ACTIVATION', 'FLT1', 'N/A')


def test_every_type_of_the_crosswalk_gives_a_known_statement_or_none():
    assert set(INTERACTION_TYPES.values()) == {*STATEMENTS, None}
    # The sign of a modification of a reaction is not a direction.
    assert INTERACTION_TYPES['INHIBITION'] == INTERACTION_TYPES['CATALYSIS']
    assert STATEMENTS[INTERACTION_TYPES['INHIBITION']][1] is None


def test_table_with_changed_columns_is_rejected():
    with pytest.raises(ValueError, match='INTERACTION_TYPE'):
        list(
            join_names(
                [{'INTERACTOR1': 'VEGFA', 'INTERACTOR2': 'FLT1'}],
                NAMES,
                INTERACTION_TYPES,
                STATEMENTS,
            )
        )


def test_empty_tables_are_rejected():
    with pytest.raises(ValueError, match='no rows'):
        list(join_names([], NAMES, INTERACTION_TYPES, STATEMENTS))
    with pytest.raises(ValueError, match='no entities'):
        hugo_names([])


def test_gmt_reader_reads_names_and_members():
    assert list(iter_gmt(table(['IGF1/2*\tna\tIGF1\tIGF2\t', '', 'RNF11\tna\tRNF11\t']))) == [
        {'name': 'IGF1/2*', 'description': 'na', 'members': ['IGF1', 'IGF2']},
        {'name': 'RNF11', 'description': 'na', 'members': ['RNF11']},
    ]
    with pytest.raises(ValueError, match='no description'):
        list(iter_gmt(table(['IGF1'])))


def test_names_written_with_and_without_a_space_are_one_entity():
    assert NAMES['IQGAP1'] == ['IQGAP1']
    assert 'IQGAP1 ' not in NAMES
