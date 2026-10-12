"""Source contract of the Baccin 2019 ligand-receptor input."""

import io
from types import SimpleNamespace

import openpyxl
import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace
from omnipath_core.source_attributes import PARTICIPANT_ROLE

from pypath.inputs_v2 import baccin2019 as mod
from pypath.inputs_v2.parsers import baccin2019 as parser


def row(*values, ligand=None, receptor=None):
    record = dict(zip(parser.COLUMNS, values, strict=True))
    record['ligand'] = ligand or record['Ligand.Mouse']
    record['receptor'] = receptor or record['Receptor.Mouse']
    return record


# Rows of SuppTable3, as the reader yields them.
RAMILOWSKI = row(
    'A2m-Lrp1', 'A2m', 'Lrp1', 'Ramilowski', 'Correct', 'Secreted', 'Other',
    '1702392,10652313, 12194978',
)
AUTHORS = row(
    'Adgre5-Cd55', 'Adgre5', 'Cd55', 'Baccin', 'Correct', 'Membrane', 'Other',
    23447688,
)
COMPLEX_RECEPTOR = row(
    'Adm-Calcrl&Ramp2', 'Adm', 'Calcrl&Ramp2', 'Baccin', 'Correct', 'Secreted',
    'Other', 10342881,
)
COMPLEX_LIGAND = row(
    'Crlf1&Clcf1-Cntfr', 'Crlf1&Clcf1', 'Cntfr', 'Baccin', 'Correct',
    'Secreted', 'Other', 16782820,
)
KEGG = row(
    'Cd226-Itgal&Itgb2', 'Cd226', 'Itgal&Itgb2', 'Baccin', 'Correct',
    'Membrane', 'Other', 'mmu04514',
)
IRRELEVANT = row(
    'C3-Cr1l', 'C3', 'Cr1l', 'Ramilowski', 'Irrelevant', 'Secreted', 'Other',
    '9988761, 10382758',
)
# The spreadsheet fused two PMIDs into one number.
FUSED_PMIDS = row(
    'Tgfb1-Acvrl1', 'Tgfb1', 'Acvrl1', 'Ramilowski', 'Correct', 'Secreted',
    'Cytokine', 1071699310716990,
)
# 7 + 7 digits, the same PMID twice.
FUSED_SAME_PMID = row(
    'Pros1-Axl', 'Pros1', 'Axl', 'Ramilowski', 'Correct', 'Secreted', 'Other',
    78670737867073,
)
# 7 + 7 digits.
FUSED_SEVEN_DIGITS = row(
    'Nrtn-Gfra1', 'Nrtn', 'Gfra1', 'Ramilowski', 'Correct', 'Secreted',
    'GrowthFactor', 91828039407096,
)
# 15 digits read as 7 + 8 and as 8 + 7.
FUSED_AMBIGUOUS = row(
    'Jag1-Notch1', 'Jag1', 'Notch1', 'Ramilowski', 'Correct', 'Membrane',
    'Other', 769772111006133,
)
# 14 digits whose 7 + 7 reading ends in a part starting with a zero.
FUSED_LEADING_ZERO = row(
    'Vegfb-Nrp1', 'Vegfb', 'Nrp1', 'Ramilowski', 'Correct', 'Secreted',
    'GrowthFactor', 40967710409677,
)
STRAY_BRACKET = row(
    'Sema3a-Nrp1', 'Sema3a', 'Nrp1', 'Ramilowski', 'Correct', 'Secreted',
    'Other', '10520994, 10934324), 15721238,10520994, 9331348',
)
NO_REFERENCE = row(
    'Hspg2-Itga2&Itgb1', 'Hspg2', 'Itga2&Itgb1', 'Baccin', 'Correct',
    'Secreted', 'Other', None,
)
INCORRECT = row(
    'App-Tnfrsf21', 'App', 'Tnfrsf21', 'Ramilowski', 'Incorrect', 'Membrane',
    'Other', 'NA',
)
NOT_CHECKED = row(
    'Bdnf-Ngfr', 'Bdnf', 'Ngfr', 'Ramilowski', 'Not Expressed', 'Secreted',
    'GrowthFactor', None,
)
SCAVENGER = row(
    'Ccl2-Ackr1', 'Ccl2', 'Ackr1', 'Baccin', 'Scavenger', 'Secreted',
    'Chemokine', 24890717,
)
LIKELY = row(
    'Tgfb2-Acvr1', 'Tgfb2', 'Acvr1', 'Baccin', 'Likely', 'Secreted',
    'Cytokine', 21102460,
)
PROCESSING = row(
    'App&Bace1-Tnfrsf21', 'App&Bace1', 'Tnfrsf21', 'Baccin', 'Correct',
    'Secreted', 'Other', 19225519,
)
STEM = row(
    'Lam-Itga10&Itgb1', 'Lam', 'Itga10&Itgb1', 'Baccin', 'Correct', 'ECM',
    'Other', 16988024,
)
ALTERNATIVES = (
    'Thbs|Thbs1|Thbs2-Itga3&Itgb1', 'Thbs|Thbs1|Thbs2', 'Itga3&Itgb1', 'Baccin',
    'Correct', 'ECM', 'Other', 16988024,
)


def extract(source_row):
    record = mod.resource.datasets()['interactions'].mapper(source_row)
    assert record is not None
    out = SilverExtractor('baccin2019', 'interactions')
    out.process_record(record, source_row, 'fixture:1', 1)
    (relation,) = [r for r in out.relations if r.predicate == 'interacts_with']
    return relation, out


def identifiers(entity):
    return [(i['ns'], i['id']) for i in entity.identifiers]


def annotations(relation, scope='relation'):
    return sorted(
        (a['term'], a['value'])
        for a in relation.annotations
        if a['scope'] == scope
    )


def members(out, complex_key):
    return sorted(
        out.entities[r.object_entity_key].identifier
        for r in out.relations
        if r.predicate == 'has_member' and r.subject_entity_key == complex_key
    )


def test_pair_is_mouse_ligand_to_receptor_with_database_and_publications():
    relation, out = extract(RAMILOWSKI)
    ligand = out.entities[relation.subject_entity_key]
    receptor = out.entities[relation.object_entity_key]
    assert identifiers(ligand) == [(Namespace.GENESYMBOL, 'A2m')]
    assert identifiers(receptor) == [(Namespace.GENESYMBOL, 'Lrp1')]
    assert {ligand.entity_type, receptor.entity_type} == {'protein'}
    # The source's own species and symbols: nothing is translated to human.
    assert {ligand.taxon, receptor.taxon} == {'10090'}
    assert relation.upstream_id == 'A2m-Lrp1'
    assert annotations(relation) == [
        ('publications', 'PMID:10652313'),
        ('publications', 'PMID:12194978'),
        ('publications', 'PMID:1702392'),
        ('supporting_data_source', 'Ramilowski2015'),
    ]


def test_authors_own_pair_names_no_supporting_database():
    relation, _ = extract(AUTHORS)
    assert annotations(relation) == [('publications', 'PMID:23447688')]


def test_pair_asserts_no_direction_role_or_ligand_class():
    relation, _ = extract(RAMILOWSKI)
    assert {a['term'] for a in relation.annotations} == {
        'in_taxon',
        'publications',
        'supporting_data_source',
        PARTICIPANT_ROLE,
    }
    assert {
        a['term'] for a in relation.annotations if a['scope'] == 'relation'
    } == {'publications', 'supporting_data_source'}
    # The ligand's compartment and class have no faithful slot.
    assert all(
        a['value'] not in {'Secreted', 'Other'} for a in relation.annotations
    )
    # The roles are annotations of the participants, not of the relation.
    assert annotations(relation, 'subject') == [
        (PARTICIPANT_ROLE, 'ligand'),
        ('in_taxon', 'NCBITaxon:10090'),
    ]
    assert annotations(relation, 'object') == [
        (PARTICIPANT_ROLE, 'receptor'),
        ('in_taxon', 'NCBITaxon:10090'),
    ]


def test_subunits_joined_by_ampersand_are_a_complex_receptor():
    relation, out = extract(COMPLEX_RECEPTOR)
    receptor = out.entities[relation.object_entity_key]
    assert receptor.entity_type == 'macromolecular_complex'
    assert identifiers(receptor) == [(Namespace.NAME, 'Calcrl&Ramp2')]
    assert receptor.taxon == '10090'
    assert members(out, receptor.entity_key) == ['Calcrl', 'Ramp2']
    assert out.entities[relation.subject_entity_key].identifier == 'Adm'
    # Membership only: no pairwise interactions among or with the subunits.
    assert sorted(r.predicate for r in out.relations) == [
        'has_member',
        'has_member',
        'interacts_with',
    ]


def roles(relation):
    return sorted(
        (a['scope'], a['value'])
        for a in relation.annotations
        if a['term'] == PARTICIPANT_ROLE
    )


def test_ligand_and_receptor_columns_give_the_participant_roles():
    relation, out = extract(RAMILOWSKI)
    assert roles(relation) == [('object', 'receptor'), ('subject', 'ligand')]
    assert out.entities[relation.subject_entity_key].identifier == 'A2m'
    assert out.entities[relation.object_entity_key].identifier == 'Lrp1'
    assert len(out.entities) == 2


@pytest.mark.parametrize(
    'source_row,subject_type,object_type',
    [
        (COMPLEX_RECEPTOR, 'protein', 'macromolecular_complex'),
        (COMPLEX_LIGAND, 'macromolecular_complex', 'protein'),
    ],
)
def test_complex_on_either_side_has_the_role_and_its_subunits_have_none(
    source_row, subject_type, object_type
):
    relation, out = extract(source_row)
    assert out.entities[relation.subject_entity_key].entity_type == subject_type
    assert out.entities[relation.object_entity_key].entity_type == object_type
    # The role is on the complex that the table lists as ligand or receptor,
    # not on its subunits (which are members of it).
    assert roles(relation) == [('object', 'receptor'), ('subject', 'ligand')]
    memberships = [r for r in out.relations if r.predicate == 'has_member']
    assert memberships
    assert all(
        a['term'] != PARTICIPANT_ROLE for m in memberships for a in m.annotations
    )
    assert all(
        a['term'] != PARTICIPANT_ROLE
        for e in out.entities.values()
        for a in e.annotations
    )


def test_complex_ligand_keeps_its_members():
    relation, out = extract(COMPLEX_LIGAND)
    ligand = out.entities[relation.subject_entity_key]
    assert ligand.entity_type == 'macromolecular_complex'
    assert members(out, ligand.entity_key) == ['Clcf1', 'Crlf1']
    assert all(
        out.entities[r.object_entity_key].taxon == '10090'
        for r in out.relations
        if r.predicate == 'has_member'
    )


def test_kegg_pathway_reference_is_a_database_not_a_publication():
    relation, _ = extract(KEGG)
    assert annotations(relation) == [('supporting_data_source', 'KEGG')]


def test_pair_judged_correct_but_irrelevant_at_homeostasis_is_emitted():
    relation, _ = extract(IRRELEVANT)
    assert annotations(relation) == [
        ('publications', 'PMID:10382758'),
        ('publications', 'PMID:9988761'),
        ('supporting_data_source', 'Ramilowski2015'),
    ]


def publications(source_row):
    relation, _ = extract(source_row)
    return [v for t, v in annotations(relation) if t == 'publications']


def test_two_pmids_run_together_are_split():
    # 8 + 8 digits (an integer cell).
    assert publications(FUSED_PMIDS) == ['PMID:10716990', 'PMID:10716993']
    # 7 + 7 digits.
    assert publications(FUSED_SEVEN_DIGITS) == ['PMID:9182803', 'PMID:9407096']
    # The same PMID twice is one publication.
    assert publications(FUSED_SAME_PMID) == ['PMID:7867073']
    # Among others in the cell, in source order.
    assert parser.pubmed_ids('1702392, 1071699310716990, 10652313') == [
        '1702392',
        '10716993',
        '10716990',
        '10652313',
    ]


@pytest.mark.parametrize(
    'token,pmids',
    [
        ('1071699310716990', ['10716993', '10716990']),
        ('91828039407096', ['9182803', '9407096']),
        # 15 digits have two readings (7 + 8 and 8 + 7).
        ('769772111006133', None),
        # 7 + 7 would end in a part starting with a zero.
        ('40967710409677', None),
        # Cut short, a single number, or fused beyond two PMIDs.
        ('58985', None),
        ('10652313', None),
        ('5898511035029', None),
        ('255887787027019005952', None),
        ('1297448112974480036308474939213342572544', None),
        ('GCID: GC01P159203', None),
    ],
)
def test_split_fused_pmids_only_when_exactly_two_pmids(token, pmids):
    assert parser.split_fused_pmids(token) == pmids


def test_ambiguous_truncated_and_oversized_numbers_are_not_publications():
    assert publications(FUSED_AMBIGUOUS) == []
    assert publications(FUSED_LEADING_ZERO) == []
    assert publications(dict(AUTHORS, **{parser.REFERENCE: 939660})) == []
    # A list fused into a number too large for an integer cell.
    assert (
        publications(dict(AUTHORS, **{parser.REFERENCE: 7.86707378670738e20}))
        == []
    )


def test_stray_bracket_and_repeated_pmids_are_cleaned():
    relation, _ = extract(STRAY_BRACKET)
    assert [v for t, v in annotations(relation) if t == 'publications'] == [
        'PMID:10520994',
        'PMID:10934324',
        'PMID:15721238',
        'PMID:9331348',
    ]


def test_pair_without_reference_has_no_publication():
    relation, _ = extract(NO_REFERENCE)
    assert annotations(relation) == []
    relation, _ = extract(dict(IRRELEVANT, **{parser.REFERENCE: 'NA'}))
    assert annotations(relation) == [
        ('supporting_data_source', 'Ramilowski2015')
    ]


@pytest.mark.parametrize(
    'source_row', [INCORRECT, NOT_CHECKED, SCAVENGER, LIKELY]
)
def test_only_pairs_the_authors_consider_correct_are_emitted(source_row):
    assert mod.map_baccin_interaction(source_row) is None


def test_precursor_with_its_proteases_is_not_a_complex():
    assert mod.map_baccin_interaction(PROCESSING) is None


def test_family_stem_identifies_no_protein():
    assert mod.map_baccin_interaction(STEM) is None
    assert mod.map_baccin_interaction(row(*ALTERNATIVES, ligand='Thbs')) is None


def test_each_alternative_is_its_own_pair():
    relation, out = extract(row(*ALTERNATIVES, ligand='Thbs1'))
    assert out.entities[relation.subject_entity_key].identifier == 'Thbs1'
    receptor = out.entities[relation.object_entity_key]
    assert members(out, receptor.entity_key) == ['Itga3', 'Itgb1']
    # The alternatives are not grouped into a complex or a family.
    assert 'Thbs2' not in {e.identifier for e in out.entities.values()}


@pytest.mark.parametrize('role', ['ligand', 'receptor'])
def test_no_relation_without_a_participant(role):
    assert mod.map_baccin_interaction(dict(RAMILOWSKI, **{role: ''})) is None


def test_unknown_verdict_or_source_is_rejected():
    with pytest.raises(ValueError, match='manual annotation'):
        mod.map_baccin_interaction(dict(RAMILOWSKI, ManualAnnotation='Maybe'))
    with pytest.raises(ValueError, match='source'):
        mod.map_baccin_interaction(dict(RAMILOWSKI, Source='CellPhoneDB'))


def workbook(*rows, header=parser.COLUMNS):
    book = openpyxl.Workbook()
    sheet = book.active
    sheet.title = parser.SHEET
    sheet.append(['SuppTable3. RNA-Magnet receptor-ligand db'])
    sheet.append([])
    sheet.append(list(header))
    for values in rows:
        sheet.append(list(values))
    content = io.BytesIO()
    book.save(content)
    return SimpleNamespace(result=io.BytesIO(content.getvalue()))


def test_reader_expands_alternatives_and_keeps_the_source_columns():
    records = list(
        parser.iter_pairs(
            workbook(
                ('Adgre5-Cd55', 'Adgre5', 'Cd55', 'Baccin', 'Correct',
                 'Membrane', 'Other', 23447688),
                ALTERNATIVES,
                ('Itgax&Itgb2-Icam1|Icam4', 'Itgax&Itgb2', 'Icam1|Icam4',
                 'Baccin', 'Correct', 'Membrane', 'Other', 16988024),
            )
        )
    )
    assert records[0] == AUTHORS
    assert [(r['ligand'], r['receptor']) for r in records[1:]] == [
        ('Thbs', 'Itga3&Itgb1'),
        ('Thbs1', 'Itga3&Itgb1'),
        ('Thbs2', 'Itga3&Itgb1'),
        ('Itgax&Itgb2', 'Icam1'),
        ('Itgax&Itgb2', 'Icam4'),
    ]
    assert {r['Ligand.Mouse'] for r in records[1:4]} == {'Thbs|Thbs1|Thbs2'}
    assert records[4]['Pair.Name'] == 'Itgax&Itgb2-Icam1|Icam4'


def test_reader_rejects_a_table_with_changed_columns():
    with pytest.raises(ValueError, match='ManualAnnotation'):
        list(
            parser.iter_pairs(
                workbook(
                    ('A2m-Lrp1', 'A2m', 'Lrp1'),
                    header=('Pair.Name', 'Ligand.Mouse', 'Receptor.Mouse'),
                )
            )
        )


def test_reader_rejects_an_empty_table():
    with pytest.raises(ValueError, match='no rows'):
        list(parser.iter_pairs(workbook()))
