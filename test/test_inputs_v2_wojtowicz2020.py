"""Source contract of the Wojtowicz 2020 cell-surface interactome input."""

import io
from types import SimpleNamespace

import openpyxl
import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import wojtowicz2020 as mod


COLUMNS = [
    'Prey_UniProt_Name', 'Bait_UniProt_Name', 'Prey_Common_Name',
    'Bait_Common_Name', 'Prey_Gene_Name', 'Bait_Gene_Name',
    'Prey-Bait_Replicate1_OD650', 'Prey-Bait_Replicate2_OD650',
    'Prey-Bait_Replicate3_OD650', 'Prey-Control_Replicate1_OD650',
    'Prey-Control_Replicate2_OD650', 'Prey-Control_Replicate3_OD650',
    'Prey-Bait_Replicate1_FoldOverBackground',
    'Prey-Bait_Replicate2_FoldOverBackground',
    'Prey-Bait_Replicate3_FoldOverBackground',
    'Prey-Bait_FoldOverBackground_Average',
    'Prey-Bait_FoldOverBackground_StdDev', 'Binding_BothOrientations',
    'PPI Status',
]  # fmt: skip
# Rows of Data S4.
NEW = dict(zip(COLUMNS, [
    '1B35', 'SIG10', 'HLA-B', 'Siglec-10', 'HLA-B', 'SIGLEC10',
    0.1343, 0.1416, 0.1386, 0.0655, 0.067, 0.0699,
    2.050381679389313, 2.1134328358208956, 1.982832618025751,
    2.04888237774532, 0.06531301670493841,
    'no', 'not found in literature or PPI databases',
]))  # fmt: skip
UNREVIEWED_BAIT = dict(zip(COLUMNS, [
    'I18BP', 'A0A0B4J1R5', 'IL18BP', 'KIR3DL3', 'IL18BP', 'KIR3DL3',
    0.2641, 0.2536, 0.245, 0.075, 0.06899999999999999, 0.0658,
    3.5213333333333336, 3.67536231884058, 3.723404255319149,
    3.6400333024976876, 0.10556641949439521,
    'no', 'not found in literature or PPI databases',
]))  # fmt: skip
NO_GENE_SYMBOL = dict(zip(COLUMNS, [
    'A2IXV5', 'CD33', 'A2IXV5', 'CD33', 'MIR cl-10', 'CD33',
    0.2875, 0.28559999999999997, 0.2891, 0.0898, 0.0788, 0.0861,
    3.201559020044543, 3.6243654822335025, 3.357723577235773,
    3.39454935983794, 0.21379529886174425,
    'no', 'not found in literature or PPI databases',
]))  # fmt: skip
HOMOPHILIC = dict(zip(COLUMNS, [
    'CD33', 'CD33', 'CD33', 'CD33', 'CD33', 'CD33',
    0.8552, 0.8389, 0.7819, 0.0868, 0.0831, 0.08039999999999999,
    9.85253456221198, 10.095066185318894, 9.725124378109456,
    9.890908375213444, 0.1879325605730674,
    'no', 'not found in literature or PPI databases',
]))  # fmt: skip


def extract(row):
    record = mod.resource.datasets()['interactions'].mapper(row)
    assert record is not None
    out = SilverExtractor('wojtowicz2020', 'interactions')
    out.process_record(record, row, 'fixture:1', 1)
    return record, out


def identifiers(entity):
    return [(i.type, i.value) for i in entity.identifiers]


def test_binding_pair_reaches_silver_extractor():
    record, out = extract(NEW)
    assert slot_name(record.predicate) == 'interacts_with'
    # The prey is the subject, the bait the object, as the table orders them.
    assert identifiers(record.subject) == [(Namespace.GENESYMBOL, 'HLA-B')]
    assert identifiers(record.object) == [
        (Namespace.GENESYMBOL, 'SIGLEC10'),
        (Namespace.NAME, 'Siglec-10'),
    ]
    assert [r.predicate for r in out.relations] == ['interacts_with']
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert {e.taxon for e in out.entities.values()} == {'9606'}
    assert {e.identifier for e in out.entities.values()} == {
        'HLA-B',
        'SIGLEC10',
    }


def test_roles_readings_and_novelty_are_not_turned_into_annotations():
    record, out = extract(NEW)
    assert {slot_name(a.term) for a in record.annotations} == {'in_taxon'}
    for entity in (record.subject, record.object):
        assert {slot_name(a.term) for a in entity.annotations} == {'in_taxon'}
    assert all(not e.annotations for e in out.entities.values())
    # The screen has no direction or effect to qualify.
    assert not any(
        'qualifier' in a['term'] for r in out.relations for a in r.annotations
    )


def test_entry_name_is_an_identifier_only_where_it_is_an_accession():
    record, _ = extract(UNREVIEWED_BAIT)
    # `I18BP` is an entry name without its species ending: not an accession.
    assert identifiers(record.subject) == [(Namespace.GENESYMBOL, 'IL18BP')]
    assert identifiers(record.object) == [
        (Namespace.UNIPROT, 'A0A0B4J1R5'),
        (Namespace.GENESYMBOL, 'KIR3DL3'),
    ]
    record, _ = extract(NO_GENE_SYMBOL)
    # The accession standing in for a common name is not repeated as a name.
    assert identifiers(record.subject) == [
        (Namespace.UNIPROT, 'A2IXV5'),
        (Namespace.GENESYMBOL, 'MIR cl-10'),
    ]


def test_homophilic_pair_is_kept_as_one_relation():
    record, out = extract(HOMOPHILIC)
    assert identifiers(record.subject) == identifiers(record.object)
    assert len(out.relations) == 1
    assert len(out.entities) == 1


def test_missing_participant_emits_no_relation():
    assert mod.map_wojtowicz_interaction(dict(NEW, Bait_Gene_Name=None)) is None
    assert mod.map_wojtowicz_interaction(dict(NEW, Prey_Gene_Name=None)) is None


def test_shared_reader_yields_the_sheet_rows_the_mapper_expects():
    book = openpyxl.Workbook()
    sheet = book.active
    sheet.title = 'Wojtowicz et al._DataS4 '
    sheet.append(COLUMNS)
    for row in (NEW, HOMOPHILIC):
        sheet.append([row[column] for column in COLUMNS])
    content = io.BytesIO()
    book.save(content)
    content.seek(0)

    dataset = mod.resource.datasets()['interactions']
    rows = list(dataset._raw_parser(SimpleNamespace(result=content)))

    assert [list(row) for row in rows] == [COLUMNS, COLUMNS]
    text = [column for column in COLUMNS if isinstance(NEW[column], str)]
    assert [[row[column] for column in text] for row in rows] == [
        [row[column] for column in text] for row in (NEW, HOMOPHILIC)
    ]
    assert all(dataset.mapper(row) is not None for row in rows)


def workbook(header, rows=()):
    book = openpyxl.Workbook()
    sheet = book.active
    sheet.title = 'Wojtowicz et al._DataS4 '
    sheet.append(header)
    for row in rows:
        sheet.append(row)
    content = io.BytesIO()
    book.save(content)
    content.seek(0)
    return SimpleNamespace(result=content)


def test_reader_rejects_a_changed_header_and_an_empty_sheet():
    dataset = mod.resource.datasets()['interactions']
    # A renamed column used to surface only as a KeyError in the mapper.
    renamed = [c if c != 'Prey_Gene_Name' else 'Prey_Gene' for c in COLUMNS]
    with pytest.raises(ValueError, match='Prey_Gene_Name'):
        list(
            dataset._raw_parser(
                workbook(renamed, [[NEW[c] for c in COLUMNS]])
            )
        )
    with pytest.raises(ValueError):
        list(dataset._raw_parser(workbook(COLUMNS)))
    assert set(mod.COLUMNS) <= set(COLUMNS)
