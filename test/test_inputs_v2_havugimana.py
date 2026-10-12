"""Havugimana 2012 complexes: reader, parsed rows and Biolink mapping."""

import io
from types import SimpleNamespace

import openpyxl
import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.havugimana import COMPLEX_NAMESPACE, map_complex, resource
from pypath.inputs_v2.parsers.havugimana import HEADER, SHEET, iter_havugimana_complexes

# Rows of Table S3A (sheet ``Predicted_complexes`` of mmc3.xls), as the sheet has them.
ROWS = [
    ['C_1', 2.0, 'O60547,P30043', 'GMDS_HUMAN,BLVRB_HUMAN', 0.0, 'O60547,P30043,', '-',
     '-', '-', '-', '-', '-', '-', '-', '-', 'Fungi'],
    ['C_62', 2.0, 'P29084,P29083', 'T2EB_HUMAN,T2EA_HUMAN', 2.0, '-', 'P29083,P29084,',
     'CORUM_152::TFIIE complex (1.89e-05)', 'REACT_2368~TFIIE [nucleoplasm] (3.11e-05)',
     'PINdb_102::TFIIE Complex- transcription (5.27e-04)', 'HPRD_2870::in vitro (2.25e-03)', '-',
     '-', 'F_334 (5.86E-06)', 'Y_664 (4.13E-05)', 'multicellularity'],
    ['C_105', 2.0, 'P46934,P11532', 'NEDD4_HUMAN,DMD_HUMAN', 0.0, 'P11532,P46934,', '-',
     '-', '-', '-', '-', '-', '-', '-', '-', 'NA'],
    ['C_133', 3.0, 'Q16637,Q9UHI6,P38432', 'SMN_HUMAN,DDX20_HUMAN,COIL_HUMAN', 2.0, 'P38432,',
     'Q16637,Q9UHI6,', 'CORUM_3284::SMN complex (GEMIN5,2,3,4, SMN) (3.21e-04)',
     'REACT_10714~SMN complex [nucleoplasm] (9.16e-04)', '-', 'HPRD_786::in vivo (7.58e-04)', '-',
     'AM000236 (1.48e-03)', '-', '-', 'vertebrates'],
]


def workbook(rows=ROWS, sheet=SHEET, header=HEADER):
    book = openpyxl.Workbook()
    table = book.active
    table.title = sheet
    table.append(['Table S3A. List of 622 Putative Complexes Identified in This Study'])
    table.append([])
    table.append(list(header))
    for row in rows:
        table.append(row)
    out = io.BytesIO()
    book.save(out)
    out.seek(0)
    return SimpleNamespace(result=out)


def parsed(rows=ROWS):
    return {row['complex_id']: row for row in iter_havugimana_complexes(workbook(rows))}


def extract(row):
    record = map_complex(row)
    assert record is not None
    out = SilverExtractor('havugimana', 'complexes')
    out.process_record(record, row, 'fixture:1', 1)
    return record, out


def test_module_imports_with_native_schema():
    assert set(resource.datasets()) == {'complexes'}
    assert resource.config.name == 'Havugimana2012'


def test_reader_gives_complexes_with_subunits_and_comparisons():
    rows = parsed()
    assert list(rows) == ['C_1', 'C_62', 'C_105', 'C_133']
    smn = rows['C_133']
    assert smn['size'] == 3 and smn['accessions'] == ['Q16637', 'Q9UHI6', 'P38432']
    assert smn['entry_names'] == ['SMN_HUMAN', 'DDX20_HUMAN', 'COIL_HUMAN']
    assert smn['novel'] == ['P38432'] and smn['annotated'] == ['Q16637', 'Q9UHI6']
    assert smn['matched_subunits'] == 2
    assert smn['best_corum'].startswith('CORUM_3284::SMN complex') and smn['best_pindb'] is None
    assert rows['C_1']['annotated'] == [] and rows['C_62']['novel'] == []


def test_reader_rejects_inconsistent_or_other_sheets():
    with pytest.raises(ValueError, match='Unexpected header'):
        list(iter_havugimana_complexes(workbook(header=HEADER[:-1] + ('Group',))))
    with pytest.raises(ValueError, match='no complexes'):
        list(iter_havugimana_complexes(workbook([])))
    for position, value, message in (
        (0, 'X_1', 'Unknown complex id'),
        (1, 3.0, 'disagree'),
        (3, 'GMDS_HUMAN,BLVRB_MOUSE', 'another species'),
        (2, 'O60547,p30043', 'malformed UniProt'),
        (5, 'O60547,', 'novel and annotated'),
    ):
        row = list(ROWS[0])
        row[position] = value
        with pytest.raises(ValueError, match=message):
            parsed([row])


def test_complex_has_its_accession_and_unmodified_subunits():
    record, out = extract(parsed()['C_133'])
    (cplx,) = [e for e in out.entities.values() if e.entity_type == 'macromolecular_complex']
    assert (cplx.namespace, cplx.identifier) == (COMPLEX_NAMESPACE, 'C_133')
    assert [i['ns'] for i in cplx.identifiers] == [COMPLEX_NAMESPACE]
    members = {e.identifier: e for e in out.entities.values() if e.entity_type == 'protein'}
    assert sorted(members) == ['P38432', 'Q16637', 'Q9UHI6']
    assert all(e.namespace == Namespace.UNIPROT and e.taxon == '9606' for e in members.values())
    assert {i['ns'] for e in members.values() for i in e.identifiers} == {
        Namespace.UNIPROT, Namespace.UNIPROT_ENTRY,
    }
    assert [r.predicate for r in out.relations] == ['has_member'] * 3


def test_comparison_columns_are_not_asserted():
    _, out = extract(parsed()['C_62'])
    (cplx,) = [e for e in out.entities.values() if e.entity_type == 'macromolecular_complex']
    # no cross-reference to the CORUM, Reactome or PINdb complex that matches best
    assert cplx.annotations == []
    assert {i['ns'] for i in cplx.identifiers} == {COMPLEX_NAMESPACE}
    assert len(out.relations) == 2


def test_equal_subunits_stay_two_complexes():
    first = parsed()['C_62']
    second = {**first, 'complex_id': 'C_9999'}
    out = SilverExtractor('havugimana', 'complexes')
    for i, row in enumerate((first, second)):
        out.process_record(map_complex(row), row, f'fixture:{i}', i)
    complexes = [e for e in out.entities.values() if e.entity_type == 'macromolecular_complex']
    assert sorted(e.identifier for e in complexes) == ['C_62', 'C_9999']
    assert len(out.relations) == 4


def test_complex_without_subunits_gives_no_record():
    assert map_complex({**parsed()['C_1'], 'accessions': [], 'entry_names': []}) is None
