"""Sequence assertions survive parser projection and per-chain mapping."""

import csv
import io

import pytest
from types import SimpleNamespace

from pypath.inputs_v2._molecular_forms import combine_forms, sequence_form
from pypath.inputs_v2.bindingdb import _target
from pypath.inputs_v2.parsers.bindingdb import (
    _iter_csv_fallback,
    _iter_duckdb_tsv,
)


def test_sequence_identity_is_exact_and_alphabet_scoped():
    first = sequence_form('m aa k')
    assert first == sequence_form('MAAK')
    assert first != sequence_form('MAAT')
    assert sequence_form('ACG', system='transcript') != sequence_form(
        'ACG', system='genomic'
    )
    assert sequence_form('sequence unavailable') is None
    assert sequence_form('ACGT', system='transcript') is None
    combined = combine_forms(
        {'isoform_identifier': {'ns': 'uniprot', 'id': 'P04637-2'}}, first
    )
    assert combined['isoform_identifier']['id'] == 'P04637-2'
    assert combined['sequence_identifiers'] == first['sequence_identifiers']


def test_bindingdb_sequence_projection_matches_csv_and_duckdb(tmp_path):
    row = {
        'BindingDB Reactant_set_id': '1',
        'Target Name': 'explicit two-chain target',
        'Number of Protein Chains in Target (>1 implies a multichain complex)': '2',
        'UniProt (SwissProt) Primary ID of Target Chain 1': 'P04637-2',
        'UniProt (SwissProt) Primary ID of Target Chain 2': 'P04637-3',
        'BindingDB Target Chain 1 Sequence': 'MAAK',
        'BindingDB Target Chain 2 Sequence': 'MAAT',
        'BindingDB Target Chain 51 Sequence': 'MAAG',
    }
    data = io.StringIO()
    writer = csv.DictWriter(data, fieldnames=list(row), delimiter='\t')
    writer.writeheader()
    writer.writerow(row)
    path = tmp_path / 'binding.tsv'
    path.write_text(data.getvalue())
    csv_rows = list(
        _iter_csv_fallback(
            SimpleNamespace(result={'data': io.StringIO(data.getvalue())})
        )
    )
    sql_rows = list(_iter_duckdb_tsv(path))
    assert csv_rows == sql_rows
    mapped = _target(sql_rows[0])
    first, second = [membership.member for membership in mapped.membership]
    assert first.molecular_form['isoform_identifier']['id'] == 'P04637-2'
    assert second.molecular_form['isoform_identifier']['id'] == 'P04637-3'
    assert (
        first.molecular_form['sequence_identifiers']
        != second.molecular_form['sequence_identifiers']
    )
    assert sql_rows[0]['BindingDB Target Chain 51 Sequence'] == 'MAAG'


def test_bindingdb_unnamed_chain_survives_with_only_explicit_sequence():
    row = {
        'Target Name': 'complex',
        'Number of Protein Chains in Target (>1 implies a multichain complex)': '2',
        'BindingDB Target Chain 1 Sequence': 'MAAK',
        'BindingDB Target Chain 2 Sequence': 'MAAT',
    }
    members = _target(row).membership
    assert len(members) == 2
    assert members[0].member.identifiers[0].type == 'protein_sequence_sha256'
    assert members[0].member.identifiers != members[1].member.identifiers


def test_bindingdb_legacy_sequence_header_is_retained(tmp_path):
    path = tmp_path / 'binding.tsv'
    path.write_text('BindingDB Target Chain  Sequence\nMAAK\n')
    row = next(_iter_duckdb_tsv(path))
    assert row['BindingDB Target Chain 1 Sequence'] == 'MAAK'


@pytest.mark.parametrize('header', ['BindingDB Target Chain Sequence 2', 'BindingDB Target Chain 2 Sequence'])
def test_bindingdb_source_sequence_numbering_is_normalized(header, tmp_path):
    path = tmp_path / 'binding.tsv'
    path.write_text(header + '\nMAAK\n')
    assert next(_iter_duckdb_tsv(path))['BindingDB Target Chain 2 Sequence'] == 'MAAK'
    row = next(_iter_csv_fallback(SimpleNamespace(result={'data': io.StringIO(path.read_text())})))
    assert row['BindingDB Target Chain 2 Sequence'] == 'MAAK'


def test_bindingdb_official_header_and_row_keep_the_reported_chain_sequence(tmp_path):
    from pathlib import Path
    fixture = Path(__file__).parent / 'data/molecular/bindingdb-202610-header-row.tsv'
    row, = _iter_duckdb_tsv(fixture)
    assert row['BindingDB Target Chain 1 Sequence'] == row['BindingDB Target Chain Sequence 1']
    assert _target(row).molecular_form['sequence_identifiers']
