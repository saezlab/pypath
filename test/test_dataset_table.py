"""Datasets with a table parser: raw() streams the SQL table's rows."""

from pypath.inputs_v2.base import Dataset

LINES = [{'id': str(i), 'keep': 'yes' if i % 3 else 'no', 'note': None if i % 2 else 'even'} for i in range(100)]


def _rows(opener=None, max_lines=None, **_kwargs):
    for line in LINES[:max_lines]:
        if line['keep'] == 'yes':
            yield {k: v for k, v in line.items() if v is not None}


def _table(db, name, opener, max_lines=None, **_kwargs):
    import pyarrow as pa

    lines = LINES[:max_lines]
    db.register('lines', pa.Table.from_pylist([{'line': i, **row} for i, row in enumerate(lines)]))
    db.execute(f"""CREATE OR REPLACE TABLE {name} AS
        SELECT row_number() OVER (ORDER BY line) - 1 AS rid, id, keep, note
        FROM lines WHERE keep = 'yes'""")
    db.unregister('lines')
    return db.execute(f'SELECT count(*) FROM {name}').fetchone()[0]


def _dataset():
    return Dataset(download=None, mapper=lambda row: row, raw_parser=_rows, raw_table=_table)


def test_table_rows_equal_row_parser_rows(tmp_path, monkeypatch):
    monkeypatch.setenv('PYPATH_TABLE_TMPDIR', str(tmp_path))
    dataset = _dataset()
    assert dataset.has_table
    assert list(dataset.raw()) == list(dataset.raw(row_parser=True)) == list(_rows())


def test_max_records_parses_a_prefix(tmp_path, monkeypatch):
    monkeypatch.setenv('PYPATH_TABLE_TMPDIR', str(tmp_path))
    rows = list(_dataset().raw(max_records=10))
    assert rows == list(_rows())[:10]
    assert list(_dataset().raw(max_records=1000)) == list(_rows())


def test_max_lines_limits_the_parse(tmp_path):
    import duckdb

    db = duckdb.connect()
    assert _dataset().table(db, 'raw', max_lines=9) == len(list(_rows(max_lines=9)))
