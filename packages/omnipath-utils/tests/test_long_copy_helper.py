"""T009 -- the long-value COPY projection.

Name keys are lowercased + whitespace-stripped (FR-002) with the original kept
as ``source_label``; structure keys (inchi/smiles) are stored verbatim, case
preserved (FR-019). Over-long values are skipped, never truncated. The DB-backed
delete-own-slice idempotency check (R6) runs only when a test DB is configured.
"""

import os

import pytest

from omnipath_utils.db._build import DatabaseBuilder

_rows = DatabaseBuilder._long_rows


class TestNameKeys:
    def test_lowercased_and_stripped(self):
        rows, skipped = _rows(
            {'  Taurine ': {'CHEBI:15891'}}, 1, 2, 3, is_name=True,
        )
        assert rows == [(1, 2, 0, 'taurine', 'Taurine', 'CHEBI:15891', 3)]
        assert skipped == 0

    def test_case_variants_collapse_to_one_key(self):
        out = {
            _rows({v: {'CHEBI:1'}}, 1, 2, 3, is_name=True)[0][0][3]
            for v in ('Taurine', 'taurine', 'TAURINE', ' TauRine  ')
        }
        assert out == {'taurine'}

    def test_label_preserves_original_case(self):
        rows, _ = _rows({'D-Glucose': {'CHEBI:17634'}}, 1, 2, 3, is_name=True)
        assert rows[0][3] == 'd-glucose'   # lookup key
        assert rows[0][4] == 'D-Glucose'   # display label


class TestStructureKeys:
    def test_verbatim_case_preserved(self):
        inchi = 'InChI=1S/C2H7NO3S/c3-1-2-7(4,5)6/h1-3H2,(H,4,5,6)'
        rows, _ = _rows({inchi: {'CHEBI:15891'}}, 1, 2, 3, is_name=False)
        assert rows[0][3] == inchi          # byte-identical, not folded
        assert rows[0][4] is None           # no display label for structures

    def test_smiles_case_sensitive(self):
        # 'C' (carbon) != 'c' (aromatic carbon): folding would corrupt it.
        smiles = 'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O'
        rows, _ = _rows({smiles: {'CHEBI:17634'}}, 1, 2, 3, is_name=False)
        assert rows[0][3] == smiles


class TestGuards:
    def test_oversized_source_skipped_not_truncated(self):
        big = 'x' * 3000
        rows, skipped = _rows({big: {'a'}}, 1, 2, 3, is_name=False, max_key=2000)
        assert rows == []
        assert skipped == 1

    def test_oversized_target_skipped(self):
        rows, skipped = _rows(
            {'aspirin': {'y' * 3000}}, 1, 2, 3, is_name=True, max_key=2000,
        )
        assert rows == []
        assert skipped == 1

    def test_empty_values_skipped(self):
        rows, skipped = _rows(
            {'  ': {'t'}, 'x': {'', '  '}}, 1, 2, 3, is_name=True,
        )
        assert rows == []
        assert skipped == 3  # empty source + 2 empty targets

    def test_limit_caps_rows(self):
        data = {f'n{i}': {'t'} for i in range(10)}
        rows, _ = _rows(data, 1, 2, 3, is_name=True, limit=3)
        assert len(rows) == 3

    def test_one_to_many_preserved(self):
        rows, _ = _rows(
            {'glucose': {'CHEBI:17234', 'CHEBI:4167'}}, 1, 2, 3, is_name=True,
        )
        assert len(rows) == 2
        assert {r[5] for r in rows} == {'CHEBI:17234', 'CHEBI:4167'}


class TestStructureBackends:
    def test_hmdb_is_a_structure_backend(self):
        # spec 011 T081: hmdb has its own precomputed inchikey column
        # (id_types.yaml's inchikey: backends), same shape as the other
        # five structure backends -- it was simply missing from this list.
        assert 'hmdb' in DatabaseBuilder._STRUCTURE_BACKENDS


@pytest.mark.skipif(
    not os.environ.get('OMNIPATH_UTILS_TEST_DB'),
    reason='set OMNIPATH_UTILS_TEST_DB to run the DB-backed idempotency check',
)
class TestSliceIdempotency:
    """delete-own-slice + COPY: a re-run leaves an equal row count (R6)."""

    def test_rerun_equal_rows(self):
        from sqlalchemy import text
        from sqlalchemy.orm import Session

        from omnipath_utils.db._connection import SCHEMA

        builder = DatabaseBuilder(db_url=os.environ['OMNIPATH_UTILS_TEST_DB'])
        builder.create_tables()
        builder.populate_id_types()
        builder.populate_backends()
        data = {'taurine': {'CHEBI:15891'}, 'glucose': {'CHEBI:17234'}}

        n1 = builder._populate_long_slice(data, 'name', 'chebi', 'chebi')
        n2 = builder._populate_long_slice(data, 'name', 'chebi', 'chebi')
        assert n1 == n2 == 2

        with Session(builder.engine) as s:
            total = s.execute(
                text(
                    f"SELECT count(*) FROM {SCHEMA}.id_mapping_long m "
                    f"JOIN {SCHEMA}.id_type st ON st.id=m.source_type_id "
                    f"JOIN {SCHEMA}.backend b ON b.id=m.backend_id "
                    "WHERE st.name='name' AND b.name='chebi'"
                )
            ).scalar()
        assert total == 2
