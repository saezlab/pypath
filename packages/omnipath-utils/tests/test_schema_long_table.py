"""T007 -- schema: the id_mapping_long long-value table and its indexes.

The long-value (name/structure) table is a *text* sibling of the database-ID
``id_mapping`` table; adding it must not change the untouched database-ID tables,
which stay ``varchar(64)`` (R2, FR-017/SC-008).
"""

from sqlalchemy import Text, String


class TestIdMappingLongSchema:
    def test_table_registered(self):
        from omnipath_utils.db._schema import Base, IdMappingLong

        assert IdMappingLong.__tablename__ == 'id_mapping_long'
        # Created by Base.metadata.create_all at schema init.
        assert 'omnipath_utils.id_mapping_long' in Base.metadata.tables

    def test_columns(self):
        from omnipath_utils.db._schema import IdMappingLong

        cols = {c.name for c in IdMappingLong.__table__.columns}
        assert {
            'source_type_id', 'target_type_id', 'ncbi_tax_id',
            'source_id', 'source_label', 'target_id', 'backend_id',
        } <= cols

    def test_value_columns_are_text(self):
        from omnipath_utils.db._schema import IdMappingLong

        c = IdMappingLong.__table__.c
        for name in ('source_id', 'source_label', 'target_id'):
            assert isinstance(c[name].type, Text), name

    def test_indexes(self):
        from omnipath_utils.db._schema import IdMappingLong

        idx = {i.name for i in IdMappingLong.__table__.indexes}
        assert 'idx_long_lookup' in idx
        assert 'idx_long_reverse' in idx

    def test_lookup_index_columns(self):
        from omnipath_utils.db._schema import IdMappingLong

        by_name = {i.name: i for i in IdMappingLong.__table__.indexes}
        lookup_cols = [c.name for c in by_name['idx_long_lookup'].columns]
        assert lookup_cols == ['source_type_id', 'target_type_id', 'source_id']
        rev_cols = [c.name for c in by_name['idx_long_reverse'].columns]
        assert rev_cols == ['target_type_id', 'source_type_id', 'target_id']

    def test_database_id_tables_untouched(self):
        """The database-ID hot-path table stays varchar(64) (SC-008)."""
        from omnipath_utils.db._schema import IdMapping

        for name in ('source_id', 'target_id'):
            t = IdMapping.__table__.c[name].type
            assert isinstance(t, String)
            assert t.length == 64
