-- Whether the served release has the MetSigDB subset; {{input}} is unused.
(SELECT count(*) > 0 FROM duckdb_tables()
 WHERE database_name = current_database() AND table_name = 'metsigdb_membership')
