-- Whether the served release has the metabolite_class subset; {{input}} is unused.
(SELECT count(*) > 0 FROM duckdb_tables()
 WHERE database_name = current_database() AND table_name = 'metabolite_class')
