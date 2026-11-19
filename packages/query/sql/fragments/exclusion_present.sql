-- Whether the served release has the exclusion metabolite subset; {{input}} is unused.
(SELECT count(*) > 0 FROM duckdb_tables()
 WHERE database_name = current_database() AND table_name = 'exclusion_metabolite')
