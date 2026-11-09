# The DuckDB database of an installation, as the Python client writes it
# (omnipath_client/database.py), which is also the database the hosted API serves a
# release from: native tables with global ids, the derived tables and indexes the
# contract names, a schema per resource, `resources`, `database_info`, and the query
# layer (Biolink names, the lookup macros and sql/*.sql). It is always written
# completely and put in place when done.

.database_file <- "omnipath.duckdb"

.record_tables <- function(record) {
  names <- vapply(record$files, function(f) f$name, "")
  all <- c(.tables(), .serving())
  all[all != .payloads() | paste0(.payloads(), ".parquet") %in% names]
}

.global_ids <- function() .contract()$global_ids

.source <- function(root, record, table, index = NULL) {
  path <- file.path(.resource_dir(root, record), paste0(table, ".parquet"))
  extra <- ""
  ids <- .global_ids()[[table]]
  if (!is.null(index) && !is.null(ids)) {
    extra <- sprintf(", (%d::BIGINT << 32) + %s AS %s", as.integer(index), ids[[1]], ids[[2]])
  }
  sprintf(
    "SELECT *, %s AS _resource, %s AS _resource_version%s FROM read_parquet(%s)",
    .quote(record$resource_id), .quote(record$version), extra, .quote(path)
  )
}

.state_rows <- function(records) {
  rows <- data.frame(
    resource_id = vapply(records, function(r) r$resource_id, ""),
    version = vapply(records, function(r) r$version, ""),
    has_evidence = vapply(records, function(r) .payloads() %in% .record_tables(r), TRUE),
    stringsAsFactors = FALSE
  )
  rows[order(rows$resource_id, rows$version), , drop = FALSE]
}

.exec <- function(con, sql, params = NULL) {
  if (is.null(params)) DBI::dbExecute(con, sql) else DBI::dbExecute(con, sql, params = params)
}

.describe <- function(con, root, records, release) {
  .exec(con, sprintf(paste0(
    "CREATE OR REPLACE TABLE database_info AS SELECT %s::VARCHAR AS data_dir, true AS materialized, ",
    "%s::VARCHAR AS query_layer, %s::VARCHAR AS release, %d::INTEGER AS database_format, ",
    "%s::TIMESTAMPTZ AS built_at, 'init'::VARCHAR AS built_by"),
    .quote(root), .quote(.query_layer()), if (is.null(release)) "NULL" else .quote(release),
    as.integer(.contract()$database_format), .quote(format(Sys.time(), "%Y-%m-%dT%H:%M:%S+00:00", tz = "UTC"))
  ))
  .exec(con, "CREATE OR REPLACE TABLE resources (resource_id VARCHAR, version VARCHAR, has_evidence BOOLEAN)")
  rows <- .state_rows(records)
  for (i in seq_len(nrow(rows))) {
    .exec(con, sprintf(
      "INSERT INTO resources VALUES (%s, %s, %s)",
      .quote(rows$resource_id[i]), .quote(rows$version[i]), if (rows$has_evidence[i]) "true" else "false"
    ))
  }
}

.resource_schema <- function(con, record) {
  schema <- .name(record$resource_id)
  .exec(con, sprintf("CREATE SCHEMA IF NOT EXISTS %s", schema))
  for (table in .record_tables(record)) {
    .exec(con, sprintf(
      "CREATE OR REPLACE VIEW %s.%s AS SELECT * FROM main.%s WHERE _resource = %s",
      schema, table, table, .quote(record$resource_id)
    ))
  }
}

.configure <- function(con, root, memory_limit) {
  .exec(con, sprintf("SET memory_limit = %s", .quote(memory_limit)))
  .exec(con, sprintf("SET temp_directory = %s", .quote(file.path(root, ".duckdb-tmp"))))
}

.native_lookups <- function() {
  vapply(.contract()$lookups, function(lookup) {
    table <- lookup[[1]]; column <- lookup[[2]]; kind <- lookup[[3]]
    global_id <- .global_ids()[[table]][[2]]
    if (kind == "at") {
      sprintf(paste0(
        "CREATE OR REPLACE MACRO op_%s_at(pairs) AS TABLE SELECT * EXCLUDE (%s) FROM %s WHERE %s IN (",
        "SELECT (i.idx << 32) + w.p.v FROM (SELECT unnest(pairs) AS p) w ",
        "JOIN op_resource_index i ON i._resource = w.p.r)"
      ), table, global_id, table, global_id)
    } else {
      sprintf(
        "CREATE OR REPLACE MACRO op_%s_with(vals) AS TABLE SELECT * EXCLUDE (%s) FROM %s WHERE %s IN (SELECT unnest(vals))",
        table, global_id, table, column
      )
    }
  }, "")
}

.query_layer_write <- function(con) {
  .exec(con, "CREATE OR REPLACE TABLE op_names (kind VARCHAR, alias VARCHAR, name VARCHAR)")
  names <- .contract()$biolink_names
  values <- vapply(names, function(n) sprintf("(%s, %s, %s)", .quote(n[[1]]), .quote(n[[2]]), .quote(n[[3]])), "")
  .exec(con, paste0("INSERT INTO op_names VALUES ", paste(values, collapse = ", ")))
  present <- DBI::dbGetQuery(con, "SELECT table_name AS n FROM duckdb_tables() WHERE schema_name = 'main' AND database_name = current_database()")$n
  for (statement in .native_lookups()) {
    table <- sub("^CREATE OR REPLACE MACRO op_(\\w+)_(at|with)\\(.*$", "\\1", statement)
    if (table %in% present) .exec(con, statement)
  }
  for (name in unlist(.contract()$macro_files)) {
    sql <- paste(readLines(file.path(.sql_dir(), name), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    .exec(con, sql)
  }
  .exec(con, sprintf("UPDATE database_info SET query_layer = %s", .quote(.query_layer())))
}

# The derived tables (each sorted) and the indexes, once the tables are written.
.derived <- function(con, present) {
  derived <- .contract()$derived_tables
  for (name in names(derived)) {
    # A derived table that replaces a written one takes its name when complete.
    replaces <- isTRUE(derived[[name]]$replaces)
    target <- if (replaces) paste0(name, "__derived") else name
    .exec(con, sprintf("CREATE TABLE %s AS SELECT * FROM (%s) ORDER BY %s",
      target, derived[[name]]$select, derived[[name]]$order))
    if (replaces) {
      .exec(con, sprintf("DROP TABLE %s", name))
      .exec(con, sprintf("ALTER TABLE %s RENAME TO %s", target, name))
    }
  }
  for (index in .contract()$indexes) {
    if (index[[2]] %in% present) {
      .exec(con, sprintf("CREATE INDEX %s ON %s(%s)", index[[1]], index[[2]], index[[3]]))
    }
  }
}

.sorted_records <- function(records) {
  records[order(vapply(records, function(r) r$resource_id, ""))]
}

.rebuild <- function(root, records, memory_limit, release = NULL) {
  staging_dir <- tempfile(".database-", tmpdir = root)
  dir.create(staging_dir)
  on.exit(unlink(staging_dir, recursive = TRUE), add = TRUE)
  staging <- file.path(staging_dir, .database_file)
  con <- DBI::dbConnect(duckdb::duckdb())
  closed <- FALSE
  on.exit(if (!closed) DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE, after = FALSE)
  # Written in the oldest supported DuckDB's storage format, so clients on another
  # DuckDB version can read it.
  .exec(con, sprintf(
    "ATTACH %s AS omnipath (STORAGE_VERSION %s)", .quote(staging), .quote(.contract()$storage_version)
  ))
  .exec(con, "USE omnipath")
  .configure(con, root, memory_limit)
  ordered <- .sorted_records(records)
  index <- stats::setNames(as.list(seq_along(ordered) - 1L), vapply(ordered, function(r) r$resource_id, ""))
  .exec(con, "CREATE TABLE op_resource_index (_resource VARCHAR, idx BIGINT)")
  for (name in names(index)) {
    .exec(con, sprintf("INSERT INTO op_resource_index VALUES (%s, %d)", .quote(name), as.integer(index[[name]])))
  }
  present <- Filter(function(t) any(vapply(ordered, function(r) t %in% .record_tables(r), TRUE)), c(.tables(), .serving()))
  for (table in present) {
    having <- Filter(function(r) table %in% .record_tables(r), ordered)
    # The columns of every resource (union by name), then each resource's rows in its
    # file's order, so rows follow their global ids.
    union <- paste(vapply(having, function(r) .source(root, r, table, index[[r$resource_id]]), ""),
      collapse = " UNION ALL BY NAME ")
    .exec(con, sprintf("CREATE TABLE %s AS SELECT * FROM (%s) LIMIT 0", table, union))
    for (record in having) {
      .exec(con, sprintf("INSERT INTO %s BY NAME %s", table, .source(root, record, table, index[[record$resource_id]])))
    }
  }
  .derived(con, present)
  for (record in ordered) .resource_schema(con, record)
  .describe(con, root, records, release)
  .query_layer_write(con)
  .exec(con, "USE memory")
  .exec(con, "DETACH omnipath")
  DBI::dbDisconnect(con, shutdown = TRUE)
  closed <- TRUE
  target <- file.path(root, .database_file)
  if (file.exists(target)) unlink(target)
  if (!file.rename(staging, target)) .fail("Could not put the database in place: ", target)
  invisible(TRUE)
}

.info <- function(root) {
  path <- file.path(root, .database_file)
  if (!file.exists(path)) return(NULL)
  con <- tryCatch(DBI::dbConnect(duckdb::duckdb(), dbdir = path, read_only = TRUE), error = function(e) NULL)
  if (is.null(con)) return(NULL)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  tryCatch({
    columns <- DBI::dbGetQuery(con, "DESCRIBE database_info")$column_name
    if (!("database_format" %in% columns)) return(NULL)
    info <- DBI::dbGetQuery(con, "SELECT query_layer, database_format FROM database_info")
    state <- DBI::dbGetQuery(con, "SELECT * FROM resources ORDER BY ALL")
    list(query_layer = info$query_layer, database_format = info$database_format, state = state)
  }, error = function(e) NULL)
}

.tables_in_sync <- function(info, installation) {
  if (is.null(info)) return(FALSE)
  wanted <- .state_rows(installation$resources)
  same_state <- nrow(info$state) == nrow(wanted) &&
    all(info$state$resource_id == wanted$resource_id) &&
    all(info$state$version == wanted$version) &&
    all(info$state$has_evidence == wanted$has_evidence)
  isTRUE(info$database_format == .contract()$database_format) && same_state
}

.in_sync <- function(root, installation) {
  info <- .info(root)
  .tables_in_sync(info, installation) && identical(info$query_layer, .query_layer())
}

.layer_only <- function(root) {
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = file.path(root, .database_file))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbBegin(con)
  .query_layer_write(con)
  DBI::dbCommit(con)
  invisible(TRUE)
}

# Changed resources mean a new database, written completely; when only the query
# layer differs (another client version wrote it), the query layer alone.
.synchronize <- function(root, installation, memory_limit = "2GB") {
  info <- .info(root)
  if (.tables_in_sync(info, installation)) {
    if (!identical(info$query_layer, .query_layer())) .layer_only(root)
    return(invisible(FALSE))
  }
  on.exit({
    temporary <- file.path(root, ".duckdb-tmp")
    if (dir.exists(temporary) && !length(list.files(temporary, all.files = TRUE, no.. = TRUE))) unlink(temporary, recursive = TRUE)
  }, add = TRUE)
  .rebuild(root, installation$resources, memory_limit, installation$release)
  invisible(TRUE)
}
