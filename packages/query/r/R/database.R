# The DuckDB database of an installation, as the Python client writes it
# (omnipath_client/database.py): views over the Parquet files or native tables with
# global ids, a schema per resource, `resources`, `database_info`, and the query
# layer (Biolink names, the installation's lookup macros and sql/*.sql).

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

.describe <- function(con, root, materialized, records) {
  .exec(con, sprintf(
    "CREATE OR REPLACE TABLE database_info AS SELECT %s::VARCHAR AS data_dir, %s::BOOLEAN AS materialized, %s::VARCHAR AS query_layer",
    .quote(root), if (materialized) "true" else "false", .quote(.query_layer())
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

.resource_schema <- function(con, root, record, materialized) {
  schema <- .name(record$resource_id)
  .exec(con, sprintf("CREATE SCHEMA IF NOT EXISTS %s", schema))
  for (table in .record_tables(record)) {
    source <- if (materialized) {
      sprintf("SELECT * FROM main.%s WHERE _resource = %s", table, .quote(record$resource_id))
    } else {
      .source(root, record, table)
    }
    .exec(con, sprintf("CREATE OR REPLACE VIEW %s.%s AS %s", schema, table, source))
  }
}

.configure <- function(con, root, memory_limit) {
  .exec(con, sprintf("SET memory_limit = %s", .quote(memory_limit)))
  .exec(con, sprintf("SET temp_directory = %s", .quote(file.path(root, ".duckdb-tmp"))))
}

.lookup_parts <- function(kind, column) {
  # The argument is evaluated once, into `wanted` (see database.py).
  if (kind == "at") {
    list(
      parameter = "pairs",
      wanted = "SELECT p.r AS r, p.v AS v FROM (SELECT unnest(pairs) AS p)",
      condition = sprintf("%s IN (SELECT v FROM wanted WHERE r = '{r}')", column)
    )
  } else {
    list(
      parameter = "vals",
      wanted = "SELECT unnest(vals) AS v",
      condition = sprintf("%s IN (SELECT v FROM wanted)", column)
    )
  }
}

.view_lookups <- function(resources) {
  vapply(.contract()$lookups, function(lookup) {
    table <- lookup[[1]]; column <- lookup[[2]]; kind <- lookup[[3]]
    parts <- .lookup_parts(kind, column)
    branches <- vapply(sort(resources), function(r) {
      sprintf("SELECT * FROM %s.%s WHERE %s", .name(r), table, gsub("{r}", r, parts$condition, fixed = TRUE))
    }, "")
    sprintf(
      "CREATE OR REPLACE MACRO op_%s_%s(%s) AS TABLE WITH wanted AS MATERIALIZED (%s) %s",
      table, kind, parts$parameter, parts$wanted, paste(branches, collapse = " UNION ALL BY NAME ")
    )
  }, "")
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

.query_layer_write <- function(con, records, materialized) {
  .exec(con, "CREATE OR REPLACE TABLE op_names (kind VARCHAR, alias VARCHAR, name VARCHAR)")
  names <- .contract()$biolink_names
  values <- vapply(names, function(n) sprintf("(%s, %s, %s)", .quote(n[[1]]), .quote(n[[2]]), .quote(n[[3]])), "")
  .exec(con, paste0("INSERT INTO op_names VALUES ", paste(values, collapse = ", ")))
  present <- c(
    DBI::dbGetQuery(con, "SELECT table_name AS n FROM duckdb_tables() WHERE schema_name = 'main' AND database_name = current_database()")$n,
    DBI::dbGetQuery(con, "SELECT view_name AS n FROM duckdb_views() WHERE schema_name = 'main' AND database_name = current_database()")$n
  )
  resources <- vapply(records, function(r) r$resource_id, "")
  statements <- if (materialized) .native_lookups() else .view_lookups(resources)
  for (statement in statements) {
    table <- sub("^CREATE OR REPLACE MACRO op_(\\w+)_(at|with)\\(.*$", "\\1", statement)
    if (table %in% present) .exec(con, statement)
  }
  for (name in unlist(.contract()$macro_files)) {
    sql <- paste(readLines(file.path(.sql_dir(), name), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    .exec(con, sql)
  }
  .exec(con, sprintf("UPDATE database_info SET query_layer = %s", .quote(.query_layer())))
}

.resource_index <- function(records, existing) {
  wanted <- vapply(records, function(r) r$resource_id, "")
  index <- existing[names(existing) %in% wanted]
  following <- if (length(existing)) max(unlist(existing)) + 1 else 0
  for (name in sort(wanted)) {
    if (is.null(index[[name]])) {
      index[[name]] <- following
      following <- following + 1
    }
  }
  index
}

.write_index <- function(con, index) {
  .exec(con, "CREATE OR REPLACE TABLE op_resource_index (_resource VARCHAR, idx BIGINT)")
  for (name in sort(names(index))) {
    .exec(con, sprintf("INSERT INTO op_resource_index VALUES (%s, %d)", .quote(name), as.integer(index[[name]])))
  }
}

.sorted_records <- function(records) {
  records[order(vapply(records, function(r) r$resource_id, ""))]
}

.rebuild <- function(root, records, materialized, memory_limit) {
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
  present <- Filter(function(t) any(vapply(ordered, function(r) t %in% .record_tables(r), TRUE)), c(.tables(), .serving()))
  if (materialized) {
    index <- .resource_index(ordered, list())
    .write_index(con, index)
    for (table in present) {
      having <- Filter(function(r) table %in% .record_tables(r), ordered)
      union <- paste(vapply(having, function(r) .source(root, r, table, index[[r$resource_id]]), ""),
        collapse = " UNION ALL BY NAME ")
      .exec(con, sprintf("CREATE TABLE %s AS SELECT * FROM (%s) LIMIT 0", table, union))
      for (record in having) {
        .exec(con, sprintf("INSERT INTO %s BY NAME %s", table, .source(root, record, table, index[[record$resource_id]])))
      }
    }
  }
  for (record in ordered) .resource_schema(con, root, record, materialized)
  if (!materialized) {
    for (table in present) {
      having <- Filter(function(r) table %in% .record_tables(r), ordered)
      union <- paste(vapply(having, function(r) sprintf("SELECT * FROM %s.%s", .name(r$resource_id), table), ""),
        collapse = " UNION ALL BY NAME ")
      .exec(con, sprintf("CREATE VIEW %s AS %s", table, union))
    }
  }
  .describe(con, root, materialized, records)
  .query_layer_write(con, records, materialized)
  .exec(con, "USE memory")
  .exec(con, "DETACH omnipath")
  DBI::dbDisconnect(con, shutdown = TRUE)
  closed <- TRUE
  target <- file.path(root, .database_file)
  if (file.exists(target)) unlink(target)
  if (!file.rename(staging, target)) .fail("Could not put the database in place: ", target)
  invisible(TRUE)
}

.update <- function(root, records, memory_limit) {
  path <- file.path(root, .database_file)
  con <- tryCatch(DBI::dbConnect(duckdb::duckdb(), dbdir = path), error = function(e) NULL)
  if (is.null(con)) return(FALSE)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  tryCatch({
    .configure(con, root, memory_limit)
    if (!isTRUE(DBI::dbGetQuery(con, "SELECT materialized FROM database_info")$materialized)) return(FALSE)
    current <- DBI::dbGetQuery(con, "SELECT * FROM resources")
    wanted <- .state_rows(records)
    key <- function(d) paste(d$resource_id, d$version, d$has_evidence, sep = "\r")
    stale <- sort(unique(c(
      current$resource_id[!key(current) %in% key(wanted)],
      wanted$resource_id[!key(wanted) %in% key(current)]
    )))
    by_name <- stats::setNames(records, vapply(records, function(r) r$resource_id, ""))
    existing <- DBI::dbGetQuery(con, "SELECT table_name AS n FROM duckdb_tables() WHERE schema_name = 'main'")$n
    rows <- DBI::dbGetQuery(con, "SELECT * FROM op_resource_index")
    index <- .resource_index(records, stats::setNames(as.list(rows$idx), rows$`_resource`))
    DBI::dbBegin(con)
    committed <- FALSE
    on.exit(if (!committed) try(DBI::dbRollback(con), silent = TRUE), add = TRUE, after = FALSE)
    .write_index(con, index)
    for (table in c(.tables(), .serving())) {
      fresh <- Filter(function(r) table %in% .record_tables(r), by_name[intersect(stale, names(by_name))])
      if (!(table %in% existing)) {
        if (length(fresh)) return(FALSE)
        next
      }
      if (length(stale)) {
        .exec(con, sprintf("DELETE FROM %s WHERE _resource IN (%s)", table, paste(.quote(stale), collapse = ", ")))
      }
      for (record in fresh) {
        .exec(con, sprintf("INSERT INTO %s BY NAME %s", table, .source(root, record, table, index[[record$resource_id]])))
      }
    }
    for (name in stale) {
      if (!is.null(by_name[[name]])) {
        .resource_schema(con, root, by_name[[name]], TRUE)
      } else {
        .exec(con, sprintf("DROP SCHEMA IF EXISTS %s CASCADE", .name(name)))
      }
    }
    .describe(con, root, TRUE, records)
    .query_layer_write(con, records, TRUE)
    DBI::dbCommit(con)
    committed <- TRUE
    TRUE
  }, error = function(e) FALSE)
}

.info <- function(root) {
  path <- file.path(root, .database_file)
  if (!file.exists(path)) return(NULL)
  con <- tryCatch(DBI::dbConnect(duckdb::duckdb(), dbdir = path, read_only = TRUE), error = function(e) NULL)
  if (is.null(con)) return(NULL)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  tryCatch({
    columns <- DBI::dbGetQuery(con, "DESCRIBE database_info")$column_name
    layer <- if ("query_layer" %in% columns) "query_layer" else "NULL::VARCHAR AS query_layer"
    info <- DBI::dbGetQuery(con, sprintf("SELECT data_dir, materialized, %s FROM database_info", layer))
    state <- DBI::dbGetQuery(con, "SELECT * FROM resources ORDER BY ALL")
    list(data_dir = info$data_dir, materialized = info$materialized, query_layer = info$query_layer, state = state)
  }, error = function(e) NULL)
}

.tables_in_sync <- function(info, root, installation) {
  if (is.null(info)) return(FALSE)
  wanted <- .state_rows(installation$resources)
  same_state <- nrow(info$state) == nrow(wanted) &&
    all(info$state$resource_id == wanted$resource_id) &&
    all(info$state$version == wanted$version) &&
    all(info$state$has_evidence == wanted$has_evidence)
  isTRUE(info$materialized == installation$materialize) &&
    (isTRUE(info$materialized) || identical(info$data_dir, root)) &&
    same_state
}

.in_sync <- function(root, installation) {
  info <- .info(root)
  .tables_in_sync(info, root, installation) && identical(info$query_layer, .query_layer())
}

.layer_only <- function(root, installation) {
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = file.path(root, .database_file))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  tryCatch({
    tables <- DBI::dbGetQuery(con, "SELECT table_name AS n FROM duckdb_tables() WHERE schema_name = 'main'")$n
    if (installation$materialize && !("op_resource_index" %in% tables)) return(FALSE)
    if (!("query_layer" %in% DBI::dbGetQuery(con, "DESCRIBE database_info")$column_name)) {
      .exec(con, "ALTER TABLE database_info ADD COLUMN query_layer VARCHAR")
    }
    DBI::dbBegin(con)
    .query_layer_write(con, installation$resources, installation$materialize)
    DBI::dbCommit(con)
    TRUE
  }, error = function(e) {
    try(DBI::dbRollback(con), silent = TRUE)
    FALSE
  })
}

.synchronize <- function(root, installation, memory_limit = "2GB") {
  info <- .info(root)
  in_sync <- .tables_in_sync(info, root, installation)
  if (in_sync && identical(info$query_layer, .query_layer())) return(invisible(FALSE))
  on.exit({
    temporary <- file.path(root, ".duckdb-tmp")
    if (dir.exists(temporary) && !length(list.files(temporary, all.files = TRUE, no.. = TRUE))) unlink(temporary, recursive = TRUE)
  }, add = TRUE)
  if (in_sync && .layer_only(root, installation)) return(invisible(TRUE))
  records <- installation$resources
  if (!(isTRUE(installation$materialize) && .update(root, records, memory_limit))) {
    .rebuild(root, records, isTRUE(installation$materialize), memory_limit)
  }
  invisible(TRUE)
}
