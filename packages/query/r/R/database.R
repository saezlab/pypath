# The DuckDB database of an installation, as the Python client writes it
# (omnipath_client/database.py), which is also the database the hosted API serves a
# release from: the resources' tables written to the build's scratch database, the
# database's own tables derived from them in build order (every entity and relation
# merged, numbered, with their identifiers, annotations, evidence and lookup tables),
# the indexes the contract names, `resources`, `database_info`, and the query layer
# (Biolink names, rule tables and sql/*.sql). It is always written completely and put in
# place when done.

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

# The release's taxon names (common name, else scientific name), empty without its
# taxonomy reference.
.taxonomy <- function(con, path = NULL) {
  if (!is.null(path) && file.exists(path)) {
    name <- "coalesce(nullif(common_name, ''), nullif(scientific_name, ''))"
    .exec(con, sprintf(paste0(
      "CREATE OR REPLACE TABLE op_taxonomy AS SELECT DISTINCT taxon_id::VARCHAR AS taxon, %s AS name ",
      "FROM read_parquet(%s) WHERE %s IS NOT NULL ORDER BY taxon"), name, .quote(path), name))
  } else {
    .exec(con, "CREATE OR REPLACE TABLE op_taxonomy (taxon VARCHAR, name VARCHAR)")
  }
}

.describe <- function(con, root, records, release, taxonomy = NULL) {
  .exec(con, sprintf(paste0(
    "CREATE OR REPLACE TABLE database_info AS SELECT %s::VARCHAR AS data_dir, true AS materialized, ",
    "%s::VARCHAR AS query_layer, %s::VARCHAR AS release, %d::INTEGER AS database_format, ",
    "%s::TIMESTAMPTZ AS built_at, 'init'::VARCHAR AS built_by, %s::VARCHAR AS taxonomy"),
    .quote(root), .quote(.query_layer()), if (is.null(release)) "NULL" else .quote(release),
    as.integer(.contract()$database_format), .quote(format(Sys.time(), "%Y-%m-%dT%H:%M:%S+00:00", tz = "UTC")),
    if (is.null(taxonomy)) "NULL" else .quote(taxonomy)
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

.configure <- function(con, root, memory_limit) {
  .exec(con, sprintf("SET memory_limit = %s", .quote(memory_limit)))
  .exec(con, sprintf("SET temp_directory = %s", .quote(file.path(root, ".duckdb-tmp"))))
}

# The Biolink names (op_names), which the shared macros read.
.names_write <- function(con) {
  .exec(con, "CREATE OR REPLACE TABLE op_names (kind VARCHAR, alias VARCHAR, name VARCHAR)")
  names <- .contract()$biolink_names
  values <- vapply(names, function(n) sprintf("(%s, %s, %s)", .quote(n[[1]]), .quote(n[[2]]), .quote(n[[3]])), "")
  .exec(con, paste0("INSERT INTO op_names VALUES ", paste(values, collapse = ", ")))
}

# The rule tables the macros read (display names, entity groups).
.query_tables_write <- function(con) {
  for (table in names(.contract()$query_tables)) {
    spec <- .contract()$query_tables[[table]]
    .exec(con, sprintf("CREATE OR REPLACE TABLE %s (%s)", table, paste(unlist(spec$columns), collapse = ", ")))
    if (length(spec$rows)) {
      values <- vapply(spec$rows, function(row) paste0("(", paste(vapply(row, .scalar_literal, ""), collapse = ", "), ")"), "")
      .exec(con, sprintf("INSERT INTO %s VALUES %s", table, paste(values, collapse = ", ")))
    }
  }
}

.query_layer_write <- function(con) {
  .names_write(con)
  .query_tables_write(con)
  for (name in unlist(.contract()$macro_files)) {
    sql <- paste(readLines(file.path(.sql_dir(), name), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    .exec(con, sql)
  }
  .exec(con, sprintf("UPDATE database_info SET query_layer = %s", .quote(.query_layer())))
}

# Where each resource table is written: the scratch database for those the build reads,
# the database itself for raw evidence, else nowhere (NULL).
.written_to <- function(table) {
  if (table %in% unlist(.contract()$work_tables)) return("work.")
  if (identical(table, .payloads())) return("")
  NULL
}

# Keys are SHA-256 hex: their first characters, then any other ("").
.slices <- c(strsplit("0123456789abcdef", "")[[1]], "")

# The derived tables (each sorted), in build order, and the indexes. The stages of the
# build ("scratch") go to the scratch database.
.derived <- function(con) {
  derived <- .contract()$derived_tables
  partition <- "\\{\\{partition ([A-Za-z0-9_.]+)\\}\\}"
  hex <- paste(sprintf("'%s'", .slices[nzchar(.slices)]), collapse = ", ")
  for (name in names(derived)) {
    target <- sprintf("%s.main.%s", if (isTRUE(derived[[name]]$scratch)) "work" else "omnipath", name)
    select <- derived[[name]]$select
    order <- derived[[name]]$order
    if (!grepl(partition, select)) {
      .exec(con, sprintf("CREATE TABLE %s AS SELECT * FROM (%s) ORDER BY %s", target, select, order))
      next
    }
    # In slices of keys by their first character, in key order: the table is sorted as
    # a whole, and no slice needs the memory of all of them.
    .exec(con, sprintf("CREATE TABLE %s AS SELECT * FROM (%s) LIMIT 0", target, gsub(partition, "false", select)))
    for (first in .slices) {
      condition <- if (nzchar(first)) sprintf("left(\\1, 1) = '%s'", first)
        else sprintf("left(coalesce(\\1, ''), 1) NOT IN (%s)", hex)
      .exec(con, sprintf("INSERT INTO %s SELECT * FROM (%s) ORDER BY %s", target, gsub(partition, condition, select), order))
    }
  }
  for (index in .contract()$indexes) {
    .exec(con, sprintf("CREATE INDEX %s ON omnipath.main.%s(%s)", index[[1]], index[[2]], index[[3]]))
  }
}

.sorted_records <- function(records) {
  records[order(vapply(records, function(r) r$resource_id, ""))]
}

.rebuild <- function(root, records, memory_limit, release = NULL, taxonomy = NULL) {
  # taxonomy: the version of the installed taxon names (references/taxonomy/<version>/)
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
  .exec(con, sprintf("ATTACH %s AS work", .quote(file.path(staging_dir, "work.duckdb"))))
  .exec(con, "USE omnipath")
  .exec(con, "SET search_path = 'omnipath.main,work.main'")
  .configure(con, root, memory_limit)
  ordered <- .sorted_records(records)
  index <- stats::setNames(as.list(seq_along(ordered) - 1L), vapply(ordered, function(r) r$resource_id, ""))
  present <- Filter(function(t) any(vapply(ordered, function(r) t %in% .record_tables(r), TRUE)), c(.tables(), .serving()))
  for (table in present) {
    where <- .written_to(table)
    if (is.null(where)) next
    having <- Filter(function(r) table %in% .record_tables(r), ordered)
    # The columns of every resource (union by name), then each resource's rows in its
    # file's order, so rows follow their global ids.
    union <- paste(vapply(having, function(r) .source(root, r, table, index[[r$resource_id]]), ""),
      collapse = " UNION ALL BY NAME ")
    .exec(con, sprintf("CREATE TABLE %s%s AS SELECT * FROM (%s) LIMIT 0", where, table, union))
    for (record in having) {
      .exec(con, sprintf("INSERT INTO %s%s BY NAME %s", where, table, .source(root, record, table, index[[record$resource_id]])))
    }
  }
  # The merged entities' summaries name taxa and use the shared macros and rules.
  .taxonomy(con, if (is.null(taxonomy)) NULL else .taxonomy_file(root, taxonomy))
  .names_write(con)
  .query_tables_write(con)
  .prepare_macros(con)
  .derived(con)
  .describe(con, root, records, release, taxonomy)
  .query_layer_write(con)
  .exec(con, "USE memory")
  .exec(con, "DETACH work")
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
    taxonomy <- if ("taxonomy" %in% columns) "taxonomy" else "NULL::VARCHAR AS taxonomy"
    info <- DBI::dbGetQuery(con, sprintf("SELECT query_layer, database_format, %s FROM database_info", taxonomy))
    state <- DBI::dbGetQuery(con, "SELECT * FROM resources ORDER BY ALL")
    list(query_layer = info$query_layer, database_format = info$database_format, state = state,
      taxonomy = if (is.na(info$taxonomy)) NULL else info$taxonomy)
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
  .tables_in_sync(info, installation) && identical(info$query_layer, .query_layer()) &&
    identical(info$taxonomy, installation$taxonomy)
}

.layer_only <- function(root) {
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = file.path(root, .database_file))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbBegin(con)
  .query_layer_write(con)
  DBI::dbCommit(con)
  invisible(TRUE)
}

# Changed resources or taxon names (the merged entities' summaries name taxa) mean a new
# database, written completely; when only the query layer differs (another client
# version wrote it), it alone.
.synchronize <- function(root, installation, memory_limit = "2GB") {
  info <- .info(root)
  if (.tables_in_sync(info, installation) && identical(info$taxonomy, installation$taxonomy)) {
    if (!identical(info$query_layer, .query_layer())) .layer_only(root)
    return(invisible(FALSE))
  }
  on.exit({
    temporary <- file.path(root, ".duckdb-tmp")
    if (dir.exists(temporary) && !length(list.files(temporary, all.files = TRUE, no.. = TRUE))) unlink(temporary, recursive = TRUE)
  }, add = TRUE)
  .rebuild(root, installation$resources, memory_limit, installation$release, installation$taxonomy)
  invisible(TRUE)
}
