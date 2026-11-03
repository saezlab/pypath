#' Connect to the installed database
#'
#' Opens the data directory's DuckDB database read-only and returns a DBI
#' connection: use it with [op_sql()], the prebuilt queries ([op_search()],
#' [op_relations()], [op_entity()]), DBI (`DBI::dbGetQuery()`) or dplyr
#' (`dplyr::tbl(con, "relation")`). Checks that every installed file is complete,
#' and writes the database again when it does not describe the installation (a moved
#' directory, an interrupted [op_init()], another client version).
#'
#' @param data_dir The data directory; see [op_data_directory()].
#' @param memory_limit DuckDB's memory limit.
#' @param threads DuckDB's threads, or `NULL` for its default.
#' @return A DBI connection; close it with [op_disconnect()].
#' @export
op_connect <- function(data_dir = NULL, memory_limit = "2GB", threads = NULL) {
  root <- op_data_directory(data_dir)
  installation <- .read_installation(root)
  if (is.null(installation) || !length(installation$resources)) {
    .fail("No OmniPath data in ", root, "; run omnipath.client::op_init() first")
  }
  for (record in installation$resources) {
    for (name in c(vapply(record$files, function(f) f$name, ""), paste0(.serving(), ".parquet"))) {
      .check_file(root, record, name)
    }
  }
  if (!.in_sync(root, installation)) {
    tryCatch(
      .synchronize(root, installation, memory_limit),
      error = function(e) .fail("The database in ", root, " is out of date; run op_init() again (", conditionMessage(e), ")")
    )
  }
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = file.path(root, .database_file), read_only = TRUE)
  DBI::dbExecute(con, sprintf("SET memory_limit = %s", .quote(memory_limit)))
  if (!is.null(threads)) DBI::dbExecute(con, sprintf("SET threads = %d", as.integer(threads)))
  for (name in names(.contract()$settings)) {
    DBI::dbExecute(con, sprintf("SET %s = %d", name, as.integer(.contract()$settings[[name]])))
  }
  con
}

.check_file <- function(root, record, name) {
  path <- file.path(.resource_dir(root, record), name)
  if (file.exists(path)) {
    real <- normalizePath(path)
    if (!startsWith(real, paste0(normalizePath(root), .Platform$file.sep))) {
      .fail("Installed file escapes the data directory")
    }
  }
  size <- NULL
  for (f in record$files) if (identical(f$name, name)) size <- f$size_bytes
  if (!file.exists(path) || (!is.null(size) && file.size(path) != size)) {
    .fail("Missing or truncated file: ", path, "; run op_init() again")
  }
}

#' Close a connection
#'
#' @param con A connection from [op_connect()].
#' @export
op_disconnect <- function(con) {
  DBI::dbDisconnect(con, shutdown = TRUE)
  invisible(NULL)
}

#' The installed resources
#'
#' @param data_dir The data directory; see [op_data_directory()].
#' @return The installed resources' catalog records.
#' @export
op_resources <- function(data_dir = NULL) {
  installation <- .read_installation(op_data_directory(data_dir))
  if (is.null(installation)) list() else installation$resources
}

#' Run SQL on the database
#'
#' `entity`, `relation` and the other tables over every installed resource (with
#' `_resource` and `_resource_version`), one resource's tables as `signor.relation`,
#' and the macros `op_search`, `op_relations` and `op_entity`. Caller-supplied local
#' SQL; not a sandbox.
#'
#' @param con A connection from [op_connect()].
#' @param query SQL.
#' @param params Values for `?` placeholders, or `NULL`.
#' @return A data frame.
#' @export
op_sql <- function(con, query, params = NULL) {
  if (is.null(params)) DBI::dbGetQuery(con, query) else DBI::dbGetQuery(con, query, params = unname(params))
}

# Prebuilt queries ---------------------------------------------------------------

# Macro parameters that take one value; every other filter takes a list, and one
# value stands for a list of one.
.scalars <- c("gene_mode", "molecular_endpoint_mode", "is_directed", "sign", "min_evidence_count", "max_results", "max_relations")
.biolink <- c(predicates = "predicate", interaction_types = "class", entity_types = "class", partner_types = "class")

.parameters <- function(con, macro, internal) {
  row <- DBI::dbGetQuery(con, sprintf(
    "SELECT parameters FROM duckdb_functions() WHERE function_name = %s LIMIT 1", .quote(macro)
  ))
  sort(setdiff(unlist(row$parameters), internal))
}

.filters <- function(con, macro, internal, filters) {
  available <- .parameters(con, macro, internal)
  unknown <- setdiff(names(filters), available)
  if (length(unknown)) {
    stop(sprintf("Unknown filters: %s; available: %s", .format_names(sort(unknown)), .format_names(available)), call. = FALSE)
  }
  sources <- unlist(filters$sources)
  installed <- DBI::dbGetQuery(con, "SELECT resource_id FROM resources ORDER BY 1")$resource_id
  missing <- setdiff(sources, installed)
  if (length(missing)) {
    stop(sprintf("Resources unavailable here: %s; installed: %s", .format_names(sort(missing)), .format_names(installed)), call. = FALSE)
  }
  # Biolink names now, not when the rows are read: the macros check them with the
  # same op_biolink, which raises an error naming an unknown one.
  for (name in intersect(names(.biolink), names(filters))) {
    if (length(unlist(filters[[name]]))) {
      tryCatch(
        DBI::dbGetQuery(con, sprintf(
          "SELECT list(name) AS names FROM op_biolink(%s, %s)", .quote(.biolink[[name]]), .literal(filters[[name]], TRUE)
        )),
        error = function(e) stop(sub("^.*?(Unknown Biolink [^\n]*).*$", "\\1", conditionMessage(e)), call. = FALSE)
      )
    }
  }
  filters
}

.arguments <- function(filters) {
  vapply(names(filters), function(name) {
    sprintf("%s := %s", name, .literal(filters[[name]], as_list = !(name %in% .scalars)))
  }, "")
}

.call <- function(con, macro, positional, filters, lazy) {
  sql <- sprintf("SELECT * FROM %s(%s)", macro, paste(c(positional, .arguments(filters)), collapse = ", "))
  if (isTRUE(lazy)) {
    if (!requireNamespace("dplyr", quietly = TRUE) || !requireNamespace("dbplyr", quietly = TRUE)) {
      stop("lazy = TRUE needs the packages dplyr and dbplyr", call. = FALSE)
    }
    return(dplyr::tbl(con, dplyr::sql(sql)))
  }
  DBI::dbGetQuery(con, sql)
}

#' Search entities as the explorer does
#'
#' Entities matching a name or identifier, merged across resources and ranked as in
#' the explorer's search. The macro `op_search`.
#'
#' @param con A connection from [op_connect()].
#' @param query A name, identifier, `ns:id` CURIE or entity key.
#' @param limit The number of entities.
#' @param ... The explorer's filters: `sources`, `entity_types`, `taxonomy_ids`,
#'   `reference_entity_keys`, `entity_ids`.
#' @param lazy Return a dplyr table instead of a data frame.
#' @return A data frame with the entity table's columns.
#' @export
op_search <- function(con, query, limit = 20, ..., lazy = FALSE) {
  filters <- .filters(con, "op_search", c("query", "max_results"), list(...))
  filters$max_results <- as.integer(limit)
  .call(con, "op_search", .literal(query), filters, lazy)
}

#' Relations as the explorer exports them
#'
#' One row per relation with its annotations, evidence and the product records its
#' evidence names; at most 100,000 relations (larger selections are SQL on
#' `relation`). The macro `op_relations`.
#'
#' @param con A connection from [op_connect()].
#' @param entities Exact names, identifiers, entity keys or group keys such as
#'   `gene:entrez:7157`, at either endpoint, or `NULL`.
#' @param limit The most relations to return.
#' @param ... The explorer's filters, such as `sources`, `predicates`,
#'   `taxonomy_ids`, `signs`, `reference_entity_keys`, `protein_entity_keys`,
#'   `isoform_identifiers`, `gene_mode`.
#' @param lazy Return a dplyr table instead of a data frame.
#' @return A data frame; annotations, evidence and product records are list columns.
#' @export
op_relations <- function(con, entities = NULL, limit = 100000, ..., lazy = FALSE) {
  most <- .contract()$relation_export_limit
  if (limit < 1 || limit > most) stop("limit must be between 1 and ", format(most, big.mark = ","), call. = FALSE)
  filters <- .filters(con, "op_relations", c("entities", "max_relations"), list(...))
  if (!is.null(entities)) filters$entities <- entities
  filters$max_relations <- as.integer(limit)
  .call(con, "op_relations", character(), filters, lazy)
}

#' One entity's details as the explorer shows them
#'
#' The entity by its exact key (as [op_search()] returns it), with its identifiers,
#' annotations, structural relationships, relation count and associations: the
#' record the API and the Python client return. The macro `op_entity`.
#'
#' @param con A connection from [op_connect()].
#' @param entity_key An entity key.
#' @param simplify Simplify the record with jsonlite (data frames for lists of
#'   records); `FALSE` keeps nested lists, as JSON has them.
#' @return A list, or `NULL` when no installed resource has the entity.
#' @export
op_entity <- function(con, entity_key, simplify = TRUE) {
  row <- DBI::dbGetQuery(con, sprintf("SELECT details FROM op_entity(%s)", .literal(as.character(entity_key))))
  if (!nrow(row)) return(NULL)
  jsonlite::fromJSON(row$details[[1]], simplifyVector = simplify)
}
