# API mode: the hosted API answers op_search(), op_relations() and op_entity() (its
# /query endpoint runs the same SQL files as local mode), as
# the Python client's ApiClient does (omnipath_client/remote.py). Answers go through
# an in-memory DuckDB database, so that they have the local answers' column types.

.is_api <- function(con) inherits(con, "omnipath_api")

.api_connection <- function(api_url, release, timeout, memory_limit, threads) {
  if (!is.character(release) || length(release) != 1 || !nzchar(trimws(release))) {
    stop("release must not be empty", call. = FALSE)
  }
  db <- DBI::dbConnect(duckdb::duckdb())
  DBI::dbExecute(db, sprintf("SET memory_limit = %s", .quote(memory_limit)))
  if (!is.null(threads)) DBI::dbExecute(db, sprintf("SET threads = %d", as.integer(threads)))
  state <- new.env(parent = emptyenv())
  state$db <- db
  state$answers <- 0L
  structure(
    list(api_url = .api_url(api_url), release = release, timeout = timeout, state = state),
    class = "omnipath_api"
  )
}

#' @export
print.omnipath_api <- function(x, ...) {
  cat("<omnipath API ", x$api_url, ", release ", x$release, ">\n", sep = "")
  invisible(x)
}

# HTTP, in one function the tests can replace.
.api_request <- function(con, method, path, query = list(), body = NULL) {
  request <- httr2::request(paste0(con$api_url, path))
  request <- httr2::req_url_query(request, !!!c(query, list(release = con$release)))
  request <- httr2::req_method(request, method)
  if (!is.null(body)) request <- httr2::req_body_json(request, body, auto_unbox = TRUE, null = "null", digits = NA)
  request <- httr2::req_options(request, connecttimeout = con$timeout)
  request <- httr2::req_error(request, is_error = function(response) FALSE)
  response <- httr2::req_perform(httr2::req_user_agent(request, "omnipath.client (R)"))
  status <- httr2::resp_status(response)
  if (status >= 400 && !(status == 404 && startsWith(path, "/entities/"))) {
    detail <- tryCatch(httr2::resp_body_json(response)$detail, error = function(e) NULL)
    detail <- if (is.null(detail)) substr(httr2::resp_body_string(response), 1, 300) else as.character(detail)
    stop(sprintf("%s %s: HTTP %d: %s", method, path, status, detail), call. = FALSE)
  }
  response
}

# The parameters of a shared query (sql/queries/<name>.sql): the filters it takes in
# both modes.
.shared_parameters <- function(query) {
  if (!file.exists(file.path(.sql_dir(), "queries", paste0(query, ".sql")))) stop("No shared query ", query, call. = FALSE)
  setdiff(names(.load_query(query)$params), "resources")
}

.api_filters <- function(macro, internal, filters) {
  available <- sort(setdiff(.shared_parameters(macro), internal))
  unknown <- setdiff(names(filters), available)
  if (length(unknown)) {
    stop(sprintf("Unknown filters: %s; available: %s", .format_names(sort(unknown)), .format_names(available)), call. = FALSE)
  }
  out <- list()
  for (name in names(filters)) {
    value <- filters[[name]]
    if (is.null(value)) next
    out[[name]] <- if (name %in% .scalars) value else as.list(unlist(value))
  }
  out
}

.api_table <- function(con) {
  con$state$answers <- con$state$answers + 1L
  sprintf("answer_%d", con$state$answers)
}

# A query by name, answered by the API from the same SQL the local client runs: a data
# frame (from the API's Parquet) or JSON text (NULL for no answer).
.api_query <- function(con, name, params, output = "table") {
  body <- if (length(params)) params else structure(list(), names = character())
  format <- if (output == "json") "json" else "parquet"
  response <- .api_request(con, "POST", paste0("/query/", utils::URLencode(name, reserved = TRUE)),
    query = list(format = format), body = body)
  if (output == "json") {
    text <- httr2::resp_body_string(response)
    return(if (identical(text, "null")) NULL else text)
  }
  path <- tempfile("omnipath-", fileext = ".parquet")
  on.exit(unlink(path), add = TRUE)
  writeBin(httr2::resp_body_raw(response), path)
  table <- .api_table(con)
  DBI::dbExecute(con$state$db, sprintf("CREATE TABLE %s AS SELECT * FROM read_parquet(%s)", table, .quote(path)))
  DBI::dbGetQuery(con$state$db, sprintf("SELECT * FROM %s", table))
}

.api_search <- function(con, query, limit, filters) {
  filters <- .api_filters("entity_search", c("query", "limit", "cursor_phase", "cursor_offset", "cursor_key"), filters)
  .api_query(con, "entity_search", c(list(query = query, limit = as.integer(limit)), filters))
}

.api_relations <- function(con, entities, limit, filters) {
  filters <- .api_filters("relations", c("entities", "max_relations"), filters)
  if (!is.null(entities)) filters$entities <- as.list(unlist(entities))
  .api_query(con, "relations", c(filters, list(max_relations = as.integer(limit))))
}

.api_entity <- function(con, entity_key, simplify) {
  text <- .api_query(con, "entity_details", list(entity_key = as.character(entity_key)), output = "json")
  if (is.null(text)) return(NULL)
  jsonlite::fromJSON(text, simplifyVector = simplify)
}
