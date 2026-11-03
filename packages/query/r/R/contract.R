# What this package shares with the Python client: inst/sql/contract.json (table
# names, the serving tables' queries, lookups, Biolink names) and the macro files,
# copied from packages/query/sql/ by
# packages/query/scripts/generate_contract.py. Do not edit them here.

.state <- new.env(parent = emptyenv())

.sql_dir <- function() {
  dir <- system.file("sql", package = "omnipath.client")
  if (!nzchar(dir)) stop("omnipath.client: the package's sql/ directory is missing")
  dir
}

.contract <- function() {
  if (is.null(.state$contract)) {
    .state$contract <- jsonlite::fromJSON(
      file.path(.sql_dir(), "contract.json"),
      simplifyVector = FALSE
    )
  }
  .state$contract
}

.tables <- function() unlist(.contract()$published_tables)
.serving <- function() unlist(.contract()$serving_tables)
.payloads <- function() .contract()$evidence_table

# The query layer's fingerprint: the first 16 hex digits of the SHA-256 of the
# contract and the macro files, as the Python client computes it, so that neither
# client writes the other's database again.
.query_layer <- function() {
  if (is.null(.state$query_layer)) {
    names <- c("contract.json", unlist(.contract()$macro_files))
    bytes <- unlist(lapply(file.path(.sql_dir(), names), function(path) {
      readBin(path, "raw", file.info(path)$size)
    }))
    con <- DBI::dbConnect(duckdb::duckdb())
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
    hash <- DBI::dbGetQuery(con, "SELECT sha256(?) AS h", params = list(rawToChar(bytes)))$h
    .state$query_layer <- substr(hash, 1, 16)
  }
  .state$query_layer
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# SQL literals ------------------------------------------------------------------

.quote <- function(x) paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
.name <- function(x) paste0('"', gsub('"', '""', x, fixed = TRUE), '"')

# One value as a SQL literal: NULL, a boolean, a number or a string.
.scalar_literal <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return("NULL")
  if (is.logical(x)) return(if (x) "true" else "false")
  if (is.numeric(x)) {
    if (x == round(x)) return(format(x, scientific = FALSE, trim = TRUE))
    return(format(x, digits = 17))
  }
  if (is.character(x) || is.factor(x)) return(.quote(as.character(x)))
  stop("Unsupported filter value of class ", class(x)[1])
}

# A filter value as a SQL literal; as_list makes one value a list of one.
.literal <- function(x, as_list = FALSE) {
  if (is.null(x)) return("NULL")
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  if (as_list) {
    return(paste0("[", paste(vapply(x, .scalar_literal, ""), collapse = ", "), "]"))
  }
  if (length(x) != 1) stop("Expected one value, got ", length(x))
  .scalar_literal(x)
}
