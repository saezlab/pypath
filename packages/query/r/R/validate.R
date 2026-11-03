# Catalog and installation records are validated once, keeping additive server
# metadata, as in the Python client (omnipath_client/contracts.py).

.fail <- function(...) stop(structure(
  class = c("omnipath_error", "error", "condition"),
  list(message = paste0(...), call = NULL)
))

.is_int <- function(x) {
  is.numeric(x) && length(x) == 1 && !is.na(x) && x == round(x)
}

.validate_url <- function(value) {
  ok <- is.character(value) && length(value) == 1 &&
    grepl("^https?://[^/@[:space:]]+(/|$)", value, ignore.case = TRUE)
  if (!ok) .fail("Expected a public HTTP(S) URL without embedded credentials")
  value
}

.validate_resource <- function(value) {
  if (!is.list(value) || is.null(names(value))) .fail("Expected a resource object")
  source <- value$resource_id
  version <- value$version
  if (!is.character(source) || length(source) != 1 || !grepl("^[a-z0-9][a-z0-9_-]*$", source)) {
    .fail("Invalid resource ID in catalog")
  }
  if (!is.character(version) || length(version) != 1 || !grepl("^[0-9]+(\\.[0-9]+)*$", version)) {
    .fail("Invalid resource version in catalog")
  }
  artifacts <- value$files
  if (!is.list(artifacts) || (length(artifacts) && !is.null(names(artifacts)))) {
    .fail("Resource files must be a list")
  }
  allowed <- paste0(.tables(), ".parquet")
  seen <- character()
  for (artifact in artifacts) {
    if (!is.list(artifact) || is.null(names(artifact))) .fail("Expected an artifact object")
    name <- artifact$name
    if (!is.character(name) || length(name) != 1 || !(name %in% allowed) || name %in% seen) {
      .fail("Unknown or duplicate artifact in catalog")
    }
    if (!.is_int(artifact$size_bytes) || artifact$size_bytes < 8) {
      .fail("Invalid artifact size in catalog")
    }
    if (!is.null(artifact$url)) .validate_url(artifact$url)
    if (!is.null(artifact$rows) && (!.is_int(artifact$rows) || artifact$rows < 0)) {
      .fail("Invalid artifact row count")
    }
    seen <- c(seen, name)
  }
  if (!all(c("entity.parquet", "relation.parquet") %in% seen)) {
    .fail("Resource must contain entity.parquet and relation.parquet")
  }
  value
}

.validate_catalog <- function(value) {
  if (!is.list(value) || (length(value) && !is.null(names(value)))) .fail("Resources must be a list")
  catalog <- list()
  for (raw in value) {
    record <- .validate_resource(raw)
    if (!is.null(catalog[[record$resource_id]])) {
      .fail("Duplicate resource in catalog: ", record$resource_id)
    }
    catalog[[record$resource_id]] <- record
  }
  catalog
}

.validate_installation <- function(value) {
  if (!is.list(value) || is.null(names(value))) .fail("Expected an installation object")
  if (!.is_int(value$format) || value$format != 1) .fail("Unsupported installation format")
  .validate_url(value$api_url)
  if (!is.character(value$release) || length(value$release) != 1 || !nzchar(trimws(value$release))) {
    .fail("Installation release must be a non-empty string")
  }
  materialize <- value$materialize %||% FALSE
  if (!is.logical(materialize) || length(materialize) != 1 || is.na(materialize)) {
    .fail("Installation materialize must be true or false")
  }
  catalog <- .validate_catalog(value$resources)
  list(
    format = 1L,
    api_url = value$api_url,
    release = value$release,
    materialize = materialize,
    resources = unname(catalog)
  )
}
