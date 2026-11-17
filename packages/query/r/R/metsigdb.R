# MetSigDB: the metabolite memberships of resource-native sets (KEGG, Reactome, WikiPathways,
# MACdb, ClassyFire). The subset is built from a whole release and served by the API, so
# only API mode answers; the Python client's ApiClient$metsigdb() is the counterpart.

.metsigdb_set_filters <- c(
  "resource", "set_type", "set_sub_type", "organism", "set", "limit", "offset", "total"
)
.metsigdb_filters <- c(
  .metsigdb_set_filters, "entity", "metabolite_structure_key", "exclude_metabolites",
  "fields", "group"
)

# The filters as the API's query parameters: a list is one comma-separated value, a
# logical `true` or `false`. The API validates the names and the values.
.metsigdb_query <- function(filters, available) {
  unknown <- setdiff(names(filters), available)
  if (length(unknown) || (length(filters) && is.null(names(filters)))) {
    stop(sprintf("Unknown filters: %s; available: %s",
      .format_names(sort(unknown)), .format_names(sort(available))), call. = FALSE)
  }
  query <- list()
  for (name in names(filters)) {
    value <- filters[[name]]
    if (is.null(value)) next
    query[[name]] <- if (is.logical(value)) {
      tolower(as.character(value))
    } else {
      paste(as.character(unlist(value)), collapse = ",")
    }
  }
  query
}

.metsigdb_get <- function(con, path, available, filters, simplify) {
  if (!.is_api(con)) {
    stop("MetSigDB memberships are answered by the API: use op_connect(mode = \"api\") ",
      "(a local installation holds resources, not the release's subsets)", call. = FALSE)
  }
  response <- .api_request(con, "GET", path, query = .metsigdb_query(filters, available))
  jsonlite::fromJSON(httr2::resp_body_string(response), simplifyVector = simplify)
}

#' Metabolite memberships of sets (MetSigDB)
#'
#' One page of metabolites in sets of KEGG, Reactome, WikiPathways, MACdb and ClassyFire,
#' from the hosted API (`/sets/metsigdb`); the API checks the filters and refuses a name
#' or value it does not publish.
#'
#' @param con An API handle from [op_connect()].
#' @param ... Filters: `resource`, `set_type`, `set_sub_type`, `organism`, `set` (the
#'   set's id), `entity` (an entity key, InChIKey, or HMDB, ChEBI, KEGG or PubChem
#'   identifier), `metabolite_structure_key`, `exclude_metabolites` (`"all"` or classes
#'   such as `"cofactors"`), `fields`, `group`, `limit`, `offset` and `total`. A value may
#'   be a vector.
#' @param simplify `TRUE` returns data frames for lists of records; `FALSE` keeps nested
#'   lists, as JSON has them.
#' @return A list: `count`, `has_more`, `total` and `rows` (or, with `group = TRUE`,
#'   `groups` of sets with their members, counted in sets).
#' @export
op_metsigdb <- function(con, ..., simplify = TRUE) {
  .metsigdb_get(con, "/sets/metsigdb", .metsigdb_filters, list(...), simplify)
}

#' The sets of MetSigDB
#'
#' The sets without their members: resource, id, label, type, sub-type, organism and size.
#'
#' @inheritParams op_metsigdb
#' @param ... Filters: `resource`, `set_type`, `set_sub_type`, `organism`, `set`,
#'   `limit`, `offset` and `total`.
#' @return A list: `count`, `has_more`, `total` and `sets`.
#' @export
op_metsigdb_sets <- function(con, ..., simplify = TRUE) {
  .metsigdb_get(con, "/sets/metsigdb/sets", .metsigdb_set_filters, list(...), simplify)
}
