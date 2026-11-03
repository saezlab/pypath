# Run by packages/query/python/tests/test_r_client.py: install the fixtures with this
# package into `r_dir`, then answer the same questions on `r_dir` and on `python_dir`
# (written by the Python client) and print them as JSON.
#
#   Rscript cross_language.R <package dir> <r_dir> <python_dir> <questions.json>

args <- commandArgs(trailingOnly = TRUE)
package <- args[[1]]
r_dir <- args[[2]]
python_dir <- args[[3]]
questions <- jsonlite::fromJSON(args[[4]], simplifyVector = FALSE)

pkgload::load_all(package, quiet = TRUE, export_all = FALSE)
fixtures <- file.path(package, "tests", "testthat", "fixtures")

install <- function(materialize) {
  testthat::local_mocked_bindings(
    .fetch_json = function(url, query = list(), timeout = 30) {
      path <- sub("^https?://[^/]+", "", url)
      if (grepl("/resources$", path)) return(jsonlite::fromJSON(file.path(fixtures, "catalog-1.json"), simplifyVector = FALSE))
      name <- sub("^.*/resources/([^/]+)/files$", "\\1", path)
      jsonlite::fromJSON(file.path(fixtures, sprintf("files-%s-1.json", name)), simplifyVector = FALSE)
    },
    .fetch_file = function(url, path, timeout = 30) {
      file.copy(file.path(fixtures, "server", sub("^https://data.example/resources/", "", url)), path, overwrite = TRUE)
      invisible(path)
    },
    .package = "omnipath.client"
  )
  op_init(c("alpha", "beta", "gamma"), data_dir = r_dir, include_evidence = TRUE,
    materialize = materialize, api_url = "https://api.example/api", progress = FALSE)
}
install(questions$materialize)

answers <- function(dir) {
  con <- op_connect(dir)
  on.exit(op_disconnect(con))
  list(
    search = lapply(questions$search, function(q) {
      rows <- do.call(op_search, c(list(con, q$query), q$filters))
      jsonlite::toJSON(rows, dataframe = "rows", na = "null", null = "null", digits = NA)
    }),
    relations = lapply(questions$relations, function(q) {
      rows <- do.call(op_relations, c(list(con, q$entities), q$filters))
      jsonlite::toJSON(rows, dataframe = "rows", na = "null", null = "null", digits = NA)
    }),
    entity = lapply(questions$entity, function(key) {
      row <- op_sql(con, sprintf("SELECT details FROM op_entity('%s')", key))
      if (nrow(row)) row$details[[1]] else "null"
    }),
    query_layer = op_sql(con, "SELECT query_layer FROM database_info")$query_layer
  )
}

cat(jsonlite::toJSON(list(r = answers(r_dir), python = answers(python_dir)), auto_unbox = TRUE, digits = NA))
