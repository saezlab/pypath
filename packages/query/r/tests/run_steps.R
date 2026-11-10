# Run by packages/query/python/tests/test_steps.py: run step-file queries with the R
# runner and print their JSON outputs (or {"error": message}) as a JSON array of strings.
#
#   Rscript run_steps.R <package dir> <database> <queries dir> <cases.json>
#
# With "package" for <queries dir>, the package's own queries, fragments and macros
# (on a release database).

args <- commandArgs(trailingOnly = TRUE)
pkgload::load_all(args[[1]], quiet = TRUE, export_all = TRUE)
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = args[[2]], read_only = TRUE)
cases <- jsonlite::fromJSON(args[[4]], simplifyVector = FALSE)
package <- identical(args[[3]], "package")
answers <- vapply(cases, function(case) {
  params <- if (length(case[[2]])) case[[2]] else list()
  answer <- tryCatch(
    if (package) {
      .run_query(con, case[[1]], params, output = "json")
    } else {
      text <- paste(readLines(file.path(args[[3]], paste0(case[[1]], ".sql")), warn = FALSE), collapse = "\n")
      .run_query(con, .parse_query(text, case[[1]]), params, output = "json",
        fragments = function(name) .fragment(name, dir = file.path(args[[3]], "fragments")), macros = FALSE)
    },
    omnipath_query_error = function(e) as.character(jsonlite::toJSON(list(error = conditionMessage(e)), auto_unbox = TRUE))
  )
  if (is.null(answer)) "null" else answer
}, "")
cat(jsonlite::toJSON(answers))
DBI::dbDisconnect(con, shutdown = TRUE)
