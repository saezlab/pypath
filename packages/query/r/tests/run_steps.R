# Run by packages/query/python/tests/test_steps.py: run step-file queries with the R
# runner and print their JSON outputs (or {"error": message}) as a JSON array of strings.
#
#   Rscript run_steps.R <package dir> <database> <queries dir> <cases.json>

args <- commandArgs(trailingOnly = TRUE)
pkgload::load_all(args[[1]], quiet = TRUE, export_all = TRUE)
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = args[[2]], read_only = TRUE)
cases <- jsonlite::fromJSON(args[[4]], simplifyVector = FALSE)
answers <- vapply(cases, function(case) {
  text <- paste(readLines(file.path(args[[3]], paste0(case[[1]], ".sql")), warn = FALSE), collapse = "\n")
  query <- .parse_query(text, case[[1]])
  params <- if (length(case[[2]])) case[[2]] else list()
  tryCatch(
    .run_query(con, query, params, output = "json"),
    omnipath_query_error = function(e) as.character(jsonlite::toJSON(list(error = conditionMessage(e)), auto_unbox = TRUE))
  )
}, "")
cat(jsonlite::toJSON(answers))
DBI::dbDisconnect(con, shutdown = TRUE)
