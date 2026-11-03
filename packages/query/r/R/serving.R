# A resource version's serving tables, derived from its published tables with the
# build's queries (contract.json), as the Python client does (serving.py).

.write_serving_tables <- function(con, directory, work, tables = .serving()) {
  contract <- .contract()
  sources <- vapply(unlist(contract$serving_sources), function(name) {
    sprintf("read_parquet(%s)", .quote(file.path(directory, paste0(name, ".parquet"))))
  }, "")
  names(sources) <- unlist(contract$serving_sources)
  options <- sprintf(
    "(FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE %d)",
    as.integer(contract$parquet_row_group_size)
  )
  for (name in tables) {
    query <- contract$serving_queries[[name]]
    select <- query$select
    for (source in names(sources)) {
      select <- gsub(paste0("{", source, "}"), sources[[source]], select, fixed = TRUE)
    }
    staging <- file.path(work, paste0(name, ".parquet"))
    .exec(con, sprintf(
      "COPY (SELECT * FROM (%s) ORDER BY %s) TO %s %s",
      select, query$order, .quote(staging), options
    ))
    if (!file.rename(staging, file.path(directory, paste0(name, ".parquet")))) {
      .fail("Could not write ", name, ".parquet in ", directory)
    }
  }
  invisible(NULL)
}
