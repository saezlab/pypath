# A catalog API and a static file host over the fixtures (see tests/make_fixtures.py),
# in place of the package's two HTTP functions. `service$downloads` lists the files
# fetched, as "alpha/1.0/entity.parquet".

fixture <- function(...) testthat::test_path("fixtures", ...)

local_service <- function(env = parent.frame()) {
  service <- new.env()
  service$catalog <- 1
  service$downloads <- character()
  service$truncated <- character()
  service$requests <- list()
  service$taxonomy <- NULL  # the release's taxon names (a Parquet file), when it has them
  testthat::local_mocked_bindings(
    .fetch_json = function(url, query = list(), timeout = 30) {
      service$requests[[length(service$requests) + 1]] <- list(url = url, query = query)
      path <- sub("^https?://[^/]+", "", url)
      if (grepl("/resources$", path)) {
        catalog <- jsonlite::fromJSON(fixture(sprintf("catalog-%d.json", service$catalog)), simplifyVector = FALSE)
        if (!is.null(service$taxonomy)) {
          catalog$references <- list(taxonomy = list(
            version = file_sha256(service$taxonomy), size_bytes = as.integer(file.size(service$taxonomy))
          ))
        }
        return(catalog)
      }
      name <- sub("^.*/resources/([^/]+)/files$", "\\1", path)
      jsonlite::fromJSON(fixture(sprintf("files-%s-%d.json", name, service$catalog)), simplifyVector = FALSE)
    },
    .fetch_file = function(url, path, timeout = 30) {
      if (grepl("/references/taxonomy", url, fixed = TRUE)) {
        service$downloads <- c(service$downloads, "taxonomy")
        file.copy(service$taxonomy, path, overwrite = TRUE)
        return(invisible(path))
      }
      relative <- sub("^https://data.example/resources/", "", url)
      service$downloads <- c(service$downloads, relative)
      file.copy(fixture("server", relative), path, overwrite = TRUE)
      if (url %in% service$truncated) {
        bytes <- readBin(path, "raw", file.size(path))
        writeBin(bytes[1:10], path)
      }
      invisible(path)
    },
    .package = "omnipath.client",
    .env = env
  )
  service$data <- file.path(withr::local_tempdir(.local_envir = env), "data")
  service$init <- function(...) {
    op_init(..., data_dir = service$data, api_url = "https://api.example/api", progress = FALSE)
  }
  service$connect <- function() op_connect(service$data)
  service
}

local_connection <- function(service, env = parent.frame()) {
  con <- service$connect()
  withr::defer(op_disconnect(con), envir = env)
  con
}

installed_names <- function(service) {
  vapply(op_resources(service$data), function(r) r$resource_id, "")
}

file_sha256 <- function(path) {
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  DBI::dbGetQuery(con, sprintf("SELECT sha256(content) AS h FROM read_blob('%s')", path))$h
}

write_taxonomy <- function(path, name) {
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  DBI::dbExecute(con, sprintf(paste0(
    "COPY (SELECT * FROM (VALUES ('9606', '9606', 'Homo sapiens', '%s'), ",
    "('10090', '10090', 'Mus musculus', 'house mouse')) ",
    "t(taxon_id, current_taxon_id, scientific_name, common_name)) TO '%s' (FORMAT PARQUET)"), name, path))
  path
}
