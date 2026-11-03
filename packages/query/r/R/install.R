# Choosing resources, downloading them and keeping the installation's manifest, as
# the Python client does (omnipath_client/installation.py). Both clients read and
# write the same data directory.

.default_api_url <- "https://omnipath-metabo-dev.schaul.click/api"
.manifest <- "installed.json"
.evidence_file <- function() paste0(.payloads(), ".parquet")
.dimensions <- c("topic", "content", "molecule", "knowledge_origin")

.api_url <- function(api_url) {
  url <- api_url %||% Sys.getenv("OMNIPATH_API_URL", "")
  if (!nzchar(url)) url <- .default_api_url
  sub("/+$", "", .validate_url(url))
}

# HTTP, in two functions the tests replace -------------------------------------

.fetch_json <- function(url, query = list(), timeout = 30) {
  request <- httr2::req_url_query(httr2::request(url), !!!query)
  request <- httr2::req_options(request, connecttimeout = timeout, low_speed_time = timeout, low_speed_limit = 1)
  response <- httr2::req_perform(httr2::req_user_agent(request, "omnipath.client (R)"))
  httr2::resp_body_json(response, simplifyVector = FALSE)
}

# A download fails when no byte arrives for `timeout` seconds, not after a total
# time: resources are up to gigabytes.
.fetch_file <- function(url, path, timeout = 30) {
  request <- httr2::req_headers(httr2::request(url), `Accept-Encoding` = "identity")
  request <- httr2::req_options(request, connecttimeout = timeout, low_speed_time = timeout, low_speed_limit = 1)
  httr2::req_perform(httr2::req_user_agent(request, "omnipath.client (R)"), path = path)
  invisible(path)
}

# The manifest -----------------------------------------------------------------

.read_installation <- function(root) {
  path <- file.path(root, .manifest)
  if (!file.exists(path)) return(NULL)
  value <- tryCatch(
    jsonlite::fromJSON(path, simplifyVector = FALSE),
    error = function(e) .fail("Invalid ", .manifest, " in ", root)
  )
  .validate_installation(value)
}

.write_installation <- function(root, installation) {
  staging <- tempfile(".installed-", tmpdir = root)
  writeLines(
    jsonlite::toJSON(installation, auto_unbox = TRUE, null = "null", digits = NA, pretty = TRUE),
    staging
  )
  if (!file.rename(staging, file.path(root, .manifest))) .fail("Could not write ", .manifest)
}

# The catalog ------------------------------------------------------------------

.catalog <- function(api_url, release, timeout) {
  response <- .fetch_json(paste0(api_url, "/resources"), list(shape = "svelte", release = release), timeout)
  if (!is.list(response) || is.null(names(response))) .fail("Expected a catalog object")
  .validate_catalog(response$resources)
}

.record_tags <- function(record, dimension) {
  tags <- Filter(function(t) is.list(t) && identical(t$dimension, dimension) && !is.null(t$id), record$tags %||% list())
  vapply(tags, function(t) t$id, "")
}

.select <- function(catalog, resources = NULL, tags = list()) {
  chosen <- catalog
  if (!is.null(resources)) {
    names <- unlist(resources)
    if (!length(names)) stop("Select at least one resource", call. = FALSE)
    unknown <- setdiff(names, names(catalog))
    if (length(unknown)) {
      stop("Unknown resources in selected release: ", .format_names(sort(unknown)), call. = FALSE)
    }
    chosen <- chosen[names(chosen) %in% names]
  }
  for (dimension in names(tags)) {
    values <- tags[[dimension]]
    if (is.null(values)) next
    wanted <- unique(unlist(values))
    known <- unique(unlist(lapply(catalog, .record_tags, dimension)))
    if (!length(wanted) || length(setdiff(wanted, known))) {
      stop(sprintf(
        "Unknown %s: %s; available: %s", dimension,
        .format_names(sort(setdiff(wanted, known))), .format_names(sort(known))
      ), call. = FALSE)
    }
    chosen <- Filter(function(r) length(intersect(.record_tags(r, dimension), wanted)) > 0, chosen)
  }
  if (!length(chosen)) stop("No resource matches the selection", call. = FALSE)
  chosen
}

.format_names <- function(x) paste0("[", paste0("'", x, "'", collapse = ", "), "]")

#' The resources of a release
#'
#' The resources a release offers, with their versions, tags, files and sizes, as
#' the catalog lists them.
#'
#' @param release A release, or "latest".
#' @param api_url The API; by default `OMNIPATH_API_URL` or the development deployment.
#' @param timeout Seconds without a response before a request fails.
#' @return A list of resource records.
#' @export
op_available <- function(release = "latest", api_url = NULL, timeout = 30) {
  unname(.catalog(.api_url(api_url), release, timeout))
}

# Downloads --------------------------------------------------------------------

.fetch_artifact <- function(artifact, path, timeout) {
  if (is.null(artifact$url)) .fail("Server does not provide a static download URL")
  .fetch_file(.validate_url(artifact$url), path, timeout)
  size <- if (file.exists(path)) file.size(path) else -1
  if (size != artifact$size_bytes || size < 8) .fail("Truncated download: ", artifact$name)
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  head <- readBin(con, "raw", 4)
  seek(con, size - 4)
  tail <- readBin(con, "raw", 4)
  magic <- charToRaw("PAR1")
  if (!identical(head, magic) || !identical(tail, magic)) .fail("Invalid Parquet download: ", artifact$name)
  invisible(path)
}

.present <- function(folder, artifacts) {
  all(vapply(artifacts, function(a) {
    path <- file.path(folder, a$name)
    file.exists(path) && file.size(path) == a$size_bytes
  }, TRUE))
}

.evidence_artifact <- function(api_url, release, record, timeout) {
  name <- record$resource_id
  listed <- .validate_resource(.fetch_json(
    paste0(api_url, "/resources/", utils::URLencode(name, reserved = TRUE), "/files"),
    list(release = release, include_evidence = "true"), timeout
  ))
  if (!identical(listed$resource_id, name) || !identical(listed$version, record$version)) {
    .fail(name, " changed during init; run init again or use a pinned release")
  }
  evidence <- Filter(function(f) identical(f$name, .evidence_file()), listed$files)
  if (!length(evidence)) .fail("No evidence artifact available for ", name)
  evidence[[1]]
}

# Make one resource version complete on disk: its files, then its serving tables. A
# version directory appears only when complete, and then replaces the resource's
# other versions. TRUE when anything was downloaded.
.install_resource <- function(con, root, record, timeout) {
  downloaded <- FALSE
  parent <- file.path(root, "resources", record$resource_id)
  target <- .resource_dir(root, record)
  published <- Filter(function(f) !identical(f$name, .evidence_file()), record$files)
  evidence <- Filter(function(f) identical(f$name, .evidence_file()), record$files)
  if (dir.exists(target) && .present(target, published)) {
    missing <- Filter(function(t) !file.exists(file.path(target, paste0(t, ".parquet"))), .serving())
    if (length(missing)) {
      work <- tempfile(".work-", tmpdir = target)
      dir.create(work)
      on.exit(unlink(work, recursive = TRUE), add = TRUE)
      .write_serving_tables(con, target, work, missing)
    }
    for (artifact in evidence) {
      if (!.present(target, list(artifact))) {
        staging <- tempfile(".download-", tmpdir = target)
        on.exit(unlink(staging), add = TRUE)
        .fetch_artifact(artifact, staging, timeout)
        if (!file.rename(staging, file.path(target, artifact$name))) .fail("Could not place ", artifact$name)
        downloaded <- TRUE
      }
    }
  } else {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    stage <- tempfile(".staging-", tmpdir = parent)
    dir.create(stage)
    on.exit(unlink(stage, recursive = TRUE), add = TRUE)
    for (artifact in record$files) .fetch_artifact(artifact, file.path(stage, artifact$name), timeout)
    work <- file.path(stage, ".work")
    dir.create(work)
    .write_serving_tables(con, stage, work)
    unlink(work, recursive = TRUE)
    if (dir.exists(target)) unlink(target, recursive = TRUE)
    if (!file.rename(stage, target)) .fail("Could not put ", record$resource_id, " in place")
    downloaded <- TRUE
  }
  for (other in list.dirs(parent, recursive = FALSE)) {
    if (!identical(normalizePath(other), normalizePath(target))) unlink(other, recursive = TRUE)
  }
  downloaded
}

#' Install resources and write the database
#'
#' Downloads resources into the data directory and prepares them for querying.
#' Selects by resource name and by the catalog's tags (`topic`, `content`,
#' `molecule`, `knowledge_origin`): several tags of one dimension are alternatives,
#' different dimensions must all match; with no selection, every resource of the
#' release. Resources already installed stay installed and move to this release's
#' versions, so running it again updates the directory. Only missing or new versions
#' are downloaded; an interrupted run keeps the resources it finished.
#'
#' Ends by writing the directory's DuckDB database (`omnipath.duckdb`): native tables
#' (the default for a new directory) or, with `materialize = FALSE`, views over the
#' files. The choice is kept on later runs unless given again. The Python client
#' `omnipath_client` reads and writes the same directory.
#'
#' @param resources Resource names, or `NULL`.
#' @param topic,content,molecule,knowledge_origin Tags to select by.
#' @param release A release, or "latest".
#' @param data_dir The data directory; see [op_data_directory()].
#' @param include_evidence Also download raw evidence payloads (kept on updates).
#' @param materialize Native tables (`TRUE`) or views (`FALSE`); `NULL` keeps the
#'   directory's choice.
#' @param api_url The API; by default `OMNIPATH_API_URL` or the development deployment.
#' @param timeout Seconds without data before a download fails.
#' @param memory_limit DuckDB's memory limit while writing.
#' @param progress Print progress messages.
#' @return The data directory, invisibly.
#' @export
op_init <- function(resources = NULL, topic = NULL, content = NULL, molecule = NULL,
                    knowledge_origin = NULL, release = "latest", data_dir = NULL,
                    include_evidence = FALSE, materialize = NULL, api_url = NULL,
                    timeout = 30, memory_limit = "2GB", progress = TRUE) {
  if (!is.character(release) || length(release) != 1 || !nzchar(trimws(release))) {
    stop("release must not be empty", call. = FALSE)
  }
  api_url <- .api_url(api_url)
  root <- op_data_directory(data_dir)
  say <- function(...) if (isTRUE(progress)) message(...)
  catalog <- .catalog(api_url, release, timeout)
  chosen <- .select(catalog, resources, list(
    topic = topic, content = content, molecule = molecule, knowledge_origin = knowledge_origin
  ))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  on.exit(unlink(file.path(root, ".duckdb-tmp"), recursive = TRUE), add = TRUE)
  .configure(con, root, memory_limit)
  installation <- .read_installation(root) %||% list(
    format = 1L, api_url = api_url, release = release, materialize = TRUE, resources = list()
  )
  materialize <- materialize %||% installation$materialize
  installed <- stats::setNames(installation$resources, vapply(installation$resources, function(r) r$resource_id, ""))
  wanted <- names(chosen)
  # Installed resources follow the release too: one directory, one release.
  targets <- catalog[names(catalog) %in% c(wanted, names(installed))]
  finish <- function() {
    if (length(installation$resources)) {
      installation$materialize <<- materialize
      .write_installation(root, installation)
      if (!.in_sync(root, installation)) {
        say("Writing the database", if (isTRUE(materialize)) " (native tables)" else "")
        .synchronize(root, installation, memory_limit)
      }
    }
  }
  tryCatch({
    number <- 0
    for (entry in targets) {
      number <- number + 1
      record <- entry
      name <- record$resource_id
      had_evidence <- !is.null(installed[[name]]) &&
        any(vapply(installed[[name]]$files, function(f) identical(f$name, .evidence_file()), TRUE))
      record$files <- Filter(function(f) !identical(f$name, .evidence_file()), record$files)
      if ((isTRUE(include_evidence) && name %in% wanted) || had_evidence) {
        record$files <- c(record$files, list(.evidence_artifact(api_url, release, record, timeout)))
      }
      size <- sum(vapply(record$files, function(f) as.numeric(f$size_bytes), 0)) / 1e6
      downloaded <- .install_resource(con, root, record, timeout)
      say(sprintf("[%d/%d] %s %s (%s)", number, length(targets), name, record$version,
        if (downloaded) sprintf("%.0f MB", size) else "up to date"))
      installed[[name]] <- record
      installation <- list(
        format = 1L, api_url = api_url, release = release, materialize = materialize,
        resources = unname(installed)
      )
      .write_installation(root, installation)
    }
  }, finally = finish())
  say(sprintf("%d resources installed in %s", length(installed), root))
  say("Database: ", file.path(root, .database_file))
  invisible(root)
}
