#' The data directory
#'
#' The directory given, else the environment variable `OMNIPATH_DATA_DIR`, else the
#' per-user default: the same directory the Python client uses (`platformdirs`'
#' user data directory for "omnipath"), so both clients share an installation.
#'
#' @param path A directory, or `NULL`.
#' @return The absolute, normalised path.
#' @export
op_data_directory <- function(path = NULL) {
  chosen <- path %||% Sys.getenv("OMNIPATH_DATA_DIR", "")
  if (!nzchar(chosen)) chosen <- .default_data_directory()
  .normpath(path.expand(chosen))
}

.default_data_directory <- function() {
  system <- Sys.info()[["sysname"]]
  if (identical(system, "Darwin")) {
    return(file.path(path.expand("~"), "Library", "Application Support", "omnipath"))
  }
  if (identical(system, "Windows")) {
    return(file.path(Sys.getenv("LOCALAPPDATA"), "omnipath", "omnipath"))
  }
  xdg <- Sys.getenv("XDG_DATA_HOME", "")
  base <- if (nzchar(xdg)) xdg else file.path(path.expand("~"), ".local", "share")
  file.path(base, "omnipath")
}

# Absolute, with "." and ".." removed lexically, like Python's os.path.normpath: the
# database records its directory, and both clients must write the same string.
.normpath <- function(path) {
  windows <- .Platform$OS.type == "windows"
  path <- gsub("\\\\", "/", path)
  absolute <- startsWith(path, "/") || grepl("^[A-Za-z]:/", path)
  if (!absolute) path <- file.path(gsub("\\\\", "/", getwd()), path)
  drive <- if (grepl("^[A-Za-z]:/", path)) substr(path, 1, 2) else ""
  if (nzchar(drive)) path <- substring(path, 3)
  parts <- strsplit(path, "/", fixed = TRUE)[[1]]
  kept <- character()
  for (part in parts) {
    if (part == "" || part == ".") next
    if (part == "..") {
      if (length(kept)) kept <- kept[-length(kept)]
    } else {
      kept <- c(kept, part)
    }
  }
  result <- paste0(drive, "/", paste(kept, collapse = "/"))
  if (windows) result <- gsub("/", "\\\\", result)
  result
}

.resource_dir <- function(root, record) {
  file.path(root, "resources", record$resource_id, record$version)
}
