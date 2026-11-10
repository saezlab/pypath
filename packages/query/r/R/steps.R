# Run a query of step files (inst/sql/queries/<name>.sql) on a release database: the
# same files and the same rules as the Python runner (omnipath_client/steps.py), where
# the file format is described. Each step is one statement; lists and values of earlier
# steps go into later steps' SQL as literals, so DuckDB reads rows through the id indexes.

.step_re <- "^--\\s*step:\\s*(\\w+)\\s+AS\\s+(table|list|value|rows)(\\s+when\\s+(.*?))?(\\s+cache)?\\s*$"
.repeat_re <- "^--\\s*repeat:\\s*(\\w+)\\s*=\\s*(.+?)\\s+next\\s+(.+?)\\s+until\\s+(\\w+)\\s+max\\s+(\\d+)\\s*$"
.placeholder_re <- "\\{\\{\\s*(#?)([A-Za-z_][A-Za-z0-9_:]*)(\\[[^]]*\\])?(\\*?)\\s*\\}\\}"
.step_types <- c("VARCHAR", "INTEGER", "BIGINT", "DOUBLE", "BOOLEAN", "JSON")

.match <- function(pattern, text) {
  m <- regmatches(text, regexec(pattern, text, perl = TRUE))[[1]]
  if (length(m) == 0) NULL else m
}

# A step file's header and steps.
.parse_query <- function(text, name = "query") {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  query <- list(name = name, params = list(), outputs = list(), examples = list(), steps = list())
  stack <- list(list())   # open step lists; the last is where steps go
  loops <- list()         # open repeat headers
  current <- NULL
  close <- function() {
    if (!is.null(current)) {
      current$sql <<- trimws(current$sql)
      n <- length(stack)
      stack[[n]][[length(stack[[n]]) + 1]] <<- current
      current <<- NULL
    }
  }
  for (line in lines) {
    stripped <- trimws(line)
    if (!is.null(m <- .match(.step_re, stripped))) {
      close()
      current <- list(type = "step", name = m[[2]], kind = tolower(m[[3]]), sql = "",
        when = if (nzchar(m[[5]])) trimws(m[[5]]) else NULL, cache = nzchar(m[[6]]))
      next
    }
    if (!is.null(m <- .match(.repeat_re, stripped))) {
      close()
      loops[[length(loops) + 1]] <- list(type = "repeat", var = m[[2]], init = m[[3]], `next` = m[[4]],
        until = m[[5]], max = as.integer(m[[6]]))
      stack[[length(stack) + 1]] <- list()
      next
    }
    if (grepl("^--\\s*end\\s*$", stripped, perl = TRUE)) {
      close()
      if (length(loops) == 0) stop(name, ": '-- end' without '-- repeat'")
      loop <- loops[[length(loops)]]
      loops[[length(loops)]] <- NULL
      loop$steps <- stack[[length(stack)]]
      stack[[length(stack)]] <- NULL
      n <- length(stack)
      stack[[n]][[length(stack[[n]]) + 1]] <- loop
      next
    }
    if (is.null(current)) {
      if (!is.null(m <- .match("^--\\s*query:\\s*(\\w+)", stripped))) {
        query$name <- m[[2]]
      } else if (!is.null(m <- .match("^--\\s*param:\\s*(\\w+)\\s+([A-Z]+(\\[\\])?)(\\s*=\\s*(.+))?$", stripped))) {
        if (!(sub("\\[\\]$", "", m[[3]]) %in% .step_types)) stop(name, ": unknown type ", m[[3]])
        default <- if (nzchar(m[[6]]) && m[[6]] != "NULL") jsonlite::fromJSON(m[[6]], simplifyVector = FALSE) else NULL
        query$params[[m[[2]]]] <- list(name = m[[2]], type = m[[3]], default = default)
      } else if (!is.null(m <- .match("^--\\s*output:\\s*(table|json)\\s+(\\w+)", stripped))) {
        query$outputs[[m[[2]]]] <- m[[3]]
      } else if (!is.null(m <- .match("^--\\s*example:\\s*(\\{.*\\})\\s*$", stripped))) {
        query$examples[[length(query$examples) + 1]] <- jsonlite::fromJSON(m[[2]], simplifyVector = FALSE)
      }
      next
    }
    current$sql <- paste0(current$sql, line, "\n")
  }
  close()
  if (length(loops) > 0) stop(name, ": '-- repeat' without '-- end'")
  if (length(query$outputs) == 0) stop(name, ": no '-- output:'")
  query$steps <- stack[[1]]
  query
}

.load_query <- function(name) {
  if (is.null(.state$queries)) .state$queries <- list()
  if (is.null(.state$queries[[name]])) {
    path <- file.path(.sql_dir(), "queries", paste0(name, ".sql"))
    if (!file.exists(path)) stop("Unknown query: ", name)
    .state$queries[[name]] <- .parse_query(paste(readLines(path, warn = FALSE), collapse = "\n"), name)
  }
  .state$queries[[name]]
}

# A SQL literal of one value, a list (`[...]`) or, for JSON, a document.
.step_literal <- function(x, type = NULL) {
  if (!is.null(type) && type == "JSON") {
    if (is.null(x)) return("NULL::JSON")
    return(paste0(.quote(as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null"))), "::JSON"))
  }
  if (!is.null(type) && grepl("\\[\\]$", type)) {
    if (is.null(x)) return(paste0("NULL::", type))
    items <- if (is.list(x)) x else as.list(x)
    if (length(items) == 0) return(paste0("[]::", type))
    return(paste0("[", paste(vapply(items, .scalar_literal, ""), collapse = ", "), "]"))
  }
  if (is.null(x)) return(if (is.null(type)) "NULL" else paste0("NULL::", type))
  if ((is.list(x) || length(x) > 1)) {
    items <- if (is.list(x)) x else as.list(x)
    if (length(items) == 0) return("[]")
    return(paste0("[", paste(vapply(items, .scalar_literal, ""), collapse = ", "), "]"))
  }
  .scalar_literal(x)
}

.coerce_param <- function(param, value) {
  if (is.null(value)) return(NULL)
  type <- param$type
  if (grepl("\\[\\]$", type)) {
    item <- list(name = param$name, type = sub("\\[\\]$", "", type))
    return(lapply(if (is.list(value)) value else as.list(value), function(v) .coerce_param(item, v)))
  }
  bad <- function() stop(param$name, " must be ", type, ", not ", deparse(value), call. = FALSE)
  if (type == "JSON") return(value)
  if (length(value) != 1) bad()
  if (type == "VARCHAR") { if (!is.character(value)) bad(); return(value) }
  if (type %in% c("INTEGER", "BIGINT")) {
    if (is.logical(value) || !(is.numeric(value) || is.character(value))) bad()
    v <- suppressWarnings(as.numeric(value)); if (is.na(v) || v != round(v)) bad(); return(v)
  }
  if (type == "DOUBLE") { if (!is.numeric(value)) bad(); return(as.numeric(value)) }
  if (type == "BOOLEAN") { if (!is.logical(value)) bad(); return(value) }
  value
}

.bind_params <- function(query, params) {
  unknown <- setdiff(names(params), names(query$params))
  if (length(unknown)) stop(query$name, ": unknown parameter(s) ", paste(sort(unknown), collapse = ", "), call. = FALSE)
  values <- list()
  for (name in names(query$params)) {
    p <- query$params[[name]]
    values[name] <- list(if (name %in% names(params)) .coerce_param(p, params[[name]]) else p$default)
  }
  values
}

.run_query <- function(con, query, params = list(), output = "table", cache = NULL, scope = NULL) {
  if (is.character(query)) query <- .load_query(query)
  if (is.null(query$outputs[[output]])) stop(query$name, " has no ", output, " output")
  env <- new.env(parent = emptyenv())
  env$values <- .bind_params(query, params)
  env$rows <- list()
  env$tables <- character()
  types <- vapply(query$params, function(p) p$type, "")

  execute <- function(sql, fetch = TRUE) {
    tryCatch(
      if (fetch) DBI::dbGetQuery(con, sql) else DBI::dbExecute(con, sql),
      error = function(e) {
        message <- conditionMessage(e)
        if (grepl("invalid: ", message, fixed = TRUE)) {
          text <- trimws(sub("[\r\n].*$", "", sub("^.*?invalid: ", "", message, perl = TRUE)))
          stop(structure(class = c("omnipath_query_error", "error", "condition"),
            list(message = text, call = NULL)))
        }
        stop(e)
      }
    )
  }
  number <- function(text) {
    text <- trimws(text)
    if (!is.null(env$values[[text]])) as.integer(env$values[[text]]) else as.integer(text)
  }
  rows_table <- function(name) {
    frame <- env$rows[[name]]
    cols <- names(frame)
    if (nrow(frame) == 0) return(paste0("(SELECT ", paste0("NULL AS ", cols, collapse = ", "), " WHERE false)"))
    values <- vapply(seq_len(nrow(frame)), function(i) {
      paste0("(", paste(vapply(cols, function(c) .scalar_literal(frame[[c]][[i]]), ""), collapse = ", "), ")")
    }, "")
    paste0("(SELECT * FROM (VALUES ", paste(values, collapse = ", "), ") AS v(", paste(cols, collapse = ", "), "))")
  }
  render <- function(sql) {
    out <- ""
    repeat {
      m <- regexpr(.placeholder_re, sql, perl = TRUE)
      if (m < 0) break
      starts <- attr(m, "capture.start"); lengths <- attr(m, "capture.length")
      part <- function(i) substr(sql, starts[i], starts[i] + lengths[i] - 1)
      count <- part(1); name <- part(2); slice <- part(3); spread <- part(4)
      text <- if (startsWith(name, "table:")) {
        rows_table(substring(name, 7))
      } else {
        if (!(name %in% names(env$values))) stop(query$name, ": unknown placeholder ", name)
        value <- env$values[[name]]
        if (nzchar(count)) {
          as.character(length(value))
        } else {
          if (nzchar(slice)) {
            bounds <- strsplit(substr(slice, 2, nchar(slice) - 1), ":", fixed = TRUE)[[1]]
            start <- number(bounds[[1]]); size <- number(bounds[[2]])
            items <- if (is.null(value)) list() else as.list(value)
            value <- items[seq_len(max(0, min(size, length(items) - start))) + start]
            spread <- "*"
          }
          if (nzchar(spread)) {
            items <- if (is.null(value)) list() else as.list(value)
            if (length(items) == 0) "NULL" else paste(vapply(items, .scalar_literal, ""), collapse = ", ")
          } else {
            .step_literal(value, if (name %in% names(types)) types[[name]] else NULL)
          }
        }
      }
      out <- paste0(out, substr(sql, 1, m - 1), text)
      sql <- substring(sql, m + attr(m, "match.length"))
    }
    paste0(out, sql)
  }
  scalar <- function(expression) {
    result <- execute(paste0("SELECT (", render(expression), ") AS v"))
    if (nrow(result) == 0) NULL else result$v[[1]]
  }
  cached <- function(step, sql, compute) {
    if (isTRUE(step$cache) && !is.null(cache)) cache(list(scope, query$name, step$name, sql), compute) else compute()
  }
  run_step <- function(step) {
    skip <- !is.null(step$when) && !isTRUE(scalar(step$when))
    sql <- render(step$sql)
    if (skip) sql <- paste0("SELECT * FROM (", sql, ") WHERE false")
    if (step$kind == "table") {
      target <- paste0("__", step$name)
      env$tables <- union(env$tables, step$name)
      if (isTRUE(step$cache) && !is.null(cache) && !skip) {
        frame <- cached(step, sql, function() execute(sql))
        duckdb::duckdb_register(con, "__cached", frame)
        on.exit(duckdb::duckdb_unregister(con, "__cached"), add = TRUE)
        execute(paste0("CREATE OR REPLACE TEMP TABLE ", target, " AS SELECT * FROM __cached"), fetch = FALSE)
      } else {
        execute(paste0("CREATE OR REPLACE TEMP TABLE ", target, " AS ", sql), fetch = FALSE)
      }
      return(invisible())
    }
    frame <- if (skip) execute(sql) else cached(step, sql, function() execute(sql))
    if (step$kind == "list") {
      env$values[step$name] <- list(if (nrow(frame) == 0) list() else as.list(frame[[1]]))
    } else if (step$kind == "value") {
      env$values[step$name] <- list(if (nrow(frame) == 0) NULL else frame[[1]][[1]])
    } else {
      env$rows[[step$name]] <- frame
      env$values[step$name] <- list(frame)
    }
  }
  run_steps <- function(items) {
    for (item in items) {
      if (identical(item$type, "repeat")) {
        env$values[item$var] <- list(scalar(item$init))
        for (round in seq_len(item$max)) {
          run_steps(item$steps)
          if (isTRUE(env$values[[item$until]]) || round == item$max) break
          env$values[item$var] <- list(scalar(item$`next`))
        }
      } else {
        run_step(item)
      }
    }
  }
  on.exit(for (name in env$tables) try(DBI::dbExecute(con, paste0("DROP TABLE IF EXISTS __", name)), silent = TRUE), add = TRUE)
  run_steps(query$steps)
  target <- query$outputs[[output]]
  if (output == "json") return(env$values[[target]])
  execute(paste0("SELECT * FROM __", target))
}
