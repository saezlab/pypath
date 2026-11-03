# Network-independent tests: a mocked catalog and file host over the fixtures, real
# DuckDB. They mirror the Python client's tests; packages/query/python/tests/test_r_client.py
# compares the two clients' answers on the same directories.

test_that("init selects by name and tag", {
  service <- local_service()
  listed <- op_available(api_url = "https://api.example/api")
  expect_equal(vapply(listed, function(r) r$resource_id, ""), c("alpha", "beta", "gamma"))
  service$init(topic = "signaling", release = "2026.09")
  expect_equal(service$requests[[2]]$query$release, "2026.09")
  expect_equal(installed_names(service), "alpha")
  expect_equal(op_resources(service$data)[[1]]$license, "CC0")
  service$init(topic = c("signaling", "metabolism"), molecule = "metabolites")
  expect_equal(sort(installed_names(service)), c("alpha", "gamma"))
  expect_error(service$init(topic = "astrology"), "Unknown topic.*available: \\['metabolism', 'signaling'\\]")
  expect_error(service$init("missing"), "Unknown resources")
  expect_error(service$init(character()), "at least one")
  expect_error(service$init("alpha", topic = "metabolism"), "No resource matches")
  service$init()
  expect_length(op_resources(service$data), 3)
})

test_that("init is incremental and follows the release", {
  service <- local_service()
  service$init("alpha")
  expect_length(service$downloads, 7)
  expect_false(any(grepl("evidence_payloads", service$downloads)))
  service$downloads <- character()
  service$init("alpha")
  expect_length(service$downloads, 0)
  service$init("beta", include_evidence = TRUE)
  expect_equal(unique(sub("/.*", "", service$downloads)), "beta")
  expect_length(service$downloads, 8)
  service$downloads <- character()
  service$init("alpha", include_evidence = TRUE)
  expect_equal(service$downloads, "alpha/1.0/evidence_payloads.parquet")
  # A new version replaces the old one, also for resources not named this time, and
  # keeps the evidence that was installed.
  service$catalog <- 2
  service$downloads <- character()
  service$init("beta")
  expect_length(service$downloads, 8)
  expect_equal(basename(list.dirs(file.path(service$data, "resources/alpha"), recursive = FALSE)), "2.0")
  con <- local_connection(service)
  expect_equal(
    op_sql(con, "SELECT * FROM resources ORDER BY ALL"),
    data.frame(resource_id = c("alpha", "beta"), version = c("2.0", "1.0"), has_evidence = c(TRUE, TRUE))
  )
  expect_equal(op_sql(con, "SELECT count(*) AS n FROM alpha.evidence_payloads")$n, 1)
})

test_that("derived serving tables match the build", {
  service <- local_service()
  service$init("alpha")
  con <- DBI::dbConnect(duckdb::duckdb())
  withr::defer(DBI::dbDisconnect(con, shutdown = TRUE))
  for (table in c("relation_endpoint", "entity_group", "entity_term", "relation_facet")) {
    rows <- function(path) {
      DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s') ORDER BY ALL", path))
    }
    derived <- rows(file.path(service$data, "resources/alpha/1.0", paste0(table, ".parquet")))
    built <- rows(fixture("server/alpha/1.0", paste0(table, ".parquet")))
    expect_gt(nrow(derived), 0)
    expect_equal(derived, built, ignore_attr = TRUE)
  }
})

queries <- function(con) {
  list(
    tables = op_sql(con, paste(
      "SELECT table_schema, table_name FROM information_schema.tables",
      "WHERE table_schema <> 'main' OR table_name NOT IN ('database_info')",
      "AND NOT starts_with(table_name, 'op_') ORDER BY ALL"
    )),
    entities = op_sql(con, "SELECT label, _resource, _resource_version FROM entity ORDER BY ALL"),
    quoted = op_sql(con, "SELECT identifier FROM entity WHERE label = ?", list("O'Brien")),
    additive = op_sql(con, "SELECT count(*) AS n FROM entity WHERE new_field IS NULL"),
    one = op_sql(con, "SELECT count(*) AS n FROM alpha.relation"),
    terms = op_sql(con, "SELECT count(*) AS n FROM beta.entity_term WHERE term = 'tp53'")
  )
}

test_that("SQL on views and native tables gives the same answers", {
  service <- local_service()
  service$init(c("alpha", "beta"), materialize = FALSE)
  con <- service$connect()
  views <- queries(con)
  expect_equal(nrow(views$entities), 6)
  expect_equal(views$quoted$identifier, c("quote", "quote"))
  expect_equal(views$additive$n, 3)
  expect_equal(op_sql(con, "SELECT DISTINCT table_type FROM information_schema.tables WHERE table_name = 'entity'")$table_type, "VIEW")
  expect_error(op_sql(con, "CREATE TABLE mine AS SELECT 1"))
  op_disconnect(con)
  service$init(c("alpha", "beta"), materialize = TRUE)
  con <- service$connect()
  expect_equal(queries(con), views)
  expect_equal(op_sql(con, "SELECT table_type FROM information_schema.tables WHERE table_name = 'entity' AND table_schema = 'main'")$table_type, "BASE TABLE")
  op_disconnect(con)
  # A new version and a new resource replace only their rows.
  service$catalog <- 2
  service$init("gamma")
  con <- local_connection(service)
  expect_equal(
    op_sql(con, "SELECT _resource, _resource_version, count(*) AS n FROM entity GROUP BY ALL ORDER BY ALL"),
    data.frame(`_resource` = c("alpha", "beta", "gamma"), `_resource_version` = c("2.0", "1.0", "1.0"), n = c(3, 3, 3), check.names = FALSE)
  )
  expect_true(jsonlite::fromJSON(file.path(service$data, "installed.json"))$materialize)
})

test_that("prebuilt queries answer like the explorer", {
  service <- local_service()
  service$init(c("alpha", "beta"), include_evidence = TRUE)
  con <- local_connection(service)
  found <- op_search(con, "TP53")
  expect_equal(names(found), c(
    "entity_key", "label", "entity_type", "namespace", "identifier", "taxon",
    "relation_count", "reference_entity_key", "gene_reference_keys"
  ))
  expect_equal(nrow(found), 1)
  expect_equal(found$label, "TP53")
  expect_equal(found$reference_entity_key, "entrez:7157")
  expect_equal(op_search(con, "11998")$label, "TP53")
  expect_equal(nrow(op_search(con, "TP53", entity_types = "small_molecule")), 0)
  expect_equal(nrow(op_search(con, "no such thing")), 0)

  relations <- op_relations(con, "TP53")
  expect_true(all(c("relation_key", "annotations", "evidence", "sources") %in% names(relations)))
  expect_equal(nrow(relations), 4)
  one <- op_relations(con, "TP53", sources = "alpha", predicates = "affects")
  expect_equal(one$subject_label, "TP53")
  expect_equal(one$object_label, "caffeine")
  expect_equal(nrow(op_relations(con, c("caffeine", "O'Brien"))), 4)
  expect_equal(nrow(op_relations(con, "gene:entrez:7157")), 4)
  expect_equal(nrow(op_relations(con, limit = 3)), 3)
  expect_error(op_relations(con, "TP53", sources = "gamma"), "unavailable")
  expect_error(op_relations(con, "TP53", source = "alpha"), "Unknown filters: \\['source'\\]")
  expect_error(op_relations(con, "TP53", predicates = "nonsense"), "Unknown Biolink predicate: nonsense")
  expect_error(op_relations(con, limit = 100001), "limit")

  details <- op_entity(con, found$entity_key)
  expect_equal(details$entity$label, "TP53")
  expect_null(op_entity(con, "TP53"))
  # The macros are SQL in the database.
  expect_equal(op_sql(con, "SELECT count(*) AS n FROM op_relations(entities := ['TP53'])")$n, 4)
})

test_that("missing data and damaged files are errors", {
  service <- local_service()
  expect_error(op_connect(file.path(service$data, "empty")), "op_init")
  service$init("alpha")
  path <- file.path(service$data, "resources/alpha/1.0/entity.parquet")
  original <- readBin(path, "raw", file.size(path))
  writeBin(original[-length(original)], path)
  expect_error(service$connect(), "Missing or truncated")
  service$init("alpha")
  expect_equal(readBin(path, "raw", file.size(path)), original)
  unlink(file.path(service$data, "omnipath.duckdb"))
  con <- local_connection(service)
  expect_equal(op_sql(con, "SELECT count(*) AS n FROM entity")$n, 3)
})

test_that("a failed download installs nothing of that resource", {
  service <- local_service()
  service$truncated <- "https://data.example/resources/beta/1.0/relation.parquet"
  expect_error(service$init(c("alpha", "beta")), "Truncated")
  expect_length(list.dirs(file.path(service$data, "resources/beta"), recursive = FALSE), 0)
  expect_equal(installed_names(service), "alpha")
  service$truncated <- character()
  service$init("beta")
  con <- local_connection(service)
  expect_equal(op_sql(con, "SELECT count(DISTINCT _resource) AS n FROM relation")$n, 2)
})

test_that("locations follow the environment and normalise", {
  withr::local_envvar(OMNIPATH_DATA_DIR = "/tmp/a/../b")
  expect_equal(op_data_directory(), "/tmp/b")
  expect_equal(op_data_directory("/x/./y/"), "/x/y")
  withr::local_envvar(OMNIPATH_DATA_DIR = "")
  expect_match(op_data_directory(), "omnipath$")
})
