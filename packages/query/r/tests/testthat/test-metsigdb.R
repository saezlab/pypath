# MetSigDB is answered by the API; these tests need no network. The Python client's tests
# (packages/query/python/tests/test_metsigdb_client.py) run the same questions against the API.

test_that("filters become the API's query parameters", {
  query <- .metsigdb_query(
    list(resource = c("Reactome", "KEGG"), total = TRUE, organism = 9606L, group = FALSE, fields = NULL),
    .metsigdb_filters
  )
  expect_equal(query$resource, "Reactome,KEGG")
  expect_equal(query$total, "true")
  expect_equal(query$group, "false")
  expect_equal(query$organism, "9606")
  expect_null(query$fields)
})

test_that("a filter nobody publishes is refused before a request", {
  expect_error(.metsigdb_query(list(hmdb = "HMDB0000538"), .metsigdb_filters), "Unknown filters.*hmdb")
  expect_error(.metsigdb_query(list("Reactome"), .metsigdb_filters), "Unknown filters")
  # The set listing takes no member-level filter.
  expect_error(.metsigdb_query(list(entity = "C00002"), .metsigdb_set_filters), "Unknown filters.*entity")
})

test_that("a local installation says where the memberships are", {
  service <- local_service()
  service$init("alpha")
  con <- local_connection(service)
  expect_error(op_metsigdb(con, resource = "Reactome"), "answered by the API")
  expect_error(op_metsigdb_sets(con), "answered by the API")
})
