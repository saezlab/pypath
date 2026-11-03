# Taxonomy labels

NCBI Taxonomy is reference data, separate from entity and relation records. The
reference contains `taxon_id`, `current_taxon_id`, `scientific_name` and
`common_name` (NCBI's GenBank common name). A merged taxon inherits its current
name while retaining its original ID for identity and filtering.

Provision the official NCBI download at
`DATA_ROOT/references/taxonomy/taxdump.tar.gz`. Download to a temporary file and
rename only on success. Source: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
The download is an explicit maintenance step; publishing and serving never
perform network requests. Replace this cached dump when updating the reference.

When publishing a release, `ReleaseStore.publish` scans the taxon columns of its
pinned entity and relation files and extracts names for those IDs from the cached
dump. It writes `references/taxonomy/<sha256>/taxonomy.parquet` atomically and pins
that checksum in `references.taxonomy` in the release manifest. The manifest also
records the source URL, source archive checksum, taxon count and missing IDs.
Unknown/deleted IDs remain visible and missing names produce a publication warning.
Chemical sentinel taxon `0` is excluded. Existing entity/relation schemas do not change.

The API caches the small release subset, with request-local selection. It returns
`taxonomyName` on entity summaries and `facetLabel` on taxonomy facets. UI formatting
is shared: common name if available, otherwise scientific name, then the original
ID. Bare IDs remain the fallback; no UI species dictionary exists.

Named releases without a reference remain readable with bare IDs. Their manifests
are not retroactively modified. `latest` uses the most recent published taxonomy
reference; newly built resources need publication to add newly encountered taxa.
Stores without a cached dump may publish legacy releases, with a warning. Once the
dump is provisioned, subsequent releases generate their references automatically.

Validation covers common/scientific names, chained merged IDs, unknown IDs,
reference checksums, entity details, both facet types, and separation of cached
facet counts from release-specific labels.
