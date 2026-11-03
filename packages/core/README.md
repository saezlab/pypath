# omnipath-core

Shared contracts for the OmniPath Parquet pipeline:

- Biolink terms, predicate rules and bundled vocabulary (`biolink.py`, `vocab/`)
- PyArrow schemas and mapped record models (`schema.py`, `silver_schema.py`)
- Identity keys and stable hashing (`keys.py`)
- Quantitative measurements and units (`measurements.py`)
- Resource metadata, versioning and layout paths (`resource_metadata.py`, `versioning.py`)

Resource manifests use integer `schema_version: 1` and `serving_schema_version: 6`.
`BuildManifest.files` maps every table file (the published tables in
`PUBLISHED_TABLES` and the serving tables in `SERVING_TABLES`) to `ManifestFile`
values with `size_bytes`, `rows` and `sha256`. `to_dict()` and `from_dict()` preserve the
published JSON shape and additional resource metadata.

`validate_build_manifest()` and `validate_release_manifest()` validate structure
without reading artifacts. PostgreSQL adds pinned-path, checksum and Parquet
schema checks; API inventory supplies its own historical naming policy through
validator callbacks. Release versions and resource versions remain independent.
Package versions and installed software digests are separate publication
provenance, recorded by `omnipath_build.provenance.runtime_provenance()`.
