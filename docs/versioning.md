# Resource versions and OmniPath releases

A **resource version** identifies one immutable set of compiled Parquet files. An
**OmniPath release** pins one version of each included resource. Creating a release
writes a small JSON manifest and reuses the resource files.

For example, OmniPath `2026.09` can contain SIGNOR `1.0.0` and ChEBI `2.1.0`.
OmniPath `2026.10` can update SIGNOR to `1.1.0` while keeping ChEBI `2.1.0`.
These are example version numbers; no releases are created automatically.

Parquet resources can update independently at any time. The API and explorer
keep Latest as their default. PostgreSQL follows a separate monthly schedule,
loading the selected resolved Parquet versions into its own fixed snapshot.
Publishing a resource or a named Parquet release does not rebuild PostgreSQL
or change its snapshot.

## Build a resource version

From the migration repository:

```bash
uv run --frozen omnipath-build build signor --version 1.0.0 --max-records 20 --output-dir data
```

Versions are explicitly assigned dotted numbers, such as `1`, `1.0.0` or
`2026.09`. They describe the compiled resource, independently of the upstream
provider's release number. Bump the version whenever data, parser behavior,
resolution rules or the output schema changes. Reusing an existing version fails.

Batch builds accept `--versions resource-versions.json`, where the file contains
`{"signor": "1.0.0", "chebi": "2.1.0"}`. Every selected source needs a version.
Alternatively, `--version` applies the same explicit number to each selected source.
The admin page has a resource-version input for individual and batch builds too.

Builds run in temporary directories. Only after successful finalization and
Parquet validation is the complete directory moved into
`data/resources/<source>/<version>/`. Failed builds remain invisible to the API.
The directory includes the eight published tables, the three serving tables,
resolution statistics and a `build_manifest.json` containing the resource version,
schema versions, timestamp, row counts and file checksums. A per-version lock prevents concurrent overwrites.

The build manifest stays at `schema_version = 1`, with filename-keyed entries
for every table file. The normalized tables use `serving_schema_version = 5`.
This changes the record contract, not the manifest layout or the independent
publication schedules. PostgreSQL loads the normalized tables directly; existing
monthly snapshots require a new load to obtain the current contract.

## Publish a pinned OmniPath release

Create a manifest, for example `release.json`:

```json
{
  "schema_version": 1,
  "version": "2026.09",
  "resources": {
    "signor": "1.0.0",
    "chebi": "2.1.0"
  }
}
```

Publish it after all referenced resource builds exist:

```bash
uv run --frozen python -m omnipath_api.releases --data-root data publish release.json
```

This validates the referenced Parquets and atomically creates
`data/releases/2026.09.json`. Published release manifests cannot be overwritten.
Latest remains the default when new releases are published.

Publishing also records `resource_manifests`: the SHA-256 of each pinned
resource's `build_manifest.json`. The PostgreSQL loader refuses a resource whose
build manifest no longer matches its pin, so a release identity always names
exact content. Sample builds (a `max_records` cap or a `datasets` subset) are
refused by both the publisher and the loader unless the release explicitly sets
`"partial_resources": true`.

The equivalent admin API operation is `POST /api/admin/releases` (manifest as
JSON body), using the
existing admin authentication. Admin deletion refuses resource versions pinned
by any release. Manage publication and deletion through these interfaces; direct
filesystem changes bypass their protections.

## Select a release

The dropdown beside the OmniPath logo offers **Latest** (the default) and each
numbered OmniPath release for the explorer,
resource catalog, evidence, statistics and downloads. Resource cards display
their pinned resource version. The selection is carried in shareable
`?release=2026.09` URLs and remembered for subsequent navigation. Switching
releases reloads the page to discard previous query results.

API clients can use `?release=2026.09` on any data endpoint, or send
`X-OmniPath-Release: 2026.09`. Query parameters take precedence over headers, then
the explorer cookie, then Latest. `GET /api/releases` lists
release manifests and the default. A resource filter can narrow a numbered release but
cannot select a different version or an unpinned resource. Unknown releases or
missing pinned artifacts fail explicitly.

**Latest** (`release=latest`) selects the highest numeric version of every
resource, including newly added resources. Version components are compared
numerically, so `1.10.0` is newer than `1.9.0`; copying an older version more
recently does not promote it. Incomplete builds are skipped. For resources with
only legacy nonnumeric versions, file modification time is the fallback;
numeric versions take precedence once available. Latest updates as builds arrive,
while numbered OmniPath releases remain fixed for reproducible analysis.

Old `release=working` links and cookies are accepted as aliases for Latest.
Legacy `default_release.json` files are ignored: the default is always Latest.

Existing hash-named directories remain readable and can be explicitly pinned
without copying them. New builds require numeric versions. This implements
immutable artifact selection; exact rebuild reproducibility still requires
recording and preserving upstream inputs, parser code and resolver snapshots.
Object-storage publication remains separate from this local filesystem workflow.
