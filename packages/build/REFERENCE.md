# Identity layer and resolver

Resource builds resolve against an identity library (`omnipath-identity-v2`),
built offline from one set of hubs. Each hub is indexed on its own, then the
identity decisions (lipid-name anchors, one-step attachment, grouping, quarantine,
gene–product links) are computed once for all hubs. Rules and decisions are
described in [entity resolution](../../core_documentation/resolution.md); the
implementation spec is [docs/identity-layer-spec.md](../../docs/identity-layer-spec.md).

From the repository root, `make hubs` and `make reference` run the whole sequence
(see the [README](../../README.md#build-and-publish-resources) for their variables):

```sh
export PYPATH_DOWNLOAD_DATADIR="$PWD/data/pypath-data"
make hubs HUB_ARGS="--max-records 0"
make reference
```

`make reference` runs these commands for every identity hub (`HUBS` and `ontology`
in `omnipath_build.identity.common`):

```sh
# Per hub: the rule-independent index and its LMDB point-lookup store.
omnipath-build build-hub-index --hub chebi --hubs-dir data/reference/hubs --output-root data/reference/hub-index
omnipath-build build-hub-kv --hub chebi --hub-index-root data/reference/hub-index
# Once for all hubs: the identity decisions and their LMDB store.
omnipath-build build-identity --hub-index-root data/reference/hub-index --output-dir data/reference/identity
```

The ontology hub only supplies term labels, so it gets an index but no kv store.
The hub kv format is `omnipath-hub-kv-v3`: the rows of the gene-symbol namespaces
(`genesymbol`, `genesymbol-syn`) are keyed by `(ns, identifier, taxon)` instead of
`(ns, identifier)`, so a taxon-scoped symbol lookup reads only the records of its
taxon. The kv is built from the `by_id` Parquet files, which already carry the
record's taxon, so an index is not rebuilt for it; `build-hub-kv` replaces a v2
store only when the v3 one is complete. The runtime reads both formats.
`build-identity` also writes the identity kv store; `build-identity-kv --force`
rebuilds it for an existing `data/reference/identity/<fingerprint>`.

A hub index is rebuilt only when its hub export changes. The identity library is
named by a fingerprint of its hub indexes and rules code; a build passes it as
`--library-dir` (or `OMNIPATH_LIBRARY_DIR`) and records the fingerprint in its
provenance. There is no candidate limit: every identifier keeps all candidates,
and matching abstains when more than one survives.

Entity records are assembled on first use and cached per fingerprint in SQLite
(`OMNIPATH_IDENTITY_CACHE`); a cached record equals a freshly built one.

`omnipath-build resolution-regression` compares two identity libraries on the same
extracted observations (`extract`, then `resolve` once per library, then `diff`).
