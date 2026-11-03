# Runtime matching policy

`policy.py` selects a policy per entity type (gene_protein, chemical, reaction,
or not matched) and the namespaces that vote. `match.py` (`LibraryMatcher`)
turns an observation into votes, resolves them through the runtime the library
manifest names, and records resolution statistics.

For an `omnipath-identity-v2` library the runtime is `IdentityRuntime`
(`../identity_runtime.py`). It derives each batch's candidates from the hub and
identity LMDB stores, applies the precedence and chemical fallback rules, and
decides with the Rust kernel without inferring a protein from a gene. Entity
records are assembled on first use and cached per library fingerprint. No
Parquet scan or join occurs during resolution; a missing or damaged store fails
the build rather than silently falling back.

Rules and outcomes are described in
[entity resolution](../../../../core_documentation/resolution.md). Building the
library is described in [the identity layer guide](../../../build/REFERENCE.md).

Runtime matching stays in `omnipath_resolver`: it imports neither parsers nor
library builders. `omnipath_build.identity` owns construction.
