# Entity resolver

`omnipath_resolver` owns complete runtime matching: the `EntityResolver` facade,
`LibraryMatcher`, identifier normalization, biological decision policy, identity
library reads (`IdentityRuntime`), entity records, and the shared namespace/key
encoding. It depends on core
and declared matching libraries; it does not import `omnipath_build`, source
parsers, replay commands, or reference compilers.

```python
from omnipath_resolver import EntityResolver, RawEntityObservation

resolver = EntityResolver(library_dir="data/reference/identity/<fingerprint>")
try:
    observation = RawEntityObservation(
        entity_key="input", entity_type="protein", namespace="uniprot",
        identifier="P04637", taxon="9606",
    )
    targets = resolver.resolve_entity_targets({"input": observation})
finally:
    resolver.close()
```

`RawEntityObservation`, `ResolvedEntityInfo`, and `ResolvedEntityTarget` are public
input/result contracts. `resolve_entities` provides scalar matching;
`resolve_entity_targets` keeps gene reference and asserted product identity separate. Readers pin one
immutable identity library, named by its fingerprint. Without a library, observations retain
their native identities. Runtime resolution reads LMDB point-lookup stores and does not scan
or create Parquet references.

Full InChIKeys take precedence for chemicals. Protein accessions remain separate.
Gene identifiers and symbols identify genes without selecting catalogue proteins.
Primary UniProt product identity uses explicit source evidence and supported gene
links; reviewed status does not disambiguate products. Taxon-qualified symbols
retain their organism scope. Conflicting evidence remains unresolved.

`omnipath_build.hubs` and `omnipath_build.identity` own hub export and the identity
library; see the [identity layer guide](../build/REFERENCE.md). The kv writers reuse
the resolver's key encoding and store layout (`observations.key`, `identity_kv`).
The Rust decision kernel lives in `rust/reference` (crate
`anchor-resolution-serving`) and is linked into the Python extension, whose
bindings remain in `src/`; Python uses the flat `omnipath_resolver/` package.

```sh
uv sync --frozen --all-packages
cargo test --locked --manifest-path packages/resolver/Cargo.toml
cargo test --locked --manifest-path packages/resolver/rust/reference/Cargo.toml
uv run pytest packages/resolver/tests packages/build/tests/test_resolver.py packages/build/tests/identity
```
