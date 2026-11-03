# Reference resolver

> **Superseded.** This describes the previous compact reference (`build-library`).
> Resource builds now resolve against the identity layer; see
> [entity resolution](../core_documentation/resolution.md) and
> [the identity layer guide](../packages/build/REFERENCE.md). Its code has been
> removed.

The reference resolver maps already-normalized identifiers to canonical chemical
and protein entities. Chemical entities are anchored primarily by full InChIKey;
protein entities are anchored primarily by primary UniProt accession. Additional
identity and relationship rules are applied while the reference is built.

The new runtime stores the result of that work in two immutable disk indexes. It
does not reconstruct mappings or enrichment for each incoming batch.

## System overview

```text
BUILD TIME
normalized hubs
  -> entity assignment around InChIKey / primary UniProt anchors
  -> admitted identifier-to-entity mappings
  -> candidate-count policy
  -> identifier index + entity index + ambiguous.parquet

RUNTIME
normalized input identifiers
  -> direct identifier-index reads
  -> Rust evidence intersection and resolution policy
  -> direct entity-index reads
  -> resolved rows with complete reference identifiers
```

The assigned Parquet reference and normalized hubs remain build, audit and
reconstruction inputs. They are not needed by the resolver at runtime.
`ambiguous.parquet` is also an audit artifact and is not read during resolution.

## Runtime indexes

The identifier index maps this scoped key:

```text
(target kind, route, namespace, identifier, optional taxon)
    -> candidate entity records + precomputed gene-evidence flags
```

One identifier may map to several entities. That is ordinary data, for example
when a gene identifier maps to several proteins. Candidate lists are stored in
full and are passed to the existing Rust decision kernel, which intersects the
available evidence and applies the chemical and gene-to-protein policies.

The entity index maps:

```text
entity ID
    -> kind, anchor, taxon, display label and all admitted identifiers
```

This second read is the complete enrichment step. It requires no joins or
aggregation. Both indexes are partitioned LMDB stores. Records use positional
MessagePack and selective Zstandard compression with one 16 KiB dictionary per
store. LMDB and the operating system provide disk indexing and page caching;
there is no application result cache.

## Ambiguous identifiers

A scoped lookup with more than 10 candidate entities is removed as a whole. Its
candidate list is never truncated, because truncation could manufacture a false
unique match. The complete excluded mapping is written to `ambiguous.parquet`
with its namespace, identifier, taxon, target, route, candidate IDs and evidence
flags.

The limit is evaluated independently for global and taxon-specific mappings. An
identifier can therefore be excluded globally and remain usable in a taxon where
it has at most 10 candidates. Corresponding aliases are removed from entity
records so forward enrichment cannot reintroduce an excluded mapping. A literal
`-` has no special rule; it follows the same candidate-count policy.

In the complete build, 218,173 of 1,891,367,762 lookup entries were excluded
(0.0115%). The 38.8 MiB ambiguity archive retains all 6,717,190 candidate links.

## Comparison with the previous resolver

| Area | Previous Parquet path | New two-index path |
| --- | --- | --- |
| Identifier lookup | Prepare batch-specific lookup data from the reference access tables | Read the complete precomputed candidate list directly |
| Decision logic | Rust policy kernel | Same Rust policy kernel |
| Enrichment | Query and aggregate reference aliases for accepted entities | Read the precomputed entity record directly |
| Runtime files | Assigned Parquet catalogue and its lookup tables | Two LMDB indexes, manifest and two small dictionaries |
| Cache | Prepared/intermediate work plus storage caches | OS-backed LMDB pages only |
| High-fanout identifiers | Large lists enter runtime and cost proportionally | Removed at build time and preserved in the ambiguity archive |
| Runtime failure boundary | Reference preparation and joins are part of each batch | Missing or corrupt index data fails directly; no fallback scan |

Under four CPU cores and a 4 GiB memory limit, the measured comparison was:

| Workload | Previous path | New path, cold-cache advice | New path, warm |
| --- | ---: | ---: | ---: |
| 11,000 real rows | 515.745 s | 4.899 s | 0.651–0.699 s |
| 1,000 broad reference rows | 18.604 s | 0.886 s | — |

All compared candidate sets, decisions and complete entity records matched under
the 10-candidate policy. A separate cold run resolved one million broadly sampled
rows in 290.241 seconds (3,445 rows/s). The 11,000-row cold run read 177 MiB from
disk, compared with 28.1 GiB for the previous path.

These timings start with normalized rows and end with complete resolved records.
They include index opening and exclude source extraction, identifier
normalization, process startup and final resource-file writing.

## Storage and build status

The complete runtime contains 288,430,295 entities and 1,891,149,589 admitted
lookup entries. Its files occupy 173.83 GiB, plus the 38.8 MiB ambiguity archive.
The compact representation reduced an equivalent JSON encoding of these two
indexes from 455.02 GiB to 170.75 GiB before policy edits, a 62.47% reduction.
That is a representation comparison, not a claim about the size of the previous
Parquet runtime.

Every converted record was decoded and compared with its source. Before applying
the candidate limit, all 14,535 saved full-reference lookup checks and all 7,850
saved complete entity records matched the previous path. After applying the
policy, every archived key was verified absent, with surviving lookups and
corrected entity records sampled across all partitions.

The full compact generation is installed at the canonical server library path,
and the compact worker image has resolved representative protein, chemical and
unmatched inputs against it. The assigned Parquet runtime and the experimental
build copies have been removed. The direct hubs-to-compact compiler is
implemented and tested end to end on fixtures; this full generation was created
by converting the completed reference build, so a fresh full-scale run of the
direct compiler is not yet measured.

The original detailed benchmark and fanout reports are not included in this
repository snapshot. The figures above are historical observations, not measurements
of a newly compiled reference. Current operations are recorded in [the build log](build-log.md).
