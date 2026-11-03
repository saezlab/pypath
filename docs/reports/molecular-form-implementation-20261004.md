# Molecular forms: implementation and validation

5 October 2026 · isolated pilot release `2026.10.4.1`

This report describes the earlier prototype checkout and its isolated pilot.
The later [integration into `parquet-migration`](molecular-parquet-migration-integration-20261005.md)
preserves that branch's schema-1 manifest layout and adds PostgreSQL/subset
compatibility. The producer revisions, contracts and measurements below refer
to the pilot; they are not measurements of a fresh migration-branch release.

The implementation follows the [decision record](../entity-resolution-and-molecular-forms.md).
SIGNOR and IntAct were rebuilt on nicesrv, validated through Parquet, API and
Python client, and published to an isolated preview. Production deployments
and release pointers are unchanged. Actual browser checks and the final
automatic Group-toggle amendment passed.

## What changed

- `inputs_v2` accepts a structured molecular form and preserves reported
  product, isoform, sequence, modification and variant specificity before
  identifier normalization. Source entity types stay at the main level.
- The resolver assigns a supported NCBI Gene reference and resolves product
  identity separately. It uses explicit, taxon-consistent gene–product links.
  Gene evidence never expands onto proteins; a preferred protein does not
  become an asserted participant. Native and ambiguous identities survive.
- The three output tables remain. Typed entity rows share a gene reference;
  reusable proteins and transcripts have separate product records. Paired
  subject/object forms stay on each relation evidence occurrence. Standalone
  forms stay on entity evidence. No entity is created for each PTM/variant
  combination. Unknown molecular information remains unspecified.
- The API, client, web explorer and exports retain those forms and referenced
  products. Exact product/isoform selection trims evidence to the same endpoint
  and occurrence. Reference and product views are explicit. Catalogue products
  remain distinct from forms observed by sources.
- One Group toggle applies gene-reference grouping and chemical-connectivity
  grouping together. Users do not choose a grouping strategy. Source types
  remain visible, and unresolved or conflicting groups remain separate.

Pypath changes are on the existing **`parquet-migration`** branch at
`9d61c6e6cf1c9c0355d1287758ade6df2d477455`. Both build dependency declarations
and `uv.lock` pin this commit. The implementation stays on the local root
branch `codex/gene-reference-molecular-forms`.

The contracts use resource schema `3.0`, serving schema `4` and disposable
serving projections `.serving/v2`. Older evidence without product pointers
remains unspecified. This continues the Parquet/DuckDB architecture; no
PostgreSQL migration was introduced.

## Actual pilot results

| Resource | Entities | Relations | Source payloads | Build time |
| --- | ---: | ---: | ---: | ---: |
| SIGNOR | 19,245 | 35,630 | 45,765 | 83.28 s |
| IntAct | 91,765 | 730,839 | 1,240,522 | 1,175.23 s |

The [build record](molecular-pilot-builds-20261004.json) contains resolution
statistics and file sizes. These are complete resource builds, not samples.

The [final output validation](molecular-pilot-serving-final-20261005.json)
ran explicit raw and projected passes with matching validator/client hashes:

| Mode | SIGNOR checks | IntAct checks | Failures |
| --- | ---: | ---: | ---: |
| Raw Parquet | 136 passed / 3 unavailable | 115 passed / 3 unavailable | 0 |
| Projected serving | 136 passed / 3 unavailable | 115 passed / 3 unavailable | 0 |

The unavailable cases are absent observed source examples: SIGNOR variants,
IntAct standalone forms, and two opposite-endpoint constraint examples in each
resource. The actual Python client was exercised through hardlinked snapshots;
it was not skipped. Validation checks exact export key multisets, populated
unique pages, paired evidence, standalone provenance, referenced-product closure,
source types, and same-occurrence product/isoform selection.

Combined-resource checks passed Gene11052 grouping with five typed members and
456 relations. Primary DUT protein P33316 has 51 relations; its explicitly
reported isoform P33316-2 has two, one per resource. Raw/projected identities,
counts and paired evidence agree. A relation key shared by both resources was
not available within the selected evidence budget; that case remains unavailable.

The isolated release manifest pins both resources to `2026.10.4.1`. Publication
completed at `2026-10-05T04:59:35Z` after all required output, reference,
combined-resource and benchmark gates passed. The cached taxonomy archive was
checksum-verified; 841 pilot taxa were included, with source taxon `3470597`
absent from that snapshot. Its ID is retained without inventing a label.

## Performance

| Actual query | Raw new data | Projected new data | Old snapshot |
| --- | ---: | ---: | ---: |
| SIGNOR relation browsing | 102 ms | 42 ms | 80 ms |
| IntAct relation browsing | 670 ms | 180 ms | 412 ms |
| Combined gene context, populated page | 427 ms | 572 ms | — |
| Combined exact protein, populated page | 37.88 s | 1.77 s | — |
| Combined isoform, populated page | 53.70 s | 4.66 s | — |

Medians use three requests with application result caches cleared. OS caches
were not flushed. Source snapshots and identity semantics differ, so the old/new
comparison does not isolate schema overhead. Populated navigation used one CPU,
a page limit of one and a 200-item evidence budget; returned page identities and
work matched. Projections improve molecular selection, not every endpoint.
They add about 4.5 MiB for SIGNOR and 96.1 MiB for IntAct. Exact isoform navigation
still takes several seconds and deserves further tuning for a broader release.
See the [populated-page measurements](molecular-populated-navigation-20261005.json).

The Python client initially spent up to three minutes rescanning and converting
nested evidence. The reviewed repair reads typed fields and expands selected
occurrences once with a resource-scoped semi join. Actual IntAct isoform product
record requests now took 2.42–2.62 s; the any/both examples took 2.44/2.39 s.
These are validation-example timings, not a general throughput guarantee.

The [corrected reference measurements](molecular-reference-runtime-20261005.json)
show a 14.868 s first reader (including checksumming 512 component files), a
0.027 s subsequent reader, 0.015–0.025 ms warm direct lookups and a 0.296 ms
five-observation matcher batch. Small warm examples do not forecast full builds.

## Reference and correctness review

The fresh reference contains 304,009,112 entity records and 1,912,480,193
identifier lookups, with exact partition write/read-back proofs. Its gene-role
component adds 44,285,162 lookups and 18,864,597 overlays, about 12.9 GiB (6.8%
of the base indexes). All 512 component partitions have exact proofs; 230 keys
exceed the ten-gene ambiguity cap and keep empty postings.

The [graph audit](molecular-form-reference-qc-20261004.json) confirms gene–product
links do not join identity components. The final indexed audit passed 25 cases,
reported one transcript coverage gap and found zero remaining defects. Supported
gene aliases have an independent lookup component, preventing irrelevant product
postings from removing valid gene symbols or labels.

`refseq_protein:NP_000537.3` maps to three primary proteins (Q53GA5, K7PPA8,
P04637), all supporting Gene7157. The gene resolves, while the product remains
ambiguous. Explicit UniProt evidence can disambiguate it; GeneID or reviewed
status alone cannot. The reported RefSeq version remains attached to evidence.

Missing transcript mappings in the existing consolidated hub snapshot cannot be
recovered by recompilation. Native transcript identities survive; improving that
coverage requires regenerating the relevant hubs.

Local checks include 580 Python passes and one skip in the earlier broad/focused
run, 131 pypath compatibility checks, 55 core checks, 109 resolver/writer checks,
nullable Arrow round trips, both Rust kernels and RefSeq ambiguity regressions.
The latest proof/client repair passed 51 validator tests, 43 client tests and
14 publication-gate scenarios. The web passed all 32 tests, Svelte reported zero
errors/warnings, and the production build passed after readable coordinate
labels replaced raw JSON. The automatic-grouping amendment passed 54 focused
API checks and real raw/projected mixed-result checks; the CURIE facet repair
passed 56 focused checks and independent real raw/projected counts. Later
bounded amendments receive their own checks;
these counts are not one aggregate suite at one revision.

Earlier failures and source/helper revisions remain archived in the isolated
server reports: nullable Arrow handling and gene lookup fixes (`attempt1-failed`),
validator memory/hardlink repair (`attempt2-validator-failed`), the superseded
successful first pass (`attempt3-streaming-validation-passed`), and the validation
mount-environment failure (`attempt4-readonly-unit-environment`). Immutable pilot
outputs were not rebuilt for serving-only repairs.

## Preview and resource coverage

Preview: http://127.0.0.1:4185/explore

Server directory: `/root/projects/omnipath-releases/20261004-molecular-forms`.
The loopback-only, read-only API uses port 18896, reached locally through 8896.
The local web uses 4185. Production API/web/data ports and pointers are unchanged.

Actual browser checks confirmed TP53 Gene7157 has 1,155 relations; selecting
P04637 yields 1,153, and returning restores 1,155. DUT isoform P33316-2 narrows
51 relations to two and displays the paired CDK1/DUT participant forms, sequence
identity and residue-11 phosphorylation. IntAct retains DUT P33316-2 paired
with NEK10 Q6ZWH5; returning restores all 51 relations. CPSF6 Gene11052 shows
five typed members and 456 relations, with one standalone section. Filtering
by RNA selects one typed member; disabling Group restores the five separate
rows. AT-7519 groups by connectivity OVPNQJVDAFNBDN and its relation tab
shows the rebuilt causal records. The final page has no grouping dropdown.

The [completion audit](molecular-completion-audit-20261005.json) records these
checks, source revisions and artifacts. [Automatic grouping](molecular-auto-grouping-20261005.json)
and [CURIE facet counts](molecular-reference-facets-20261005.json) also passed
on actual raw/projected pilots. Facet counts retain existing resource-row
semantics; grouped member counts deduplicate typed entity keys across resources.
Source/display amendments are pinned to `e6cb197` and `f07cc3e`; their guards
proved reference files, immutable pilot outputs, release manifest, validator
and client were unchanged. The earlier full-output validation and benchmark
reports identify their producer revisions; focused actual-data checks cover
these subsequent serving amendments.

All **45 inputs_v2 resource modules were audited** in the
[resource coverage report](molecular-form-input-coverage-20261004.md). SIGNOR and
IntAct have full rebuilt data validation. BindingDB and Reactome have targeted
preservation fixes and fixtures, not full rebuilds. Remaining work includes
Reactome controllers/transcripts, Recon3D selectors, BindingDB sequences, ChEMBL
assay variants, UniProt feature catalogue, miRBase sequence/coordinates and
BRENDA state sections. Broader resource rebuilds, global PTM/variant search and
production deployment remain outside this isolated pilot.
