# Parquet migration: derivations and products

This milestone builds the remaining shared derivations and product adapters from
a fixed, already resolved release in PostgreSQL.

## Implemented

- Ontology edges, source-scoped and release-wide shortest-path closure, term
  inspection and release-visible hierarchy counts. Subclass and part-of paths
  remain separate; cycles do not produce self ancestry.
- Source counts, resource overlap and source-event reaction contexts with
  explicit participant roles, coefficients, verified compartments and diagnostics.
- MetSigDB membership extraction for all five current product rules.
- MetaLinksDB, LIANA and Reactions registry presets and streaming query adapters.
- COSMOS metabolic edges, connectors, reversible halves and label diagnostics.
- One transaction for selected subset rebuilds, with exact release stamps and
  rollback of earlier product changes if a later product fails.

Authoritative per-resource records, aliases, evidence, qualifiers and original
payload text remain available. The canonical `statement` view includes ontology;
`relation` selects graph statements, keeping graph incidence counts distinct.

## Agreed behavior and compatibility limits

MetSigDB uses published chemical entities, including fallback identities.
The existing Parquet schema does not retain the legacy matched flag, and
resolver-labelled aliases occur on both matched and fallback entities. The
chosen policy changes no Parquet files and repeats no resolution.

Network results preserve published qualified statement keys rather than merging
all equal endpoints. Paging applies to source records before grouping, and
contributors may span pages. LIANA roles require matching explicit payload IDs;
transport requires catalytic and compartment evidence. Missing resources and
these limits are reported with query results.

COSMOS reaction indexes identify separate source events. They do not reproduce
the legacy merge of equal participant multisets across sources. Unknown direction
does not justify reverse edges; GPR associations do not become catalysts.

## Verification

Combined PostgreSQL, subset and core tests: **183 passed**.
Lint, formatting and standalone package checks: passed. Built PostgreSQL and subset wheels imported from an isolated environment without pypath, the builder or resolver installed; all **6 orchestration tests passed** against those installed wheels.

Fixtures cover protein interactions, chemical memberships, ontology cycles and
cross-resource paths, reversible and orphan reactions, typed stoichiometry,
ambiguous aliases, namespace collisions, conflicting payload rollback and
deterministic repeat builds. Source fixtures contain at most 20 original records;
derived memberships and reaction spokes are not truncated again. No upstream
resource downloads or full builds ran for this milestone.

These checks establish the implementation on representative bounded data.
They do not establish full-release parity, production performance or compatibility
with every existing web query. The next step is to exercise representative
consumer queries and compare expected product output before backend benchmarking.
