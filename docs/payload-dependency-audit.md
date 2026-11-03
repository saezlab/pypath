> Historical audit of the earlier record-layout pipeline. Original paths below identify that version; those readers have since been removed. See [current architecture](architecture.md).

# Evidence payload dependency audit

Reviewed 2026-09-30 on the `parquet-migration` branch. This checklist covers the
active packages and scripts in this repository. The retired `legacy/postgres`
implementation is not part of the migrated execution path.

The migration publishes the missing information from `inputs_v2` as
annotations, retains it in the existing PostgreSQL annotation/evidence tables,
and removes raw-record dependence from PostgreSQL products. Complete raw records
remain in the versioned Parquet artifacts.

## Product computations

Before this change, three active areas interpreted 18 named raw fields
to compute PostgreSQL product information. These consumers now read annotations. This is a list of reader paths, not a claim that
every scientifically useful field in every input resource is already mapped.

### Reaction derivation

Reader: `packages/omnipath_postgres/src/omnipath_postgres/reactions.py`.
COSMOS and reaction/transport network views consume its derived tables; they do
not independently parse source payloads.

| Raw fields read | Former purpose | Annotation replacement |
| --- | --- | --- |
| `direction`, `conversion_direction` | Original direction assertion, effective orientation, reversibility and diagnostics | Source-event direction attributes on individual evidence occurrences |
| `participant_role`, `participant_chebi`, `participant_compartment` | Verify membership position, role and identifier before attaching a compartment | Correctly paired participant compartment annotation; preserve source member identity and occurrence attribution |
| `reactants`, `products` | Packed Recon3D/Human-GEM member identifiers, compartments and coefficients | Existing participant annotations, with missing fields mapped upstream |
| `reactant_kegg_id`, `product_kegg_id`, `reactant_stoichiometry`, `product_stoichiometry` | KEGG membership verification and original coefficient text | Existing member identifiers and coefficient annotations; retain symbolic values |

Relevant input modules: `rhea.py`, `recon3d.py`, `metatlas.py` (Human-GEM),
`kegg.py` and `reactome.py`. `kegg_metabolic.py` is an alias, not a separate
mapping implementation. The reaction reader is generic: its coverage must also
include any other source that publishes a molecular-activity subject with
`has_input`, `has_output` or `enabled_by`.

| Input module | Datasets to cover |
| --- | --- |
| `rhea.py` | `reactions`, `metabolic_reactions`, `transport_reactions` |
| `recon3d.py` | `reactions`, `metabolic_reactions`, `transport_reactions` |
| `metatlas.py` / Human-GEM | `reactions`, `metabolic_reactions`, `transport_reactions` |
| `kegg.py` | `reactions`, including factory-created resources |
| `reactome.py` | `reactions` |

Do not add every field again. Current mappings already provide compartments and
stoichiometry in Recon3D/Human-GEM/Reactome and coefficients in KEGG. Existing
artifact versions can predate these mappings. Rhea's current parser does not
provide an explicit conversion direction; preserve unknown rather than infer it
from an equation.

Compatibility requirements:

- Keep conversion direction distinct from gene-effect qualifiers, graph
  `is_directed` and sign. Adding these attributes must not change statement keys.
- Preserve KEGG's existing treatment of right-to-left source assertions: the
  parser has already oriented its inputs and outputs.
- Preserve variable coefficients such as `n` and `(n+1)` as original values.
- Keep each source event separate by resource, version, source, dataset and row;
  do not merge conflicting direction or compartment assertions globally.
- Direction must reach participant-relation evidence, not only the reaction
  entity's annotations: compiled entity annotations merge by entity key without
  preserving the individual source row identity.
- Verify missing/invalid member identifiers and aligned compartment arrays in
  the input mapper. A positional mapping must not assign the next member's
  compartment to an unresolved member.
- Retain meaningful conflict/unknown diagnostics. The old implementation also
  checks full raw-text conflicts and non-object payload shapes; decide their
  validation treatment explicitly rather than silently losing those checks.

### MACdb MetSigDB

Reader: `packages/omnipath_subsets/src/omnipath_subsets/metsigdb/build.py`.
The only interpreted raw field is `Trait_Type`, from entity-owned trait records.
It supplies `metsigdb_membership.set_sub_type`.

Change the `trait` dataset in `inputs_v2/macdb.py` to annotate the trait entity, then query that
annotation in MetSigDB. Preserve the existing resource/version/entity scope and
the current blank-value and multiple-value selection behavior.

### LIANA ligand/receptor roles

Reader: `packages/omnipath_subsets/src/omnipath_subsets/network_views/_query.py`.
The interpreted fields are `Ligand ENSEMBL ID`, `Ligand Species ID`,
`Ligand Symbols`, and the three corresponding `Receptor` fields.

Change the `interactions` dataset in `inputs_v2/connectomedb.py`, which supplies ConnectomeDB2025, to preserve
roles on each interaction's evidence. The final values must identify the correct
published endpoints after canonicalization. A protein's global entity annotation
cannot express its role in a particular interaction. Endpoint-scoped role
attributes must follow canonical endpoint flips. Do not assume canonical
subject means ligand, combine complementary fields from different evidence rows,
or accept ambiguous/contradictory role assignments.

## Additional PostgreSQL payload dependencies removed

These did not introduce further semantic source attributes, but also required
changes to omit full raw records:

- Binary network queries always join `payloads` and include the results in
  their output. This affects both LIANA and MetaLinksDB. ChEMBL mechanisms also
  filter the carried raw records to the selected evidence occurrences.
- Reaction derivation copies a full event record into
  `reaction_context.payload_json`. Reaction network output includes it through
  `to_jsonb(context)`; removing the base payload table alone is insufficient.
- PostgreSQL schema creation, projection, COPY columns, manifest-count
  comparisons, indexes, result counts, documentation and tests currently assume
  a populated payload table. Adjust this plumbing while preserving strict
  verification of the original three Parquet artifacts and their checksums.

## Raw Parquet readers that remain

- `omnipath_api/queries/evidence.py` exposes original records for inspection.
- API resource listing/download/ZIP paths can include the raw payload file;
  inventory and admin/web pages inspect its metadata.
- `omnipath_build/reference/replay_resources.py` replays original payloads
  through source mappers to prepare reference inputs. It needs complete records,
  not just the attributes used by PostgreSQL products.
- Build extraction, writing, consolidation, complex-key remapping and migration
  smoke checks preserve/verify the raw artifact. They are not PostgreSQL product
  queries and should retain that artifact.

The filtered network export service reads compiled relations, not source
payloads. No additional attribute extraction is needed there.

## Verification checks

1. Test each affected mapper with at most 20 original source records. Check
   annotation values and exact source/participant attribution, including missing
   identifiers, opposite roles, contradictory directions and symbolic amounts.
2. Verify annotation survival through Parquet and PostgreSQL, including existing
   measurements, units, comparators, qualifiers and duplicate evidence.
3. Load small annotation-complete fixtures without raw payload storage, make the
   input Parquets unavailable after loading, and rebuild/query reactions,
   COSMOS, MetSigDB, LIANA, MetaLinksDB and reaction/transport views.
4. Compare product meaning and provenance with the current behavior, allowing
   only explicitly documented removal of raw-record output.
5. Recheck all active product code for payload-table joins, raw-field parsing and
   indirect full-record output. Keep raw API inspection/download and reference
   replay tests passing, and keep independent Parquet/API versus monthly
   PostgreSQL release selection unchanged.

The code changes are implemented and reviewed. The native bounded replay on
nicesrv rebuilt seven resources successfully, consuming 20 original records per
resource (140 total). It used the existing immutable resolution reference and
private output paths, without source downloads. The new PostgreSQL import and
PG-only product checks passed; counts and limits are recorded in the
[validation report](source-attributes-validation.md).

The exact rebuild scope for a future full release is Rhea, Recon3D, MetAtlas
(Human-GEM), KEGG, Reactome, MACdb and ConnectomeDB2025. The first three expose
`reactions`, `metabolic_reactions` and `transport_reactions`; KEGG and Reactome
expose `reactions`. Retain every previously published dataset when rebuilding a
complete immutable resource version, including unaffected datasets. The capped
versions are private validation samples, not complete production replacements.

Reaction evidence also carries ordinary `prov:wasDerivedFrom=urn:sha256:...`
and `dcterms:type` attributes. The SHA covers the exact original stored source
text. Import compares matching raw-row hashes/types with published evidence
before discarding the raw text; rebuilding checks source-event hash consistency.
Thus conflicts and unsupported shapes remain visible without a new context table.
Legacy reaction versions lacking these attributes fail with rebuild guidance.

No additional currently discoverable source emits eligible molecular-activity
membership edges. Custom or organism-specific KEGG factory artifacts must also
be rebuilt if their resolved content contains `molecular_activity` subjects with
`has_input`, `has_output` or `enabled_by` and missing source-record SHA metadata.

These checks cover all former product payload readers. They do not establish
full-resource parity or claim every scientifically useful raw field has a mapped
attribute. API raw inspection/download and reference replay retain raw Parquets.
