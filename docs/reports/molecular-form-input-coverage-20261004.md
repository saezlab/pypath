# Molecular-form input coverage

Static audit · 4 October 2026 · updated for current `parquet-migration`

## Scope and interpretation

Audited all **45 resource modules** declaring `resource = Resource(...)` in
pypath `inputs_v2`, initially at
`51aedf4a0f37ef66dce31274d614866a0f04ebb2`, then reconciled with the current
remote `parquet-migration` head
`e27d72e095f53b1a0fd76d63bbe059c10f8a48d8`, including their discoverable parser
modules, builder fields and query projections. Changes live in the isolated
`/private/tmp/omnipath-molecular-pypath` checkout; the original sibling checkout
is unchanged. The molecular commits were reapplied without conflicts onto that
existing branch: pilot commit `545b7f5ed44a5414bd01271d0bdb90dd8de1dbd8`, then
additional mapper commit/final head
`9d61c6e6cf1c9c0355d1287758ade6df2d477455` on **`parquet-migration`**. The earlier
isolated feature branch is only a local recovery snapshot and is not a
publication target. Current review patches are based on `e27d72e09`; the source
diff contains eight Python files, excluding tests/data/caches/Git metadata.
Review artifacts:
`/private/tmp/omnipath-molecular-pypath-parquet-migration.patch` and
`/private/tmp/omnipath-molecular-pypath-parquet-migration-source.patch`.

The newer remote changes are preserved: portable shared-core dependency
`>=0.1.0,<0.2`, source context for conversion direction/participant role/trait
type, and cellular-location evidence in Reactome, Recon3D, Human-GEM and Rhea.
Those attributes do not alter the form findings below. Integration requires the
matching `omnipath_core.source_attributes` definitions alongside the new form
API. That compatibility module has been restored using the exact deployed
definitions supplied by the parent task, preserving the upstream annotation
terms without importing other architecture.

This is a code-coverage audit, **not a download audit** of every upstream release.
“No detailed field selected” means the current parser/mapper does not select an
explicit feature; it does not claim that the upstream database has no such data.
The generic build helper can preserve an asserted primary UniProt isoform/chain,
RefSeq or Ensembl product identifier when it survives parsing. It cannot recover
discarded suffixes, infer a transcript from a gene symbol, or reconstruct a PTM
or variant from arbitrary names, comments or catalogue cross-references.

Source types remain main-level types. Bare UniProt accessions do not assert the
canonical sequence. Empty forms mean unspecified. Catalogue feature alternatives
must not be converted into a combined observed molecular form.

## Implemented source coverage

| Resource / dataset | Before | Current structured retention and limit | Code evidence |
| --- | --- | --- | --- |
| SIGNOR / interactions | `_normalize_signor_identifier` removes `-PRO_…`; feature column is description text. | Capture primary isoform/chain IDs before normalization; exact source PTM labels and MOD accessions become individual modifications; explicit mutation feature labels become variants. Exact ranges retain positions, fuzzy ranges retain descriptions. Alternative identifier lists do not become participant assertions. | `pypath/inputs_v2/signor.py::_participant_builder`; `pypath/inputs_v2/_molecular_forms.py::mitab_participant_form` |
| SIGNOR / complexes, protein families | Member isoform accessions survive identifiers but have no occurrence form. | `MembersFromList` now accepts an indexed molecular-form callback; each asserted member isoform retains its own form. | `signor.py::_group_member_form`, `_group_members`; `pypath/internals/tabular_builder.py::MembersFromList` |
| IntAct / interactions | `_parse_identifier_pairs` removes `-PRO_…`; no structured participant feature mapping. | Same primary participant capture as SIGNOR. Source feature text is retained as description even where no approved structured feature mapping applies. Source type controls capture; gene endpoints do not acquire forms from UniProt alternative cross-references. | `pypath/inputs_v2/intact.py::interactor_a_builder`, `interactor_b_builder`; `_molecular_forms.py` |
| BindingDB / interactions | The UniProt regex rejects explicit isoform and chain suffixes; multichain targets are already represented as separate members. | Accept explicit `-n` and `-PRO_…` suffixes and retain each chain's form through `_target`. No mutation is inferred from `Target Name`. Target sequence columns are still omitted by the narrow parser projection. | `pypath/inputs_v2/bindingdb.py::f`, `_target_chain_form`, `_target`; `pypath/inputs_v2/parsers/bindingdb.py::_BINDINGDB_COLUMNS` |
| Reactome / reaction participants | BioPAX `modificationType` labels and `sequencePosition` values are flattened into `participant_modification`; feature membership and positional status are lost. | Capture individual modification features before flattening, transport the nested form as escaped JSON, and attach it to the matching membership occurrence. Keep exact `EQUAL` sites/interval boundaries; other statuses and missing statuses produce null exact positions with source ranges/status retained in description. Explicit source gene types now remain gene. Parser cache version advances 11 → 12. | `pypath/inputs_v2/parsers/reactome.py::_participant_molecular_form`, `_extract_participant_data`, `_flatten_participants`; `pypath/inputs_v2/reactome.py::_participant_molecular_form`, `reactions_schema` |

MITAB classification uses an explicit source-label crosswalk, including the
spellings in the local SIGNOR causalTab export. Unknown feature labels remain
source descriptions; substring similarity is not used to assign biology.
Feature parsing follows the
[HUPO-PSI MITAB 2.7 specification](https://github.com/HUPO-PSI/miTab/blob/master/PSI-MITAB27Format.md).

## Actionable remaining detail

1. **Reactome controllers and family members.** `_iterate_controls` builds
   controller/member dictionaries independently of `_extract_participant_data`.
   Their BioPAX feature lists do not yet reach forms. Extend that path while
   preserving each member's own physical-entity URI and state. Do not union
   distinct alternative member forms into a controller form. `_extract_xrefs_from_props`
   currently keeps UniProt but lacks explicit RefSeq/Ensembl transcript handling;
   add validated source namespace mappings. Reaction parsing also excludes
   transcription/translation in `_is_transcription_or_translation`, so transcript
   knowledge requires a source-model expansion. Non-`modificationType` features,
   variant sequence changes, and `notFeature` assertions need distinct modeling;
   they must not be guessed from labels or treated as ordinary modifications.

2. **Recon3D source product selectors.** `parsers/recon3d.py::_strip_isoform`
   removes `_ATn`-style source suffixes. `_parse_gene_rule::atom`, `_parse_genes`,
   `_parse_catalysis` and `_parse_enzyme_complexes` project/deduplicate these IDs.
   Reactions retain the original `gene_reaction_rule` in raw payloads; standalone
   gene rows retain only the stripped `entrez_id` and name. Keep original source
   selector IDs alongside gene projection and preserve separate GPR alternatives
   before deduplication. Establish what each source `_ATn` refers to and an
   explicit mapping before populating protein/transcript keys or claiming that
   `_AT1` is UniProt isoform 1. Current `recon3d.py` intentionally emits gene
   members and gene clauses; source selector retention must not change those to
   protein assertions without evidence. No guessed reconstruction was added.

3. **BindingDB sequence identity.** `_BINDINGDB_COLUMNS` projects primary chain
   accessions/names but omits target-chain sequences. Preserve supported sequence
   identifiers or establish a versioned/content identity for explicit chain
   sequences, then attach any source-reported mutation to the correct chain.
   Both DuckDB and CSV parser paths must retain identical information.

4. **ChEMBL assay-specific variants.** `parsers/chembl.py` selects
   `component_sequences` accession/type/description and aggregates component
   accessions with `GROUP_CONCAT`. Its current SQL does not select explicit
   variant records or sequence-version/coordinate references. Verify the archive
   schema and carry source variant/assay links as per-assay observations before
   changing forms. A target's component catalogue must not assert that every
   assay used every listed transcript/isoform or combine distinct assay variants.

5. **UniProt catalogue features.** `uniprot.py::_protein_description_pairs`
   deliberately extracts only `/note` text from Mutagenesis/Transmembrane,
   stripping feature coordinates from serving annotations. The selected TSV
   fields remain in parsed raw payloads. Add individually addressable reference
   sequence/features with versioned coordinates and provenance; do not populate
   one occurrence form with all annotated mutagenesis/PTM alternatives.

6. **miRBase products.** `mirbase.py::_parse_entry` retains precursor `MI…`,
   mature `MIMAT…` and explicit maturation links, but ignores feature positions
   and sequence lines. Preserve source sequences and precursor coordinate ranges
   if selected for this model. Mature products already have distinct source
   namespaces; do not replace them with inferred protein or transcript mappings.

7. **BRENDA.** `parsers/brenda.py::RECORDS_ENABLED` selects enzyme/protein,
   activator/inhibitor/cofactor, kinetic and reference records; other source
   sections are dropped. Its protein→EC mapper does not attach mutation/PTM
   statements or per-assay state. Add source sections with verified semantics
   and explicit protein-record linkage before attaching variants to enzyme data.

8. **Other resources.** Validate primary accession fixtures across the inventory
   below, including versioned ENSP/ENST/RefSeq and isoforms within complex members.
   Expanding a current source projection requires a matching resource fixture,
   paired evidence test and pilot build. Do not infer molecular forms from gene
   symbols, expression datasets, broad disease mutation predicates, or labels.

## Complete resource inventory

Paths below are relative to the audited pypath checkout. `I` means identifiers
remain available to generic capture; `D` means detailed occurrence forms are
implemented; `P` means a parser/projection loses or flattens potentially relevant
information; `N` means the selected representation has no molecular participant
form to construct. These describe the current code, not every upstream export.

| Resource module | Coverage / exact evidence | Next action |
| --- | --- | --- |
| `bindingdb.py` | D/P: primary chain columns and multichain members; narrow parser omits sequence columns. | Retain explicit sequence identities and validated per-chain mutations. |
| `brenda.py` | I/P: protein `UniProt`→EC associations; `parsers/brenda.py::RECORDS_ENABLED` filters sections. | Add verified state sections with protein-record linkage. |
| `cellchat.py` | I: gene symbols, protein subunits/cofactors; `parsers/cellchat.py` normalizes gene lists. | Keep protein typing without inventing isoforms. |
| `cellinker.py` | I: selected protein/metabolite ligand-receptor, enzyme/transporter participants. | Check primary product IDs when supplied; no detailed form field selected. |
| `cellphonedb.py` | I: partner IDs and complex `interactors`; primary protein IDs survive. | Verify complex member isoforms; no detailed feature mapping. |
| `chebi.py` | N: chemical identities, ontology links and identifier translations. | No protein form inference from chemical modifications/names. |
| `chembl.py` | I/P: target component accessions/types; `parsers/chembl.py` projects and aggregates components. | Audit assay-linked variant tables and coordinate references. |
| `chemont.py` | N: ontology terms. | Preserve ontology semantics. |
| `connectomedb.py` | I: Ensembl IDs and symbols for protein-typed ligand/receptor participants. | Keep gene-vs-product namespace distinction; no detailed feature selected. |
| `corum.py` | I: member UniProt regex retains `-n`; `parsers/corum.py` yields TSV rows. | Validate member occurrence forms; source comment fields remain raw. |
| `drugcentral.py` | I: target `accession`, `swissprot`, gene symbol; `ACT_COMMENT` is source context. | Do not infer variants from comments. |
| `foodb.py` | N: foods, compounds, content measurements. | No macromolecular form asserted. |
| `go.py` | N: ontology terms. | Catalogue semantics only. |
| `guidetopharma.py` | I: explicit ligand/target types and UniProt/Ensembl/Gene IDs. | Preserve only asserted product identifiers; no mutation/PTM field selected. |
| `hmdb.py` | N: chemical entities and mappings; parser contains metabolite annotations. | Protein cross-references are not participant states. |
| `hpo.py` | I: gene/HPO annotations with `ncbi_gene_id`. | Keep gene-level assertions; no protein expansion. |
| `icellnet.py` | I: gene-symbol protein participants and multigene complexes. | No sequence specificity is supplied by selected identifiers. |
| `imm1415.py` | N: metabolites and chemical cross-references. | No macromolecular form asserted. |
| `intact.py` | D: paired participant forms from primary MITAB IDs/features. | Extend only approved feature terms; preserve unrecognized features raw. |
| `lipidmaps.py` | N: lipid chemical identities and molecular descriptors. | Lipid nomenclature is chemical identity, not a protein PTM assertion. |
| `macdb.py` | N: metabolite–trait associations. | No product form inferred from trait names. |
| `mebocost.py` | I: protein/enzyme/sensor participants represented by gene names. | No sequence/PTM field selected. |
| `metatlas.py` | I: Human-GEM ENSG GPR clauses; parser documentation confirms no isoform suffix. | Preserve gene rules; do not infer ENST/ENSP products. |
| `mirbase.py` | I/P: distinct precursor/mature accessions and maturation links; coordinates/sequence omitted. | Retain mature sequence and precursor coordinate context explicitly. |
| `mondo.py` | I: HGNC gene targets and mutation-class disease predicates. | Broad mutation predicates do not identify a particular variant. |
| `mrclinksdb.py` | I: protein receptors/transporters and complex members from IDs/symbols. | Validate primary accession specificity; no detailed feature selected. |
| `neuronchat.py` | I: contributor and receptor-subunit genes in parser records. | No transcript isoform inferred from transcriptomics descriptions. |
| `nichenet.py` | I: ligand/receptor symbols and source/database fields. | No detailed occurrence form selected. |
| `omnipath_ontology.py` | N: ontology vocabulary. | Metadata only. |
| `pfocr.py` | I: figure-derived genes/chemicals selected by parser identity lists. | No specific mutation inferred from figure/text context. |
| `phenol_explorer.py` | N: foods, chemicals, contents and citations. | No macromolecular participant form selected. |
| `psi_mi.py` | N: PSI-MI ontology terms. | Feature concepts are vocabulary, not observed modified proteins. |
| `ptfi.py` | N: foods and chemical membership. | No macromolecular participant form selected. |
| `rampdb.py` | I/N: identifier translations and gene/chemical/pathway records. | Catalogue mappings do not assert occurrence isoforms. |
| `reactome.py` | D/P: reaction participant features fixed; controls use a separate parser path. | Extend controller/member state and transcript reference paths explicitly. |
| `recon3d.py` | I/P: gene-only GPR assertions; `_ATn` stripped/deduplicated by parser. | Preserve source selectors before projection; validate selector-to-product mappings. |
| `refmet.py` | N: chemical reference identities/descriptors. | No macromolecular form selected. |
| `rhea.py` | I: reactions, chemical participants and explicit UniProt enzyme memberships. | Reference reactions do not assert a particular PTM/variant. |
| `signor.py` | D: paired causalTab features and nested isoform member forms. | Extend approved feature vocabulary with real source fixtures. |
| `slctables.py` | I: protein-typed gene names and transporter ontology links. | Gene identifiers do not imply a particular protein sequence. |
| `stitch.py` | I: protein STRING IDs and chemical interactions. | Keep native protein fallback; no exact sequence/PTM field selected. |
| `swisslipids.py` | N: lipid chemical ontology/identifier hierarchy. | No protein PTM inference from lipid classes. |
| `tcdb.py` | I: `uniprot` and `transporter_uniprot` accessions. | Validate specific protein IDs; no detailed state field selected. |
| `uniprot.py` | I/P: reference accessions, biological sequence, feature narratives. | Model catalogue features individually; do not assert combined states. |
| `wikipathways.py` | I: RDF endpoint IDs; `parsers/wikipathways.py::_IDENTIFIER_PATTERNS` retains full identifier segments. | Verify explicit RDF state constructs before adding feature extraction. |

## Validation and deployment boundary

Focused tests are `test/test_molecular_forms.py`, `test/test_signor_identifiers.py`
and `test/test_tabular_execution.py` in the isolated pypath checkout, plus the
shared core molecular-form/Parquet tests. Fixtures exercise paired features,
standalone evidence, source type preservation, ambiguous sequence references,
unknown/fuzzy positions, explicit cache dependencies, nested SIGNOR members,
BindingDB complex chains, and Reactome feature transport. These tests make no
network calls and perform no full resource rebuild.
After reconciliation, all **131 focused pypath/source-context/measurement tests**
pass, including the existing remote branch's statement-identity and Parquet
source-context regressions. The shared core suite also passes **55 tests** with
Ruff passing. Pypath's new MITAB helper passes its configured Ruff
checks; existing unrelated unused imports/variables and lambda style findings
remain in older pypath files.

SIGNOR and IntAct have passed full isolated rebuilds, raw/projected output
validation and actual API/client checks; the
[implementation report](molecular-form-implementation-20261004.md) records results.
BindingDB/Reactome additional changes have focused fixture coverage only; full
resource build, output comparison and performance/storage validation are pending.
The final pypath commits were published as a fast-forward on the existing
`parquet-migration` branch at the user’s request. Root dependency declarations
and `uv.lock` pin `9d61c6e6cf1c9c0355d1287758ade6df2d477455`.
