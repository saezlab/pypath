# Chemical and food input migration

The thirteen resource modules now choose native Biolink classes, slots and qualified predicates. Identifiers use `Namespace`. Following the explicit user preference, source names are `Namespace.NAME` identifiers and synonyms/abbreviations are `Namespace.SYNONYM` identifiers, not annotations. They are appended after existing accessions so the source identifier ordering is preserved. Formulas remain native chemical-formula attributes. Download/preparation, source parsing and biological mapping remain separate. Declarative builders use the shared implicit caches; there are no new source-specific mapping caches.

## Quantitative observations

BindingDB, ChEMBL and DrugCentral retain **Ki, Kd, IC50, EC50, kon and koff** using the exact [BAO endpoint definitions](https://github.com/BioAssayOntology/BAO/blob/master/vocabulary/bao_vocabulary_result.owl): respectively `BAO:0000192`, `BAO:0000034`, `BAO:0000190`, `BAO:0000188`, `BAO:0000480` and `BAO:0000479`. A native Biolink `QuantityValue`, wrapped with source-field provenance, carries numeric value, units and comparison operator together. Thus `IC50 <=10 nM` remains a numeric, unit-bearing bound, not a narrative or an inhibition claim. The original operator also preserves `<=`, `>=` and approximate comparisons beyond Biolink's three binary-relation enum values.

Food concentrations use [`PATO:0000033`, concentration of](https://www.ebi.ac.uk/ols4/ontologies/pato/classes?iri=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FPATO_0000033). The original per-member field distinguishes reported value, mean, median, minimum, maximum and standard deviation; each retains its corresponding unit and member. This is source-statistic provenance, not an inferred conversion between statistics. Unit strings are retained without numerical conversion; they are not claimed to be globally standardized UCUM values. STITCH confidence values use `has_confidence_score` with separate source fields for combined and action scores; ChEMBL's assay score uses `chembl_confidence_score`.

BindingDB and ChEMBL `pchembl_value`, BindingDB `pH` and `Temp (C)`, and Human-GEM `lower_bound` / `upper_bound` use native `has_quantitative_value` with exact source-field provenance. Heterogeneous pChEMBL values are not mislabeled as pXC50. pH and pChEMBL have no asserted unit; Celsius uses UCUM `Cel`. Human-GEM preserves explicit `flux_units` when present; the audited YAML does not declare flux units, so these are left unspecified rather than assumed. The parser now retains the two numeric bounds instead of only their direction projection.

Malformed numbers, non-finite values and unsupported comparison operators omit **only that quantitative annotation**, retaining the record and source evidence. Other assay endpoints, classifications, flags and experimental context without a reviewed representation remain in parsed source evidence; they are not dumped into `description`. Consequently, their presence in source evidence does not imply a corresponding typed serving attribute. Only genuine source prose uses `description`. Existing reviewed `chemrof:mass`, `chemrof:monoisotopic_mass` and `chemrof:charge` properties retain chemical properties; formulas use `has_chemical_formula`.

The companion serving change adds typed annotation quantities and entity-linked payloads for standalone entity records. Resource rows with no identifiable relation endpoint remain exclusions, rather than fabricated graph edges. These changes require rebuilt resource versions; old releases are not silently reinterpreted.

## Resource decisions

| Module | Representation and corrections |
|---|---|
| BindingDB | Measured chemical–target association; potency alone does not assert inhibition. PubChem substances and compounds use distinct namespaces. Explicit multichain targets become complexes with `has_part` protein members, preserving chains 2–50 in the parser instead of assigning the entire measurement to chain 1. Unknown target-organism text does not become a guessed taxon. Bare numbers in ligand names do not become ChEBI aliases. |
| ChEMBL | Source-declared molecule/target classes; explicit reviewed action codes become `affects` with direction, explicit binding agents become physical interactions, and other reported associations remain broad. Measurement types never become predicates. Ensembl gene/transcript/protein accessions retain their correct namespaces. A selectivity group does not become an inferred protein family. |
| DrugCentral | Structure-reference preparation moved out of the mapper. Multiple target accessions form an unresolved group identified by its source name and retaining all protein members; words such as “complex” in its name do not classify an assembly. Exact reviewed action codes provide direction. |
| FooDB | Food `has_part` chemical, with measurements attached to each membership. Compound identifiers are projected separately and shared by FooDB accession. Nonpositive taxon placeholders are not NCBITaxon identifiers. |
| HMDB | Chemical entities and explicit ChemOnt associations. Localization and disease text does not create inferred biological edges. Translation namespace matches the graph namespace. The snapshot audit below deliberately did not download/join ChemOnt; association coverage is fixture-tested, not claimed from that audit. |
| LIPID MAPS | Chemical entities with real structure identifiers, formula and monoisotopic mass. Classification labels are retained in source evidence. Translation uses the same namespace as graph input. |
| Human-GEM / metatlas | Model reaction `has_input`/`has_output` edges, source stoichiometry, and associated genes identified by ENSG. EC codes are topics, not reaction aliases. Model gene-rule AND groups are unresolved gene groups, not inferred physical protein complexes; a flattened gene list does not assert that every gene independently catalyses the reaction. Bounds are numeric quantities; compartments remain source context. |
| Phenol-Explorer | Food `has_part` chemical with aligned quantitative observations. Numeric food and compound IDs have different namespaces, preventing cross-domain collisions. Member publication lists are split into valid PMID CURIEs. |
| PTFI | Food/specimen composition, formula and aligned quantities. Preserve explicit source analyte IDs when no PTF accession is embedded in the name. Formula parsing requires the documented formula/accession name pattern; ordinary names no longer produce fabricated partial formulas. |
| RaMP-DB | A pathway association connects separate analyte and pathway nodes; a pathway ID is no longer an alias of its analyte. Gene records use Gene. Release discovery is deferred until download, so importing a mapper never requests the network. |
| RefMet | All records use ChemicalEntity without substring-based lipid inference. Corrected `pubchem_cid` and stripped source-header whitespace so RefMet primary identifiers survive parsing. |
| STITCH | Preserve explicit acting endpoint orientation, including protein→chemical. Direction qualifiers require source orientation; unoriented statements do not acquire a causal direction. Ensembl protein namespace is explicit. Source modes that do not justify a more specific predicate remain associations. |
| SwissLipids | Chemical entities; source parent hierarchy becomes `subclass_of`. Placeholder InChI values are excluded. Translation uses the graph namespace. |

## Verification and scope

33 offline tests exercise every module through `SilverExtractor`, including unknown actions, endpoint direction, misleading complex names, multiple members, taxonomy/publications, absence of invented localization, distinct pathway identity, formula negatives, typed IC50 bounds, assay context quantities, Human-GEM bound parsing and preservation, malformed numeric values and per-member measurement alignment. Name-identifier tests additionally check SwissLipids aliases, accession ordering, nested chemical names, direct target-group names, and exclusion of a whole-complex name from unnamed chains. All thirteen modules import successfully. Polars was installed in the existing runtime for the FooDB parser.

Available local snapshots were read without downloading or publishing. These snapshot counts precede the final name-as-identifier adjustment; that adjustment was verified with the offline tests, without repeating large audits. Counts below are parsed rows processed through mapping and extraction, **not final globally consolidated entity/relation counts**. A cap checks a prefix, not a representative random sample or full release validation.

| Snapshot / dataset | Audited rows | Result / limit |
|---|---:|---|
| RefMet | 203,854, full | All map/extract; no errors. |
| DrugCentral interactions | 19,378, full | All map/extract. One source comparison operator `-` excludes its quantitative annotation only. |
| ChEMBL molecules | 1,000 / 1,422,531 | All map/extract. |
| ChEMBL activities | 1,000 / 3,459,458 | All map/extract. |
| ChEMBL mechanisms | 7,568, full | 6,990 map/extract; **578 lack target IDs** and cannot form relations. No extraction errors. |
| RaMP analyte / metabolite_class / analytehaspathway / pathway / source / chem_props | 1,000 each | All map/extract; underlying tables contain 463,257 / 572,582 / 1,355,476 / 122,936 / 1,051,927 / 289,754 rows respectively. |
| HMDB metabolites | 1,000 | All map/extract; no ChemOnt join in this audit. |
| LIPID MAPS lipids | 1,000 | All map/extract. |
| SwissLipids | 1,000 | All map/extract; 262 hierarchy observations. |
| PTFI | 223 local specimens, full directory | All map/extract; 261,488 member observations. Local directory coverage is not a claim of complete remote coverage. |
| BindingDB | 1,000 raw rows | All map/extract; parser normalization exercised without the optional potency filter. Multichain targets covered by explicit fixtures. |
| Human-GEM metabolites / reactions / AND groups | 4,165 / 12,931 / 182, full | All map/extract; 79,127 reaction-participant observations and 807 group-member observations. The local companion cross-reference table supplied the parser lookup. |
| FooDB compounds | 1,000 | All map/extract. Full food-content join not rerun; membership/quantities fixture-tested. |
| Phenol-Explorer foods | 459, full | All map/extract; 6,962 membership observations. |
| STITCH human | 1,000 parsed merged interactions | All map/extract; source parser performed link/action matching and its existing score filter. Mouse not audited. |

This verifies migration behavior and documents known exclusions; it is not a full resource rebuild, quantitative completeness audit, release publication or performance benchmark. Generated audit scripts and raw results are not retained in the repository.
