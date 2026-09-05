# Communication and pharmacology inputs: native Biolink migration

The eleven modules below retain separate download, source parsing and mapping stages. Biological records use native Biolink 4.4.4 classes/slots/enums and the shared `Namespace` registry; resource/license/update enums remain configuration metadata. Automatic mapper caches apply through the shared builder, without resource-specific cache decorators.

Names are identifiers under the user's explicit contract: source names use `Namespace.NAME`, and synonyms/abbreviations use `Namespace.SYNONYM`. They are not annotations. Actual accessions retain their position before supplemental names; accessionless source entities can use the name as unresolved identity. Chemical formulas, publications, taxonomy and measurements remain native attributes. Display names and source metadata do not establish identifier equivalence or biological effects. This user-requested boundary overrides the earlier modeling document's advice to use `slots.name`/`slots.synonym` annotations.

| Module | Modeling and corrections |
| --- | --- |
| CellChat | Explicit protein symbols and complex membership use `Protein`/`MacromolecularComplex`. Exact lookup in the source `geneInfo.Symbol` column handles older species files without symbol columns. Unknown labels remain `NamedThing`; labels no longer generate guessed complex members. Cofactor collections use `NamedThing` grouping rather than asserting physical complexes. Cofactor stimulation/inhibition becomes `affects` with increased/decreased direction; unknown effects raise an error. PMID/PMC references become publications. |
| Cellinker | Protein and chemical entities; ligand/receptor observations use `interacts_with`. Enzyme/transporter-metabolite tables use `associated_with` rather than fictitious reaction/transport event entities. Missing accessions may use an explicitly supplied protein name as unresolved identity. Supporting database names use `supporting_data_source`. PubChem CID and SID remain separate identifier systems. |
| CellPhoneDB | Protein accessions and source complexes remain distinct. The existing source convention for synthetic metabolite system names is retained: these are chemicals, not complexes of their biosynthesis/supporting proteins; chemicals have no species taxonomy. The fifth source complex subunit is now retained. Publication references use native slots; source names use generic name identifiers. |
| ConnectomeDB | Source species crosswalk yields NCBITaxon CURIEs. HGNC/Ensembl gene identifiers retain their distinct domains. Source evidence strings do not imply directness, and words within location prose do not assert localization. |
| Guide to Pharmacology | Explicit classes, exact action-to-mechanism crosswalks, native effect qualifiers and quantitative values, described below. Target accessions have a distinct namespace from ligand accessions. PubChem substance and compound accessions are separate. Human reference aliases are not attached to nonhuman assay targets. A missing reference-table record can retain its explicit target accession; missing gene aliases do not imply a gene family. |
| ICELLNET | Explicit multi-subunit participants remain complexes with constituent proteins; missing interaction participants are excluded. Source spreadsheet-specific gene repair remains local. PMID references use `publications`; generated pair labels use generic name identifiers, not invented resource accession systems. |
| iMM1415 | Chemical identifiers, names and SBO topics remain structured. KEGG drug identifiers are no longer mislabeled as KEGG compounds. Missing optional annotation dictionaries parse safely. The compartment label does not become an inferred GO localization. |
| MEBOCOST | HMDB chemicals and protein symbols use `interacts_with`, with taxonomy on proteins. Publications and supporting databases are distinct native attributes. Actual free-text evidence comments may remain descriptions. |
| MRClinksDB | Chemical identifiers/formulas, receptor proteins/complexes and publications use native slots. CID/SID parsing examines the explicit prefix and preserves both namespaces; SID-first values cannot become compound identifiers. Transporter rows become protein–chemical associations, without fabricated transport event nodes. |
| NeuronChat | Explicit small-molecule tokens use the existing source crosswalk; their chemicals do not inherit organism taxonomy. Receptor subunit groups retain membership. Descriptive interaction classes do not become arbitrary predicates or invented direction qualifiers. |
| NicheNet | Explicit protein symbol pairs with species context use `interacts_with`; missing participants are excluded. Pair labels use generic name identifiers, not invented resource accession systems. |

## Quantitative pharmacology

`pIC50`, `pEC50`, `pKi` and `pKd` use their native Biolink slots with numeric `QuantityValue` data. **All reported high, low and median values remain separate**, through `Measurement.source_field` provenance such as `pIC50 Affinity High`. The source's logarithmic values are not converted to concentrations. `pA2`/`pKB`, for which the pinned schema lacks corresponding slots, remain numeric under `has_quantitative_value` with the source measurement name in `source_field`; they are not relabeled as another affinity measure.

The original linear affinity high/low/median fields are also retained as numeric values in nM with the source IC50/Ki/etc. endpoint name and the original comparison operator. A comparator on a linear value is not blindly copied onto its logarithmic counterpart. Original assay descriptions and action comments are legitimate prose descriptions.

Exact source actions use `CausalMechanismQualifierEnum`, including agonism, partial/biased agonism, antagonism, inhibition, inverse agonism and binding. Broader source text is not interpreted by substring. Binding uses the nondeprecated `directly_physically_interacts_with`; the obsolete `binds` predicate is not emitted. Full agonism uses the broader native agonism mechanism; the source retains its more specific wording. Unspecified aspects and ambiguous net directions are omitted.

## Remaining structured-field gaps

No biological metadata is converted to a generic labeled description. Supporting resource labels have a faithful native slot; numerical affinity data has a typed representation. Some legacy attributes still require a separately reviewed model: free-text ClassyFire class hierarchies without ontology accessions, broad ligand/receptor role labels, source-specific signaling/pathway classifications, arbitrary subcellular-location prose, and Guide to Pharmacology's general approval/labeling flags. `clinical_approval_status` describes a chemical–disease assertion, so it is not used for a general approved-drug flag. These source fields remain in parsed evidence, **not claimed as migrated serving attributes**. A source-specific ontology crosswalk or an explicit structured source-attribute contract is needed to expose them without inventing meanings. This limitation remains work to review; it is not a claim of zero information loss in the serving annotation projection.

## Verification on available local snapshots

All available raw snapshots below were parsed and every mapped row passed through `SilverExtractor`, without network downloads or retained generated outputs. Counts are parsed top-level rows, not distinct serving entities/relations; complex membership creates additional observations.

| Resource/datasets | Parsed rows | Mapper exclusions / limits |
| --- | ---: | --- |
| CellChat human interactions/cofactor interactions/complexes/cofactors | 3,233 / 9,196 / 338 / 32 | No mapper exceptions or exclusions |
| CellChat zebrafish same datasets | 2,774 / 734 / 267 / 62 | No mapper exceptions or exclusions |
| CellChat mouse | 1,000 / 1,000 / 337 / 33 cached parsed rows | Raw RDA is truncated (`Compressed data ended before the end-of-stream marker`); mapper tested against available parsed caches. First two datasets are capped caches, not full raw audits. |
| Cellinker human/mouse protein interactions | 7,402 / 5,808 | No mapper exceptions or exclusions |
| Cellinker human/mouse metabolite interactions | 794 / 678 | No mapper exceptions or exclusions |
| Cellinker human/mouse enzymes | 1,962 / 1,962 | No mapper exceptions or exclusions; explicitly named but accessionless proteins remain unresolved |
| Cellinker human/mouse transporters | 244 / 244 | No mapper exceptions or exclusions |
| CellPhoneDB interactions/complexes | 2,911 / 362 | No mapper exceptions or exclusions |
| ConnectomeDB all species | 45,768 | No mapper exceptions or exclusions |
| Guide to Pharmacology | 24,207 | 256 lack any molecular target accession (including organism-level assays); 23,951 mapped. One target absent from the cached target metadata table is retained by its explicit accession. |
| ICELLNET interactions/complexes | 1,669 / 172 | No mapper exceptions or exclusions after parser participant filtering |
| iMM1415 | Representative fixtures only | Configured local `bigg/iMM1415.json` snapshot unavailable |
| MEBOCOST human/mouse | 793 / 790 | No mapper exceptions or exclusions |
| MRClinksDB human/mouse interactions | 794 / 678 | No mapper exceptions or exclusions |
| MRClinksDB human/mouse transporters | 244 / 244 | Three mouse rows lack protein accessions; their explicit enzyme names now preserve unresolved protein identity |
| NeuronChat human/mouse | 190 / 183 | No mapper exceptions or exclusions |
| NicheNet human/mouse | 4,986 / 5,668 | No mapper exceptions or exclusions |

These are mapping/extraction audits, not canonical resolver audits or newly published immutable resource versions. Identifier completeness, unsupported source classifications and the missing/truncated snapshots are not hidden behind a successful import. `test/test_inputs_v2_communication_migration.py` covers every module through representative fixtures, native identifier domains, negative biological assertions, all quantitative statistics/comparators, and participant/member completeness.

Validation: **42 regression tests passed**; all eleven owned modules pass the undefined/unused-name checks, and the changed files pass whitespace validation. Temporary audit scripts and generated results were removed; only this document and durable tests remain.
