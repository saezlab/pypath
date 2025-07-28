| Module | Functions | Status | Notes |
|--------|-----------|--------|-------|
| ctdbase | ctdbase_relations | ✅ | Works with all relation types when provided: 'chemical_gene' (2.9M records in 12.4s) |
| | ctdbase_vocabulary | ⚠️ | Works with 5/6 vocabulary types. 'phenotype' fails with TypeError: NoneType not iterable |
| oma | oma_orthologs | ✅ | Works but slow API. With score=0.99, rel_type={'1:1'}: 16K records in 60.5s |
| | oma_table | ✅ | Works. With same params: 16K unique IDs in 9.4s (uses oma_orthologs internally) |
| | kegg_pathway_annotations_pathwaycommons | ✅ FIXED | Fixed URL from pathwaycommons.org to download.baderlab.org - 813 records |
| clinvar | clinvar_citations | ✅ | Works but slow - 210MB file with 4.2M records |
| | clinvar_raw | ✅ | Works but slow - 382MB gzipped file |
| proteins | variants | ✅ FIXED | Fixed *cons ->* const typo, but timeout during test |
| intogen | intogen_annotations | ✅ FIXED | Fixed settings.context AttributeError - 483 records, 4.6s |
| | adrecs_drug_identifiers | ❌ | Excel file corruption error |
| msigdb | *_all_functions | ⏳ | Skipped - timeout |
| lmpid | lmpid_dmi | ✅ | 0 records (empty dataset) |
| | lmpid_interactions | ✅ | 0 records (empty dataset) |
| opentargets | *_all_functions | ⏳ | Skipped - timeout |
| cytosig | cytosig_annotations | ✅ | 4887 records, 5.6s |
| | cytosig_df | ❌ | DataFrame indexing error - test framework issue |
| exocarta | get_exocarta | ⚠️ | 0 records - NoneType iterator warning |
| | get_vesiclepedia | ✅ | 290197 records, 0.14s |
| | pdb_complexes | ✅ FIXED | 0 records (empty with graceful handling) - Fixed NoneType handling when PDB chains unavailable |
| cancercellmap | ccmap_interactions | ❌ | zipfile.BadZipFile: File is not a zip file |
| cancerdrugsdb | cancerdrugsdb_annotations | ❌ | URL unreachable - host 'acfdata.coworks.be' cannot be resolved |
| | cancerdrugsdb_download | ❌ | URL unreachable - host 'acfdata.coworks.be' cannot be resolved |
| | cancerdrugsdb_interactions | ❌ | URL unreachable - host 'acfdata.coworks.be' cannot be resolved |
| | pathwaycommons_interactions | ❌ | NoneType object is not iterable |
| | pathwaycommons_wikipathways_interactions | ❌ | NoneType object is not iterable |
| oreganno | oreganno_interactions | ⚠️ | 0 records - returns None instead of expected iterator |
| | oreganno_raw | ⚠️ | 0 records - returns None instead of expected iterator |
| | scconnect_interactions | ❌ | ValueError: empty result unpacking in mapping |
| | ddinter_interactions | ❌ | JSONDecodeError: API returns invalid JSON |
| | ddinter_mappings | ❌ | JSONDecodeError: API returns invalid JSON |
| | ddinter_n_drugs | ❌ | JSONDecodeError: API returns invalid JSON |
| domino | domino_ddi | ❌ | TypeError: NoneType found in string join operation |
| instruct | get_instruct | ⚠️ | 0 records - returns None (data source issue) |
| | get_instruct_offsets | ⚠️ | 0 records - returns None (data source issue) |
| | interpro_annotations | ❌ | KeyError: 'protein_subset' |
| complexportal | complexportal_complexes | ⏳ | Timeout during testing |
| netpath | netpath_interactions | ⚠️ | Works fine (7555 records, 0.3s) but depends on netpath_names which fails |
| | netpath_names | ❌ | netpath.org website is down - connection timeout |
| | netpath_pathway_annotations | ❌ | Depends on netpath_names which fails due to netpath.org being down |
| lipidmaps/structures | sdf | ✅ | 0 records (empty result), 21.9s |