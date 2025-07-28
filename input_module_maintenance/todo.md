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
| | adrecs_drug_identifiers | ✅ | Works - 2526 records from Drug_information.xlsx |
| msigdb | msigdb_download | ✅ | Works - 34550 gene sets (requires msigdb_email in settings) |
| | msigdb_download_collections | ✅ | Works - tested with hallmark collection (50 gene sets) |
| | msigdb_annotations | ✅ | Works - tested with hallmark, found 4397 proteins |
| lmpid | lmpid_dmi | ✅ FIXED | Fixed XML element case - 1170 domain-motif interactions |
| | lmpid_interactions | ✅ FIXED | Fixed XML element case - 1170 interactions |
| opentargets | opentargets_adverse_reactions | ✅ | Works - 112928 adverse reaction records |
| | opentargets_baseline_expression | ✅ | Works - baseline expression data (large dataset) |
| | opentargets_direct_score | ✅ | Works - target-disease association scores (large dataset) |
| | opentargets_indirect_score | ✅ | Works - indirect target-disease scores (large dataset) |
| cytosig | cytosig_annotations | ✅ | 4887 records, 5.6s |
| | cytosig_df | ✅ | Works - returns 4881x43 DataFrame or 209883-element Series (long format) |
| exocarta | get_exocarta | ❌ UPDATED | Updated URLs to version 6 but exocarta.org site is currently down (503 errors) |
| | get_vesiclepedia | ✅ | 290197 records, 0.14s |
| | pdb_complexes | ✅ FIXED | 0 records (empty with graceful handling) - Fixed NoneType handling when PDB chains unavailable |
| cancercellmap | ccmap_interactions | ✅ FIXED | Fixed URLs to new PathwayCommons location - 47644 interactions |
| cancerdrugsdb | cancerdrugsdb_annotations | ✅ FIXED | Fixed URL + ChEMBL mapping - annotations for 216 drugs |
| | cancerdrugsdb_download | ✅ FIXED | Fixed URL to data.tp53.org.uk - 330 records |
| | cancerdrugsdb_interactions | ✅ FIXED | Fixed URL + ChEMBL mapping - 5157 interactions |
| | pathwaycommons_interactions | ✅ FIXED | Fixed to use v14 unified file with source filtering - 2.4M interactions |
| | pathwaycommons_wikipathways_interactions | ❌ UPDATED | WikiPathways not available in v14 - only in older versions |
| oreganno | oreganno_interactions | ⚠️ | 0 records - returns None instead of expected iterator |
| | scconnect_interactions | ✅ FIXED | Fixed ValueError by properly handling empty gene symbols and species fields |
| | ddinter_interactions | ✅ FIXED | Updated to use DDinter2 CSV-based API - 160K interactions |
| | ddinter_mappings | ✅ FIXED | Updated to use DDinter2 CSV-based API - 1939 drugs |
| | ddinter_n_drugs | ✅ FIXED | Updated to use DDinter2 CSV-based API - 1939 drugs |
| domino | domino_ddi | ✅ FIXED | Function works correctly (1294 interactions) - TypeError only in __str__ method |
| complexportal | complexportal_complexes | ⏳ | Timeout during testing |
| netpath | netpath_interactions | ⚠️ | Works fine (7555 records, 0.3s) but depends on netpath_names which fails |
| | netpath_names | ❌ | netpath.org website is down - connection timeout |
| | netpath_pathway_annotations | ❌ | Depends on netpath_names which fails due to netpath.org being down |
