# PyPath Input Module Testing Status Report
Generated: 2025-07-03

## Summary
- Total modules: 190
- Total functions: 509
- Testing progress: 0/190 modules tested

## Module Status

### Legend
- ✅ Working: Module functions correctly
- ⚠️ Warning: Module works but has issues (e.g., slow, partial data)
- ❌ Failed: Module has errors
- 🔄 In Progress: Currently testing
- ⏳ Pending: Not yet tested

## Detailed Results

### Batch 1: transmir through signor (Tested: 2025-07-03)

| Module | Functions | Status | Notes |
|--------|-----------|--------|-------|
| transmir | transmir_interactions | ✅ | 2678 records, 1.1s |
| lncdisease | lncdisease_interactions | ✅ | 478 records, 0.4s |
| icellnet | icellnet_annotations | ✅ FIXED | 1194 records - Fixed column names |
| | icellnet_complexes | ✅ FIXED | 156 records - Fixed column names |
| | icellnet_interactions | ✅ FIXED | 1647 records - Fixed column names & generator test |
| ipi | ipi_uniprot | ✅ | 86754 records, 0.1s |
| ca1 | ca1_interactions | ❌ | Server returns 403 Forbidden - Science journal blocks automated downloads |
| wang | cui_interactions | ❌ | EMBOPress website blocks automated downloads - 'fname' attr fixed |
| | hsn_interactions | ✅ | 62937 records, 0.6s |
| | wang_annotations | ✅ FIXED | 1500 records - Added error handling for failed datasets |
| | wang_interactions | ✅ | 62937 records, 0.4s |
| ctdbase | ctdbase_relations | ⚠️ | Requires relation_type: 'chemical_gene', 'chemical_disease', 'disease_pathway', 'chemical_phenotype', or 'gene_disease' |
| | ctdbase_vocabulary | ⚠️ | Requires vocabulary_type: 'chemical', 'gene', 'disease', 'pathway', 'anatomy', or 'phenotype' |
| hpmr | hpmr_annotations | ✅ | 1141 records |
| | hpmr_complexes | ✅ | 0 records (empty result) |
| | hpmr_interactions | ✅ | 634 records |
| laudanna | laudanna_directions | ✅ | 77555 records, 0.7s |
| | laudanna_effects | ✅ | 63438 records, 0.6s |
| compleat | compleat_complexes | ✅ | 9692 records, 2.0s |
| | compleat_raw | ✅ | 9704 records |
| signor | signor_complexes | ✅ | 4889 records |
| | signor_enzyme_substrate | ✅ FIXED | 12669 records - Fixed with hardcoded AA mapping |
| | signor_interactions | ✅ | 93928 records, 0.2s |
| | signor_pathway_annotations | ✅ | 775 records, 37.7s (slow!) |
| | signor_pathways | ✅ | 2 records |
| | signor_protein_families | ✅ | 91 records |
| dorothea | dorothea_full_raw | ❌ | KeyError: 0 - pandas indexing issue |
| | dorothea_interactions | ✅ | 309009 records, 0.0s |
| | dorothea_interactions_old | ✅ | 0 records (empty result) |
| | dorothea_old_csv | ✅ | 0 records (empty result) |
| | dorothea_rda_raw | ❌ | KeyError: 0 - pandas indexing issue |
| | tfregulons_interactions | ✅ | 309009 records, 0.0s |
| | tfregulons_interactions_old | ✅ | 0 records (empty result) |

| uniprot | get_uniprot_sec | ✅ | 72333 records, 0.0s |
| | query_builder | ✅ | 1 records (empty query builder) |
| | swissprot_deleted | ✅ | 326 records, 0.3s |
| | uniprot_data | ✅ | 20420 records, 8.9s |
| | uniprot_families | ✅ | 14485 records, 9.9s |
| | uniprot_keywords | ✅ | 20420 records, 11.0s |
| | uniprot_locations | ✅ | 17249 records, 12.2s |
| | uniprot_ncbi_taxids_2 | ✅ | 27781 records, 0.0s |
| | uniprot_taxonomy | ✅ | 558215 records, 0.7s |
| | uniprot_tissues | ✅ | 10169 records, 11.5s |
| | uniprot_topology | ✅ | 5244 records, 50.1s |
| | uniprot_ncbi_taxids | ❌ | settings.context missing |
| | uniprot_query | ❌ | empty query causes unpack error |
| | trembl_deleted | ⚠️ | Requires user confirmation (5GB memory) |
| | uniprot_deleted | ⚠️ | Requires user confirmation (5GB memory) |
| | *_requires_args | ⚠️ | Functions requiring specific arguments |

| phosphatome | phosphatome_annotations | ❌ | Download blocked - likely 403 Forbidden like ca1 |

| progeny | progeny_annotations | ✅ FIXED | 18581 records, 23.3s - Fixed R data parsing with pyreadr |
| | progeny_raw | ✅ FIXED | DataFrame returned - Test framework sampling issue only |

| expasy | expasy_enzyme_classes | ✅ | 354 records, 0.2s |
| | expasy_enzymes | ✅ | 8405 records, 2.6s |

| biogrid | biogrid_all_interactions | ✅ | 8721 records, 80.5s (large dataset) |
| | biogrid_interactions | ✅ | 7409 records, 0.6s |

| celltypist | celltypist_annotations | ✅ | 460 records, 1.4s |