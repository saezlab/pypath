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

### Summary for Batch 1:
- **Total functions tested**: 26
- **Successful**: 17 (65.4%)
- **Failed**: 8 (30.8%)
- **Warnings**: 1 (3.8%)

### Common Issues Found:
1. **Missing required arguments**: ctdbase functions need specific parameters
2. **File format changes**: icellnet expects 'PubMed ID' column
3. **Science download issues**: ca1 and wang_annotations have file handling errors
4. **Attribute errors**: Various modules have changed APIs
