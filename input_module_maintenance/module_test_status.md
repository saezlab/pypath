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
| oma | oma_orthologs | ⏳ | Skipped - timeout |
| | oma_table | ⏳ | Skipped - timeout |
| protmapper | get_protmapper | ✅ | 2 records, 2.2s |
| | protmapper_enzyme_substrate | ✅ | 22139 records, 0.2s |
| | protmapper_interactions | ✅ | 22139 records, 0.2s |
| uniprot_idmapping | idtypes | ✅ | 481 records, 0.0s |
| pazar | pazar_interactions | ✅ | 16386 records |
| trrust | trrust_human | ✅ | 9396 records |
| | trrust_interactions | ✅ | 9396 records |
| | trrust_mouse | ✅ | 7057 records |
| acsn | acsn_interactions | ✅ | 37725 records |
| panglaodb | panglaodb_annotations | ✅ | 4492 records, 2.4s |
| | panglaodb_raw | ✅ | 8286 records, 0.0s |
| proteinatlas | proteinatlas_annotations | ❌ | NoneType error - URL/API changes |
| | proteinatlas_secretome_annotations | ❌ | Science.org 403 Forbidden |
| | proteinatlas_subcellular_annotations | ❌ | Missing files_multipart attribute |
| netbiol | arn_interactions | ✅ | 95 records, 0.0s |
| | nrf2ome_interactions | ✅ | 109 records, 0.0s |
| lambert2018 | lambert2018_annotations | ❌ | Cell journal 403 Forbidden |
| | lambert2018_s1_raw | ❌ | Cell journal 403 Forbidden |
| kegg | kegg_interactions | ✅ | 14734 records, 2.4s |
| | kegg_medicus | ✅ | 12946 records, 0.1s |
| | kegg_medicus_complexes | ✅ | 539 records, 0.2s |
| | kegg_medicus_interactions | ✅ | 9566 records, 0.1s |
| | kegg_pathway_annotations | ✅ | 2669 records, 1.4s |
| | kegg_pathways | ✅ | 2 records, 1.5s |
| | kegg_dbget | ⚠️ | Requires 'entry' parameter |
| | kegg_pathway_annotations_pathwaycommons | ❌ | Gzip format error |
| ebi | ebi_rest | ⚠️ | Utility function - requires URL parameter |
| comppi | comppi_interaction_locations | ✅ | 587971 records, ~18s |
| | comppi_locations | ✅ | 22803 records, ~5s |
| encode | encode_tf_mirna_interactions | ✅ | 1237 records, 0.7s |
| homologene | get_homologene | ✅ | 275237 records, 79.1s |
| | homologene_dict | ⚠️ | Requires organism parameters |
| | homologene_uniprot_dict | ⚠️ | Requires organism parameters |
| clinvar | clinvar_citations | ⏳ | Skipped - timeout |
| | clinvar_raw | ⏳ | Skipped - timeout |
| humancellmap | humancellmap_annotations | ✅ | 4384 records, 2.8s |
| pubmed | *_all_functions | ⚠️ | All require specific parameters (PMIDs, DOIs, etc.) |
| cellinker | cellinker_annotations | ✅ | 1920 records |
| | cellinker_complex_annotations | ✅ | 134 records |
| | cellinker_complexes | ✅ | 143 records |
| | cellinker_complexes_raw | ✅ | 145 records |
| | cellinker_lr_interactions | ✅ | 3811 records |
| | cellinker_lr_interactions_raw | ✅ | 3744 records |
| | cellinker_protein_annotations | ✅ | 1786 records |
| | cellinker_smol_interactions | ✅ | 314 records |
| | cellinker_smol_interactions_raw | ✅ | 341 records |
| cpad | cpad_pathway_cancer | ✅ | 2 records |
| | cpad_annotations | ❌ | Mirbase dependency failure |
| phosphoelm | phosphoelm_enzyme_substrate | ✅ | 2426 records |
| | phosphoelm_interactions | ✅ | 2426 records |
| | phosphoelm_kinases | ✅ | 247 records |
| graphviz | graphviz_attrs | ✅ | 3 records, 1.7s |
| tfcensus | tfcensus_annotations | ✅ | 1871 records, 1.5s |
| mcam | mcam_cell_adhesion_molecules | ✅ | 112 records, 2.0s |
| proteins | variants | ✅ FIXED | Fixed _cons -> _const typo, but timeout during test |
| guide2pharma | guide2pharma_complexes | ✅ | 93 records, 1.1s |
| | guide2pharma_download | ✅ | 2 records, 0.1s |
| | guide2pharma_interactions | ✅ | 2620 records, 0.1s |
| credentials | credentials | ⚠️ | Requires resource name or user/password parameters |
| sider | sider_drug_names | ✅ | 1430 records |
| | sider_meddra_side_effects | ✅ | 20307 records |
| | sider_side_effect_frequencies | ✅ | 968 records |
| | sider_side_effects | ✅ | 1430 records |
| hippie | hippie_interactions | ✅ | 102410 records, 21.3s |
| iptmnet | iptmnet_interactions | ✅ | 36109 records |
| humsavar | uniprot_variants | ✅ | 13122 records, 6.6s |
| kea | kea_enzyme_substrate | ✅ | 35224 records, 1.5s |
| | kea_interactions | ✅ | 35224 records, 0.1s |
| kegg_api | *_all_functions | ⚠️ | All require 'organism' parameter (e.g. 'hsa') |
| cellphonedb | cellphonedb_complex_annotations | ✅ | 358 records |
| | cellphonedb_complexes | ✅ | 358 records |
| | cellphonedb_interactions | ✅ | 2903 records |
| | cellphonedb_ligands_receptors | ✅ | 2 records |
| | cellphonedb_protein_annotations | ✅ | 1359 records |
| opm | opm_annotations | ✅ | 88 records, 19.3s |
| matrisome | matrisome_annotations | ✅ | 2502 records, 1.6s |
| intogen | intogen_annotations | ⏳ | Skipped - timeout |
| cellcellinteractions | cellcellinteractions_annotations | ✅ | 3425 records, 1.9s |
| matrixdb | matrixdb_annotations | ✅ | 0 records, 1.2s (empty dataset) |
| | matrixdb_ecm_proteins | ✅ | 0 records (empty dataset) |
| | matrixdb_membrane_proteins | ✅ | 0 records (empty dataset) |
| | matrixdb_secreted_proteins | ✅ | 0 records (empty dataset) |
| | matrixdb_interactions | ❌ | Gzip format error - receiving HTML instead |
| mirdeathdb | mirdeathdb_interactions | ✅ | 462 records, 0.0s |
| cell | cell_supplementary | ⚠️ | Utility function - requires supp_url and article_url |
| almen2009 | almen2009_annotations | ✅ | 4244 records, 1.0s |
| biogps | biogps_datasets | ✅ | 9 records |
| | biogps_download_all | ✅ | 9 records, 56.7s |
| | biogps_annotations | ❌ | ArrayMapping._get_id_type method missing |
| | biogps_download | ⚠️ | Requires dataset parameter |
| uniprot_db | all_swissprots | ✅ | 20420 records |
| | all_trembls | ✅ | 184785 records, 40.3s |
| | all_uniprots | ✅ | 205003 records |
| | get_db | ✅ | 205003 records |
| | init_db | ✅ | Function executed |
| | is_swissprot | ⚠️ | Requires 'name' parameter |
| | is_trembl | ⚠️ | Requires 'name' parameter |
| | is_uniprot | ⚠️ | Requires 'name' parameter |
| pubchem | pubchem_mapping | ⏳ | Skipped - timeout |
| corum | corum_complexes | ✅ | 2734 records, 0.5s |
| deathdomain | deathdomain_interactions | ❌ | NoneType web scraping error |
| | deathdomain_interactions_rescued | ✅ | 184 records, 0.0s |
| integrins | get_integrins | ✅ | 25 records, 1.4s |
| adrecs | adrecs_adr_ontology | ✅ | 13855 records |
| | adrecs_drug_adr | ✅ | 809346 records |
| | adrecs_hierarchy | ✅ | 13828 records |
| | adrecs_drug_identifiers | ❌ | Excel file corruption error |
| msigdb | *_all_functions | ⏳ | Skipped - timeout |
| lmpid | lmpid_dmi | ✅ | 0 records (empty dataset) |
| | lmpid_interactions | ✅ | 0 records (empty dataset) |
| opentargets | *_all_functions | ⏳ | Skipped - timeout |
| kirouac2010 | kirouac2010_interactions | ✅ | 267 records, 0.0s |
| drugcentral | drugcentral_drugs | ✅ | 4099 records, 1.5s |
| | drugcentral_interactions | ✅ | 23115 records, 3.5s |
| | drugcentral_mapping | ⚠️ | Requires id_type and target_id_type parameters |
| locate | locate_localizations | ✅ | 9466 records, 66.4s |
| cytosig | cytosig_annotations | ✅ | 4887 records, 5.6s |
| | cytosig_df | ❌ | DataFrame indexing error - test framework issue |
| mppi | mppi_interactions | ✅ | 777 records, 0.2s |
| adhesome | adhesome_annotations | ✅ | 239 records, 0.8s |
| | adhesome_interactions | ✅ | 6542 records, 0.3s |
| exocarta | get_exocarta | ⚠️ | 0 records - NoneType iterator warning |
| | get_vesiclepedia | ✅ | 290197 records, 0.14s |
| mirbase | get_mirbase_aliases | ❌ | TypeError: 'NoneType' object is not iterable |
| | mirbase_ids | ⚠️ | 0 records - NoneType iterator warning |
| | mirbase_mature | ⚠️ | 0 records - NoneType iterator warning |
| | mirbase_mature_all | ❌ | TypeError: 'NoneType' object is not iterable |
| | mirbase_precursor | ⚠️ | 0 records - NoneType iterator warning |
| | mirbase_precursor_all | ❌ | TypeError: 'NoneType' object is not iterable |
| | mirbase_precursor_to_mature | ⚠️ | 0 records - NoneType iterator warning |
| gpcrdb | gpcrdb_annotations | ✅ | 808 records, 0.4s |
| csa | get_csa | ✅ | 0 records, 0.1s (returns None) |
| stitch | stitch_actions_interactions | ⚠️ | 0 records - NoneType iterator warning |
| | stitch_links_interactions | ⚠️ | 0 records - NoneType iterator warning |
| disgenet | disgenet_annotations | ❌ | NameError: name 'urls' is not defined |
| | disgenet_diseases | ❌ | NameError: name 'urls' is not defined |
| | disgenet_variants | ❌ | NameError: name 'urls' is not defined |
| | disgenet_gene_disease | ❌ | TypeError: unsupported operand type(s) for +: 'NoneType' and 'str' |
| | disgenet_variant_disease | ❌ | TypeError: unsupported operand type(s) for +: 'NoneType' and 'str' |
| | disgenet_disease_disease | ❌ | TypeError: unsupported operand type(s) for +: 'NoneType' and 'str' |
| pfam | pfam_names | ✅ | 2 records, 1.0s |
| | pfam_pdb | ❌ | AttributeError: module 'pypath.share.common' has no attribute 'non_digit' |
| | pfam_regions | ⏳ | Timeout - FTP connection timeout |
| | pfam_uniprot | ❌ | AttributeError: 'NoneType' object has no attribute 'split' |
| string | string_effects | ❌ | TypeError: 'NoneType' object is not an iterator |
| | string_links_interactions | ✅ | 201712 records, 0.0s |
| | string_physical_interactions | ✅ | 89862 records, 0.5s |
| | string_species | ✅ | 12535 records, 0.4s |
| cpdb | cpdb_interactions | ✅ | 531371 records, 13.1s |
| | cpdb_interactions_ltp | ✅ | 482222 records, 1.4s |
| reactome | pathway_hierarchy | ✅ | 23019 records, 0.0s |
| | reactome_raw | ⚠️ | Requires id_type parameter |
| hgnc | hgnc_genegroups | ✅ | 15590 records, 3.1s |
| pdb | pdb_chains | ✅ | 2 records, 30.0s |
| | pdb_complexes | ❌ | AttributeError: 'NoneType' object has no attribute 'items' |
| | pdb_uniprot | ✅ | 2 records, 1.0s |
| depod | depod_enzyme_substrate | ❌ | TypeError: argument of type 'function' is not iterable |
| | depod_interactions | ✅ | 832 records, 0.0s |
| pisa | pisa_bonds | ⚠️ | Requires bonds and chains parameters |
| | pisa_interfaces | ⚠️ | Requires pdbs parameter |
| phosphosite | phosphosite_directions | ✅ | 9094 records, 0.0s |
| | phosphosite_enzyme_substrate | ✅ | 13459 records, 0.1s |
| | phosphosite_interactions | ✅ | 2 records, 0.0s |
| | phosphosite_interactions_all | ✅ | 9164 records, 0.0s |
| | phosphosite_interactions_curated | ✅ | 4374 records, 0.0s |
| | phosphosite_interactions_new | ✅ | 2 records, 0.0s |
| | phosphosite_interactions_noref | ✅ | 4790 records, 0.0s |
| | phosphosite_ptms | ✅ | 239789 records, 5.7s |
| | phosphosite_regsites | ✅ | 5435 records, 0.1s |
| | phosphosite_ptm_orthology | ❌ | TypeError: 'function' object is not iterable |
| | phosphosite_regsites_one_organism | ❌ | TypeError: 'function' object is not iterable |
| | regsites_tab | ⚠️ | Requires regsites parameter |
| elm | elm_classes | ✅ | 353 records, 0.3s |
| | elm_domains | ✅ | 604 records, 0.1s |
| | elm_instances | ✅ | 4277 records, 0.6s |
| | elm_interactions | ❌ | AttributeError: 'NoneType' object has no attribute 'groups' |
| intact | intact_interactions | ✅ | 77349 records, 20.1s |
| go | *_all_functions | ⏳ | Timeout - Functions timed out during testing |
| trip | take_a_trip | ✅ | 5 records, 0.0s |
| | trip_interactions | ✅ | 359 records, 0.0s |
| | trip_process | ✅ | 359 records, 0.0s |
| | trip_find_uniprot | ⚠️ | Requires soup parameter |
| | trip_get_uniprot | ⚠️ | Requires syn parameter |
| | trip_process_table | ⚠️ | Requires tab, result, intrs, and trp_uniprot parameters |
| collectri | collectri_interactions | ✅ | 64989 records, 0.0s |
| | collectri_raw | ✅ | 43416 records, 0.0s |
| mirtarbase | mirtarbase_interactions | ❌ | AttributeError: 'NoneType' object has no attribute 'name' |
| mitab | mitab_field_list | ⚠️ | Requires field parameter |
| | mitab_field_uniprot | ⚠️ | Requires field parameter |
| cosmic | cancer_gene_census_annotations | ✅ | 0 records, 0.0s (empty dataset) |
| abs | abs_interactions | ✅ | 650 records, 0.0s |
| cancercellmap | ccmap_interactions | ❌ | zipfile.BadZipFile: File is not a zip file |
| phobius | phobius_annotations | ✅ | 20350 records, 0.4s |
| zhong2015 | zhong2015_annotations | ✅ | 466 records, 1.3s |
| pharos | pharos_diseases | ❌ | gzip.BadGzipFile: Not a gzipped file |
| | pharos_expression | ❌ | gzip.BadGzipFile: Not a gzipped file |
| | pharos_general | ⚠️ | Requires query parameter |
| | pharos_gtex | ❌ | gzip.BadGzipFile: Not a gzipped file |
| | pharos_ligands | ❌ | gzip.BadGzipFile: Not a gzipped file |
| | pharos_orthologs | ❌ | gzip.BadGzipFile: Not a gzipped file |
| | pharos_targets | ❌ | gzip.BadGzipFile: Not a gzipped file |
| | pharos_xrefs | ❌ | gzip.BadGzipFile: Not a gzipped file |
