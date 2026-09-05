# Legacy adapter resource inventory

Snapshot: 2026-09-05, current working tree (including existing uncommitted changes). This is a code inventory, not an availability or biological validation report. No legacy downloads were executed.

Companion: [mapping design](inputs-v2-legacy-adapter-mapping.md).

## Scope and counting

“Missing” means that the resource has no corresponding resource module in `pypath/inputs_v2`, after resolving the aliases below. Existing v2 resources are excluded from the main backlog even when a particular legacy dataset has not been reproduced. A resource can occur in multiple categories. Network labels and versioned datasets are retained rather than counted as independent upstream databases.

The inventory combines: all complex input registrations in `core/complex.py`; all enzyme–substrate registrations in `resources/data/resources.json`; every literal `NetworkInput` definition in `resources/data_formats.py`; and a recursive scan of legacy functions, supplemented with nonstandard structural-interaction entry points. Private parsers, annotations-only functions, ontology/ID mapping and general disease-association APIs are not treated as ordinary molecular network resources. Exceptions and related APIs appear below so they do not disappear behind this boundary.

Static references establish intended support, not that a function imports or an upstream endpoint works. Generated registry variants inherit the base definitions; the tables group these by input module and preserve distinct literal labels/entry points. This does not claim exhaustive runtime expansion of every configuration combination.

## 1. Complexes

Nine implemented source/version candidates across seven modules, plus two blocked/empty candidates. HPMR and 3DComplex are listed for completeness, not counted as productive candidates. scConnect needs an iterable-of-complexes normalizer; the others use dictionary values.

| Resource / version | Legacy entry point | Contract / qualification |
| --- | --- | --- |
| Compleat | [compleat.compleat_complexes](../pypath/inputs/compleat.py#L57) | Dictionary of Complex; `predicted` controls inclusion. |
| ComplexPortal | [complexportal.complexportal_complexes](../pypath/inputs/complexportal.py#L31) | Dictionary of Complex; use `return_details=False`; separately preserve details if requested. |
| Havugimana2012 | [havugimana.havugimana_complexes](../pypath/inputs/havugimana.py#L51) | Dictionary of Complex. |
| hu.MAP1 | [humap.humap1_complexes](../pypath/inputs/humap.py#L28) | Dictionary of Complex; source IDs may be absent. |
| hu.MAP2 | [humap.humap2and3_complexes](../pypath/inputs/humap.py#L54) | Pin `humap_version=2`; record `min_confidence`. |
| hu.MAP3 | [humap.humap2and3_complexes](../pypath/inputs/humap.py#L54) | Pin `humap_version=3`; record `min_confidence`. |
| PDB | [pdb.pdb_complexes](../pypath/inputs/pdb.py#L189) | Dictionary of Complex; preserve PDB identifiers and member counts. |
| SPIKE | [spike.spike_complexes](../pypath/inputs/spike.py#L195) | Dictionary of Complex; record `min_confidence`. |
| scConnect | [scconnect.scconnect_complexes](../pypath/inputs/scconnect.py#L137) | Set of Complex, not dictionary; specify organism. Not in core complex default list. |
| HPMR | [hpmr.hpmr_complexes](../pypath/inputs/hpmr.py#L336) | Function documents an empty result / ambiguous complex evidence; do not enable as productive complex source without review. |
| 3DComplex | [threedcomplex.threedcomplex_complexes](../pypath/inputs/threedcomplex.py#L34) | Stub: raises `NotImplementedError`; cannot be fixed by an adapter. |

Embedded `Complex` participants also occur in Baccin2019 and CollecTRI interaction outputs; these require the same participant converter but are not separate complex-download resources. Complex grouping must not be inferred from names or arbitrary collections.

## 2. Enzyme–substrate

Ten missing resources have explicit enzyme–substrate registrations. SIGNOR is the eleventh registered source and already has a v2 module. These registered outputs use dictionaries; configuration must retain the enzyme/substrate identifier domains rather than assume UniProt everywhere.

| Resource | Legacy entry point | Enzyme ID | Substrate ID / field | Species contract |
| --- | --- | --- | --- | --- |
| DEPOD | [depod.depod_enzyme_substrate](../pypath/inputs/depod.py#L77) | `uniprot` | `"uniprot"` | Accepts species selection |
| HPRD | [hprd.hprd_enzyme_substrate](../pypath/inputs/hprd.py#L69) | `genesymbol` | `[["refseqp", "substrate_refseqp"]]` | Human-only registry configuration |
| KEA | [kea.kea_enzyme_substrate](../pypath/inputs/kea.py#L93) | `uniprot` | `"uniprot"` | Human-only registry configuration |
| Li2012 | [li2012.li2012_enzyme_substrate](../pypath/inputs/li2012.py#L81) | `genesymbol` | `"genesymbol"` | Human-only registry configuration |
| MIMP | [mimp.mimp_enzyme_substrate](../pypath/inputs/mimp.py#L26) | `genesymbol` | `[["genesymbol", "substrate"], ["refseqn", "substrate_refseq"]]` | Human-only registry configuration |
| PhosphoNetworks | [phosphonetworks.phosphonetworks_enzyme_substrate](../pypath/inputs/phosphonetworks.py#L26) | `genesymbol` | `"genesymbol"` | Human-only registry configuration |
| PhosphoSite | [phosphosite.phosphosite_enzyme_substrate](../pypath/inputs/phosphosite.py#L46) | `uniprot` | `"uniprot"` | Accepts species selection |
| ProtMapper | [protmapper.protmapper_enzyme_substrate](../pypath/inputs/protmapper.py#L61) | `uniprot` | `"uniprot"` | Human-only registry configuration |
| dbPTM | [dbptm.dbptm_enzyme_substrate](../pypath/inputs/dbptm.py#L30) | `genesymbol` | `"uniprot"` | Accepts species selection |
| phosphoELM | [phosphoelm.phosphoelm_enzyme_substrate](../pypath/inputs/phosphoelm.py#L30) | `uniprot` | `"uniprot"` | Accepts species selection |

Additional candidate: [iptmnet.iptmnet_interactions](../pypath/inputs/iptmnet.py#L50) returns `IptmnetInteraction` namedtuples with `enzyme`, `substrate`, both isoforms, `ptm_type`, residue, position, score and references. It is registered as a network, not an enzyme–substrate input; an explicit field rename makes it eligible for the same adapter. Thus the enzyme–substrate backlog is **10 registered resources + iPTMnet**.

Important variants: PhosphoSite must use `raw=True` (`raw=False` constructs legacy DomainMotif objects); `strict=True` avoids the function's cross-species inference path. ProtMapper must use `interactions=False` to retain sites. `dbptm_enzyme_substrate_old` is a historical alternative, not another source. PhosphoSite PTM-only, regulatory-site and orthology tables are related datasets, not interchangeable enzyme–substrate records.

## 3. Interactions / networks

### Registered input definitions

72 missing legacy module groups occur in `NetworkInput` definitions below. Multiple definitions can intentionally use different filters, arguments, confidence levels or reference policies. The same endpoint pair in different datasets does not justify discarding its evidence. **Unresolved** indicates that the referenced function was not found by the static definition scan and must be checked/repaired before adaptation.

| Resource labels | Legacy module | Registered entry points | Registry evidence |
| --- | --- | --- | --- |
| ABS | abs | [abs.abs_interactions](../pypath/inputs/abs.py#L24) | [definitions](../pypath/resources/data_formats.py#L2133) |
| ACSN | acsn | [acsn.acsn_interactions](../pypath/inputs/acsn.py#L27) | [definitions](../pypath/resources/data_formats.py#L236) |
| Adhesome | adhesome | [adhesome.adhesome_interactions](../pypath/inputs/adhesome.py#L29) | [definitions](../pypath/resources/data_formats.py#L653) |
| AlzPathway | alzpathway | [alzpathway.alzpathway_interactions](../pypath/inputs/alzpathway.py#L38) | [definitions](../pypath/resources/data_formats.py#L1242) |
| Baccin2019 | baccin2019 | [baccin2019.baccin2019_interactions](../pypath/inputs/baccin2019.py#L32) | [definitions](../pypath/resources/data_formats.py#L3058) |
| BioGRID | biogrid | [biogrid.biogrid_interactions](../pypath/inputs/biogrid.py#L29) | [definitions](../pypath/resources/data_formats.py#L1055) |
| CA1 | ca1 | [ca1.ca1_interactions](../pypath/inputs/ca1.py#L27) | [definitions](../pypath/resources/data_formats.py#L449) |
| CancerCellMap | cancercellmap | [cancercellmap.ccmap_interactions](../pypath/inputs/cancercellmap.py#L27) | [definitions](../pypath/resources/data_formats.py#L1075) |
| CancerDrugsDB | cancerdrugsdb | [cancerdrugsdb.cancerdrugsdb_interactions](../pypath/inputs/cancerdrugsdb.py#L60) | [definitions](../pypath/resources/data_formats.py#L3259) |
| CellCall | cellcall | [cellcall.cellcall_interactions](../pypath/inputs/cellcall.py#L109) | [definitions](../pypath/resources/data_formats.py#L3159) |
| CellTalkDB | celltalkdb | [celltalkdb.celltalkdb_interactions](../pypath/inputs/celltalkdb.py#L72) | [definitions](../pypath/resources/data_formats.py#L694) |
| CollecTRI | collectri | [collectri.collectri_interactions](../pypath/inputs/collectri.py#L155) | [definitions](../pypath/resources/data_formats.py#L2481) |
| CollecTRI2 | collectri2 | [collectri2.collectri2_interactions](../pypath/inputs/collectri2.py#L234) | [definitions](../pypath/resources/data_formats.py#L2514) |
| CPDB | cpdb | [cpdb.cpdb_interactions](../pypath/inputs/cpdb.py#L25) | [definitions](../pypath/resources/data_formats.py#L1843) |
| dbPTM | dbptm | [dbptm.dbptm_interactions](../pypath/inputs/dbptm.py#L153) | [definitions](../pypath/resources/data_formats.py#L1430) |
| DeathDomain | deathdomain | [deathdomain.deathdomain_interactions_rescued](../pypath/inputs/deathdomain.py#L80) | [definitions](../pypath/resources/data_formats.py#L545) |
| DEPOD | depod | [depod.depod_interactions](../pypath/inputs/depod.py#L30) | [definitions](../pypath/resources/data_formats.py#L1308) |
| DIP | dip | [dip.dip_interactions](../pypath/inputs/dip.py#L26) | [definitions](../pypath/resources/data_formats.py#L1114) |
| DOMINO | domino | [domino.domino_interactions](../pypath/inputs/domino.py#L216) | [definitions](../pypath/resources/data_formats.py#L1395) |
| DoRothEA | dorothea | `dorothea.dorothea_interactions_old` **unresolved**<br>[dorothea.dorothea_interactions](../pypath/inputs/dorothea.py#L350) | [definitions](../pypath/resources/data_formats.py#L2408) |
| ELM | elm | [elm.elm_interactions](../pypath/inputs/elm.py#L127) | [definitions](../pypath/resources/data_formats.py#L1372) |
| EMBRACE | embrace | [embrace.embrace_interactions](../pypath/inputs/embrace.py#L119) | [definitions](../pypath/resources/data_formats.py#L3086) |
| ENCODE_tf-mirna | encode | [encode.encode_tf_mirna_interactions](../pypath/inputs/encode.py#L26) | [definitions](../pypath/resources/data_formats.py#L2792) |
| HINT | hint | [hint.hint_interactions](../pypath/inputs/hint.py#L64) | [definitions](../pypath/resources/data_formats.py#L1220) |
| HIPPIE | hippie | [hippie.hippie_interactions](../pypath/inputs/hippie.py#L30) | [definitions](../pypath/resources/data_formats.py#L1674) |
| HPMR | hpmr | [hpmr.hpmr_interactions](../pypath/inputs/hpmr.py#L356) | [definitions](../pypath/resources/data_formats.py#L2983) |
| HPRD-phos, HPRD | hprd | [hprd.hprd_interactions](../pypath/inputs/hprd.py#L47)<br>[hprd.hprd_interactions_htp](../pypath/inputs/hprd.py#L55) | [definitions](../pypath/resources/data_formats.py#L1451) |
| HTRIdb | htri | [htri.htri_interactions](../pypath/inputs/htri.py#L26) | [definitions](../pypath/resources/data_formats.py#L2217) |
| Lit-BM-17, HI-II, Lit-BM-13, HuRI, HI-union, Yu2011, Yang2016, Vidal HI-III, HI-III | huri | [huri.lit_bm_17_interactions](../pypath/inputs/huri.py#L187)<br>[huri.rolland_hi_ii_14](../pypath/inputs/huri.py#L31)<br>[huri.lit_bm_13_interactions](../pypath/inputs/huri.py#L153)<br>[huri.huri_interactions](../pypath/inputs/huri.py#L230)<br>[huri.hi_union_interactions](../pypath/inputs/huri.py#L240)<br>[huri.yu2011_interactions](../pypath/inputs/huri.py#L235)<br>[huri.yang2016_interactions](../pypath/inputs/huri.py#L245)<br>`huri.vidal_hi_iii` **unresolved**<br>`huri.hi_iii` **unresolved** | [definitions](../pypath/resources/data_formats.py#L1198) |
| InnateDB | innatedb | [innatedb.innatedb_interactions](../pypath/inputs/innatedb.py#L27) | [definitions](../pypath/resources/data_formats.py#L1159) |
| iPTMnet | iptmnet | [iptmnet.iptmnet_interactions](../pypath/inputs/iptmnet.py#L50) | [definitions](../pypath/resources/data_formats.py#L1495) |
| iTALK | italk | [italk.italk_interactions](../pypath/inputs/italk.py#L45) | [definitions](../pypath/resources/data_formats.py#L3103) |
| KEA | kea | [kea.kea_interactions](../pypath/inputs/kea.py#L39) | [definitions](../pypath/resources/data_formats.py#L1473) |
| Kirouac2010 | kirouac2010 | [kirouac2010.kirouac2010_interactions](../pypath/inputs/kirouac2010.py#L33) | [definitions](../pypath/resources/data_formats.py#L2961) |
| Laudanna-effects, Laudanna-directions | laudanna | [laudanna.laudanna_effects](../pypath/inputs/laudanna.py#L58)<br>[laudanna.laudanna_directions](../pypath/inputs/laudanna.py#L26) | [definitions](../pypath/resources/data_formats.py#L1005) |
| Li2012 | li2012 | [li2012.li2012_interactions](../pypath/inputs/li2012.py#L50) | [definitions](../pypath/resources/data_formats.py#L1610) |
| LMPID | lmpid | [lmpid.lmpid_interactions](../pypath/inputs/lmpid.py#L74) | [definitions](../pypath/resources/data_formats.py#L1327) |
| LncRNADisease | lncdisease | [lncdisease.lncdisease_interactions](../pypath/inputs/lncdisease.py#L27) | [definitions](../pypath/resources/data_formats.py#L2849) |
| lncrnadb | lncrnadb | [lncrnadb.lncrnadb_interactions](../pypath/inputs/lncrnadb.py#L30) | [definitions](../pypath/resources/data_formats.py#L2873) |
| LRdb | lrdb | [lrdb.lrdb_interactions](../pypath/inputs/lrdb.py#L40) | [definitions](../pypath/resources/data_formats.py#L3036) |
| Macrophage | macrophage | [macrophage.macrophage_interactions](../pypath/inputs/macrophage.py#L28) | [definitions](../pypath/resources/data_formats.py#L523) |
| MatrixDB | matrixdb | [matrixdb.matrixdb_interactions](../pypath/inputs/matrixdb.py#L28) | [definitions](../pypath/resources/data_formats.py#L1261) |
| MIMP | mimp | [mimp.mimp_interactions](../pypath/inputs/mimp.py#L139) | [definitions](../pypath/resources/data_formats.py#L1589) |
| miR2Disease | mir2disease | [mir2disease.mir2disease_interactions](../pypath/inputs/mir2disease.py#L27) | [definitions](../pypath/resources/data_formats.py#L2587) |
| miRDeathDB | mirdeathdb | [mirdeathdb.mirdeathdb_interactions](../pypath/inputs/mirdeathdb.py#L26) | [definitions](../pypath/resources/data_formats.py#L2608) |
| miRecords | mirecords | [mirecords.mirecords_interactions](../pypath/inputs/mirecords.py#L28) | [definitions](../pypath/resources/data_formats.py#L2667) |
| miRTarBase | mirtarbase | [mirtarbase.mirtarbase_interactions](../pypath/inputs/mirtarbase.py#L30) | [definitions](../pypath/resources/data_formats.py#L2691) |
| MPPI | mppi | [mppi.mppi_interactions](../pypath/inputs/mppi.py#L26) | [definitions](../pypath/resources/data_formats.py#L1093) |
| ncRDeathDB | ncrdeathdb | [ncrdeathdb.ncrdeathdb_interactions](../pypath/inputs/ncrdeathdb.py#L27) | [definitions](../pypath/resources/data_formats.py#L2634) |
| ARN, NRF2ome | netbiol | [netbiol.arn_interactions](../pypath/inputs/netbiol.py#L60)<br>[netbiol.nrf2ome_interactions](../pypath/inputs/netbiol.py#L69) | [definitions](../pypath/resources/data_formats.py#L477) |
| NetPath | netpath | [netpath.netpath_interactions](../pypath/inputs/netpath.py#L39) | [definitions](../pypath/resources/data_formats.py#L1136) |
| ORegAnno | oreganno | [oreganno.oreganno_interactions](../pypath/inputs/oreganno.py#L50) | [definitions](../pypath/resources/data_formats.py#L2238) |
| PathwayCommons | pathwaycommons | [pathwaycommons.pathwaycommons_interactions](../pypath/inputs/pathwaycommons.py#L79) | [definitions](../pypath/resources/data_formats.py#L75) |
| PAZAR | pazar | [pazar.pazar_interactions](../pypath/inputs/pazar.py#L26) | [definitions](../pypath/resources/data_formats.py#L2196) |
| PDZBase | pdzbase | [pdzbase.pdzbase_interactions](../pypath/inputs/pdzbase.py#L30) | [definitions](../pypath/resources/data_formats.py#L566) |
| phosphoELM | phosphoelm | [phosphoelm.phosphoelm_interactions](../pypath/inputs/phosphoelm.py#L98) | [definitions](../pypath/resources/data_formats.py#L1347) |
| PhosphoNetworks | phosphonetworks | [phosphonetworks.phosphonetworks_interactions](../pypath/inputs/phosphonetworks.py#L77) | [definitions](../pypath/resources/data_formats.py#L1569) |
| PhosphoPoint | phosphopoint | [phosphopoint.phosphopoint_interactions](../pypath/inputs/phosphopoint.py#L26) | [definitions](../pypath/resources/data_formats.py#L1547) |
| PhosphoSite, PhosphoSite_noref | phosphosite | [phosphosite.phosphosite_interactions_curated](../pypath/inputs/phosphosite.py#L1123)<br>[phosphosite.phosphosite_interactions_noref](../pypath/inputs/phosphosite.py#L1140) | [definitions](../pypath/resources/data_formats.py#L1287) |
| ProtMapper | protmapper | [protmapper.protmapper_interactions](../pypath/inputs/protmapper.py#L149) | [definitions](../pypath/resources/data_formats.py#L1632) |
| Ramilowski2015 | ramilowski2015 | [ramilowski2015.ramilowski_interactions](../pypath/inputs/ramilowski2015.py#L30) | [definitions](../pypath/resources/data_formats.py#L2934) |
| reaction resources, NCI-PID | reaction | [reaction.get_reactions](../pypath/inputs/reaction.py#L1165)<br>[reaction.pid_interactions](../pypath/inputs/reaction.py#L771) | [definitions](../pypath/resources/data_formats.py#L50) |
| scConnect | scconnect | [scconnect.scconnect_interactions](../pypath/inputs/scconnect.py#L157) | [definitions](../pypath/resources/data_formats.py#L813) |
| SignaLink3 | signalink | [signalink.signalink_interactions](../pypath/inputs/signalink.py#L33) | [definitions](../pypath/resources/data_formats.py#L409) |
| SLC | slc | [slc.slc_interactions](../pypath/inputs/slc.py#L114) | [definitions](../pypath/resources/data_formats.py#L3299) |
| SPIKE, SPIKE_LC | spike | [spike.spike_interactions](../pypath/inputs/spike.py#L33) | [definitions](../pypath/resources/data_formats.py#L387) |
| talklr | talklr | [talklr.talklr_interactions](../pypath/inputs/talklr.py#L47) | [definitions](../pypath/resources/data_formats.py#L773) |
| TransmiR | transmir | [transmir.transmir_interactions](../pypath/inputs/transmir.py#L27) | [definitions](../pypath/resources/data_formats.py#L2771) |
| TRIP | trip | [trip.trip_interactions](../pypath/inputs/trip.py#L302) | [definitions](../pypath/resources/data_formats.py#L365) |
| TRRUST | trrust | [trrust.trrust_interactions](../pypath/inputs/trrust.py#L38) | [definitions](../pypath/resources/data_formats.py#L2363) |
| Wang, Cui2007, HumanSignalingNetwork | wang | [wang.wang_interactions](../pypath/inputs/wang.py#L121)<br>[wang.cui_interactions](../pypath/inputs/wang.py#L135)<br>[wang.hsn_interactions](../pypath/inputs/wang.py#L51) | [definitions](../pypath/resources/data_formats.py#L902) |
| Wojtowicz2020 | wojtowicz2020 | [wojtowicz2020.wojtowicz2020_interactions](../pypath/inputs/wojtowicz2020.py#L80) | [definitions](../pypath/resources/data_formats.py#L3120) |

`reaction.get_reactions` aggregates upstream datasets; its source names must not be flattened into a new primary database. NCI-PID remains missing, whereas the Reactome wrapper in the same module is excluded because Reactome already exists in v2. PathwayCommons likewise needs source/via provenance and per-dataset identity. The SLC interaction tables are retained as a candidate: `inputs_v2/slctables.py` exports transporter reference data and is not evidence that the SLC substrate interaction table is covered.

### Additional datasets within registered module groups

These entry points are not all separately registered above. In particular, PANTHER is a missing upstream resource exposed through the shared `reaction` module.

| Resource / module | Additional entry points | Qualification |
| --- | --- | --- |
| PANTHER / reaction | [reaction.panther_interactions](../pypath/inputs/reaction.py#L776) | Source-specific BioPAX-derived interactions; no PANTHER v2 resource. |
| ACSN / reaction | [reaction.acsn_interactions_2](../pypath/inputs/reaction.py#L766) | Alternative ACSN parser, not another database. |
| BioGRID | [biogrid.biogrid_all_interactions](../pypath/inputs/biogrid.py#L98) | Broader interaction dataset; inspect evidence types. |
| CPDB | [cpdb.cpdb_interactions_ltp](../pypath/inputs/cpdb.py#L65) | Filtered low-throughput view. |
| DeathDomain | [deathdomain.deathdomain_interactions](../pypath/inputs/deathdomain.py#L28) | Alternative to registered rescued input. |
| HuRI | `hi_ii_interactions`, `hi_i_interactions`, `lit_bm_interactions` in [huri.py](../pypath/inputs/huri.py) | Additional versioned datasets; do not replace stale names without checking the intended release. |
| PhosphoSite | `phosphosite_interactions`, `phosphosite_interactions_new`, `phosphosite_interactions_all` in [phosphosite.py](../pypath/inputs/phosphosite.py) | Different return shapes and combined views; choose an explicit contract. |
| DOMINO, Li2012, LMPID | [domino_ddi](../pypath/inputs/domino.py#L237), [li2012_dmi](../pypath/inputs/li2012.py#L144), [lmpid_dmi](../pypath/inputs/lmpid.py#L85) | Domain/domain-motif object outputs need the structural mapping route, not ordinary row adaptation. |

PathwayCommons dynamically generates per-source views for NCI-PID, KEGG, CORUM, BIND, HPRD, WikiPathways, INOH, BioGRID, NetPath, Reactome, DIP, IntAct, PANTHER and PhosphoSite, including transcriptional variants ([source list and factory](../pypath/resources/data_formats.py#L174)). **BIND and INOH are additional missing upstream resources accessible through this aggregator**, rather than standalone legacy modules. NCI-PID and PANTHER also have the reaction entry points listed above. Already-native sources in that list are not additional missing resources, but their via-PathwayCommons observations remain distinct provenance.

### Additional interaction functions outside those registered module groups

These have interaction-oriented entry points but no literal `NetworkInput` registration in the scanned file. They need an explicit schema/contract audit before enabling; a function name alone does not establish adapter compatibility. Functions in a module already listed above (for example additional HuRI releases or alternate PhosphoSite views) are variants of that source, not additional resource rows.

| Legacy module / source | Entry points |
| --- | --- |
| biogrid_new | [biogrid_new.biogrid_interactions](../pypath/inputs/biogrid_new.py#L43)<br>[biogrid_new.biogrid_all_interactions](../pypath/inputs/biogrid_new.py#L104) |
| comppi | [comppi.comppi_interaction_locations](../pypath/inputs/comppi.py#L34) |
| ddinter | [ddinter.ddinter_drug_interactions](../pypath/inputs/ddinter.py#L142)<br>[ddinter.ddinter_interactions](../pypath/inputs/ddinter.py#L192) |
| dgidb | [dgidb.dgidb_interactions](../pypath/inputs/dgidb.py#L32) |
| drugbank | [drugbank.drugbank_raw_interactions](../pypath/inputs/drugbank.py#L97)<br>[drugbank.drugbank_interactions](../pypath/inputs/drugbank.py#L168) |
| negatome | [negatome.negatome_interactions](../pypath/inputs/negatome.py#L26) |
| pathophenodb | [pathophenodb.disease_pathogen_interactions](../pypath/inputs/pathophenodb.py#L46) |
| pepcyber | [pepcyber.pepcyber_interactions](../pypath/inputs/pepcyber.py#L38) |
| proteinatlas | [proteinatlas.proteinatlas_interactions](../pypath/inputs/proteinatlas.py#L458) |
| string | [string.string_links_interactions](../pypath/inputs/string.py#L102)<br>[string.string_physical_interactions](../pypath/inputs/string.py#L191) |
| twosides | [twosides.twosides_interactions](../pypath/inputs/twosides.py#L31) |

`biogrid_new` is another BioGRID implementation, not a second database; choose one maintained implementation. `pathophenodb` is disease–pathogen evidence and belongs in a separate association mapping if included; `comppi` adds compartment-specific interaction information rather than a generic two-column edge. `twosides` contains drug-pair adverse-event associations; do not flatten its event context into an ordinary physical interaction.

### Structural and nonstandard interaction APIs

| Source | Entry points | Adapter boundary |
| --- | --- | --- |
| InWeb / InWeb_IM | [imweb.get_inweb](../pypath/inputs/imweb.py#L150); [imweb.get_imweb](../pypath/inputs/imweb.py#L90) | Resource-specific table / access flow; no NetworkInput contract. |
| Interactome3D | [i3d.get_i3d](../pypath/inputs/i3d.py#L29) | Structural interaction records; preserve chain/domain evidence. |
| INstruct | [instruct.get_instruct](../pypath/inputs/instruct.py#L26) | Structural records and separate offsets. |
| 3did | [threedid.get_3did_ddi](../pypath/inputs/threedid.py#L41); [threedid.get_3did_dmi](../pypath/inputs/threedid.py#L451) | Domain interactions; nested records/optional interface returns. |
| 3DComplex | [threedcomplex.threedcomplex_ddi](../pypath/inputs/threedcomplex.py#L43); [threedcomplex.threedcomplex_contacts](../pypath/inputs/threedcomplex.py#L169) | Implemented contact/DDI APIs despite complex-download stub. |
| PISA | [pisa.pisa_interfaces](../pypath/inputs/pisa.py#L100) | Requires PDB inputs; interface objects, optional unmapped return. |
| iELM | [ielm.get_ielm](../pypath/inputs/ielm.py#L129); [ielm.get_ielm_huge](../pypath/inputs/ielm.py#L41) | Requires input network; derived analysis, not standalone source acquisition. |
| Switches.ELM | [switches_elm.get_switches_elm](../pypath/inputs/switches_elm.py#L38) | PTM-dependent interaction regulation; dedicated site/mechanism mapping. |

These are not promised coverage of the initial row-based network adapter. `mitab.mitab_interactions` is a shared parser, not a resource; `reaction.process_complex` is a processing helper. Disease/phenotype/clinical association APIs such as CTDbase, DisGeNET, OpenTargets, DISEASES, ClinVar, ADReCS, SIDER and OFFSIDES require a separate association inventory if that scope is added. CTDbase also has chemical–gene relations, which are worth a later dedicated audit; they are not silently classified as compatible network rows here.

### Stale / non-module NetworkInput references

| Resource label | Configured input | Action |
| --- | --- | --- |
| HI-III | `/home/denes/Documents/pw/data/hi3-2.3.tsv` | [Review definition](../pypath/resources/data_formats.py#L1869); direct file/URL input, outside legacy-function adapter |
| ENCODE-distal | `http://encodenets.gersteinlab.org/enets3.Distal.txt` | [Review definition](../pypath/resources/data_formats.py#L2154); direct file/URL input, outside legacy-function adapter |
| ENCODE-proximal | `http://encodenets.gersteinlab.org/enets2.Proximal_filtered.txt` | [Review definition](../pypath/resources/data_formats.py#L2175); direct file/URL input, outside legacy-function adapter |
| ORegAnno | `get_oreganno_old` | [Review definition](../pypath/resources/data_formats.py#L2561); resolve historical name to an existing function or remove stale registration |
| Negatome | `negatome_pairs` | [Review definition](../pypath/resources/data_formats.py#L3378); resolve historical name to an existing function or remove stale registration |
| GDSC | `gdsc.sif` | [Review definition](../pypath/resources/data_formats.py#L3399); resolve historical name to an existing function or remove stale registration |

The stale-reference table includes `gdsc.sif` (no legacy module found). The module table also exposes older HuRI function names and `dorothea.dorothea_interactions_old`. `negatome.negatome_interactions` exists even though its registry uses the historical bare `negatome_pairs` name. Negative interaction evidence requires explicit negation support; never emit it as a positive interaction.

## Alias resolution and resources already present

| Legacy module(s) | Existing v2 module | Qualification |
| --- | --- | --- |
| guide2pharma | guidetopharma | Same Guide to Pharmacology resource. |
| cellchatdb | cellchat | Same CellChat resource. |
| connectomedb2025 | connectomedb | Verified v2 config name is ConnectomeDB2025 and downloads Current-Release. |
| connectomedb2020 | connectomedb | Resource family exists; the 2020 release is not claimed reproduced. Ramilowski2015 remains a separately labeled legacy dataset in the network backlog. |
| new_intact | intact | Alternate implementation of existing resource. |
| new_signor | signor | Alternate implementation of existing resource. |
| new_stitch | stitch | Alternate implementation of existing resource. |
| reactome_old; reaction.reactome_interactions | reactome | Same upstream resource; representations differ. |
| kegg_api | kegg | API/helper family of existing resource; not proof every endpoint is migrated. |

Same-name v2 modules exclude, among others, CORUM, CellPhoneDB, Cellinker, ICELLNET, KEGG, SIGNOR, IntAct, STITCH, BindingDB, DrugCentral, NicheNet and MRClinksDB. Presence must not be confused with dataset parity: separately check SIGNOR enzyme–substrate/complex outputs and KEGG Medicus complexes before claiming those legacy paths are replaced. Avoid creating parallel legacy wrappers that duplicate a native v2 dataset.

## Repeating the inventory

1. Enumerate `pypath/inputs/**/*.py` and `pypath/inputs_v2/*.py` without importing modules or downloading data.
2. Parse the complex class `input_method` definitions, the JSON `inputs.enzyme_substrate` entries, and all `NetworkInput(...)` calls; retain labels, functions, arguments and source locations.
3. Compare by resource identity using the alias table, not only filename equality.
4. Supplement with public complex, enzyme–substrate and interaction functions and the structural APIs listed above; separate generic helpers, derived analyses and incomplete functions.
5. Verify entry-point existence and return-shape fixtures before moving a row from inventory to implementation. Validate runtime-generated variants when selecting an actual resource configuration.

This inventory deliberately makes no claim that legacy datasets are operational, complete, licensed for a particular deployment, or semantically identical to their native v2 counterparts.
