# Legacy parity of the inputs_v2 ports

14 October 2026 · pypath `15de653` · build cache on nicesrv (`/root/data/pypath-cache`) ·
harness [`scripts/legacy_parity.py`](../../../scripts/legacy_parity.py)

## What was compared

For each resource, the legacy pypath input (`pypath.inputs`) and its inputs_v2
port read **the same cached file**, and both are reduced to pairs in the
source's own identifiers: interaction partners, an entity and its annotation
(class, cell type, phenotype, pathway), a member and its group, a metabolite
and its cross-references; complexes are compared by their member sets. Signs
are compared where both sides give one.

- Nothing is downloaded. inputs_v2 opens the cached files; a legacy download
  (through `Curl`, the download manager, `requests` or `urllib`) is served the
  inputs_v2 file of the same request, or fails with `NotCached`. Where legacy
  writes the URL of the same file differently (a GitHub spelling, OmniPath's
  rescued copy), the spec says so.
- Legacy identifier translation (`pypath.utils.mapping`) is replaced by
  identity, so legacy keeps the source's identifiers as inputs_v2 does.
- A pair found only the other way round (an undirected interaction listed in
  one orientation) counts as found.
- HMDB is compared on the first 20,000 metabolites of the file on both sides;
  everything else in full.

The JSON file of each resource has the counts, up to 25 examples of each kind
of difference and the files each side read, with their checksums.

## Results

| Resource | Legacy keys | Found in v2 | Only legacy | Only v2 | Signs agree | Verdict |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| adhesome | 3,002 | 2,794 (93.1 %) | 208 | 177 | 149 / 155 | v2 takes the 2007 network from the paper (decision 9) |
| arn | 95 | 95 (100.0 %) | 0 | 0 | 28 / 28 | identical |
| ca1 | 1,259 | 1,056 (83.9 %) | 203 | 0 | 551 / 553 | small molecules left out (kept, as decided) |
| cancercellmap | – | – | – | – | – | legacy needs a file v2 does not read |
| cancerdrugsdb | 6,290 | 6,290 (100.0 %) | 0 | 23 | – | identical; v2 has 23 more |
| cellcall | 1,141 | 1,141 (100.0 %) | 0 | 306 | – | identical; v2 adds the extended set |
| cellcellinteractions | 3,877 | 3,877 (100.0 %) | 0 | 0 | – | identical |
| cellchat | 3,233 | 3,229 (99.9 %) | 4 | 0 | – | **fixed**: metabolite ligands |
| cellphonedb | 2,903 | 2,903 (100.0 %) | 0 | 4 | – | **fixed**: metabolite ChEBI terms, HLA partners |
| cellphonedb_roles | 2,533 | 2,533 | 0 | 0 | – | **added**: protein and complex roles |
| celltalkdb | 3,398 | 3,398 (100.0 %) | 0 | 0 | – | identical |
| celltypist | 1,163 | 1,163 (100.0 %) | 0 | 0 | – | identical |
| collectri | 43,536 | 43,416 (99.7 %) | 120 | 120 | 19,587 / 19,587 | miRNA targets named differently (decision 1) |
| compleat | 9,702 | 9,702 | 0 | 1 | – | identical |
| complexportal | – | – | – | – | – | not the same file (PSI-MI XML vs ComplexTAB) |
| corum | 2,734 | 2,733 | 1 | 1 | – | v2 cleans a source typo |
| cpad | 2,225 | 2,225 (100.0 %) | 0 | 8,183 | – | identical; v2 has more pairs |
| cpdb | 531,371 | 522,618 (98.4 %) | 8,753 | 5,198 | – | legacy keeps alternatives joined with "." |
| cytosig | 209,883 | 209,883 (100.0 %) | 0 | 0 | – | identical |
| dbptm | 1,107 | 1,100 (99.4 %) | 7 | 1 | – | v2 splits "MAPK1 OR MAPK3" |
| dip | 2,281 | 2,280 (100.0 %) | 1 | 5,509 | – | identical (legacy on the full file) |
| domino | 4,549 | 4,538 (99.8 %) | 11 | 745 | – | nearly identical |
| dorothea | 1,019,144 | 1,019,144 (100.0 %) | 0 | 0 | 5,702 / 5,702 | identical |
| elm | 2,128 | 2,128 (100.0 %) | 0 | 0 | – | identical |
| encode | 1,237 | 1,237 (100.0 %) | 0 | 0 | – | identical (TF-miRNA table) |
| gpcrdb | 813 | 811 (99.8 %) | 2 | 753 | – | v2 strips HTML from two names |
| guidetopharma | 742 | 380 (51.2 %) | 362 | 23,757 | 289 / 289 | **fixed**: peptide gene symbols; legacy expands subunits (decision 13) |
| havugimana | 622 | 622 | 0 | 0 | – | identical |
| hippie | 546,755 | 538,622 (98.5 %) | 8,133 | 0 | – | identical at legacy's threshold (decision 2) |
| hmdb | 52,929 | 52,928 (100.0 %) | 1 | 1 | – | **fixed**: KEGG IDs (first 20,000 metabolites) |
| hpo | 274,064 | 274,064 (100.0 %) | 0 | 0 | – | identical |
| hprd | 4,264 | 4,256 (99.8 %) | 8 | 5 | – | nearly identical (legacy with in vitro records) |
| htri | 18,160 | 18,160 (100.0 %) | 0 | 0 | – | identical |
| humap | 15,326 | 15,326 | 0 | 6,244 | – | identical; v2 adds huMAP2 |
| icellnet | 1,669 | 1,657 (99.3 %) | 12 | 12 | – | legacy mangles 12 DKK names |
| imm1415_metabolites | 7,478 | 7,478 (100.0 %) | 0 | 499 | – | identical |
| instruct | 6,585 | 6,585 (100.0 %) | 0 | 0 | – | identical |
| italk | 2,575 | 2,575 (100.0 %) | 0 | 0 | – | identical |
| kea | 7,181 | 7,181 (100.0 %) | 0 | 21 | – | identical |
| kinasedotcom | 1,266 | 1,260 (99.5 %) | 6 | 309 | – | identical (legacy reads the header as a kinase) |
| lambert2018 | 4,404 | 3,162 (71.8 %) | 1,242 | 94 | – | v2 lists the TFs only (decision 4) |
| laudanna_directions | 75,393 | 75,393 (100.0 %) | 0 | 0 | – | identical |
| laudanna_effects | 62,050 | 62,050 (100.0 %) | 0 | 0 | 40,752 / 40,752 | identical |
| li2012 | 185 | 185 (100.0 %) | 0 | 0 | – | identical |
| lipidmaps | – | – | – | – | – | legacy fails on the current file |
| lncdisease | 313 | 311 (99.4 %) | 2 | 0 | – | nearly identical |
| lncrnadb | 301 | 301 (100.0 %) | 0 | 4 | – | identical |
| lrdb | 3,239 | 3,239 (100.0 %) | 0 | 0 | – | identical |
| macrophage | 1,702 | 1,669 (98.1 %) | 33 | 506 | – | **fixed**: partners read as legacy (decision 8) |
| matrisome | – | – | – | – | – | not the same file (rescued XLS vs current sheet) |
| mcam | 386 | 386 (100.0 %) | 0 | 0 | – | identical |
| mir2disease | – | – | – | – | – | legacy fails on the file |
| mirdeathdb | 395 | 395 (100.0 %) | 0 | 4 | – | **fixed**: each miRBase product |
| mirecords | 2,611 | 2,387 (91.4 %) | 224 | 199 | – | v2 drops "UNKNOWN" targets, cleans names |
| mppi | 465 | 430 (92.5 %) | 35 | 767 | – | legacy junk ("-", isoforms "_2") |
| ncrdeathdb | 5,427 | 5,399 (99.5 %) | 28 | 107 | – | **fixed**: each miRBase product |
| netpath | 4,216 | 2,383 (56.5 %) | 1,833 | 3 | – | v2 keeps complexes (decision 3) |
| nrf2ome | 109 | 109 (100.0 %) | 0 | 0 | 41 / 41 | identical |
| oreganno | 1,488 | 1,456 (97.8 %) | 32 | 9,290 | – | **fixed**: multi-name cells; v2 has all species (decision 7) |
| panglaodb | 8,286 | 8,253 (99.6 %) | 33 | 0 | – | nearly identical (renamed genes) |
| pazar | 10,045 | 10,043 (100.0 %) | 2 | 2 | – | identical |
| pdzbase | 276 | 276 (100.0 %) | 0 | 0 | – | identical |
| phobius | 9,622 | 9,622 (100.0 %) | 0 | 0 | – | identical |
| phosphatome | 863 | 863 (100.0 %) | 0 | 511 | – | identical |
| phosphonetworks | 1,821 | 1,821 (100.0 %) | 0 | 0 | – | identical |
| phosphopoint | 9,269 | 9,269 (100.0 %) | 0 | 0 | – | identical |
| phosphosite | 7,184 | 7,183 (100.0 %) | 1 | 70 | – | one legacy RefSeq bug; v2 adds hamster/frog substrates |
| progeny | 274,143 | 274,143 (100.0 %) | 0 | 0 | – | identical |
| psi_mi | 1,665 | 1,665 (100.0 %) | 0 | 0 | – | identical |
| ramilowski2015 | 1,894 | 1,894 (100.0 %) | 0 | 0 | – | identical |
| recon3d_metabolites | 8,388 | 8,388 (100.0 %) | 0 | 0 | – | identical |
| recon3d_reactions | 33,727 | 33,727 (100.0 %) | 0 | 17,438 | – | identical |
| signalink | 21,259 | 21,259 (100.0 %) | 0 | 0 | 1,669 / 1,669 | identical |
| signor | – | – | – | – | – | not the same file (CSV export vs causalTab) |
| signor_complexes | 681 IDs | 681 | 0 | 9 | – | v2 keeps complexes legacy overwrote |
| spike | 8,327 | 7,652 (91.9 %) | 675 | 36,352 | – | legacy expands complexes to members (decision 12) |
| spike_complexes | 154 | 153 | 1 | 425 | – | identical |
| swisslipids | 1,325,241 | 1,140,644 (86.1 %) | 184,597 | 5 | – | identical on IDs; the rest are legacy "NONE" placeholders |
| talklr | 2,422 | 2,422 (100.0 %) | 0 | 0 | – | identical |
| tcdb | 24,747 | 24,747 (100.0 %) | 0 | 129 | – | identical |
| tfcensus | 3,821 | 3,084 (80.7 %) | 737 | 0 | – | v2 leaves out class x (decision 5) |
| topdb | 9,341 | 8,321 (89.1 %) | 1,020 | 0 | – | v2 leaves out PDB-chain entries (decision 6) |
| trrust | 8,427 | 8,427 (100.0 %) | 0 | 0 | 4,859 / 4,859 | identical |
| wang | 62,937 | 62,812 (99.8 %) | 125 | 0 | 41,259 / 41,259 | nearly identical |
| zhong2015 | 467 | 467 (100.0 %) | 0 | 0 | – | identical |

"Legacy keys" are distinct legacy pairs (or complexes); "found" includes
pairs found the other way round; "only v2" counts v2 relations with no legacy
pair, often because v2 imports more of the file.

## Fixed in inputs_v2 (rebuilt on nicesrv as 2026.10.14, not yet released)

- **CellChat**: 988 human (995 mouse) interactions had the enzyme, or an
  enzyme-transporter complex, as the ligand of a metabolite. The ligand is
  now the metabolite with its ChEBI term; the proxy genes are associated with
  it in new `*_metabolite_genes` datasets.
- **CellPhoneDB**: the protein and complex roles (receptor, secreted,
  transmembrane, peripheral, integrin, with the receptor and secreted classes)
  were not imported; they are, and agree with legacy's reading. Metabolite
  ligands carry their ChEBI terms (48 of 52). HLA-A, -B, -C are proteins.
- **Macrophage**: the partner cells are read as legacy read them (complexes
  split into their genes, forms in brackets removed): 1,658 of the 1,702
  legacy pairs found, 174 before.
- **ORegAnno**: a factor or gene cell of several names (`STAT1/STAT3`) gives
  its first name as the symbol, as legacy.
- **miRDeathDB, ncRDeathDB**: each miRBase accession listed for a microRNA
  gives a relation, as legacy; before, a symbol with several products had
  none.
- **HMDB**: no metabolite had its KEGG compound ID (a pattern with two groups);
  now they do.
- **Guide to Pharmacology**: peptide ligands have the gene symbol the
  interactions file gives for them, as legacy named them.

## For decision

Places where v2 does not do what legacy did, or adds what the source does not
state. Each is how the module is now; nothing was changed for these.

1. **CollecTRI**: legacy renames miRNA targets to the mature miRNA
   (`hsa-miR-223`); v2 keeps the gene symbol (`MIR223`). 120 pairs.
2. **HIPPIE**: legacy keeps interactions of score 0.75 or more; v2 keeps all
   (1.1 million relations) with the score. At 0.75 the two agree.
3. **NetPath**: an interaction of three or more proteins is a complex in v2;
   legacy made every pair of its members, and also each protein with itself.
4. **Lambert 2018**: legacy annotated all 2,765 candidate genes with their DBD
   family and an `is_tf` flag; v2 has the 1,639 genes judged TFs.
5. **TF census**: v2 leaves out the loci of class x (judged not TFs); legacy
   kept them with their class.
6. **TopDB**: v2 has the proteins named by UniProt entry names; legacy also
   kept about 1,020 entries named by PDB chain (`2K3C_A`).
7. **ORegAnno**: v2 imports all species; legacy took one organism (human by
   default). v2 also keeps the dash of a one-name cell (`HLA-A`), which legacy
   removed.
8. **Macrophage**: v2 removes each bracket of a cell; legacy removed from the
   first to the last bracket and lost the members between
   (`NFKB1(p50):RELA(p65)` was NFKB1 alone). v2 also keeps the rows with a gene
   as partner, which legacy left out; in 38 of them the gene is partner 1 of an
   activation by a complex, which reads reversed but is how the sheet states it.
9. **Adhesome**: for the 2007 network v2 takes the interactions and their
   direction from the publication, as the website lists every interaction in
   both orientations; 158 legacy pairs from the website are not among them.
10. **CellChat and CellPhoneDB metabolite ChEBI terms**: the CellChat table of
    ChEBI terms, and the CellPhoneDB terms for metabolites the source gives
    without one, are identifiers the source does not state (each looked up by
    name in ChEBI). By the rule "do not add what the source does not provide",
    these would be names only.
11. **GO**: the v2 module imports the ontology only; legacy's GO inputs read
    the GOA gene annotations, which no v2 module imports.
12. **SPIKE**: legacy wrote a complex partner as its members' Entrez IDs; v2
    refers to the SPIKE complex (with its members in the complex dataset).
13. **Guide to Pharmacology**: a peptide ligand of several subunits (activins,
    inhibins: `INHBA|INHBB`) is one ligand in v2, with its subunits; legacy
    paired each subunit gene with the target. Legacy also keeps human targets
    only. (Single-gene peptide ligands now have their gene symbol, as legacy
    named them.)

v2 also cleans source errors that legacy kept, without a decision needed:
CORUM `Q9NZL9-`, GPCRdb HTML in names, dbPTM `X OR Y` kinases, miRecords
`UNKNOWN` targets, SwissLipids and HMDB malformed IDs, PhosphoSite RefSeq
substrates, SIGNOR complexes with the same members, ICELLNET DKK names.

## Not comparable on the same file

| Resource | Why |
| --- | --- |
| SIGNOR interactions, ComplexPortal, Matrisome | the two generations read different exports |
| HPMR, phospho.ELM, MIMP, iPTMnet, PepCyber, Reactome, Rhea, CancerCellMap, CancerSEA, MetAtlas | legacy reads files v2 does not |
| Baccin 2019, LMPID, LOCATE, scConnect, 3DID, I3D | legacy needs UniProt or Pfam tables |
| SwitchesELM, DisGeNET | legacy calls a web service (OLS, the DisGeNET API) |
| Almen 2009 | legacy keys the proteins by IPI, v2 by UniProt |
| LIPID MAPS, miR2Disease, PDB | the legacy function fails on the current file |
| GO | see decision 11 |
| MEBOCOST | no legacy input |

These are left for the served-level comparison against the legacy web service.

## Notes on the harness

- The legacy package does not import cleanly any more: module-level lookups
  call UniChem, miRBase and Ensembl. The harness stubs them, and serves
  organism names for human, mouse and rat from a small table.
- `pypath.inputs.rdata` (legacy) patches the `rdata` converter on import,
  after which the v2 CellChat parser fails; each resource runs in its own
  process, v2 before legacy.
- Run on nicesrv, from the build checkout:
  `PYPATH_DOWNLOAD_DATADIR=/root/data/pypath-cache uv run --frozen --with pycurl python scripts/legacy_parity.py --out <dir>`.
