# The ChEBI role walk and the audit of `chebi_role` on release 2026.11.27

2026-11-27

This report covers step 6 of stage two of the plan
`saezverse/projects/omnipath/plans/2026-11-metabolite-classification-logic.md`,
and the ChEBI part of step 7. Step 6 compares the role walk from
`CHEBI:25212` with the list of 45 roles that the plan used as a stand-in.
Step 7 audits a sample of the bucket `chebi_role`. The earlier report
`metabolite-class-stage-two-measurements-20261127.md` measured release
2026.11.19, before the ChEBI parser kept the role branch.

**Result: the walk reaches 60 roles and 16,253 molecules. The list of 45 roles
reaches 16,136, all of them inside the walk.** With all three signals on,
`chebi_role` holds 11,912 molecules, 9,182 of them `entity` rows. In the
audit sample of these rows, 39 of 50 are specific metabolites. Only 7 of 50
carry a human or mammalian role. The rest are natural products of plants,
fungi, bacteria and marine organisms. In a second sample, restricted to human
and mammalian roles, 21 of 50 are classes or placeholders, and all 21 have no
InChIKey.

## What this report measures

| Item | Value |
|---|---|
| Release | `2026.11.27`, read from `/workspace/instances/metsig1/data-2026.11.19` |
| Resources that changed | `chebi` and `hmdb` at `2026.11.27`. Every other resource, KEGG, Rhea, Recon3D and metatlas included, stays at `2026.11.19`. |
| Code | `pypath` at `0f4978696`. The subset code last changed in `70e6006e4`. |
| Build | `parquet_build.build_tables` with build id `audit` |
| Scripts | In `/tmp/claude-1029/-workspace-instances-metsig1/97c1ab55-a53d-4d90-8da3-8c1f34a8d599/scratchpad/stage2-chebi/`, listed at the end |

The ChEBI build of 2026.11.27 also reads a newer upstream ChEBI. It holds
about 13,300 new ChEBI identifiers. So a change between the two releases can
come from the parser fix or from new ChEBI content. The list of 45 roles
separates the two: it reaches 16,063 molecules on 2026.11.19 and 16,136 on
2026.11.27. The new content adds 73 carriers to the list, and it removes
none. So most of the growth of the bucket comes from the walk.

**Carrier.** A carrier of a role is the subject of a ChEBI
`has_chemical_role` edge to that role, with the type `chemical_entity`. This
is the rule of `claim_chebi_role.sql`.

**The list of 45 roles.** It is the set of ChEBI entities whose label ends in
"metabolite", without "antimetabolite", and that are the object of at least
one `has_chemical_role` edge. ChEBI has 68 entities whose label ends in
"metabolite". The 23 others have no carrier. 21 of them are molecules, not
roles, for example "Valsartan metabolite" and "JWH 018 N-pentanoic acid
metabolite". The list is the same 45 roles on both releases.

## Step 6: the walk and the list

| Figure | Walk from `CHEBI:25212` | List of 45 roles |
|---|---|---|
| Roles | 60, depth 0 to 5 | 45 |
| Roles with at least one carrier | 57 | 45 |
| Carriers | 16,253 | 16,136 |
| Carriers of both | 16,136 | 16,136 |
| Carriers of one only | 117 | 0 |

Every role of the list is in the walk. By depth, the walk holds 1, 10, 18,
19, 6 and 6 roles, from depth 0 to depth 5.

**Roles in the walk and not in the list.** Their labels do not end in
"metabolite". The carrier count is that of the role itself.

| Role | Label | Depth | Carriers |
|---|---|---|---|
| `CHEBI:27311` | volatile oil component | 3 | 159 |
| `CHEBI:25442` | mycotoxin | 2 | 68 |
| `CHEBI:27026` | toxin | 1 | 51 |
| `CHEBI:38231` | phytotoxin | 2 | 37 |
| `CHEBI:26115` | phytoalexin | 3 | 22 |
| `CHEBI:62215` | allelochemical | 1 | 13 |
| `CHEBI:64584` | uremic toxin | 2 | 10 |
| `CHEBI:72316` | virulence factor | 2 | 9 |
| `CHEBI:88048` | cyanotoxin | 2 | 8 |
| `CHEBI:61015` | nephrotoxin | 2 | 6 |
| `CHEBI:78505` | venom | 2 | 4 |
| `CHEBI:64477` | monolignol | 3 | 3 |
| `NIERWTDCCGPRCQ-ZFVALMAISA-N` | caudatoside | 3 | 0 |

Three terms of the walk have no carrier of their own: prokaryotic
metabolite (`CHEBI:75787`), Daphnia metabolite (`CHEBI:83057`) and
caudatoside. The first two have carriers below them. Caudatoside is a
molecule, not a role. ChEBI gives it a `subclass_of` edge to phytotoxin. It
has no carrier, so it changes no count.

**Carriers that only the walk reaches.** The 117 molecules come through
these roles: mycotoxin (35), volatile oil component (31), toxin (17),
phytoalexin (13) and phytotoxin (12). Cyanotoxin and allelochemical give 3
each, monolignol and virulence factor 2 each, and nephrotoxin 1. A molecule can count under more
than one role. 100 of the 117 have an InChIKey. Examples: fumitremorgin C,
cyclopiazonic acid, notoamide D and botrydial (mycotoxins), thymol,
beta-eudesmol and 2,4-xylenol (volatile oil components), erdasporine A and B
(toxins). These are natural products of fungi and plants. No carrier is only
in the list.

**The root role carries a third of the walk.** 5,317 molecules carry
`CHEBI:25212` itself. 5,022 of them carry no other role of the walk. ChEBI
gives these molecules no organism.

### The main sub-roles

"Roles below" counts the role and every role under it. "Carriers" counts the
molecules with a role of that subtree. "In the bucket" counts the carriers
that end as `entity` rows of `chebi_role` in the full build. These are not in
HMDB and not in the connectivity group of an HMDB molecule. The
subtrees overlap.

| Role | Roles below | Carriers | In the bucket | Of these, without an InChIKey |
|---|---|---|---|---|
| eukaryotic metabolite (`CHEBI:75763`) | 38 | 9,454 | 4,342 | 614 |
| animal metabolite (`CHEBI:75767`) | 19 | 4,971 | 1,685 | 515 |
| mammalian metabolite (`CHEBI:75768`) | 9 | 4,374 | 1,198 | 491 |
| human metabolite (`CHEBI:77746`) | 5 | 3,013 | 641 | 157 |
| mouse metabolite (`CHEBI:75771`) | 1 | 1,641 | 444 | 231 |
| rat metabolite (`CHEBI:86264`) | 1 | 186 | 49 | 7 |
| Papio hamadryas metabolite (`CHEBI:137684`) | 1 | 158 | 149 | 143 |
| human xenobiotic metabolite (`CHEBI:76967`) | 2 | 665 | 234 | 4 |
| drug metabolite (`CHEBI:49103`) | 1 | 192 | 63 | 1 |
| plant metabolite (`CHEBI:76924`) | 10 | 3,574 | 1,787 | 65 |
| fungal metabolite (`CHEBI:76946`) | 7 | 1,789 | 1,000 | 89 |
| Saccharomyces cerevisiae metabolite (`CHEBI:75772`) | 1 | 763 | 178 | 67 |
| prokaryotic metabolite (`CHEBI:75787`) | 8 | 1,952 | 1,120 | 111 |
| bacterial metabolite (`CHEBI:76969`) | 6 | 1,949 | 1,119 | 110 |
| Escherichia coli metabolite (`CHEBI:76971`) | 1 | 678 | 158 | 68 |
| marine metabolite (`CHEBI:76507`) | 3 | 505 | 358 | 15 |
| algal metabolite (`CHEBI:84735`) | 1 | 143 | 36 | 9 |
| xenobiotic metabolite (`CHEBI:76206`) | 6 | 965 | 406 | 18 |
| fundamental metabolite (`CHEBI:78675`) | 1 | 155 | 9 | 5 |

The subtree of human metabolite holds human urinary metabolite, human blood
serum metabolite, human xenobiotic metabolite and drug metabolite. The
subtree of mammalian metabolite adds mouse, rat and Papio hamadryas
metabolite. So the two roots of the second sample give the same set: the
subtree of `CHEBI:75768`, 9 roles and 4,374 carriers.

Two parts of this subtree are not endogenous metabolites. 184 of its 1,198
rows in the bucket carry only drug metabolite or human xenobiotic metabolite.
146 carry only Papio hamadryas metabolite, and 138 of these are lipid species
with a Goslin name and no InChIKey.

## Task 2: the buckets with all three signals

The rule file is the packaged one with the stage-two rows uncommented:
`rules_full.yaml`. The build takes about 4 seconds.

| Bucket | Molecules | `entity` rows | `structure_key` rows | `no_key` | `shared_key` groups |
|---|---|---|---|---|---|
| `hmdb` | 272,373 | 217,916 | 54,457 | 21 | 34,639 |
| `chebi_role` | 11,912 | 9,182 | 2,730 | 947 | 1,848 |
| `metabolic_network` | 15,979 | 12,894 | 3,085 | 5,881 | 2,117 |

The sources claim 217,916 molecules (HMDB), 16,253 (the ChEBI walk) and
23,766 (the reaction resources). Of the 16,253 molecules that the walk
claims, 4,012 are HMDB molecules and 3,059 are in the group of one. 9,182
stay in `chebi_role`.

**Compared with 2026.11.19.** The left column gives the figures of the
earlier report.

| Figure | 2026.11.19 | 2026.11.27 |
|---|---|---|
| Molecules that the ChEBI walk claims | 5,317 | 16,253 |
| `chebi_role` molecules | 4,608 | 11,912 |
| `chebi_role` through a group | 768 | 2,730 |
| `chebi_role` `entity` rows | 3,840 | 9,182 |
| `chebi_role` without an InChIKey | 257 | 947 |
| `hmdb` molecules | 272,030 | 272,373 |
| `metabolic_network` molecules | 18,494 | 15,979 |
| `metabolic_network` `entity` rows | 14,537 | 12,894 |

`chebi_role` grows to 2.6 times its size. Signal 2 comes before signal 3, so
the wider walk also takes molecules from `metabolic_network`. The build with
the root `CHEBI:75768` below shows the size of that effect: 2,034 molecules. The bucket
`hmdb` grows by 343 molecules, all `structure_key` rows. The HMDB build of
2026.11.27 and the new ChEBI content can both cause that. The HMDB claims
stay at 217,916.

**The overlap report.** It counts the molecules that each combination of
signals claims, before the collapse to one class. The right column gives the
figures of the plan, from release 2026.11.19 and the list of 45 roles.

| Signals | 2026.11.27, walk | Plan, 2026.11.19, list |
|---|---|---|
| `hmdb` | 212,054 | 212,071 |
| `metabolic_network` | 17,155 | 17,182 |
| `chebi_role` | 9,707 | 9,577 |
| `chebi_role`, `metabolic_network` | 2,534 | 2,507 |
| `chebi_role`, `hmdb`, `metabolic_network` | 2,227 | 2,211 |
| `hmdb`, `metabolic_network` | 1,850 | 1,866 |
| `chebi_role`, `hmdb` | 1,785 | 1,768 |

The walk now gives about the figures of the trial of the plan. So the
gate figures of the plan, 1,419,750 relations for all three buckets, are
again a fair guide. This report does not count the gates again.

**A build with the root `CHEBI:75768`.** The rule file `rules_mammal.yaml`
changes only the root. `chebi_role` then holds 1,689 molecules: 1,198
`entity` rows, 491 `structure_key` rows, and 491 `no_key` rows.
`metabolic_network` grows back to 18,013.

## Step 7: the audit samples

Each sample holds 50 `entity` rows of `chebi_role` of the full build, in the
order of `hash(entity_key)`.

- (b) Every `entity` row. The bucket holds 9,182.
- (h) The `entity` rows whose molecule carries a role below
  `CHEBI:77746` or `CHEBI:75768`. The bucket holds 1,198. The first 7 rows of
  (h) are the 7 rows of (b) with such a role.

### How the audit classifies a molecule

The script `audit.py` gives each molecule its label, identifier, InChIKey,
`supported_by`, roles, cross-references and ChEBI definition. It also looks
for an HMDB molecule with the same name or synonym, or one that shares a
ChEBI, KEGG or CAS identifier. `probe_hmdb.py` then looks up suspected
duplicates in HMDB by name. The classes are those of the earlier report, in
the same order:

| Class | Rule |
|---|---|
| Placeholder or class | A protein, a group, a polymer, or a molecule bound to a carrier protein. Also a lipid species without double-bond positions, and a name that leaves the ceramide or acyl open |
| Duplicate of an HMDB molecule | The same compound as an HMDB molecule: a charge form, a salt, or the same molecule under another structure |
| Drug or xenobiotic | A drug, a metabolite of a drug, a solvent, a pesticide or a breakdown product of one |
| Specific metabolite | A defined molecule that occurs in some organism, and none of the above |
| Unclear | A salt that HMDB does not list as its parent, a radical, or a name that does not fix the structure |

This audit puts a salt in the duplicates when HMDB lists the parent
molecule. The earlier report found no such salt, so the change does not alter
its figures.

### The counts

| Class | (b) `chebi_role` | (h) human or mammalian roles |
|---|---|---|
| Specific metabolite | 39 | 15 |
| Placeholder or class | 4 | 21 |
| Duplicate of an HMDB molecule | 1 | 7 |
| Drug or xenobiotic | 3 | 6 |
| Unclear | 3 | 1 |
| Rows with a human or mammalian role | 7 | 50 |
| Rows with the root role only | 22 | 0 |
| Rows without an InChIKey | 5 | 21 |

**Sample (b).** The 39 specific metabolites are mostly natural products. By
role, definition or name, about 17 come from plants, 10 from bacteria and 6
from fungi. Two are mouse lipids. Daphnia, a coral and a dinoflagellate give
one each, and one has no organism. 22 rows carry only the root role "metabolite".
Hordatine A is again the one duplicate. The four placeholders and classes are
ferrocytochrome c, phlorotannin and two lipid species. 3-bromo-1-propanol and
sulfadiazine hydroxylamine carry a mammalian role and are xenobiotic. The
figures agree with the sample of the earlier report: 45 specific of 50 then,
39 now. The wider walk adds species roles, and with them more lipid species
and xenobiotics.

**Sample (h).** All 21 placeholders and classes have no InChIKey. 12 of them
are lipid species, such as `PC 18:1/20:4` and `DG 40:4`. 8 of these 12 carry
the Papio hamadryas role. The others are two proteins and
decanoyl-[acyl-carrier protein]. Also the methyl group, two gangliosides with
an unnamed ceramide, dolichol of any chain length, a lysophosphatidylcholine
with an unnamed acyl, and alkylphosphonate. Of the 29 rows with an InChIKey,
15 are specific metabolites, 7 duplicates and 6 drugs or xenobiotics. The
last one is unclear: the heparan sulfate record of item 4 of the next
section. The duplicates are three
salts (oleoylcarnitine, 3-hydroxybutyrate and glycochenodeoxycholate),
methoxide, and three molecules where HMDB has a wrong structure.

**HMDB has a wrong structure in three cases of (h).** HMDB gives 14-HDoHE
(`HMDB0060044`) a cyclopropene ring. It gives 3-mercapto-2-methyl-1-butanol
(`HMDB0032373`) the SMILES `CCC(C)C(O)S`. It gives ganglioside GD1a
(d18:1/24:1) (`HMDB0011803`) an acetonyl group in place of the
N-acetyl group of GalNAc. Each structure has another connectivity than the
ChEBI molecule, so the group cannot join them. ChEBI cross-references the
first one to its HMDB accession.

## What looks wrong

1. **ChEBI gives proteins and groups a metabolite role.** Ferrocytochrome c
   and decanoyl-[acyl-carrier protein] carry human metabolite and mouse
   metabolite. Reduced ferredoxin and the methyl group carry mouse
   metabolite. A root below `CHEBI:25212` does not remove them.
2. **The human subtree holds drug metabolites.** Drug metabolite
   (`CHEBI:49103`) and human xenobiotic metabolite (`CHEBI:76967`) are below
   human metabolite. 184 rows of the restricted bucket carry only these two.
3. **A molecule is a subclass of a role.** ChEBI gives caudatoside
   (`NIERWTDCCGPRCQ-ZFVALMAISA-N`) a `subclass_of` edge to phytotoxin. The
   walk takes it as a role. It has no carrier, so no count changes. A filter
   of the walk to role terms would make the walk safe from such edges.
4. **One entity mixes a polymer and a disaccharide.** `HP_dp02_0003` has an
   InChIKey, but it also carries the ChEBI identifiers 14392 and 17174 of
   heparan sulfate, and ChEBI's definition of the polymer.
5. **Wrong HMDB structures hide duplicates.** The three cases of the section
   above. The other agent measures the HMDB status, so this report only
   names them. HMDB also gives the label "Cedrene-beta" to
   `HMDB0304654`, which has the structure of a ganglioside.

## Suggestion for the human decision

**This is a suggestion. The figures above inform the decision but do not
take it.**

1. **Do not switch on `chebi_role` as it is for a human gate.** The bucket
   adds 9,182 molecules to `hmdb`, and 78 per cent of the sample are real
   metabolites. But only 13 per cent of the rows carry a human or mammalian
   role. The rest are natural products of other organisms, and 40 per cent
   carry the generic role "metabolite" alone. For a gate that does not care
   about the organism, the bucket is fair as it is: few placeholders and few
   duplicates.
2. **For a human gate, restrict the root to `CHEBI:75768` and require an
   InChIKey.** The root change needs no code: one string in the rule file.
   It gives 1,198 `entity` rows. The InChIKey rule removes every placeholder
   and class of sample (h), and no other row. That leaves about 707 rows.
   The rule needs a new option of the ChEBI source, for example
   `require_structure_key: true`, a condition in `claim_chebi_role.sql`, and
   a test.
3. **Then decide on the drug metabolites.** About a quarter of the 707 rows
   are drug or xenobiotic metabolites. To remove them, the walk needs a list
   of roots to leave out, here `CHEBI:49103` and `CHEBI:76967`. That is
   another option of the rule file, a change of the walk and a test.
4. **Expect duplicates in the restricted bucket.** 7 of 29 keyed rows of (h)
   are HMDB molecules under another identity. A link of salts to their parent
   would catch three. A correction of the wrong HMDB structures is a matter
   for the HMDB input, not for this subset.
5. **Keep the walk, not the list.** The walk reaches every carrier of the
   list, and 117 more through toxin and volatile oil roles. It needs no
   label rule.

## How to repeat the measurement

Run the scripts from
`/tmp/claude-1029/-workspace-instances-metsig1/97c1ab55-a53d-4d90-8da3-8c1f34a8d599/scratchpad/stage2-chebi/`
with the Python of `../venv/`, in this order:

| Script | What it does | Time |
|---|---|---|
| `walk_vs_list.py` | Runs the walk and the list, counts carriers and sub-roles. Writes `walk_vs_list.json`, `walk_roles.parquet`, `hm_roles.parquet` and `hcr.parquet`. | 2 s |
| `explore2.py` | Compares the list on 2026.11.19 and 2026.11.27, and lists the subtrees | 5 s |
| `build.py` | Writes `rules_full.yaml`, `rules_mammal.yaml` and `rules_hmdb_chebi.yaml`, and builds each. Writes `build_*.json` and `mc_<rules>.parquet`. | 12 s |
| `profile.py` | Counts the bucket rows for each sub-role | 3 s |
| `audit.py` | Draws the samples (b) and (h). Writes `audit_b.json` and `audit_h.json`. | 8 s |
| `probe_hmdb.py` | Looks up HMDB records by name | 2 s |
| `tables.py` | Holds the classes read by hand, and prints the counts and the sample tables | 1 s |

## Appendix: the sample tables

The columns are the label, the namespace and identifier of the entity, the
InChIKey and `supported_by`. Then come the roles of the walk that the
molecule carries, whether one of them is human or mammalian, and the class. A
dash means no value.

<!-- vale off -->

### (b) `chebi_role`, `entity` rows

| # | Label | Identifier | InChIKey | Supported by | Roles | Mammalian | Class |
|---|---|---|---|---|---|---|---|
| 1 | ferrocytochrome c | chebi `CHEBI:16928` | - | chebi, kegg, metatlas | mouse metabolite, S. cerevisiae metabolite, human metabolite | yes | Placeholder (protein) |
| 2 | 1-(2-oxindole-3-acetyl)-beta-D-glucose | inchikey | `MLVJTFNIARPUAO-LCNORDJOSA-N` | chebi | Arabidopsis thaliana metabolite | - | Specific |
| 3 | staphyloferrin B | inchikey | `SIAZVTIHOHTZDD-PWJLMRLQSA-N` | chebi | bacterial metabolite | - | Specific |
| 4 | mitchellene B | inchikey | `BFNUVKSKPNBQRT-IGLSZRPSSA-N` | chebi | plant metabolite | - | Specific |
| 5 | (24S)-24,25-epoxycucurbitadienol | inchikey | `FKLMMFZHGLDEMF-CUGODZHQSA-N` | chebi, rhea | plant metabolite | - | Specific |
| 6 | Stemofuran P | inchikey | `LVVYOWMHRJEABT-UHFFFAOYSA-N` | chebi | metabolite | - | Specific |
| 7 | chaetomanone | inchikey | `ZVLQJGNLXSDFDN-IYFFJJEWSA-N` | chebi | metabolite | - | Specific |
| 8 | O-(hydroxypropionyl)carnitine | chebi `CHEBI:86072` | - | chebi | metabolite | - | Unclear (position not given; HMDB holds 3-hydroxypropionylcarnitine) |
| 9 | 2''-O-(4-sinapoylarabinosyl)isoscoparin | inchikey | `YOFXNYRWIIFBRE-QOLNKJPBSA-N` | chebi | metabolite | - | Specific |
| 10 | PC 18:2/20:1 | goslin `sn_position:PC 18:2/20:1` | - | chebi | Papio hamadryas metabolite | yes | Class (lipid species) |
| 11 | aureothin | inchikey | `GQKXCBCSVYJUMI-WACKOAQBSA-N` | chebi, kegg, rhea | bacterial metabolite | - | Specific |
| 12 | 7,4'-dihydroxy-8,3'-dimethoxyflavanone | inchikey | `UKCAIYJHVFWSRN-UHFFFAOYSA-N` | chebi | metabolite | - | Specific |
| 13 | GlcCer(1) 18:0;3OH/22:0 | inchikey | `BCRCIDJHDMYLDJ-WGYZNKDKSA-N` | chebi | mouse metabolite | yes | Specific |
| 14 | 2-epi-5-epi-valiolone 7-phosphate(2-) | inchikey | `WBHJALJLPRKCBI-MVIOUDGNSA-L` | chebi, kegg, rhea | bacterial metabolite | - | Specific |
| 15 | 2-(2-hydroxyphenyl)ethanol | inchikey | `ABFCOJLLBHXNOU-UHFFFAOYSA-N` | chebi | metabolite | - | Specific |
| 16 | 6a,7-didehydro-1,9,10-trihydroxy-2-methoxy-6-methylaporphinium trifluoroacetate | inchikey | `SBKNQDPIZMRYQD-UHFFFAOYSA-N` | chebi | metabolite | - | Unclear (salt) |
| 17 | 3'-L-tyrosyl-AMP | inchikey | `IEDNIXUXZOSRER-URQYDQELSA-N` | chebi | Mycoplasma genitalium metabolite | - | Specific |
| 18 | phlorotannin | chebi `CHEBI:71222` | - | chebi | algal metabolite | - | Class (oligomers) |
| 19 | aspyridone A | inchikey | `LIBKJCZUBOETPB-RYUDHWBXSA-N` | chebi | fungal metabolite | - | Specific |
| 20 | pentalenolactone D | inchikey | `MRLXXQBBRNRWDA-LIEMUPCESA-N` | chebi | bacterial metabolite | - | Specific |
| 21 | 14-hydroxyplatensic acid methyl ester | inchikey | `FBBBYTFYYNYGMH-MKJSHFIQSA-N` | chebi | metabolite | - | Specific |
| 22 | (21S)-21,22-epoxypolypoda-8(26)-13,17-trien-3beta-ol | inchikey | `YOXYJGZFQYXPHZ-AEBUGLRVSA-N` | chebi, kegg, rhea | plant metabolite | - | Specific |
| 23 | 3-bromo-1-propanol | inchikey | `RQFUZUMFPRMVDX-UHFFFAOYSA-N` | chebi | mammalian metabolite | yes | Xenobiotic |
| 24 | DG 40:4 | goslin `species:DG 40:4` | - | chebi | mouse metabolite | yes | Class (lipid species) |
| 25 | (3S,4E)-3-methyldec-4-en-1-yl sulfate | inchikey | `QNUMXENJPAUAPZ-WSKFYRRCSA-M` | chebi | Daphnia pulex metabolite | - | Specific |
| 26 | venturamide B | inchikey | `BJJVHBJYOXGRRO-XEQNPHJVSA-N` | chebi | metabolite | - | Specific |
| 27 | coptisine | inchikey | `XYHOBCMEDLZUMP-UHFFFAOYSA-N` | chebi | metabolite | - | Specific |
| 28 | Decarbamoylgonyautoxin-4 | inchikey | `AXKNFGCPJPRVCA-RGDLXGNYSA-N` | chebi | marine metabolite | - | Specific |
| 29 | hordatine A | inchikey | `KVYNYRIOAYQBFK-AIIPJEMGSA-N` | chebi, kegg | metabolite | - | Duplicate (`HMDB0030461`) |
| 30 | chetracin A | inchikey | `PFPUKJNFELCZAQ-AAZDRVSSSA-N` | chebi | metabolite | - | Specific |
| 31 | hydroxyversicolorone(1-) | inchikey | `JGXCLZAVTLWCBF-LXGOIASLSA-M` | chebi, kegg, rhea | metabolite | - | Specific |
| 32 | 3-methylindolyl precursor of nosiheptide | inchikey | `MMGRHNBDIMLRQF-AKJSKULZSA-N` | chebi | bacterial metabolite | - | Specific |
| 33 | Turnagainolide B | inchikey | `SMWGFJCKWZCJNQ-QWGWBJRGSA-N` | chebi | metabolite | - | Specific |
| 34 | (2S)-5-hydroxy-7-methoxy-6,8-dimethylflavanone | inchikey | `KZNAGLGLKRUIPD-AWEZNQCLSA-N` | chebi | plant metabolite | - | Specific |
| 35 | Euscaphic acid C | inchikey | `FJBWUBWZSSBDEN-XVNGACRRSA-N` | chebi | metabolite | - | Specific |
| 36 | Maoesin E, (rel)- | inchikey | `YLAQTEHVFWEICH-MQAPUGLSSA-N` | chebi | metabolite | - | Specific |
| 37 | 8-hydroxy-12-oxoabieta-9(11),13-dien-20-oic 8,20-lactone | inchikey | `UFXFRAKIDXOMGU-CWFSZBLJSA-N` | chebi | plant metabolite | - | Specific |
| 38 | sulfadiazine hydroxylamine | inchikey | `BHNNTIGNLIYDOR-UHFFFAOYSA-N` | chebi | drug metabolite | yes | Drug (metabolite of sulfadiazine) |
| 39 | mitchellene D | inchikey | `KRUUDBOQHFKKKU-SNCMFLLNSA-N` | chebi | plant metabolite | - | Specific |
| 40 | (2E,6R)-8-hydroxy-2,6-dimethyl-2-octenoic acid | inchikey | `CDAJACJYYZHNJA-CBVZUMJXSA-N` | chebi | metabolite | - | Specific |
| 41 | aureothricin | inchikey | `UGZYFXMSMFMTSM-UHFFFAOYSA-N` | chebi | bacterial metabolite | - | Specific |
| 42 | eucalyptal C | inchikey | `LLCSUGKJCVVJBS-LQEBHLGXSA-N` | chebi | metabolite | - | Specific |
| 43 | dimethenamid OXA | inchikey | `HOYCASTVMCEOTP-UHFFFAOYSA-N` | chebi | marine xenobiotic metabolite | - | Xenobiotic (herbicide breakdown) |
| 44 | kurarinol | inchikey | `XMUPAAIHKAIUSU-QRQCRPRQSA-N` | chebi | plant metabolite | - | Specific |
| 45 | Athrolide B | inchikey | `BHEXSLGHYNIAGL-WDCRBBJTSA-N` | chebi | metabolite | - | Specific |
| 46 | 22alpha,23alpha-epoxy-solanida-1,4-dien-3-one, (rel)- | inchikey | `UGAJCJIZYPKGTJ-WUOUAPBYSA-N` | chebi | metabolite | - | Specific |
| 47 | Conoideocrellide B | inchikey | `XUGFNOLBANLPEQ-JNVOQPTMSA-N` | chebi | metabolite | - | Specific |
| 48 | 1,2-bis-[(9Z,11Z)-octadecadienoyl]-sn-glycero-3-phosphoethanolamine | inchikey | `CMCPQFKZYMRAKT-UHSBSJHCSA-N` | chebi | mouse metabolite | yes | Specific |
| 49 | glycyl radical | inchikey | `DXMZNKNJOBUNRO-UHFFFAOYSA-N` | chebi | fundamental metabolite | - | Unclear (radical) |
| 50 | Brialalepolide C | inchikey | `CEVYVRCWJDJLEW-LXIKAAAESA-N` | chebi | metabolite | - | Specific |

### (h) `chebi_role`, `entity` rows with a human or mammalian role

| # | Label | Identifier | InChIKey | Supported by | Roles | Mammalian | Class |
|---|---|---|---|---|---|---|---|
| 1 | ferrocytochrome c | chebi `CHEBI:16928` | - | chebi, kegg, metatlas | mouse metabolite, S. cerevisiae metabolite, human metabolite | yes | Placeholder (protein) |
| 2 | PC 18:2/20:1 | goslin `sn_position:PC 18:2/20:1` | - | chebi | Papio hamadryas metabolite | yes | Class (lipid species) |
| 3 | GlcCer(1) 18:0;3OH/22:0 | inchikey | `BCRCIDJHDMYLDJ-WGYZNKDKSA-N` | chebi | mouse metabolite | yes | Specific |
| 4 | 3-bromo-1-propanol | inchikey | `RQFUZUMFPRMVDX-UHFFFAOYSA-N` | chebi | mammalian metabolite | yes | Xenobiotic |
| 5 | DG 40:4 | goslin `species:DG 40:4` | - | chebi | mouse metabolite | yes | Class (lipid species) |
| 6 | sulfadiazine hydroxylamine | inchikey | `BHNNTIGNLIYDOR-UHFFFAOYSA-N` | chebi | drug metabolite | yes | Drug (metabolite of sulfadiazine) |
| 7 | 1,2-bis-[(9Z,11Z)-octadecadienoyl]-sn-glycero-3-phosphoethanolamine | inchikey | `CMCPQFKZYMRAKT-UHSBSJHCSA-N` | chebi | mouse metabolite | yes | Specific |
| 8 | PE 16:1/18:0 | goslin `sn_position:PE 16:1/18:0` | - | chebi | Papio hamadryas metabolite | yes | Class (lipid species) |
| 9 | PC 18:1/20:4 | goslin `sn_position:PC 18:1/20:4` | - | chebi | Papio hamadryas metabolite | yes | Class (lipid species) |
| 10 | crassicauline A | inchikey | `GAZDXIGXYWVWQX-QVAFJCLZSA-N` | chebi | plant metabolite, human urinary metabolite | yes | Xenobiotic (plant alkaloid) |
| 11 | GlcCer(1) 18:0;3OH/16:0 | inchikey | `BLGKYYVFGMKTEZ-UKDORTLVSA-N` | chebi, rhea | mouse metabolite | yes | Specific |
| 12 | 7-(2-amino-2-carboxyethyl)-5-hydroxy-3,4-dihydro-2H-1,4-benzothiazin-3-one | inchikey | `IRQPIGVFKBQNEM-UHFFFAOYSA-N` | chebi | human metabolite | yes | Specific |
| 13 | (13R)-S-glutathionyl-(14S)-hydroxy-(4Z,7Z,9E,11E,16Z,19Z)-docosahexaenoic acid | inchikey | `NWYKIASJHFGEAN-CPMAONNBSA-N` | chebi | human metabolite | yes | Specific |
| 14 | (alpha-D-mannosyl)(5)-beta-D-mannosyldiacetylchitobiosyldiphosphodolichol | chebi `CHEBI:37634` | - | chebi | human metabolite | yes | Class (dolichol of any length) |
| 15 | alpha-Neu5Ac-(2->3)-beta-D-Gal-(1->3)-beta-D-GalNAc-(1->4)-beta-D-Gal-(1->4)-beta-D-Glc-(1<->1')-Cer | chebi `CHEBI:62370` | - | chebi | mouse metabolite | yes | Placeholder (unnamed ceramide) |
| 16 | 7,16,17-trihydroxy-(14E,19Z)-docosa-8,10,12,14,19-pentaenoic acid | inchikey | `LDVCZYSUILUYKN-OSLFFFSGSA-N` | chebi | human xenobiotic metabolite | yes | Specific |
| 17 | (R)-oleoylcarnitine hydrochloride | inchikey | `JPDKUMDAMZLOJR-KFQPROOXSA-N` | chebi | human metabolite | yes | Duplicate (salt, `HMDB0005065`) |
| 18 | 1,4-dimethylcyclohexane | inchikey | `QRMPKOFEUHIBNM-UHFFFAOYSA-N` | chebi | plant metabolite, human metabolite | yes | Specific |
| 19 | PC 22:6/22:6 | goslin `sn_position:PC 22:6/22:6` | - | chebi | Papio hamadryas metabolite | yes | Class (lipid species) |
| 20 | PC 14:0/18:2 | goslin `sn_position:PC 14:0/18:2` | - | chebi | Papio hamadryas metabolite | yes | Class (lipid species) |
| 21 | 3''-hydroxycannabidiol | inchikey | `OEOXRPUFTSADIH-UQJFVLDMSA-N` | chebi | human xenobiotic metabolite | yes | Drug (metabolite of cannabidiol) |
| 22 | methoxide | inchikey | `NBTOZLQBSIZIKS-UHFFFAOYSA-N` | chebi | human metabolite | yes | Duplicate (charge form, `HMDB0001875`) |
| 23 | (7Z,10Z,13Z,15E,17S,19Z)-17-hydroperoxydocosapentaenoic acid | inchikey | `PQYYUXHETXYOCV-DKMPAXASSA-N` | chebi | human xenobiotic metabolite | yes | Specific |
| 24 | Hex3HexNAcNeuAc2Cer 18:1;O2/24:1 | inchikey | `LZIGFMROBFWALA-GDLFOJDGSA-N` | chebi | mouse metabolite | yes | Duplicate (`HMDB0011803`, wrong HMDB structure) |
| 25 | 1-octadecenoyl-2-octadecadienoyl-sn-glycero-3-phospho-(1'-sn-glycerol) | goslin `sn_position:PG 18:1/18:2` | - | chebi | mouse metabolite | yes | Class (lipid species) |
| 26 | (13R)-S-cysteinyl-(14S)-hydroxy-(4Z,7Z,9E,11E,16Z,19Z)-docosahexaenoic acid | inchikey | `GIIVKOKEBBFSDI-AIYHDIOKSA-N` | chebi | human metabolite | yes | Specific |
| 27 | 1-hexadecyl-2-formyl-sn-glycero-3-phosphocholine | inchikey | `BBASNTQPDXPDFB-RUZDIDTESA-N` | chebi | mouse metabolite | yes | Specific |
| 28 | decanoyl-[acyl-carrier protein] | chebi `CHEBI:4349` | - | chebi, kegg, metatlas | mouse metabolite, E. coli metabolite, human metabolite | yes | Placeholder (carrier protein) |
| 29 | HP_dp02_0003 | inchikey | `KNKWELGWBDSUOG-PNRHKHKDSA-N` | chebi, kegg | mouse metabolite | yes | Unclear (polymer record with a disaccharide name) |
| 30 | sodium 3-hydroxybutyrate | inchikey | `NBPUSGBJDWCHKC-UHFFFAOYSA-M` | chebi | human metabolite | yes | Duplicate (salt, `HMDB0000011`) |
| 31 | Cer 18:1(4E);1OH,3OH/16:0;2OH | inchikey | `YWZGQPJRQXHVQX-HOXYDHSRSA-N` | chebi | mouse metabolite | yes | Specific |
| 32 | GM2 | chebi `CHEBI:51013` | - | chebi, kegg, metatlas | mouse metabolite | yes | Placeholder (unnamed ceramide) |
| 33 | (+/-)-14-HDoHE | inchikey | `ZNEBXONKCYFJAF-BGKMTWLOSA-N` | chebi | human xenobiotic metabolite | yes | Duplicate (`HMDB0060044`, wrong HMDB structure) |
| 34 | bromobenzene | inchikey | `QARVLSVVCXYDNA-UHFFFAOYSA-N` | chebi, kegg, metatlas | mouse metabolite | yes | Xenobiotic |
| 35 | 1,2-dioctadecadienoyl-sn-glycero-3-phospho-(1'-sn-glycerol) | chebi `CHEBI:85001` | - | chebi | mouse metabolite | yes | Class (lipid species) |
| 36 | methyl group | chebi `CHEBI:32875` | - | chebi, rhea | mouse metabolite | yes | Placeholder (group) |
| 37 | PC 15:1/18:1 | goslin `sn_position:PC 15:1/18:1` | - | chebi | Papio hamadryas metabolite | yes | Class (lipid species) |
| 38 | PC 15:0/18:2 | goslin `sn_position:PC 15:0/18:2` | - | chebi | Papio hamadryas metabolite | yes | Class (lipid species) |
| 39 | 2-tetradecanoyl-sn-glycero-3-phosphocholine | inchikey | `IGZPHNNYPPAPLA-OAQYLSRUSA-N` | chebi, rhea | plant metabolite, human metabolite | yes | Specific |
| 40 | sodium glycochenodeoxycholate | inchikey | `AAYACJGHNRIFCT-YRJJIGPTSA-M` | chebi | human metabolite | yes | Duplicate (salt, `HMDB0000637`) |
| 41 | 2-acyl-sn-glycero-3-phosphocholine(1+) | chebi `CHEBI:16728` | - | chebi, kegg | mouse metabolite | yes | Class (unnamed acyl) |
| 42 | 6-phospho-2-dehydro-D-gluconic acid | inchikey | `ZKUSPPOKDDRMIU-JJYYJPOSSA-N` | chebi, kegg | mouse metabolite | yes | Specific |
| 43 | 2-methyl-3-sulfanylbutan-1-ol | inchikey | `RFMHFOPFUZZBAD-UHFFFAOYSA-N` | chebi, rhea | plant metabolite, human metabolite | yes | Duplicate (`HMDB0032373`, wrong HMDB structure) |
| 44 | CHEBI:151275 | inchikey | `NSXKHTGRJJVIIM-UHHJCXHESA-N` | chebi | human metabolite | yes | Specific |
| 45 | deoxycholic acid 24-O-(beta-D-glucuronide) | inchikey | `MZVCLCMWDVOQAH-GYSLYLNXSA-N` | chebi | human urinary metabolite | yes | Specific |
| 46 | reduced ferredoxin | chebi `CHEBI:17513` | - | chebi, kegg, metatlas | mouse metabolite | yes | Placeholder (protein) |
| 47 | alkylphosphonate | chebi `CHEBI:60983` | - | chebi | mouse metabolite | yes | Class |
| 48 | cyclosporin A sulfate | inchikey | `CQDXGGWUGMQJHS-CGLBZJNRSA-M` | chebi, metatlas, recon3d | drug metabolite, human xenobiotic metabolite | yes | Drug (metabolite of cyclosporin A) |
| 49 | PC 18:1/20:0 | goslin `sn_position:PC 18:1/20:0` | - | chebi | Papio hamadryas metabolite | yes | Class (lipid species) |
| 50 | PC 36:4 | goslin `species:PC 36:4` | - | chebi | human metabolite | yes | Class (lipid species) |

<!-- vale on -->
