# Stage two of the metabolite classes: measurements

2026-11-27

This report gives the measurements of stage two of the plan
`saezverse/projects/omnipath/plans/2026-11-metabolite-classification-logic.md`.
It covers step 8, the contribution of each reaction resource. It also covers
the part of step 7 that the current release allows, the audit samples.

**Result: the four reaction resources add 53,937 MetaLinksDB-shaped relations
to the gate `[hmdb, chebi_role]`, a gain of 4.0 per cent.** In the audit
sample of `metabolic_network`, 17 of 50 molecules are generic placeholders. All
17 have no InChIKey. The variants that the groups pull into `hmdb` are real
duplicates of HMDB molecules in 38 of 50 cases, and drugs in 10.

## What this report measures

| Item | Value |
|---|---|
| Release | `2026.11.19`, read from `/workspace/instances/metsig1/data-2026.11.19` |
| Code | `pypath` at `051a2639b`. The subset code last changed in `70e6006e4`. |
| Build | `parquet_build.build_tables` with all three signals of the rule file |
| Scripts | In `/tmp/claude-1029/-workspace-instances-metsig1/97c1ab55-a53d-4d90-8da3-8c1f34a8d599/scratchpad/stage2/`, listed at the end |

The rule file of each build is the packaged one with the stage-two rows
uncommented. Signal 3 changes from build to build: no resource, all four, each
resource alone, and each resource removed. That makes ten builds. Each takes
about 2.5 seconds.

**The MetaLinksDB-shaped relation set.** A relation is in the set when its
predicate is `interacts_with`, `affects` or `associated_with`, one end is a
chemical, and the other end has a type in `_PARTNER_TYPES` of
`network_views/_definitions.py`. The type of an entity is the one that
`release_entity` gives it: the first resource in name order that has one. The
set holds 6,261,809 rows of `release_relation`. The rows are relations of each
resource, not distinct pairs. The figures of the plan count the same rows.

**The gate.** A gate with buckets B admits a relation when one of its ends
has a row in `metabolite_class` with a bucket in B.

The set reads every resource. BindingDB is in it, although the `metalinksdb`
preset leaves BindingDB out.

## The buckets

| Bucket | Molecules | Through a group | Without an InChIKey |
|---|---|---|---|
| `hmdb` | 272,030 | 54,114 | 21 |
| `chebi_role` | 4,608 | 768 | 257 |
| `metabolic_network` | 18,494 | 3,957 | 6,071 |

The ChEBI role walk claims 5,317 molecules. Today it reaches the direct
carriers of `CHEBI:25212` only, so the figures of `chebi_role` are provisional.

## The gates, measured again

| Gate | Relations | Figure of the plan |
|---|---|---|
| `[hmdb]` | 1,349,547 | 1,349,547 |
| `[hmdb, chebi_role]` | 1,354,756 | |
| `[chebi_role]` | 5,212 | |
| `[metabolic_network]` | 53,991 | 38,491 |
| `[hmdb, chebi_role, metabolic_network]` | 1,408,693 | 1,419,750 |

The figure for `[hmdb]` holds. The other figures of the plan came from a trial
with a list of 45 ChEBI roles. That list reached more molecules than the role
walk of today. So in the trial more molecules moved to `chebi_role`, and fewer
stayed in `metabolic_network`.

The gate `[hmdb, chebi_role]` admits the same 1,354,756 relations with and
without signal 3. This agrees with the build: signal 3 cannot move a molecule
into an earlier bucket.

## What each reaction resource contributes

**Molecules.** The four resources claim 23,766 molecules. Of these, 14,537 end
in `metabolic_network` as `entity` rows. The others are in an earlier bucket.

| Resource | Claims | Claims alone among the four | Ends in `metabolic_network` | Of these, alone | Of these, without an InChIKey |
|---|---|---|---|---|---|
| `kegg` | 9,533 | 5,373 | 5,273 | 3,515 | 1,892 |
| `metatlas` | 4,119 | 1,544 | 2,195 | 1,123 | 1,607 |
| `recon3d` | 2,755 | 705 | 1,355 | 506 | 851 |
| `rhea` | 14,131 | 10,326 | 8,404 | 6,769 | 2,542 |

Of the 14,537 molecules in `metabolic_network`, 6,769 come from Rhea alone and
3,515 from KEGG alone. Metatlas and Recon3D share 780. Recon3D adds 506
molecules that no other resource claims. Metatlas and Recon3D also give the
largest share of molecules without an InChIKey: 73 and 63 per cent.

**Relations.** The gain is the count for the gate
`[hmdb, chebi_role, metabolic_network]` minus 1,354,756, the count for
`[hmdb, chebi_role]`. "Alone" builds signal 3 from one resource. "Unique" is the
loss when the resource leaves the set of four.

| Resource | Gate with the resource alone | Gain alone | Gate without the resource | Unique gain |
|---|---|---|---|---|
| `kegg` | 1,384,131 | 29,375 | 1,400,954 | 7,739 |
| `metatlas` | 1,376,628 | 21,872 | 1,398,580 | 10,113 |
| `recon3d` | 1,361,808 | 7,052 | 1,407,712 | 981 |
| `rhea` | 1,379,483 | 24,727 | 1,401,010 | 7,683 |
| all four | 1,408,693 | 53,937 | 1,354,756 | 53,937 |

The unique gains add up to 26,516. The other 27,421 relations need a molecule
that two or more resources claim, so no single removal loses them.

## Which resources the added relations come from

The resource here is the one that states the relation, not the reaction
resource that claims the chemical.

| Resource of the relation | Added relations, all four |
|---|---|
| `pathwaycommons` | 26,328 |
| `metatlas` | 10,073 |
| `stitch` | 3,579 |
| `brenda` | 2,711 |
| `tcdb` | 2,658 |
| `rhea` | 1,772 |
| `intact` | 1,494 |
| `recon3d` | 1,331 |
| `pmidb` | 1,205 |
| `chembl` | 1,097 |
| 12 others | 1,689 |

The unique gain of each reaction resource, by the resource of the relation:

| Reaction resource | Largest parts of the unique gain |
|---|---|
| `kegg` | `pathwaycommons` 4,051, `pmidb` 1,205, `tcdb` 708, `brenda` 659, `stitch` 568 |
| `metatlas` | `metatlas` 5,349, `pathwaycommons` 4,344, `stitch` 295 |
| `recon3d` | `recon3d` 582, `imm1415` 249, `pathwaycommons` 143 |
| `rhea` | `pathwaycommons` 3,526, `rhea` 1,355, `stitch` 1,123, `brenda` 764 |

The reaction resources also state relations of their own: the transport
relations that the build projects from their reactions. Metatlas, Rhea and
Recon3D give 13,176 of the added relations this way. Half the unique gain of
metatlas is its own relations.

**A few molecules carry most of the gain.** The 53,937 added relations have
2,944 molecules of `metabolic_network` at their chemical end.

| Molecule | Added relations | Claimed by | Note |
|---|---|---|---|
| N(alpha)-methyl-L-tryptophan | 6,893 | `kegg`, `rhea` | 6,886 from Pathway Commons, which labels this InChIKey "Irosustat" |
| arsenite(3-) | 4,218 | `metatlas` | 3,977 from Pathway Commons |
| 2-methylcholine | 1,666 | `kegg`, `rhea` | |
| selenium atom | 1,607 | `kegg` | |
| lipid (`CHEBI:18059`) | 1,317 | `kegg` | A class. 1,205 from PMIDB |
| 3,3',5,5'-tetrabromobisphenol A | 1,134 | `kegg` | A xenobiotic |
| H+ | 847 | `metatlas` | A fallback without an InChIKey. HMDB holds the proton |
| protein polypeptide chain | 793 | `kegg`, `metatlas` | A placeholder |

**The InChIKey splits the gain.** Molecules without an InChIKey carry 21,297 of
the added relations, through 887 molecules. Among the molecules that one
resource alone claims, the relations split so:

| Claimed by | Relations, with an InChIKey | Relations, without an InChIKey | Molecules |
|---|---|---|---|
| `kegg` alone | 5,611 | 3,056 | 618 |
| `metatlas` alone | 4,574 | 5,616 | 183 |
| `recon3d` alone | 101 | 896 | 117 |
| `rhea` alone | 2,836 | 3,133 | 454 |
| two or more | 15,436 | 8,603 | 733 |

A relation counts in more than one row when both of its ends are molecules of
`metabolic_network`.

## The audit samples

Each group holds 50 molecules of the full build, in the order of
`hash(entity_key)`. The groups:

- (a) `metabolic_network`, `entity` rows. The bucket holds 14,537 of them.
- (b) `chebi_role`, `entity` rows. The bucket holds 3,840 of them. This sample
  is provisional, because the ChEBI parser fix widens the role walk.
- (c) `hmdb`, `structure_key` rows: the variants that a group pulls in. The
  bucket holds 54,114 of them.

### How the audit classifies a molecule

The script `audit.py` gives each molecule its label, identifier, InChIKey,
`supported_by`, the resources that hold it, and its cross-references. It also
looks for an HMDB molecule with the same label, or one that shares a ChEBI,
KEGG or CAS identifier. The class then comes from these rules and from a
reading of each row:

| Class | Rule |
|---|---|
| Placeholder or class | The label names a protein, a residue or group of a macromolecule, a polymer, a molecule bound to a carrier protein (`[ACP]`), or a lipid class without defined chains |
| Duplicate of an HMDB molecule | The same compound as an HMDB molecule. In (a) and (b), the same label or a shared cross-reference, then a check of the structure. In (c), a charge form, a stereo form, or the stereo-defined form of an HMDB entry without stereo |
| Drug or xenobiotic | A drug, a drug candidate, a radiotracer, a pesticide or a breakdown product of one |
| Specific metabolite | A defined molecule that occurs in some organism, and none of the above |
| Unclear | A salt with an isolation counterion, or a name that does not fix the structure |

A row takes the first class of the table that applies. A drug that is also a
variant of an HMDB drug counts as a drug.

The search for duplicates finds only exact labels and shared identifiers. So
the count of duplicates in (a) and (b) is a lower limit.

### The counts

| Class | (a) `metabolic_network` | (b) `chebi_role` | (c) `hmdb` variants |
|---|---|---|---|
| Specific metabolite | 28 | 45 | 2 |
| Placeholder or class | 17 | 2 | 0 |
| Duplicate of an HMDB molecule | 2 | 1 | 38 |
| Drug or xenobiotic | 3 | 0 | 10 |
| Unclear | 0 | 2 | 0 |

**Group (a).** The 17 placeholders are proteins, residues of RNA, DNA and
peptides, and polymers such as glycogen. The others are an acyl carrier
protein and glycosphingolipids with an unnamed ceramide. All 17 have no InChIKey. The
33 other molecules include 8 without an InChIKey. Two of those are fallbacks
of molecules that HMDB holds: L-dopachrome and indoxyl. Most specific
metabolites are bacterial, plant or fungal, and often a charge form that HMDB
does not list.

**Group (b).** Nearly every carrier is a natural product of a plant, a fungus,
a bacterium or a marine organism. ChEBI gives these molecules only the generic
role "metabolite". Hordatine A is in HMDB under another connectivity, with the
same name and KEGG identifier. Two lipid names give a class, not one molecule.

**Group (c).** Most of the 38 duplicates are lipids with another stereo
descriptor. The others are charge forms, such as L-tryptophanium, and
stereo-defined forms of HMDB entries without stereo. The 10 drugs are
variants of drugs that HMDB lists, for example etrasimod, 17-AAG and
cefmetazole. Two molecules are other stereoisomers that are metabolites of
their own: sn-glycerol 1-phosphate, the enantiomer of HMDB's glycerol
3-phosphate, and (-)-neoisodihydrocarveol.

## What looks wrong

1. **Fallbacks of HMDB molecules enter `metabolic_network`.** A fallback has no
   structure key, so the groups cannot place it with its HMDB molecule.
   Examples are L-dopachrome (`CHEBI:15772`), indoxyl (`MAM03702`), H+ and K+
   of metatlas. H+ alone adds 847 relations.
2. **Placeholders enter `metabolic_network`.** A reaction states proteins,
   residues and classes as participants, and signal 3 takes every participant.
   A rough label match finds 1,973 placeholder labels among the 14,537 rows.
   The sample suggests some 34 per cent, so the label match finds too few.
3. **One wrong identity carries 13 per cent of the gain.** Pathway Commons gives
   the InChIKey of N(alpha)-methyl-L-tryptophan the name "Irosustat". Irosustat
   is a drug of another structure. The 6,886 relations come from that input,
   not from the subset.
4. **`via_entity_key` can name a less useful HMDB entry.** L-tryptophanium
   points to (±)-tryptophan, `HMDB0030396`, not L-tryptophan, because the rule
   takes the lowest key. That is correct by the rule, but a reader can mistake
   it.
5. **A cross-reference of HMDB is wrong.** HMDB gives "Prenyl glucoside" the
   KEGG identifier `C09808` of phellamurin. This does not change the subset,
   because signal 1 reads no cross-reference.

## Suggestions for the decisions of stage two

These are suggestions for a human decision. The figures above inform the
decision but do not take it.

1. **Keep Rhea and KEGG.** They give the most molecules that no other resource
   claims, 6,769 and 3,515. They give unique gains of 7,683 and 7,739
   relations. Of the relations of the molecules that only KEGG claims, 65 per
   cent have an InChIKey. For Rhea, the figure is 48 per cent.
2. **Keep metatlas only with a fallback filter.** It has the largest unique
   gain, 10,113 relations. Half are its own transport relations. The molecules
   that only metatlas claims carry 10,190 added relations, and 5,616 of these
   go to molecules without an InChIKey. These include ion fallbacks that
   duplicate HMDB molecules.
3. **Recon3D adds little.** Its unique gain is 981 relations, and 582 of these
   are its own. Metatlas holds most of its molecules. Leaving it out costs
   0.07 per cent of the gate.
4. **Add a placeholder filter to signal 3.** A requirement for an InChIKey
   removes every placeholder of the sample. It also removes 8 of 33 other
   molecules, 2 of them duplicates. On the bucket, it removes 6,071 rows and
   21,297 relations. A filter on ChEBI classes such as macromolecule, residue
   and group keeps more real molecules, but needs the ChEBI parser fix first.
5. **Link the fallbacks to HMDB.** A fallback with a ChEBI or KEGG identifier
   that an HMDB molecule also carries can take the bucket `hmdb`. That needs a
   rule of its own and a test.
6. **Decide `chebi_role` after the parser fix.** Today it adds 5,209 relations
   to `[hmdb]`, and its carriers are mostly natural products of other
   organisms. Draw sample (b) again from the wider walk.
7. **Leave the groups of `hmdb` as they are.** The variants are duplicates or
   drugs of HMDB, and no variant in the sample is a placeholder. The drugs come
   from the breadth of HMDB. The status filter of step 9 is the place to act
   on them.

## How to repeat the measurement

The scripts are in
`/tmp/claude-1029/-workspace-instances-metsig1/97c1ab55-a53d-4d90-8da3-8c1f34a8d599/scratchpad/stage2/`.
Run them from that directory with `/workspace/instances/metsig1/pypath/.venv/bin/python`,
in this order:

| Script | What it does | Time |
|---|---|---|
| `measure.py` | Writes the ten rule files `rules_*.yaml` and builds each variant. Builds the relation set and counts the gates. Writes `results.json`, `mc_full.parquet`, `claim_full.parquet` and `added.parquet`. | 31 s |
| `audit.py` | Draws the three samples from `mc_full.parquet`. Writes `audit_<group>.json`. | a few minutes |
| `probe_top.py` | Shows where the relations of the top molecules come from | about 1 min |
| `probe_split.py` | Splits the gain by the claiming resources and by InChIKey | 1 s |

## Appendix: the sample tables

The columns are the label, the namespace and identifier of the entity, the
InChIKey, `supported_by` and the class. A dash means no value.

<!-- vale off -->

### (a) `metabolic_network`, `entity` rows

| # | Label | Identifier | InChIKey | Supported by | Class |
|---|---|---|---|---|---|
| 1 | ferrocytochrome c | chebi `CHEBI:16928` | - | kegg, metatlas | Placeholder (protein) |
| 2 | Phospho-beta-adrenergic receptor | metanetx `MNXM93809` | - | kegg, metatlas | Placeholder (protein) |
| 3 | 5-Methyl-5,6,7,8-tetrahydromethanopterin | kegg `C04488` | - | kegg | Specific |
| 4 | Glycogen, structure 4 (glycogenin-1,6-{2[1,4-Glc], [1,4-Glc]}) | metanetx `MNXM11716` | - | metatlas, recon3d | Placeholder (polymer) |
| 5 | UDP-2-acetamido-2,6-dideoxy-beta-L-talose(2-) | inchikey | `INJACODUUNZJCO-GOGQOHCTSA-L` | kegg, rhea | Specific |
| 6 | L-dopachrome | chebi `CHEBI:15772` | - | metatlas | Duplicate (`HMDB0001430`, fallback) |
| 7 | Trans,Cis-Dodeca-2,5-Dienoyl Coenzyme A | human_gem_metabolite `MAM03333` | - | metatlas | Specific |
| 8 | Keratan sulfate I, degradation product 20 | metanetx `MNXM8791` | - | metatlas, recon3d | Placeholder (polymer) |
| 9 | (24S)-24,25-epoxycucurbitadienol | inchikey | `FKLMMFZHGLDEMF-CUGODZHQSA-N` | rhea | Specific |
| 10 | S-methyl-trans-thiourocanate(1-) | inchikey | `ONIAGRIJGZWVQU-NSCUHMNNSA-M` | rhea | Specific |
| 11 | (E)-3-diazocoumarate | inchikey | `HIPZKJMIBCPNSL-DUXPYHPUSA-N` | rhea | Specific |
| 12 | N(alpha)-citryl-N(epsilon)-acetyl-N(epsilon)-hydroxylysine | inchikey | `HPOPJZNNZCIWFX-CUVJYRNJSA-N` | kegg | Specific |
| 13 | 3-hydroxypentadecanoyl-[ACP] | human_gem_metabolite `MAM00794` | - | metatlas | Placeholder (carrier protein) |
| 14 | ADP-2'-imine-ribofurano[1',2':4,5]oxazolidine(2-) residue | chebi `CHEBI:184300` | - | rhea | Placeholder (residue) |
| 15 | (3S,4S)-3-hydroxytetradecane-1,3,4-tricarboxylate | inchikey | `QFOFNCNFUGQWTO-DYVFJYSZSA-K` | rhea | Specific |
| 16 | herquline A(1+) | inchikey | `UFKNDVKQCSBIQE-OCANJJRCSA-O` | rhea | Specific |
| 17 | 2-Amino-4,5-dihydroxy-6-oxo-7-(phosphooxy)heptanoate | inchikey | `OABFYXXSGQYCAM-ZNRXLLLOSA-N` | kegg | Specific |
| 18 | (4S)- 2,3-Dehydroflavan-3,4-diol | inchikey | `QNAVJRDCFCKNDZ-ZDUSSCGKSA-N` | kegg | Specific |
| 19 | aureothin | inchikey | `GQKXCBCSVYJUMI-WACKOAQBSA-N` | kegg, rhea | Specific |
| 20 | [Methylmalonyl-CoA:pyruvate carboxytransferase] | metanetx `MNXM18527` | - | kegg | Placeholder (protein) |
| 21 | (Gal)4 (Glc)1 (GlcNAc)3 (LFuc)3 (Cer)1 | metanetx `MNXM7942` | - | metatlas, recon3d | Placeholder (unnamed ceramide) |
| 22 | 4-guanidinobutanal | inchikey | `VCOFTLCIPLEZKE-UHFFFAOYSA-N` | kegg | Specific |
| 23 | 4-fluorobenzoyl-CoA | inchikey | `ORILDBXADZMSHG-TYHXJLICSA-N` | kegg | Xenobiotic |
| 24 | 3'-end 2'-deoxyribonucleotide 3'-(L-tyrosyl-phosphate)(2-) residue | chebi `CHEBI:157696` | - | rhea | Placeholder (residue) |
| 25 | 2-aminophenylglyoxylate | inchikey | `MQMWPBBDMIYYMI-UHFFFAOYSA-M` | kegg, rhea | Specific |
| 26 | Alanyl-O-citryl-aminoehanol | kegg `C23100` | - | kegg | Specific |
| 27 | trans-2-carboxybenzylidenepyruvate(2-) | inchikey | `APKXMKWCGDBYNV-AATRIKPKSA-L` | kegg | Specific |
| 28 | 2-epi-5-epi-valiolone 7-phosphate(2-) | inchikey | `WBHJALJLPRKCBI-MVIOUDGNSA-L` | kegg, rhea | Specific |
| 29 | 3-[3-(trifluoromethyl)phenyl]-1-methylurea | inchikey | `SPOVEJAPXRPMBI-UHFFFAOYSA-N` | rhea | Xenobiotic (herbicide breakdown) |
| 30 | 5,7-dimethyl-2-methylamino-4-(3-pyridylmethyl)-1,3-benzothiazol-6-yl sulfate(1-) | inchikey | `LPZDTQCWCMPFND-UHFFFAOYSA-M` | rhea | Drug (metabolite of a candidate) |
| 31 | 5-ureidoimidazole-4-carboxylic acid | inchikey | `PQSKWIHIUOGLHX-UHFFFAOYSA-N` | kegg | Specific |
| 32 | 5-Bromo-L-tryptophan | kegg `C23128` | - | kegg | Specific |
| 33 | Lc4Cer | chebi `CHEBI:90800` | - | rhea | Placeholder (unnamed ceramide) |
| 34 | C-terminal Gly-Gly-AMP(1-) group | chebi `CHEBI:90618` | - | rhea | Placeholder (group) |
| 35 | poly(RNA)-3'-uridine ribonucleotide polyanion | chebi `CHEBI:173116` | - | rhea | Placeholder (polymer) |
| 36 | 1-heptadecanoyl-2-stearoyl-sn-glycero-3-phosphate(2-) | inchikey | `VEIHRUNGGHIRIO-PSXMRANNSA-L` | rhea | Specific |
| 37 | Reduced flavodoxin | kegg `C02745` | - | kegg | Placeholder (protein) |
| 38 | (25R)-3-Oxocholest-4-en-26-oyl-CoA | kegg `C23011` | - | kegg | Specific |
| 39 | 3-O-(beta-D-glucopyranosyl)-isomotiol | inchikey | `LVHPZNOIRZPLHW-SKFYXBMVSA-N` | rhea | Specific |
| 40 | (13R,14S)-dihydroxy-(4Z,7Z,9E,11E,16Z,19Z)-docosahexaenoate | inchikey | `ALWYOLKNLLFCAY-GYROAIJESA-M` | rhea | Specific |
| 41 | (4aR,6S,7R)-2-(adenosyl(5')diphospho(5)ribosyl)-6-hydroxy-7-acetamidosulfomoyl-4-carbamoyl-4a,5,6,7-tetrahydro-2H-cyclopenta[c]pyridine-7-carboxylate | inchikey | `HLBMZIBPRKXPSY-LBFKCWOPSA-K` | rhea | Specific |
| 42 | (21S)-21,22-epoxypolypoda-8(26)-13,17-trien-3beta-ol | inchikey | `YOXYJGZFQYXPHZ-AEBUGLRVSA-N` | kegg, rhea | Specific |
| 43 | 3'-end uridine-uridine 2',3'-cyclophosphate(3-) residue | chebi `CHEBI:173080` | - | rhea | Placeholder (residue) |
| 44 | agmatidine 5'-phosphate(1+) residue | chebi `CHEBI:83545` | - | rhea | Placeholder (residue) |
| 45 | 2-Methylthio-N6-L-threonylcarbamoyladenine in tRNA | metanetx `MNXM147735` | - | kegg | Placeholder (residue) |
| 46 | O-(S-decanoylpantetheine-4-phosphoryl)serine(1-) residue | chebi `CHEBI:78468` | - | rhea | Placeholder (residue) |
| 47 | Indoxyl | human_gem_metabolite `MAM03702` | - | metatlas | Duplicate (`HMDB0004094`, fallback) |
| 48 | alpha-(hydroxymethyl)serine | inchikey | `ZRPDXDBGEYHEBJ-UHFFFAOYSA-N` | rhea | Specific |
| 49 | 1'-carboxy-chondrochloren B(1-) | inchikey | `GOTORXHCOMNKGN-HKTLKUQXSA-M` | rhea | Specific |
| 50 | cholesterol-ester-13,16,19-doco | human_gem_metabolite `MAM01460` | - | metatlas | Specific |

### (b) `chebi_role`, `entity` rows

| # | Label | Identifier | InChIKey | Supported by | Class |
|---|---|---|---|---|---|
| 1 | Stemofuran P | inchikey | `LVVYOWMHRJEABT-UHFFFAOYSA-N` | chebi | Specific |
| 2 | chaetomanone | inchikey | `ZVLQJGNLXSDFDN-IYFFJJEWSA-N` | chebi | Specific |
| 3 | O-(hydroxypropionyl)carnitine | chebi `CHEBI:86072` | - | chebi | Unclear (position not given) |
| 4 | 2''-O-(4-sinapoylarabinosyl)isoscoparin | inchikey | `YOFXNYRWIIFBRE-QOLNKJPBSA-N` | chebi | Specific |
| 5 | 7,4'-dihydroxy-8,3'-dimethoxyflavanone | inchikey | `UKCAIYJHVFWSRN-UHFFFAOYSA-N` | chebi | Specific |
| 6 | 2-(2-hydroxyphenyl)ethanol | inchikey | `ABFCOJLLBHXNOU-UHFFFAOYSA-N` | chebi | Specific |
| 7 | 6a,7-didehydro-1,9,10-trihydroxy-2-methoxy-6-methylaporphinium trifluoroacetate | inchikey | `SBKNQDPIZMRYQD-UHFFFAOYSA-N` | chebi | Unclear (salt) |
| 8 | 14-hydroxyplatensic acid methyl ester | inchikey | `FBBBYTFYYNYGMH-MKJSHFIQSA-N` | chebi | Specific |
| 9 | venturamide B | inchikey | `BJJVHBJYOXGRRO-XEQNPHJVSA-N` | chebi | Specific |
| 10 | coptisine | inchikey | `XYHOBCMEDLZUMP-UHFFFAOYSA-N` | chebi | Specific |
| 11 | hordatine A | inchikey | `KVYNYRIOAYQBFK-AIIPJEMGSA-N` | chebi, kegg | Duplicate (`HMDB0030461`) |
| 12 | chetracin A | inchikey | `PFPUKJNFELCZAQ-AAZDRVSSSA-N` | chebi | Specific |
| 13 | hydroxyversicolorone(1-) | inchikey | `JGXCLZAVTLWCBF-LXGOIASLSA-M` | chebi, kegg, rhea | Specific |
| 14 | Turnagainolide B | inchikey | `SMWGFJCKWZCJNQ-QWGWBJRGSA-N` | chebi | Specific |
| 15 | Euscaphic acid C | inchikey | `FJBWUBWZSSBDEN-XVNGACRRSA-N` | chebi | Specific |
| 16 | Maoesin E, (rel)- | inchikey | `YLAQTEHVFWEICH-MQAPUGLSSA-N` | chebi | Specific |
| 17 | (2E,6R)-8-hydroxy-2,6-dimethyl-2-octenoic acid | inchikey | `CDAJACJYYZHNJA-CBVZUMJXSA-N` | chebi | Specific |
| 18 | eucalyptal C | inchikey | `LLCSUGKJCVVJBS-LQEBHLGXSA-N` | chebi | Specific |
| 19 | Athrolide B | inchikey | `BHEXSLGHYNIAGL-WDCRBBJTSA-N` | chebi | Specific |
| 20 | 22alpha,23alpha-epoxy-solanida-1,4-dien-3-one, (rel)- | inchikey | `UGAJCJIZYPKGTJ-WUOUAPBYSA-N` | chebi | Specific |
| 21 | Conoideocrellide B | inchikey | `XUGFNOLBANLPEQ-JNVOQPTMSA-N` | chebi | Specific |
| 22 | Brialalepolide C | inchikey | `CEVYVRCWJDJLEW-LXIKAAAESA-N` | chebi | Specific |
| 23 | chaetoglocin C | inchikey | `BKRGUIJDEOVXCR-ZCFIWIBFSA-N` | chebi | Specific |
| 24 | stachsterol | inchikey | `GEYCTESILRJHOA-RBNQBTAKSA-N` | chebi | Specific |
| 25 | 3,4,5-tribromo-2-(2,4-dibromophenoxy)phenol | inchikey | `LNZHBUPVHNJGJG-UHFFFAOYSA-N` | chebi | Specific |
| 26 | seco-plakortolide O | inchikey | `BBUUFCRHKRIKHW-XRXSLGMHSA-N` | chebi | Specific |
| 27 | Pseudopalmatine | inchikey | `CLFBXKHKECKSQM-UHFFFAOYSA-N` | chebi | Specific |
| 28 | 15-norlabda-8(20),12E-diene-14-carboxaldehyde-19-oic acid | inchikey | `NBGCWDAYASHSEK-WCBDFPTQSA-N` | chebi | Specific |
| 29 | kaempferol 7-O-glucuronide | inchikey | `HKMMHEQNCQNUNN-JENRNSKYSA-N` | chebi | Specific |
| 30 | 6-O-Nicotinoyl-7-O-acetylscutebarbatine G, (rel) | inchikey | `MCZQKBCMIDJMDT-RJGRQJGLSA-N` | chebi | Specific |
| 31 | vaccihein A | inchikey | `HTYXVCATPMWMSN-UHFFFAOYSA-N` | chebi | Specific |
| 32 | eurysterol B sulfonic acid | inchikey | `QPVBWEUWNJWRLD-XPVUCKRVSA-N` | chebi | Specific |
| 33 | phellamurin | inchikey | `GRDZTDZJQRPNCN-YIANMRPHSA-N` | chebi | Specific |
| 34 | trichagmalin D | inchikey | `MYVFKKIXMOXVPS-KLCWGGGHSA-N` | chebi | Specific |
| 35 | PC O-44:2 | goslin `species:PC O-44:2` | - | chebi | Class (lipid species) |
| 36 | Scutebarbatine J(rel) | inchikey | `QKISVSHUNJJKNY-CZIZJKDCSA-N` | chebi | Specific |
| 37 | 3''O-acetylpolygalacin D | inchikey | `PVYHFAATIMELDB-WGMQSPGRSA-N` | chebi | Specific |
| 38 | FA 20:2;oxo | goslin `structure_defined:FA 20:2;oxo` | - | chebi | Class (lipid species) |
| 39 | Nocazine A | inchikey | `PIZYGUXDBZYKHO-AXPXABNXSA-N` | chebi | Specific |
| 40 | Xylopimarane | inchikey | `GQECNZJTMLRIAM-DWRPHZCTSA-N` | chebi | Specific |
| 41 | mycocyclosin | inchikey | `PYDUENRHWXNYLP-KBPBESRZSA-N` | chebi, kegg, rhea | Specific |
| 42 | isoscoparin | inchikey | `KOMUHHCFAXYRPO-DGHBBABESA-N` | chebi, kegg | Specific |
| 43 | methyl kulonate | inchikey | `ZLOYHDGVJFCYJK-RLIKWOGQSA-N` | chebi | Specific |
| 44 | kabiramide K | inchikey | `BURLCNUJBHPOOO-RUWUMMDESA-N` | chebi | Specific |
| 45 | hydrangenol 4'-O-beta-D-apiofuranosyl-(1->6)-beta-D-glucopyranoside | inchikey | `NCBLEDBCMRNNAV-NFZGRANLSA-N` | chebi | Specific |
| 46 | Wewakazole | inchikey | `LRFQMOZAQWOOPZ-IBFZZIGPSA-N` | chebi | Specific |
| 47 | rel-16-Angeloyloxy-13alpha-isobutanoyloxy-4beta,9alpha,20-trihydroxytiglia-1,5-diene-3,7-dione | inchikey | `YYLJXVKJJBZWBK-VFLXVHJUSA-N` | chebi | Specific |
| 48 | (cyclo(N-methyl-L-Phe-L-Val-D-Ile-L-Leu-L-Pro)) | inchikey | `NSBVGLQWOAQURC-ZFWJRURNSA-N` | chebi | Specific |
| 49 | maoecrystal Z | inchikey | `WAOUILKFOLDIAX-NJOKITKISA-N` | chebi | Specific |
| 50 | chaetocin C | inchikey | `VSNZIVJJANWYJU-IRQSJMDVSA-N` | chebi | Specific |

### (c) `hmdb`, `structure_key` rows

Every row of this group has `supported_by` empty, except where the table
names resources. The last column but one names the HMDB molecule of
`via_entity_key`.

| # | Label | InChIKey | Supported by | Via | Class |
|---|---|---|---|---|---|
| 1 | CL 18:2(9Z,12Z)/16:1(9Z)/18:2(9Z,12Z)/18:1(11Z) | `CSJTXGPVMRYBBU-UZNPABKHSA-L` | - | `HMDB0058742`, same name | Duplicate |
| 2 | PE(20:2(11Z,14Z)/P-18:1(9Z)) | `IBIDTCQKTHZZPM-HUESYALOSA-N` | - | `HMDB0009315`, same lipid | Duplicate |
| 3 | 1-(15Z-tetracosenoyl)-2-(6Z,9Z,12Z-octadecatrienoyl)-sn-glycero-3-phospho-L-serine | `PXDXCPYWDQRHSF-NNHJAWLDSA-M` | - | `HMDB0112919`, PS 24:1/18:3 | Duplicate |
| 4 | CID:69446282 | `KUFRQPKVAWMTJO-SBNHGUPZSA-N` | - | `HMDB0244766`, 17-Dimethylaminogeldanamycin | Drug |
| 5 | 3alpha,17alpha-dihydroxy-5beta-pregnan-20-one | `LKQDFQLSEHWIRK-UKBVDAKRSA-N` | - | `HMDB0244771`, no stereo | Duplicate |
| 6 | etrasimod(1-) | `MVGWUTBTXDYMND-QGZVFWFLSA-M` | - | `HMDB0252121`, Etrasimod | Drug |
| 7 | PTEROCARPANQUINONE | `BUEHOQICFGQBNY-ACJLOTCBSA-N` | - | `HMDB0254184`, Pterocarpanquinone | Drug |
| 8 | ALTROPANE I-123 | `GTQLIPQFXVKRKJ-HYRAKNMCSA-N` | - | `HMDB0248262`, Altropane | Drug (radiotracer) |
| 9 | (2S)-2-[[[4-[[(2R)-2-amino-3-mercaptopropyl]amino]-2-(1-naphthalenyl)phenyl]-oxomethyl]amino]-4-methylpentanoic acid methyl ester | `XVWPFYDMUFBHBF-CLOONOSVSA-N` | - | `HMDB0252710`, no stereo | Drug |
| 10 | PMX-205 | `VATFHFJULBPYLM-ILOBPARPSA-N` | - | `HMDB0256653` | Drug |
| 11 | Acid thyroliberin | `ITYONPBTNRIEBA-SRVKXCTJSA-N` | - | `HMDB0256375`, Pglu-his-pro | Duplicate |
| 12 | 1-(9Z,12Z,15Z-octadecatrienoyl)-2-(11Z,14Z-eicosadienoyl)-3-(6Z,9Z,12Z,15Z-octadecatetraenoyl)-sn-glycerol | `YPKAPKOHLOHKML-UEZATWCYSA-N` | - | `HMDB0055268`, same TG | Duplicate |
| 13 | 1-(4Z,7Z,10Z,13Z,16Z,19Z-docosahexaenoyl)-2-(9Z,12Z-octadecadienoyl)-sn-glycero-3-phospho-(1'-sn-glycerol-3'-phosphate) | `WIACNPDTUUWXQR-CFWJWXGNSA-K` | - | `HMDB0116475`, charge form | Duplicate |
| 14 | 3-hydroxyoctadecanoyl-CoA | `WZMAIEGYXCOYSH-AJEXUEQGSA-J` | - | `HMDB0012715`, same name | Duplicate |
| 15 | 1-(11E-octadecenoyl)-2-(9Z-tetradecenoyl)-3-(7Z,10Z,13Z,16Z-docosatetraenoyl)-sn-glycerol | `QIBYPZKSFQCMRD-ULNYJJSBSA-N` | - | `HMDB0049231`, same TG | Duplicate |
| 16 | sn-glycerol 1-phosphate | `AWUCVROLDVIAJX-VKHMYHEASA-N` | kegg | `HMDB0000126`, Glycerol 3-phosphate | Specific (enantiomer) |
| 17 | 5,6-dihydroxyindole-2-carboxylate | `YFTGOBNOJKXZJC-UHFFFAOYSA-M` | recon3d, rhea | `HMDB0001253`, charge form | Duplicate |
| 18 | 1-(11Z-octadecenoyl)-2-(7Z,10Z,13Z,16Z-docosatetraenoyl)-sn-glycero-3-phospho-(1'-sn-glycerol) | `QFTIPVNWXHDAQV-ADJDFPOESA-M` | - | `HMDB0010626`, charge form | Duplicate |
| 19 | 19alpha-hydroxyasiatic acid-28-O-beta-D-glucopyrannoside | `WKKBYJLXSKPKSC-JVJIQXRHSA-N` | - | `HMDB0035691`, Quercilicoside A, no stereo | Duplicate |
| 20 | CID:90660380 | `AYUNIORJHRXIBJ-ZIKJQUQBSA-N` | - | `HMDB0244762`, 17-AAG | Drug |
| 21 | 1-(15Z-tetracosenoyl)-2-(4Z,7Z,10Z,13Z,16Z,19Z-docosahexaenoyl)-3-(7Z,10Z,13Z,16Z,19Z-docosapentaenoyl)-sn-glycerol | `VHWRAFCCARCRSA-VGMDJCHVSA-N` | - | `HMDB0052378`, same TG | Duplicate |
| 22 | CID:6437840 | `JOSHOGBFUULHNI-LYKFAKFTSA-N` | - | `HMDB0253266`, Hydroxyethyl retinamide | Drug |
| 23 | PC(22:0/18:1(11Z)) | `ZXUGYTRXWUMLQE-YACUFSJGSA-N` | - | `HMDB0008531`, same lipid | Duplicate |
| 24 | (2S)-2-hydroxytetradecanoic acid | `JYZJYKOZGGEXSX-ZDUSSCGKSA-N` | - | `HMDB0002261`, no stereo | Duplicate |
| 25 | (-)-neoisodihydrocarveol | `KRCZYMFUWVJCLI-UTLUCORTSA-N` | kegg, rhea | `HMDB0302494`, another stereoisomer | Specific |
| 26 | TG 16:0/20:1(11Z)/20:1(11Z) | `OJBGBKUEIBXOEG-YYEKOYDESA-N` | - | `HMDB0005386`, same name | Duplicate |
| 27 | 1-(9Z-tetradecenoyl)-2-(7Z,10Z,13Z,16Z,19Z-docosapentaenoyl)-3-(11Z-eicosenoyl)-sn-glycerol | `ZNNDXFQWHSFNNR-GDCVHKBMSA-N` | - | `HMDB0048349`, same TG | Duplicate |
| 28 | Myricetin 3-neohesperidoside | `LUITYVJNCILVDI-NVMXCJJVSA-N` | - | `HMDB0038287`, same name | Duplicate |
| 29 | 1,4-dideoxy-1,4-[(S)-[(2S,3S)-2, 4-dihydroxy-3-(sulfooxy)butyl]episulfoniumylidene]-D-arabinitol | `SOWRVDSZMRPKRG-XESHOGEWSA-N` | - | `HMDB0257455`, Salacinol | Duplicate |
| 30 | CHEBI:156835 | `AAKHDGZVKNZPFG-XERMAKKKSA-N` | - | `HMDB0248833`, B Lewis b | Duplicate |
| 31 | CID:2626 | `SNBUBQHDYVFSQF-UHFFFAOYSA-N` | - | `HMDB0014419`, Cefmetazole | Drug |
| 32 | HC 5:1(3);oxo | `LABTWGUMFABVFG-UHFFFAOYSA-N` | - | `HMDB0001184`, (3E)-pent-3-en-2-one | Duplicate |
| 33 | 1-(7Z,10Z,13Z,16Z-docosatetraenoyl)-2-(4Z,7Z,10Z,13Z,16Z,19Z-docosahexaenoyl)-3-(9Z,12Z,15Z-octadecatrienoyl)-sn-glycerol | `MIQPWXAIHLQOQB-JFTKSENDSA-N` | - | `HMDB0054894`, same TG | Duplicate |
| 34 | PI(22:5(4Z,7Z,10Z,13Z,16Z)/16:0) | `RQOVEBMFWCAUPR-OXRHWTJRSA-N` | - | `HMDB0009916`, same lipid | Duplicate |
| 35 | 1-(5Z,8Z,11Z,14Z-eicosatetraenoyl)-2-(9Z-tetradecenoyl)-sn-glycero-3-phospho-L-serine | `GHAIZGVIVLPOCT-PFJWXCNISA-M` | - | `HMDB0012430`, charge form | Duplicate |
| 36 | CL 18:0/22:5(4Z,7Z,10Z,13Z,16Z)/18:1(9Z)/18:0 | `IASIQNHCHIMJSQ-KOIJUYCZSA-L` | - | `HMDB0057288`, same name | Duplicate |
| 37 | TG 20:0/18:3(9Z,12Z,15Z)/18:2(9Z,12Z) | `ZQURXMDEQOSXBZ-ISARBORQSA-N` | - | `HMDB0046086`, same TG | Duplicate |
| 38 | Armillyl everninate | `KIBURFBDQVFUHO-ZXSIQIRYSA-N` | - | `HMDB0035849`, same name | Duplicate |
| 39 | 1-eicosanoyl-2-(4Z,7Z,10Z,13Z,16Z,19Z-docosahexaenoyl)-3-(6Z,9Z,12Z-octadecatrienoyl)-sn-glycerol | `JASRGXRLHMKFMV-GQLITWBWSA-N` | - | `HMDB0046217`, same TG | Duplicate |
| 40 | ST 28:1;O;Hex | `FWNZEKQVBDXWKA-NYENXMMQSA-N` | rhea | `HMDB0303802`, Campesterol glucoside | Duplicate |
| 41 | 1-(8Z,11Z,14Z-eicosatrienoyl)-2-eicosanoyl-3-(11E-octadecenoyl)-sn-glycerol | `HYPRJUDAMZSDSL-FLWPIAGDSA-N` | - | `HMDB0049165`, same TG | Duplicate |
| 42 | PIP(22:2(13Z,16Z)/16:0) | `IRUVBUYWSSKKPM-CSDMGBHBSA-N` | - | `HMDB0010023`, same lipid | Duplicate |
| 43 | 1-pentadecanoyl-2-(13Z-docosenoyl)-3-(8Z,11Z,14Z-eicosatrienoyl)-sn-glycerol | `JPQBIVGPVQLOCO-JOWIJVFESA-N` | - | `HMDB0043354`, same TG | Duplicate |
| 44 | L-tryptophanium | `QIVBCDIJIAJPQS-VIFPVBQESA-O` | - | `HMDB0030396`, (±)-Tryptophan | Duplicate |
| 45 | 1-Lyso-2-arachidonoyl-phosphatidate | `BDCFJMBXZCIVRH-CGRWFSSPSA-N` | - | `HMDB0012496`, same name | Duplicate |
| 46 | 1-hexadecanoyl-2-(13Z,16Z-docosadienoyl)-3-(8Z,11Z,14Z-eicosatrienoyl)-sn-glycerol | `ZNIJKWJHVJGSJB-IVKAUKIJSA-N` | - | `HMDB0044405`, same TG | Duplicate |
| 47 | 1-pentadecanoyl-2-(5Z,8Z,11Z,14Z,17Z-eicosapentaenoyl)-3-(13Z-docosenoyl)-sn-glycerol | `BBYNELOLAOBLQM-KMFBYXSPSA-N` | - | `HMDB0043726`, same TG | Duplicate |
| 48 | (S)-methcathinone(1+) | `LPLLVINFLBSFRP-QMMMGPOBSA-O` | rhea | `HMDB0041927`, Methcathinone | Drug |
| 49 | (-)-Menthyl acetate | `XHXUANMFYXWVNG-ADEWGFFLSA-N` | kegg | `HMDB0303918`, DL-Menthyl acetate | Duplicate |
| 50 | 1-pentadecanoyl-2-(8Z,11Z,14Z,17Z-eicosatetraenoyl)-3-(13Z,16Z-docosadienoyl)-sn-glycerol | `OTUYWLCDLGLHJW-AERSBBFDSA-N` | - | `HMDB0043704`, same TG | Duplicate |

<!-- vale on -->

A "same TG" row names the same fatty acids in another order of the sn
positions, or with another stereo descriptor. A swap of sn-1 and sn-3 changes
only the stereo block of the InChIKey.
