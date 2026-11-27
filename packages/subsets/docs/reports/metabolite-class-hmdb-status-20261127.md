# HMDB status and the metabolite-class gate on release 2026.11.27

2026-11-27

This report gives the figures for step 9 of stage two of the plan
`saezverse/projects/omnipath/plans/2026-11-metabolite-classification-logic.md`:
does the bucket `hmdb` narrow by status? It also measures the gates of step 8
again on the new release. The report takes no decision. The section
"Options for the decision" states the choices and their costs.

**Result: option (b), detected and quantified only, removes 223,742 of 272,373
rows from `hmdb`. It removes 322,126 of 1,349,547 relations from the gate
`[hmdb]`.** A narrowing that removes only the predicted entries, option (c),
removes 95,977 rows but only 339 relations. Only 27 of the 95,354 predicted
HMDB molecules have a MetaLinksDB-shaped relation. In the audit sample of
expected molecules, 17 of 50 are drugs and 16 are food or plant compounds.
Four are endogenous human metabolites, and acetyl-CoA is one of them.

## What this report measures

| Item | Value |
|---|---|
| Release | `2026.11.27`, read from `/workspace/instances/metsig1/data-2026.11.19`. The comparison reads `2026.11.19` from the same root. |
| Resource versions | `chebi` and `hmdb` at `2026.11.27`, both built from `pypath` at `0f4978696`. The other 167 resources are the same as in `2026.11.19`. |
| Code | `pypath` at `0f4978696`. The subset code last changed in `70e6006e4`. |
| Build | `parquet_build.build_tables(con, build_id = 'audit', rules = read_rules(path))` |
| Scripts | In `/tmp/claude-1029/-workspace-instances-metsig1/97c1ab55-a53d-4d90-8da3-8c1f34a8d599/scratchpad/stage2-hmdb/`, listed at the end |

**The rule files.** "Stage one" is the packaged rule file, with the HMDB
signal only. "Full" is the packaged file with the two stage-two rows
uncommented. All signals are on in "full".

**The narrowings.** The build has no status condition. The script `build.py`
replaces the text of `claim_resource.sql` for one build at a time, in memory.
The replacement adds one condition: the HMDB record carries an
`hmdb:status` value in the admitted set. The repository stays unchanged. The
rest of the cascade, the groups included, runs as the packaged code runs it.

| Option | Admitted statuses |
|---|---|
| (a) | all four, the build of today |
| (b) | `quantified`, `detected` |
| (c) | `quantified`, `detected`, `expected` |
| (d) | every status except `predicted` |

**The relation set and the gate** are the ones of the report
`metabolite-class-stage-two-measurements-20261127.md`. A relation is in the set
when its predicate is `interacts_with`, `affects` or `associated_with`. One end
must be a chemical. The other end must have a type in `_PARTNER_TYPES` of
`network_views/_definitions.py`. The type of an entity is the one that
`release_entity` gives it. A gate with buckets B admits a relation when one end
has a row in `metabolite_class` with a bucket in B.

The set holds 6,261,809 rows on both releases. The same scripts on
`2026.11.19` give the bucket sizes and the five gate counts of the earlier
report, to the row. So the method of that report and this method agree.

## The status in the release

The `hmdb` resource of `2026.11.27` holds 217,916 chemicals and 3,350
`ontology_class` entities. Every chemical carries `hmdb:status`. No
`ontology_class` entity carries it.

| Status | HMDB chemicals |
|---|---|
| `quantified` | 3,385 |
| `detected` | 20,921 |
| `expected` | 98,256 |
| `predicted` | 95,354 |
| total | 217,916 |

One chemical, (-)-parthenolide (`KTEXNACQROZXEV-UHFFFAOYSA-N`), has two
statuses, `predicted` and `detected`. The tables count it as `detected`, its
best status. So `predicted` counts 95,354 molecules and 95,355 annotation rows.

## Status per bucket member

The bucket `hmdb` is the same with the stage-one rule file and with the full
one: 272,373 rows. The other signals cannot claim a molecule that HMDB claims.

**`entity` rows.** All 217,916 HMDB chemicals are `entity` rows of `hmdb`. The
status of a row is the status of its own record.

| Status | `entity` rows | With a relation in the set | Without an InChIKey |
|---|---|---|---|
| `quantified` | 3,385 | 1,373 | 0 |
| `detected` | 20,921 | 7,087 | 4 |
| `expected` | 98,256 | 4,930 | 0 |
| `predicted` | 95,354 | 27 | 17 |
| total | 217,916 | 13,417 | 21 |

**`structure_key` rows.** A variant has no HMDB record of its own. Its status
here is the best status among the HMDB members of its connectivity group.

| Best status in the group | `structure_key` rows | Groups | With a relation in the set | Claimed by another signal |
|---|---|---|---|---|
| `quantified` | 6,646 | 2,719 | 1,807 | 1,536 |
| `detected` | 16,371 | 9,382 | 7,498 | 1,739 |
| `expected` | 30,737 | 22,054 | 4,242 | 3,230 |
| `predicted` | 703 | 484 | 36 | 147 |
| total | 54,457 | 34,639 | 13,583 | 6,652 |

"Claimed by another signal" counts the variants with a non-empty
`supported_by` in the full build.

The relations follow the same split. Of the 1,349,547 relations of `[hmdb]`,
the best end of 472,296 is quantified. The best end of 479,387 is detected,
and of 396,987 expected. 877 relations reach `hmdb` only through a
predicted molecule. Here an HMDB molecule takes the status of its own record,
and a variant the best status of its group.

## What each narrowing does to the buckets

All builds use the full rule file. A group stays in `hmdb` when at least one of
its HMDB members has an admitted status.

| Bucket | (a) | (b) | (c) | (d) |
|---|---|---|---|---|
| `hmdb`, `entity` rows | 217,916 | 24,306 | 122,562 | 122,562 |
| `hmdb`, `structure_key` rows | 54,457 | 24,325 | 53,834 | 53,834 |
| `hmdb`, total | 272,373 | 48,631 | 176,396 | 176,396 |
| `chebi_role` | 11,912 | 18,970 | 11,954 | 11,954 |
| `metabolic_network` | 15,979 | 20,574 | 16,302 | 16,302 |
| all rows | 300,264 | 88,175 | 204,652 | 204,652 |

**Option (d) is option (c).** Every HMDB chemical carries one of the four
statuses. The two builds give the same rows, with no difference in any column.
The rest of this report gives (c) only.

**Option (b).** The narrowing excludes 193,610 HMDB records. 1,308 of them stay
in `hmdb` as `structure_key` rows, because their group holds a detected or
quantified molecule. 1,302 of these are expected records, for example
acetyl-CoA. The `structure_key` rows of `hmdb` include these 1,308.

**Option (c).** The narrowing excludes the 95,354 predicted records. 80 of them
stay in `hmdb` through their group.

**Where the excluded rows go.** The table counts the rows of `hmdb` under (a)
that leave `hmdb`. "`entity` row under (a)" is an HMDB record. The brackets
give how many of the rows are `entity` rows after the narrowing.

| Destination | (b), `entity` row under (a) | (b), `structure_key` row under (a) | (c), `entity` row under (a) | (c), `structure_key` row under (a) |
|---|---|---|---|---|
| `chebi_role` | 3,060 (1,701 `entity`) | 3,998 (1,358 `entity`) | 17 (3 `entity`) | 25 (15 `entity`) |
| `metabolic_network` | 1,914 (914 `entity`) | 2,681 (1,338 `entity`) | 109 (2 `entity`) | 214 (127 `entity`) |
| no row | 187,328 | 24,761 | 95,148 | 464 |
| total | 192,302 | 31,440 | 95,274 | 703 |

Under (b), 97 per cent of the excluded HMDB records get no row. The predicted
records leave almost whole under both options: 95,168 of 95,354 get no row
under (b).

**The molecules with relations.** Under (b), 4,254 HMDB records and 4,278
variants with a relation leave `hmdb`. Of the 4,254 records, 1,267 go to
`chebi_role`, 642 to `metabolic_network` and 2,345 get no row. Of the same
4,254 records, 1,244 have a DrugBank cross-reference in HMDB. 874 have a
DrugCentral record, 2,136 a FooDB record and 1,885 a ChEMBL record. Under (c),
15 HMDB records and 36 variants with a relation leave `hmdb`.

The lost molecules with the most relations under (b):

| Molecule | Relations | Under (b) |
|---|---|---|
| Benzo[a]pyrene | 16,347 | `chebi_role` |
| Hydrogen Ion | 10,875 | `metabolic_network` |
| Copper sulfate | 6,309 | no row |
| Silica | 4,766 | no row |
| Ethyl carbamate | 4,397 | `chebi_role` |
| Edetic Acid | 4,331 | `metabolic_network` |
| Indomethacin | 3,873 | `chebi_role` |
| Vorinostat | 3,680 | no row |
| Benzo[a]pyrene-7,8-dihydrodiol-9,10-oxide | 3,466 | `metabolic_network` |
| Temozolomide | 3,253 | no row |
| Acrylamide | 3,161 | `metabolic_network` |
| H2 | 3,148 | `chebi_role` |

Under (c), the lost molecules with the most relations are acylcarnitine and
acyl-CoA variants. The largest carries 57 relations.

## The gates

**Today's rules on both releases.** Option (a), all signals on.

| Gate | 2026.11.19 | 2026.11.27 | Change |
|---|---|---|---|
| `[hmdb]` | 1,349,547 | 1,349,547 | 0 |
| `[hmdb, chebi_role]` | 1,354,756 | 1,381,437 | +26,681 |
| `[hmdb, chebi_role, metabolic_network]` | 1,408,693 | 1,419,874 | +11,181 |
| `[chebi_role]` | 5,212 | 31,929 | +26,717 |
| `[metabolic_network]` | 53,991 | 38,464 | -15,527 |

The new ChEBI build gives the role walk the full role tree. Signal 2 claims
16,253 molecules, against 5,317 on `2026.11.19`. 2,527 molecules move from
`metabolic_network` to `chebi_role`. The bucket `hmdb` gains 347
`structure_key` rows. All 347 are ChEBI chemicals that the old ChEBI build did
not hold, and none has a relation in the set. So `[hmdb]` stays the same.
The report of the parallel ChEBI measurement gives the role walk in detail.

**The narrowings on 2026.11.27.**

| Gate | (a) | (b) | (c) |
|---|---|---|---|
| `[hmdb]` | 1,349,547 | 1,027,421 | 1,349,208 |
| `[hmdb, chebi_role]` | 1,381,437 | 1,158,367 | 1,381,123 |
| `[hmdb, chebi_role, metabolic_network]` | 1,419,874 | 1,252,231 | 1,419,844 |
| `[chebi_role]` | 31,929 | 131,064 | 31,955 |
| `[metabolic_network]` | 38,464 | 93,898 | 38,748 |

The stage-one rule file gives the same `[hmdb]` counts for (a), (b) and (c).

Under (b), `[hmdb]` loses 322,126 relations, 23.9 per cent. The gate with all
three buckets loses 167,643, 11.8 per cent. So signals 2 and 3 take back
154,483 of the lost relations. Under (c), `[hmdb]` loses 339 relations and the
gate with all three buckets loses 30.

**The resources that lose most.** The resource is the one that states the
relation. The last column counts the relations that the gate with all three
buckets also loses under (b).

| Resource | `[hmdb]` under (a) | Lost under (b) | Share | Also lost from all three buckets |
|---|---|---|---|---|
| `pathwaycommons` | 482,118 | 124,280 | 25.8% | 57,885 |
| `stitch` | 465,521 | 103,653 | 22.3% | 54,657 |
| `chembl` | 101,760 | 29,002 | 28.5% | 22,502 |
| `brenda` | 144,888 | 24,665 | 17.0% | 11,329 |
| `bindingdb` | 36,301 | 9,673 | 26.6% | 7,878 |
| `rhea` | 25,485 | 9,167 | 36.0% | 0 |
| `drugcentral` | 17,470 | 8,836 | 50.6% | 8,029 |
| `tcdb` | 11,320 | 2,144 | 18.9% | 235 |
| `cancerdrugsdb` | 4,622 | 2,057 | 44.5% | 1,621 |
| `guidetopharma` | 7,931 | 2,037 | 25.7% | 1,698 |
| `metatlas` | 16,001 | 1,802 | 11.3% | 0 |
| `dgidb` | 6,141 | 1,349 | 22.0% | 733 |
| 18 others | | 3,461 | | |

The drug resources lose the largest share: DrugCentral loses half its
relations, CancerDrugsDB 45 per cent. Rhea, metatlas and Recon3D lose
relations from `[hmdb]` but none from the gate with all three buckets. Their
molecules go to `metabolic_network`.

Under (c), the loss is small everywhere: metatlas 196, Recon3D 65, Pathway
Commons 22, STITCH 18, ChEMBL 17, BindingDB 12, BRENDA 4, iMM1415 4 and
SIGNOR 1.

## The audit of predicted and expected molecules

**The samples.** Each sample takes HMDB `entity` rows of the full build with
one status only and at least one relation in the set. The order is
`hash(entity_key)` and the limit is 50.

- Predicted: 27 molecules are eligible, so the sample holds all 27.
- Expected: 4,930 molecules are eligible. The sample holds 50.

**How the audit classifies a molecule.** The script `audit.py` gives each
molecule its label, own HMDB accession, InChIKey and relations by resource. It
also gives the resources that hold it, its DrugBank and ChEMBL
cross-references and the HMDB description. `groupmates.py` adds the other HMDB members of its group with
their status. The class comes from a reading of each row. `classify.py` holds
the classes.

| Class | Rule |
|---|---|
| Endogenous human metabolite | A defined molecule of human metabolism: an intermediate of a pathway, a conjugate, or a second record or charge form of such a molecule |
| Food or plant compound | A compound of a food, a plant or a spice, or a flavor or food additive |
| Drug or drug metabolite | A drug, a drug candidate, or a product of the metabolism of a drug |
| Microbial metabolite | A product of a bacterium or a fungus that humans do not make |
| Lipid species (generic) | A member of a series that HMDB lists by acyl chain, without evidence for the member itself, or a lipid name that does not fix the structure |
| Other or unclear | An element, an inorganic salt, an environmental xenobiotic, or a compound of unclear origin |

A row takes the first class of the table that applies. The classes come from
one reading. A second reader can move borderline rows, for example the N-acyl
amino acids between "endogenous" and "generic".

### The counts

The relations column counts the relations of the set at the sample molecules.
The last column counts the molecules that stay in `hmdb` through their group:
under (c) for the predicted sample, under (b) for the expected sample.

**Predicted, 27 molecules:**

| Class | Molecules | Relations | Stay in `hmdb` under (c) |
|---|---|---|---|
| Endogenous human metabolite | 12 | 532 | 11 |
| Food or plant compound | 0 | 0 | 0 |
| Drug or drug metabolite | 0 | 0 | 0 |
| Microbial metabolite | 1 | 12 | 1 |
| Lipid species (generic) | 14 | 32 | 0 |
| Other or unclear | 0 | 0 | 0 |

**Expected, 50 molecules:**

| Class | Molecules | Relations | Stay in `hmdb` under (b) |
|---|---|---|---|
| Endogenous human metabolite | 4 | 883 | 2 |
| Food or plant compound | 16 | 230 | 1 |
| Drug or drug metabolite | 17 | 759 | 0 |
| Microbial metabolite | 3 | 20 | 0 |
| Lipid species (generic) | 1 | 35 | 0 |
| Other or unclear | 9 | 226 | 1 |

**The predicted sample.** The 12 endogenous molecules are acyl-CoAs of amino
acid and fatty acid breakdown, and bile acid conjugates. Five are second HMDB
records of quantified molecules, for example N-arachidonoylglycine and
taurochenodeoxycholate. The other seven are acyl-CoAs, for example
acetoacetyl-CoA, which HMDB also lists in an expected record. Under (c), 11
of the 12 stay in `hmdb` through their group. The 14 generic lipids are N-acyl
amino acids, acylcarnitines and one oxidized diacylglycerol. Most have one
relation, from BindingDB.

**The expected sample.** The 17 drugs carry 759 relations. They include
valsartan, valdecoxib, rimonabant and methylphenidate. The 4 endogenous
molecules carry 883 relations. Acetyl-CoA alone carries 878 of them. Its
record `HMDB0001206` is expected. The stereo-free record `HMDB0247926` is
detected, so under (b) acetyl-CoA stays in `hmdb` through the group. Two
endogenous molecules leave `hmdb` under (b): (2E,4Z)-decadienoyl-CoA and
1-octadecylglycerone 3-phosphate. Both go to `metabolic_network`. Of the 46
molecules that leave `hmdb` under (b), 16 go to `chebi_role` or
`metabolic_network`, and 30 get no row.

## What looks wrong

1. **HMDB marks central metabolites as expected.** Acetyl-CoA
   (`HMDB0001206`), the proton (`Hydrogen Ion`) and H2 are expected. The
   proton carries 10,875 relations and moves to `metabolic_network` under (b).
   A status narrowing measures what HMDB has detected, not what is a
   metabolite.
2. **HMDB holds second records of one molecule with another status.**
   N-arachidonoylglycine is quantified in `HMDB0005096` and predicted in
   `HMDB0241918`. Taurochenodeoxycholate, the glycine conjugates of bile acids
   and N-oleoylglycine show the same pattern. The groups repair this when the
   records share the connectivity.
3. **The `hmdb` column can name another record.** It takes the lowest HMDB
   accession of the entity, and the resolver adds accessions of other records.
   For 663 of 217,916 HMDB molecules it names another record. For example,
   N(tele)-methylhistaminium shows `HMDB0000898`, the record of
   1-methylhistamine, but its own record is `HMDB0062574`. A join of the status
   on this column gives the wrong status for these molecules.
4. **An HMDB description mixes two compounds.** The description of N-oleoyl
   glycine (`HMDB0241916`) gives it the names MK-0518 and raltegravir. The
   HMDB description of 5-fluoro-1H-indole-2-carboxylic acid says that
   Lactobacillus makes it. This does not change the subset.
5. **The count of 95,355 predicted entries is a count of annotation rows.**
   The molecules count 95,354, because parthenolide carries two statuses.

## Options for the decision

These are the options for a human decision. The figures inform the decision
but do not take it. All figures are for release `2026.11.27` with all signals
on.

| Option | `hmdb` rows | `[hmdb]` relations | All three buckets, relations |
|---|---|---|---|
| (a) all statuses | 272,373 | 1,349,547 | 1,419,874 |
| (b) quantified and detected | 48,631 | 1,027,421 (-23.9%) | 1,252,231 (-11.8%) |
| (c) all but predicted, the same as (d) | 176,396 | 1,349,208 (-0.03%) | 1,419,844 (-0.002%) |

1. **Keep all statuses, option (a).** The gate keeps every HMDB relation. The
   bucket keeps 95,354 predicted molecules, and 27 of them have a relation.
   The bucket also keeps the drugs and food compounds of the expected records.
   The plan lists this breadth as a known limit.
2. **Remove the predicted entries, option (c).** The bucket loses 95,977 rows,
   35 per cent. The gate `[hmdb]` loses 339 relations. In the sample, the lost
   predicted molecules are mostly generic lipids with one relation each. Most
   endogenous predicted molecules stay in `hmdb` through a group. The change
   makes the bucket smaller and more exact, but barely changes a gated network.
3. **Keep detected and quantified only, option (b).** The bucket loses 82 per
   cent of its rows. The gate `[hmdb]` loses 322,126 relations. With all three
   buckets, the loss is 167,643 relations. The lost molecules include many
   drugs: 17 of 50 in the expected sample, and DrugCentral loses half its
   relations. They also include food compounds, xenobiotics and central
   metabolites. Of the four endogenous molecules of the sample, acetyl-CoA
   stays only through its group, and two leave `hmdb`. The proton leaves too.
   So (b) removes drugs, but it also removes metabolites that HMDB has not
   measured.
4. **Narrow by status in the gate, not in the bucket.** The subset can keep
   all statuses and publish the best status of each row as a column. A dataset
   then picks the statuses it wants. This keeps one table for all datasets.
   It needs a column and a gate key that the plan does not have. This report
   measured no more than the figures above for it.

Any narrowing needs a code change: a status condition on the resource signal,
with a test. The status of a `structure_key` row then comes from its group,
as in this report.

## How to repeat the measurement

The scripts are in
`/tmp/claude-1029/-workspace-instances-metsig1/97c1ab55-a53d-4d90-8da3-8c1f34a8d599/scratchpad/stage2-hmdb/`.
Run them from that directory with
`/tmp/claude-1029/-workspace-instances-metsig1/97c1ab55-a53d-4d90-8da3-8c1f34a8d599/scratchpad/venv/bin/python`,
in this order:

| Script | What it does | Time |
|---|---|---|
| `build.py [release]` | Writes `rules_s1.yaml` and `rules_full.yaml`. Builds the stage-one and full variants with the status narrowings (a) to (d). Writes `mc_*`, `claim_*`, `assigned_*`, `grouped_*` and `hmdb_status.parquet`. With `2026.11.19`, builds (a) only. | 25 s |
| `relset.py <release>` | Builds the MetaLinksDB-shaped relation set once and writes `relset_<release>.parquet` | 10 s |
| `status.py` | Tasks 1 and 2: the status per bucket member and the destinations under each narrowing | 2 s |
| `gates.py` | Task 3: the gates, the losses by resource and the relations by status | 1 min |
| `lost_profile.py` | The molecules with relations that leave `hmdb` under (b) and (c) | 30 s |
| `diff_releases.py` | The rows that differ between the two releases | 20 s |
| `audit.py`, then `groupmates.py`, then `own_accession.py` | Task 4: draws the samples to `audit_<status>.json`, adds the group members and the own HMDB accession | 1 min |
| `classify.py` | Holds the audit classes, prints the counts and the sample tables | 1 s |

The outputs are `status.out`, `gates.out` and `classify.out` in the same
directory.

## Appendix: the sample tables

The columns are the label, the own HMDB accession, the InChIKey, the relations
in the set, the resource with most of them and the class. The last column gives
the row of the molecule under the narrowing: (c) for the predicted sample, (b)
for the expected sample.

<!-- vale off -->

### Predicted, all 27 eligible molecules

| # | Label | HMDB | InChIKey | Relations | Top partner resource | Class | Note | Under (c) |
|---|---|---|---|---|---|---|---|---|
| 1 | N-Oleoyl Histidine | `HMDB0241960` | `QJFAABXOMYWNGP-UHFFFAOYSA-N` | 1 | bindingdb (1) | Lipid species (generic) | N-acyl amino acid series | no row |
| 2 | 3-methylbut-2-enoyl-CoA | `HMDB0301419` | `BXIPALATIYNHJN-UHFFFAOYSA-N` | 21 | stitch (21) | Endogenous human metabolite | 3-methylcrotonyl-CoA, leucine breakdown | hmdb structure_key |
| 3 | N-Oleoyl Glutamic acid | `HMDB0241959` | `UUFVSGRELDGPGL-UHFFFAOYSA-N` | 1 | bindingdb (1) | Lipid species (generic) | N-acyl amino acid series | metabolic_network structure_key |
| 4 | (2E)-3-methylpent-2-enedioyl-CoA | `HMDB0301420` | `GXKSHRDAHFLWPN-UHFFFAOYSA-N` | 6 | stitch (6) | Endogenous human metabolite | 3-methylglutaconyl-CoA, leucine breakdown | chebi_role structure_key |
| 5 | N-Oleoyl Alanine | `HMDB0241953` | `CHQZBIVWDKADIA-UHFFFAOYSA-N` | 1 | bindingdb (1) | Lipid species (generic) | N-acyl amino acid series | no row |
| 6 | N-Oleoyl Aspartic acid | `HMDB0241956` | `XFBKKGFUMFHHHK-UHFFFAOYSA-N` | 1 | bindingdb (1) | Lipid species (generic) | N-acyl amino acid series | no row |
| 7 | N-Oleoyl Glycine | `HMDB0241916` | `HPFXACZRFJDURI-UHFFFAOYSA-N` | 1 | bindingdb (1) | Endogenous human metabolite | second record of quantified N-oleoylglycine | hmdb structure_key |
| 8 | 3-hydroxyhexanedioyl-CoA | `HMDB0301432` | `OTEACGAEDCIMBS-UHFFFAOYSA-N` | 12 | stitch (12) | Microbial metabolite | 3-hydroxyadipyl-CoA, bacterial adipate breakdown | hmdb structure_key |
| 9 | N-Oleoyl Arginine | `HMDB0241954` | `YZYYKAWTXPUEGT-UHFFFAOYSA-N` | 1 | bindingdb (1) | Lipid species (generic) | N-acyl amino acid series | no row |
| 10 | N-Oleoyl Asparagine | `HMDB0241955` | `MZSJQKKCFOSWBW-UHFFFAOYSA-N` | 1 | bindingdb (1) | Lipid species (generic) | N-acyl amino acid series | metabolic_network structure_key |
| 11 | 3-oxobutanoyl-CoA | `HMDB0301412` | `OJFDKHTZOUZBOS-UHFFFAOYSA-N` | 48 | stitch (48) | Endogenous human metabolite | acetoacetyl-CoA | hmdb structure_key |
| 12 | N-Stearoyl Serine | `HMDB0241949` | `CLHUCXCVFSEJRR-UHFFFAOYSA-N` | 1 | chembl (1) | Lipid species (generic) | N-acyl amino acid series | no row |
| 13 | Henicosanoylcarnitine | `HMDB0241596` | `DZUIUBMZQFPGGK-UHFFFAOYSA-N` | 5 | recon3d (5) | Lipid species (generic) | acylcarnitine series | metabolic_network entity |
| 14 | N-Myristoyl Glutamic acid | `HMDB0242051` | `MTJZWYHTZFVEGI-UHFFFAOYSA-N` | 2 | chembl (1) | Lipid species (generic) | N-acyl amino acid series | no row |
| 15 | Chenodeoxycholyltaurine | `HMDB0242411` | `BHTRKEVKTKCXOH-UHFFFAOYSA-N` | 5 | bindingdb (5) | Endogenous human metabolite | second record of quantified TCDCA | hmdb structure_key |
| 16 | 2-(1-oxoeicosa-5,8,11,14-tetraenylamino)acetic acid | `HMDB0241918` | `YLEARPUNMCCKMP-UHFFFAOYSA-N` | 375 | stitch (372) | Endogenous human metabolite | second record of quantified N-arachidonoylglycine | hmdb structure_key |
| 17 | N-Oleoyl Leucine | `HMDB0241962` | `UMOAAMQGRRCHPA-UHFFFAOYSA-N` | 1 | bindingdb (1) | Lipid species (generic) | N-acyl amino acid series | metabolic_network structure_key |
| 18 | N-Arachidonoyl Alanine | `HMDB0241989` | `ZECSOKFEQQDUCP-UHFFFAOYSA-N` | 1 | bindingdb (1) | Lipid species (generic) | N-acyl amino acid series | chebi_role structure_key |
| 19 | DG(8:0/18:1(12Z)-O(9S,10R)/0:0) | `HMDB0297175` | `SUPCGUICHYILCS-ZIMUHACTSA-N` | 1 | pathwaycommons (1) | Lipid species (generic) | oxidized diacylglycerol series | no row |
| 20 | N-Docosahexaenoyl Glycine | `HMDB0242014` | `BEYWKSIFICEQGH-UHFFFAOYSA-N` | 1 | bindingdb (1) | Lipid species (generic) | N-acyl amino acid series | no row |
| 21 | (2S,3R)-3-Hydroxy-2-methylbutanoyl-CoA | `HMDB0301410` | `PEKYNTFSOBAABV-UHFFFAOYSA-N` | 18 | stitch (18) | Endogenous human metabolite | isoleucine breakdown | hmdb structure_key |
| 22 | CoA 20:3(8,11,14) | `HMDB0301353` | `FJWJALRUNNZIBB-UHFFFAOYSA-N` | 16 | stitch (16) | Endogenous human metabolite | dihomo-gamma-linolenoyl-CoA | hmdb structure_key |
| 23 | (9Z,12Z,15Z,18Z,21Z)-Tetracosa-9,12,15,18,21-pentaenoylcarnitine | `HMDB0241635` | `AQYQJMYFFSCGBI-UHFFFAOYSA-N` | 14 | metatlas (12) | Lipid species (generic) | acylcarnitine series | metabolic_network entity |
| 24 | 2-Decenoyl-CoA | `HMDB0300871` | `MGNBGCRQQFMNBM-UHFFFAOYSA-N` | 24 | stitch (24) | Endogenous human metabolite | beta-oxidation intermediate | hmdb structure_key |
| 25 | (2E,6Z)-Dodeca-2,6-dienoyl-CoA | `HMDB0301034` | `BNPQDIKRZDRREL-UHFFFAOYSA-N` | 12 | stitch (12) | Endogenous human metabolite | linoleate beta-oxidation intermediate | hmdb structure_key |
| 26 | Chenodeoxycholylglycine | `HMDB0242392` | `GHCZAUBVMUEKKP-UHFFFAOYSA-N` | 4 | bindingdb (4) | Endogenous human metabolite | second record of quantified GCDCA | hmdb structure_key |
| 27 | Deoxycholylglycine | `HMDB0242414` | `WVULKSPCQVQLCU-UHFFFAOYSA-N` | 2 | stitch (1) | Endogenous human metabolite | second record of quantified GDCA | hmdb structure_key |

### Expected, 50 of 4,930 eligible molecules

| # | Label | HMDB | InChIKey | Relations | Top partner resource | Class | Note | Under (b) |
|---|---|---|---|---|---|---|---|---|
| 1 | 5-fluoro-1H-indole-2-carboxylic acid | `HMDB0303392` | `WTXBRZCVLDTWLP-UHFFFAOYSA-N` | 13 | chembl (6) | Other or unclear | fluorinated; HMDB says isolated from Lactobacillus | no row |
| 2 | 1-Deoxynojirimycin | `HMDB0035359` | `LXBIFEVIBLOUGU-UHFFFAOYSA-N` | 59 | stitch (33) | Food or plant compound | mulberry alkaloid | chebi_role structure_key |
| 3 | Ripretinib | `HMDB0304864` | `CEFJVGZHQAGLHS-UHFFFAOYSA-N` | 18 | drugcentral (6) | Drug or drug metabolite | kinase inhibitor | no row |
| 4 | 2-hydroxy-2,2-diphenylacetic acid (1,1-dimethyl-3-piperidin-1-iumyl) ester | `HMDB0015591` | `GKNPSSNBBWDAGH-UHFFFAOYSA-N` | 29 | stitch (16) | Drug or drug metabolite | mepenzolate | no row |
| 5 | Hexylcaine | `HMDB0014616` | `DKLKMKYDWHYZTD-UHFFFAOYSA-N` | 12 | stitch (6) | Drug or drug metabolite | local anaesthetic | no row |
| 6 | Glycitin | `HMDB0002219` | `OZBAVEKZGSOMOJ-UHFFFAOYSA-N` | 1 | chembl (1) | Food or plant compound | soy isoflavone | chebi_role structure_key |
| 7 | Antrafenine | `HMDB0015488` | `NWGGKKGAFZIVBJ-UHFFFAOYSA-N` | 7 | stitch (4) | Drug or drug metabolite | analgesic | no row |
| 8 | Valsartan | `HMDB0014323` | `ACWBQPMHZXGDFX-QFIPXVFZSA-N` | 129 | pathwaycommons (74) | Drug or drug metabolite | angiotensin receptor blocker | no row |
| 9 | Methylphenidate | `HMDB0014566` | `DUGOZIWVEXMGBE-UHFFFAOYSA-N` | 91 | stitch (65) | Drug or drug metabolite | stimulant | no row |
| 10 | Acetyl-CoA | `HMDB0001206` | `ZSLZBFCDCINBPY-ZSJPKINUSA-N` | 878 | stitch (425) | Endogenous human metabolite | stereo-free record HMDB0247926 is detected | hmdb structure_key |
| 11 | 8-O-4'-Diferulic acid | `HMDB0304640` | `AHJKRLASYNVKDZ-UHFFFAOYSA-N` | 38 | pathwaycommons (38) | Food or plant compound | plant cell wall | chebi_role entity |
| 12 | Potassium nitrate (KNO3) | `HMDB0303531` | `FGIUAXJPYTZDNR-UHFFFAOYSA-N` | 40 | brenda (28) | Other or unclear | inorganic salt, food additive | no row |
| 13 | 1-Phenyl-1-propanol | `HMDB0031627` | `DYUQAZSOFZSPHD-UHFFFAOYSA-N` | 1 | chembl (1) | Food or plant compound | flavouring | no row |
| 14 | Dodecadienoic acid | `HMDB0302621` | `HQSBWLQFLLMPKC-UHFFFAOYSA-N` | 35 | stitch (35) | Lipid species (generic) | double-bond positions not given | no row |
| 15 | Rimonabant | `HMDB0015623` | `JZCPYUJPEARBJL-UHFFFAOYSA-N` | 103 | pathwaycommons (38) | Drug or drug metabolite | CB1 inverse agonist | no row |
| 16 | Morusinol | `HMDB0030620` | `AFOKZNPZDXHDHD-UHFFFAOYSA-N` | 8 | chembl (3) | Food or plant compound | mulberry | no row |
| 17 | Aldophosphamide | `HMDB0060433` | `QMGUSPDJTPDFSF-UHFFFAOYSA-N` | 59 | stitch (59) | Drug or drug metabolite | cyclophosphamide metabolite | metabolic_network entity |
| 18 | Valdecoxib | `HMDB0005033` | `LNPDTQAFDNKSHK-UHFFFAOYSA-N` | 167 | chembl (40) | Drug or drug metabolite | NSAID | no row |
| 19 | Fumitremorgin C | `HMDB0038642` | `DBEYVIGIPJSTOR-FHWLQOOXSA-N` | 16 | chembl (7) | Microbial metabolite | Aspergillus | chebi_role entity |
| 20 | Amidosulfonic acid | `HMDB0034830` | `IIACRCGMVDHOTQ-UHFFFAOYSA-N` | 61 | brenda (25) | Other or unclear | sulfamic acid, packaging contaminant | no row |
| 21 | Ganciclovir | `HMDB0015139` | `IRSCQMHQWWYFCW-UHFFFAOYSA-N` | 43 | chembl (16) | Drug or drug metabolite | antiviral | no row |
| 22 | CoA 10:2(2E,4Z) | `HMDB0304510` | `FASAKYLWSRDQOH-IMVFQKDNSA-J` | 2 | pathwaycommons (2) | Endogenous human metabolite | linoleate beta-oxidation intermediate | metabolic_network entity |
| 23 | 1-octadecylglycerone 3-phosphate | `HMDB0011142` | `DSLZSTBYLRGPIK-UHFFFAOYSA-N` | 2 | pathwaycommons (2) | Endogenous human metabolite | ether lipid synthesis | metabolic_network structure_key |
| 24 | Lochnocarpol A | `HMDB0302553` | `HCNLDGTUMBOHKT-UHFFFAOYSA-N` | 3 | chembl (2) | Food or plant compound | prenylated flavanone | chebi_role structure_key |
| 25 | Liquiritigenin | `HMDB0029519` | `FURUXTVZLHCCNA-AWEZNQCLSA-N` | 34 | pathwaycommons (16) | Food or plant compound | liquorice flavanone | chebi_role entity |
| 26 | 4-Carboxynevirapine | `HMDB0060759` | `WDXCMIDQWGYHIN-UHFFFAOYSA-N` | 1 | pathwaycommons (1) | Drug or drug metabolite | nevirapine metabolite | no row |
| 27 | 1,2-Cyclohexanedione | `HMDB0031344` | `OILAIQUEIWYQPH-UHFFFAOYSA-N` | 20 | brenda (18) | Food or plant compound | flavouring | metabolic_network entity |
| 28 | N-hydroxy-4-aminobiphenyl | `HMDB0062723` | `MYVLYOJYVMLSFA-UHFFFAOYSA-N` | 19 | stitch (14) | Other or unclear | xenobiotic, carcinogen metabolite | chebi_role entity |
| 29 | N(tele)-methylhistaminium | `HMDB0062574` | `FHQDWPCFSJMNCT-UHFFFAOYSA-O` | 1 | pathwaycommons (1) | Endogenous human metabolite | charge form of quantified 1-methylhistamine | hmdb structure_key |
| 30 | Curlone | `HMDB0037068` | `JIJQKFPGBBEJNF-UHFFFAOYSA-N` | 1 | brenda (1) | Food or plant compound | turmeric | no row |
| 31 | Lutetium | `HMDB0302224` | `PSDMOPINLDTFSZ-UHFFFAOYSA-N` | 1 | brenda (1) | Other or unclear | element | no row |
| 32 | Clorazepate | `HMDB0014766` | `XDDJGVMJFWAHJX-UHFFFAOYSA-N` | 35 | pathwaycommons (16) | Drug or drug metabolite | benzodiazepine | no row |
| 33 | Loracarbef | `HMDB0014590` | `JAPHQRWPEGVNBT-UTUOFQBUSA-N` | 2 | drugcentral (2) | Drug or drug metabolite | antibiotic | no row |
| 34 | Bevantolol | `HMDB0015409` | `HXLAFSUPPDYFEO-UHFFFAOYSA-N` | 23 | stitch (11) | Drug or drug metabolite | beta blocker | no row |
| 35 | 4-Hydroxy-4-(methylnitrosoamino)-1-(3-pyridinyl)-1-butanone | `HMDB0062404` | `UIPLZNGNFRQXQE-UHFFFAOYSA-N` | 71 | stitch (71) | Other or unclear | xenobiotic, tobacco nitrosamine metabolite | metabolic_network entity |
| 36 | (+)-gallocatechin gallate | `HMDB0041598` | `WMBWREPUVVBILR-GHTZIAJQSA-N` | 6 | brenda (4) | Food or plant compound | tea catechin | hmdb structure_key |
| 37 | S-Allylmercaptocysteine | `HMDB0302210` | `WYQZZUUUOXNSCS-YFKPBYRVSA-N` | 2 | stitch (2) | Food or plant compound | garlic | no row |
| 38 | Dihydrocapsaicin | `HMDB0038457` | `XJQPQKLURWNAAH-UHFFFAOYSA-N` | 27 | pathwaycommons (19) | Food or plant compound | chili pepper | no row |
| 39 | 2-Undecyl-4(1H)-quinolinone | `HMDB0032996` | `KKSIGOUTTVXAKQ-UHFFFAOYSA-N` | 3 | bindingdb (2) | Food or plant compound | rue alkaloid | no row |
| 40 | Cassitoroside | `HMDB0041289` | `PGDDDJBSORSPAG-UHFFFAOYSA-N` | 1 | bindingdb (1) | Food or plant compound | Cassia seed | no row |
| 41 | (S)-carnitinium | `HMDB0062634` | `PHIQHXFUZVPYII-LURJTMIESA-O` | 11 | brenda (11) | Other or unclear | D-carnitine cation, enantiomer of the endogenous form | hmdb structure_key |
| 42 | Vaporole | `HMDB0001382` | `OWFXIOWLTKNBAP-UHFFFAOYSA-N` | 1 | chembl (1) | Drug or drug metabolite | isoamyl nitrite | no row |
| 43 | 2,5-dihydroxybenzaldehyde | `HMDB0004062` | `CLFRCXCBWIQVRN-UHFFFAOYSA-N` | 8 | stitch (5) | Other or unclear | gentisaldehyde, origin unclear | chebi_role entity |
| 44 | Adenosyl cobinamide | `HMDB0006903` | `KQXSPGAEBZWHMC-QMUWONGRSA-M` | 1 | pathwaycommons (1) | Microbial metabolite | cobamide synthesis intermediate | chebi_role entity |
| 45 | 1,3-diphenylurea | `HMDB0032066` | `GWEHVDNNLFDJLR-UHFFFAOYSA-N` | 24 | chembl (10) | Food or plant compound | HMDB: coconut milk | chebi_role entity |
| 46 | Phenyl salicylate | `HMDB0032018` | `ZQBAKBUEJOMQEX-UHFFFAOYSA-N` | 5 | pathwaycommons (4) | Drug or drug metabolite | salol | no row |
| 47 | Sulfathiazole | `HMDB0015619` | `JNMRHUJNCSQMMB-UHFFFAOYSA-N` | 34 | brenda (16) | Drug or drug metabolite | sulfa drug | no row |
| 48 | Polyporic acid | `HMDB0302338` | `HZKFHDXTSAYOSN-UHFFFAOYSA-N` | 3 | stitch (2) | Microbial metabolite | fungal pigment | metabolic_network structure_key |
| 49 | Diethyl sulfide | `HMDB0031666` | `LJSQFQKUNVCTIA-UHFFFAOYSA-N` | 2 | pathwaycommons (2) | Food or plant compound | flavour | no row |
| 50 | 2-Nitro-p-cresol | `HMDB0059811` | `SYDNSSSQVSOXTN-UHFFFAOYSA-N` | 2 | brenda (2) | Other or unclear | xenobiotic, nitrophenol | no row |

<!-- vale on -->
