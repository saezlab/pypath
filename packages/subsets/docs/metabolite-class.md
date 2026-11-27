# Metabolite classes

> **Status: stage one.** The subset fills the bucket `hmdb`. Buckets 2 and 3 are
> in the code and switch on after the audit of stage two of the plan
> `saezverse/projects/omnipath/plans/2026-11-metabolite-classification-logic.md`.
> The code is in `packages/subsets/omnipath_subsets/metabolite_class/`.

The subset `metabolite_class` answers one question: does this molecule count as
a metabolite? A network dataset uses the answer to keep the relations of
metabolites. It leaves out the relations of other chemicals, such as drug
candidates and screening compounds.

## The idea in plain terms

The release holds 3.6 million chemicals. Most come from drug and screening
catalogs, and only a small share are metabolites. No single resource lists every
metabolite, so the subset asks a few resources in a fixed order and records the
first one that says yes.

| Order | Bucket | The molecule is in the bucket when… |
|---|---|---|
| 1 | `hmdb` | the Human Metabolome Database (HMDB) has a record of it |
| 2 | `chebi_role` | ChEBI says it has the role "metabolite", or a narrower role such as "human metabolite" (stage two) |
| 3 | `metabolic_network` | a reaction in KEGG, Rhea, Recon3D or Human-GEM uses it or makes it (stage two) |
| — | `none` | nothing says yes |

A molecule gets one bucket, the first that applies. A later source never
overrides an earlier one, and no source takes a molecule out. Buckets 2 and 3
switch on after an audit. Until then the subset fills `hmdb` alone.

## Variants of one molecule

The same molecule often appears more than once in the release: as a charged
form, a salt, or another stereoisomer. Each form has its own InChIKey, so each
is its own entity. Resolution keeps them apart on purpose (see
[resolution](../../resolver/docs/resolution.md)).

The subset groups them. The first 14 characters of an InChIKey describe the
skeleton of the molecule, without charge or stereochemistry. Every entity with
the same first 14 characters joins one group, and every member of the group
gets the best bucket that any member has.

For example, HMDB records cortisol under one InChIKey. KEGG and FooDB record a
stereoisomer of it under another. Both keys start with `JYGXADMDTFJGBT`, so both
are in the bucket `hmdb`.

A row shows how it got its bucket:

- `match_kind = entity`: a source said yes to this molecule itself.
- `match_kind = structure_key`: a source said yes to another member of the
  group. `via_entity_key` names that member.

`exclusion_metabolites` uses the same two values and the same column.

Some entities have no InChIKey. Examples are lipids known only by name, and
records the resolver could not match to a chemical database (native
fallbacks). They join no group and keep their own bucket. The subset keeps
them, so a dataset that uses their native identifier can still find them.

## What the tables hold

**`metabolite_class`**: one row per molecule in a bucket. A molecule with no row
is in `none`.

| Column | Content |
|---|---|
| `entity_key` | The entity in the release |
| `label` | A readable name: the HMDB name, else the ChEBI name, else the release label |
| `inchikey` | The InChIKey, or empty when the molecule has none |
| `hmdb` | The HMDB accession, for example `HMDB0000063`, or empty |
| `metabolite_class` | The bucket |
| `match_kind` | `entity` or `structure_key`, as above |
| `via_entity_key` | For a `structure_key` row, the group member that set the bucket |
| `supported_by` | The resources that say yes to the molecule itself, for example `[hmdb, chebi, rhea]` |
| `build_id` | The digest of the release manifest the build read |

**`metabolite_class_report`**: the figures of each build.

Each report row has a `kind`, a bucket, a count `n` and the details as JSON text.

- `class_size`: the number of molecules in each bucket.
- `overlap`: how many molecules each combination of sources claims.
- `shared_key`: the groups that pulled in variants, and the variants.
- `no_key`: the molecules a source claimed that have no InChIKey, so the
  subset cannot add variants for them.

On the 2026.11.19 release, with all three signals switched on:

| Bucket | Molecules | Added through a group | Without an InChIKey |
|---|---|---|---|
| `hmdb` | 272,030 | 54,114 | 21 |
| `chebi_role` | 4,608 | 768 | 257 |
| `metabolic_network` | 18,494 | 3,957 | 6,073 |

Stage one fills `hmdb` alone, with the same 272,030 molecules: the other signals
do not change which molecules HMDB or its groups claim. The ChEBI role reaches
few molecules until stage two repairs the ChEBI parser, which leaves out the
role tree today. The build of the release takes some ten seconds.

## Using it in your subset

The build makes the subset from a published release, like the other Parquet
subsets. It sits next to them:

```
<data_root>/subsets/<release>/metabolite_class/metabolite_class.parquet
```

**Keep only metabolites.** Join your chemical end to the table and pick the
buckets you want:

```sql
SELECT r.*
FROM my_relations r
JOIN metabolite_class m ON m.entity_key = r.chemical_entity_key
WHERE m.metabolite_class IN ('hmdb');
```

**In the query layer,** every relation query takes the filter
`metabolite_classes`: `relation_page`, `relations`, `relation_counts`,
`relation_facets` and `relation_filter_options`. A relation passes when one of
its ends is in a named bucket:

```python
client.relations("NR3C1", metabolite_classes=["hmdb"])
```

The API takes the same filter as `metaboliteClasses`. An unknown bucket name
gives an empty result, not the whole network. A release built without this
subset refuses the filter with an error, and answers as before without it.

**Choose the buckets on purpose.** A bucket is the first source that says yes.
So `metabolic_network` does not mean every molecule in a reaction: a molecule
that HMDB also lists is in `hmdb`. To get every molecule a source supports,
read `supported_by`.

**Leave out the variants** if you want only the molecules a source named, by
adding `AND m.match_kind = 'entity'`.

**Drop entities without an InChIKey** if your subset joins other tables on the
InChIKey, by adding `AND m.inchikey IS NOT NULL`.

## How the gate decides

Think of a network dataset as a guest list. The release holds millions of
candidate pairs: a chemical and the molecule it acts on, usually a protein. The
gate is a doorman with a checklist. A pair gets in only when every check passes.
The first "no" turns it away.

| # | Question | Gate key | Turns away, for example |
|---|---|---|---|
| 1 | Is one side a chemical? | `entity_types` | A protein–protein pair |
| 2 | Is the other side a molecule, such as a protein or gene? | `partner_types` | A chemical linked to an organism or an ontology term |
| 3 | Is the link an interaction (`interacts_with`, `affects`, `associated_with`)? | `predicates` | "This chemical has the role antibacterial", a label and not an interaction |
| 4 | Is the chemical a metabolite? | `metabolite_classes` | A drug candidate from a screening catalog |

Checks 1 to 3 are the existing gate of the MetaLinksDB view. Check 4 reads this
subset.

**How check 4 decides.** The doorman holds a list, the table
`metabolite_class`. Each metabolite on it carries a badge, its bucket, which
names the source that vouched for it. The dataset names the badges it accepts,
for example `["hmdb"]`. A pair gets in when its chemical is on the list with an
accepted badge. A protein never carries a badge, so the check reads the chemical
side only.

### Four example pairs

These cases follow a build on the 2026.11.19 release.

**Cortisol and its receptor NR3C1.** Checks 1 to 3 pass. Cortisol is on the
list with the badge `hmdb`. The pair gets in.

**A stereoisomer of cortisol and the same receptor.** KEGG records cortisol
with a different stereo notation, so the release holds a separate entity. HMDB
does not list that entity. Its InChIKey starts with the same 14 characters as
the cortisol of HMDB, so the subset gives it the badge `hmdb`. Its row has
`match_kind = structure_key`, and `via_entity_key` names the cortisol of HMDB.
The pair gets in.

Without the groups, a notation difference would turn this pair away. On the
2026.11.19 release, the groups add some 279,000 pairs to the gate `["hmdb"]`.

**A screening compound of ChEMBL and a kinase.** Checks 1 to 3 pass. No source
vouches for the compound as a metabolite, so it is not on the list. The doorman
turns the pair away. Most pairs of ChEMBL and BindingDB end here. The gate
`["hmdb"]` keeps 2.5 per cent of the ChEMBL relations and 1.9 per cent of the
BindingDB relations.

**A role of cortisol from Guide to Pharmacology.** The link says that cortisol
has an antibacterial role. Check 3 fails, because a role is a label. The doorman
turns the pair away before check 4.

### What the accepted badges give

The badge names the first source that vouched, in the order HMDB, ChEBI,
reactions.

| The dataset accepts | It gets | Pairs |
|---|---|---|
| `["hmdb"]` | Every molecule that HMDB lists, with its stereo and charge variants | 1,349,547 |
| `["hmdb", "chebi_role", "metabolic_network"]` | Every metabolite that any source vouches for | 1,419,750 |
| `["metabolic_network"]` | The metabolites that reactions know and HMDB and ChEBI do not | 38,491 |

The pairs are MetaLinksDB-shaped relations of the 2026.11.19 release. The figure
for `["hmdb"]` holds for stage one. The other two come from a trial that read
45 ChEBI roles by name, and stage two measures them again.

The third row does not hold every molecule of a reaction. A molecule that HMDB
also lists carries the badge `hmdb`.

### A misspelt badge admits nobody

A dataset that asks for a badge that does not exist, for example `["hmbd"]`,
gets an empty result. The other choice is to admit every pair, and then a typo
switches the filter off without a sign. An empty result shows the mistake at
once.

**In one line:** the gate keeps a pair when the link is an interaction and the
chemical is on the list with an accepted badge. The variants of a listed
metabolite share its badge.

## Tests

The tests of the subset are `packages/subsets/tests/test_metabolite_class.py`,
next to the tests of `exclusion_metabolites`. They build the subset over the
small release of `metabolite_class_fixture.py`, written straight to Parquet. The release holds cortisol and its
stereoisomer, a charged form, a ChEMBL compound, a native fallback, a lipid, two
reaction molecules and a receptor. The tests need no database.

Each rule on this page has one test, and the name of the test states the rule.
For example, a test shows that the stereoisomer of cortisol gets the badge
`hmdb` through its group. The tests of the filter are
`packages/query/python/tests/test_metabolite_classes.py` and
`packages/api/tests/test_metabolite_classes_filter.py`. One of them shows that a
misspelt badge gives an empty result.

One more test compares the finished tables with the frozen rows of
`golden/metabolite_class.json`. A
change to the subset that alters its output fails that test. Review the
difference first. When the change is deliberate, rewrite the frozen rows:

```sh
UPDATE_GOLDEN=1 uv run pytest packages/subsets/tests/test_metabolite_class.py
```

A new signal comes with a test of its own. Change this page in the same commit.

## What to keep in mind

- The bucket `hmdb` includes HMDB's predicted and expected entries, not only
  measured ones.
- The chemicals in the release carry no species. The buckets do not say whether
  a molecule occurs in humans.
- The release does not record which entities are native fallbacks, so the
  table cannot mark them.
- `inchikey IS NULL` finds every entity without a structure key, fallbacks
  included.
- A drug can be a metabolite. Cortisol, adenosine and lactate are both, and
  they are in `hmdb`.
