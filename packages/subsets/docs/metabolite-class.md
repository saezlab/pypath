# Metabolite classes

> **Status: planned.** This page describes the subset as designed in the plan
> `saezverse/projects/omnipath/plans/2026-11-metabolite-classification-logic.md`.
> The code does not exist yet. The figures come from a trial run on the
> 2026.11.19 release.

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
| 2 | `chebi_role` | ChEBI says it has the role "metabolite", or a narrower role such as "human metabolite" (planned) |
| 3 | `metabolic_network` | a reaction in KEGG, Rhea, Recon3D or Human-GEM uses it or makes it (planned) |
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

**`metabolite_class_report`**: the figures of each build.

- `class_size`: the number of molecules in each bucket.
- `overlap`: how many molecules each combination of sources claims.
- `shared_key`: the groups that pulled in variants, and the variants.
- `no_key`: the molecules a source claimed that have no InChIKey, so the
  subset cannot add variants for them.

On the 2026.11.19 trial run:

| Bucket | Molecules | Added through a group | Without an InChIKey |
|---|---|---|---|
| `hmdb` | 272,030 | 54,114 | 21 |
| `chebi_role` | 11,767 | 2,698 | 930 |
| `metabolic_network` | 15,991 | 3,078 | 5,882 |

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

**In a network view,** name the buckets in the gate. The query macros apply the
filter to the chemical end of each relation:

```python
"metabolite_classes": ["hmdb"],
```

An unknown bucket name gives an empty result, not the whole network.

**Choose the buckets on purpose.** A bucket is the first source that says yes.
So `metabolic_network` does not mean every molecule in a reaction: a molecule
that HMDB also lists is in `hmdb`. To get every molecule a source supports,
read `supported_by`.

**Leave out the variants** if you want only the molecules a source named, by
adding `AND m.match_kind = 'entity'`.

**Drop entities without an InChIKey** if your subset joins other tables on the
InChIKey, by adding `AND m.inchikey IS NOT NULL`.

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
