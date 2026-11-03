# pypath `master` changes missing from `parquet-migration`

Checked 2026-10-25 against `origin/master`.

- **Branch point:** `6176be278` (2026-07-05)
- **`master` tip:** `a9553acd9` (2026-09-25)
- **`parquet-migration` tip:** `e2de1c22b` (2026-10-19)

`master` has 59 commits since the branch point: 55 normal commits and 4 merges, all in the `pypath` repo. None of the 55 is on the branch as a cherry-pick. One change, the KEGG fix, was redone separately on the branch.

**Don't merge or cherry-pick `master` into the branch.** `master` code uses the old builder vocabulary (`EntityTypeCv`, `IdentifierNamespaceCv`, `INTERACTION` with `MembershipBuilder`). The branch uses biolink classes (`ChemicalEntity`, `Protein`), `omnipath_core.naming.Namespace` and `RelationBuilder` with `slots.*` predicates. Each item below has to be redone by hand. The template for interaction resources is `pypath/inputs_v2/mrclinksdb.py` on the branch.

## Status

MetaChatDB and PMI-DB were ported on 2026-10-25, and MetaChatDB is also in the MetalinksDB preset.

Hormone2Cell was ported too, but kept off this branch for now. The module and its tests are on the local pypath branch `hormone2cell-biolink`. In a validation build its 20 chemical hormones, which have names only, did not resolve; it needs a reviewed ChEBI crosswalk.

HSM and the fixes in section 2 are not ported.

What each port leaves out:

- **MetaChatDB:** the rows that are copies of MEBOCOST, MRClinksDB and MetalinksDB, as Nico's version did. Rows whose sensor is both a receptor and a transporter become two statements.
- **PMI-DB:** the tested non-interactions. The protein table's GO terms are a copy of UniProt's and aren't emitted. The metabolite and gene tables are joined into each pair, so they aren't separate datasets.
- **Hormone2Cell:** receptors that aren't "established", group entries and prohormones. Also not ported: the expression atlas (S3 and S6), and the producing gene sets in S2D except for naming a protein hormone's gene.

## 1. New resources (only on `master`)

| Resource | Commits | Status on `master` | What it holds |
|---|---|---|---|
| **MetaChatDB** | `079c9ae75` `2e2b45a5c` `28c1b16ca` `6c6dfd05c` `a9553acd9` | Complete | Metabolite–sensor pairs (receptors, nuclear receptors, transporters), human and mouse, with PMIDs. The parser keeps only rows MetaChatDB itself reports, which drops rows it imported from MEBOCOST, MRClinksDB and MetalinksDB. |
| **PMI-DB** | `2094b8038` `dd72c2419` `27d9d5d18` `523a54352` `5c81903ac` `530440e26` `48bd12fd6` `9954fb7ca` `0a1031580` | Complete | Protein–metabolite interactions from experimental assays (e.g. LiP-SMap), with protein info (GO) and metabolite info. Only rows with `interaction == TRUE` are kept. |
| **Hormone2Cell** | `768e00e6c` `e27387202` `d9ac30cda` `fb9bb04d6` … `c9ad989b3` `0c818a41a` (26 commits) | Complete | Hormone → receptor → cell type, from the Hormone Cell Atlas supplementary tables. |
| **HSM** (Human Serum Metabolome) | `184fd67b4` `6df041938` `f2cd8113a` `bb0ec1139` | **Unfinished**: download and parser only; schema and `Resource` are commented out | Serum metabolite records (XML), similar to HMDB. |

### What porting one resource takes

1. **pypath:** add the `ResourceCv` entry to `pypath/internals/cv_terms/resources.py`. Nico's accession numbers are all still free on the branch: OM:1182 MetaChatDB, OM:1169 PMI-DB, OM:1168 Hormone2Cell, OM:1199 HSM.
2. **pypath:** rewrite `pypath/inputs_v2/<name>.py` for the branch:
   - `ResourceConfig` takes `topics`, `content`, `knowledge_origin` and `summary`, replacing `primary_category`.
   - Pairs become a `RelationBuilder`:
     - binding uses `slots.interacts_with`;
     - transport uses `slots.affects` plus `TRANSPORT_QUALIFIERS`, as in the MRClinksDB transporter datasets.
   - Species uses `slots.in_taxon` with `NCBITaxon:<id>`, PMIDs use `slots.publications` with `PMID:<id>`, and identifiers use `Namespace.*`.
   - Compound classes use the ChemOnt association pattern from `pypath/inputs_v2/_class_terms.py`, which replaces Nico's ClassyFire class → `MoleculeSubtypeCv` map.
3. **pypath:** add a test, e.g. port `tests/test_pmidb.py` from `master`.
4. **omnipath-build:**
   - Regenerate `packages/core/omnipath_core/vocab/resource_metadata.yaml` with `scripts/generate_resource_metadata.py`.
   - Add a `Namespace` entry in `packages/core/omnipath_core/naming.py` only if the resource issues its own identifiers. MetaChatDB doesn't.
5. **omnipath-build, to reach MetalinksDB:** add the resource to `_CURATED_SOURCES` in `packages/subsets/omnipath_subsets/network_views/_definitions.py`.
6. **omnipath-build, to reach COSMOS:** COSMOS currently takes only reactions from `REACTION_SOURCES = ("kegg", "rhea", "metatlas", "recon3d")` (`packages/subsets/omnipath_subsets/cosmos/build.py`). Metabolite–receptor and metabolite–transporter pairs don't go into it, from MetaChatDB or anywhere else. Including them needs a new COSMOS component, not just a list entry. **This is an open decision.**
7. Rebuild on nicesrv, then check row counts and identifier resolution.

### Points to decide per resource

- **MetaChatDB:**
  - `Sensor.Type` gives receptor, nuclear receptor or transporter. Transporter rows should probably become transport relations rather than `interacts_with`.
  - Sensors are identified only by gene symbol and protein name. There is no UniProt ID, so resolution depends on the symbol plus `in_taxon`.
  - `Long.Range.Channel` (Blood, Urine, …) and the WikiPathways IDs need a place in the biolink model, or should be dropped.
- **PMI-DB:**
  - Metabolites come as KEGG and HMDB IDs.
  - The `formula` column is tagged `SMILES` on `master`. Check what it actually holds.
  - The yeast taxon `4932` is marked "presumably" in Nico's code.
- **Hormone2Cell:**
  - `master` adds new CV terms `HORMONE` (OM:0051) and `RECEPTOR` (OM:0052). On the branch these are biolink classes or roles, not entity types.
  - It also uses `ABBREVIATED_NAME`.
- **HSM:** needs its schema written first. It overlaps heavily with HMDB, so it may not be worth porting now.

## 2. Fixes to code that already exists on the branch

| Commit | Change | On branch? | What to do |
|---|---|---|---|
| `d2cd08e38` `826754d72` | KEGG `conv/pubchem` IDs are PubChem **substance** IDs, not compound IDs | Yes, done differently (`Namespace.PUBCHEM_SUBSTANCE` in `kegg.py`) | Nothing. Optionally port the test `tests/test_kegg_pubchem_namespace.py`. |
| `4308c722b` | `inputs.uniprot`: `organism_id:None` leaked into UniProt queries | No (`pypath/inputs/uniprot.py:493`) | Apply the one-line fix (`and val is not None`). |
| `19b30006e` | `reactome_old`: row generator referred to an undefined `fields` | No (`pypath/inputs/reactome_old.py:107`) | Apply the one-line fix. |
| `292e872bb` | HMDB: parse and emit the FooDB cross-reference | No | Add `foodb_id` to `parsers/hmdb.py` and a `Namespace.FOODB` identifier to `hmdb.py` (`FOODB` already exists in `naming.py`). |
| `941901da5` | HMDB: a failed ChemOnt download no longer aborts the whole HMDB load | No (`hmdb.py:265`) | Wrap `chemont_name_map(...)` in `try`/`except`, log a warning, and fall back to `{}`. |
| `d4078d5b4` | `ResourceConfig.mints`: which identifier namespaces a resource is the authority for (ChEBI mints `chebi`, …) | No equivalent on the branch | Check first whether the resolver still needs this, since source authority may already be handled elsewhere. If it does, port it with `Namespace` values. |
| `33f37fbaa` | New identifier type `LIPID_NAME` (OM:0209, Goslin lipid nomenclature) | No | Port it if the metabolic-domain lipid work should emit standardized lipid names; otherwise skip it. |

## Suggested order

1. The one-line fixes and the HMDB fixes. They're small, separate, and need no decisions.
2. MetaChatDB, which is what the colleague asked for. It goes into MetalinksDB via `_CURATED_SOURCES`, and COSMOS is decided separately.
3. PMI-DB, then Hormone2Cell.
4. HSM, which has to be finished on `master` or written straight for the branch.
