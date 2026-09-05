# Modeling inputs with Biolink

Use the installed, pinned `biolink-model` package (currently 4.4.4) as the biological
schema. `signor.py`, `chebi.py`, and `uniprot.py` are the rewritten examples.
Rewrite remaining inputs individually; do not add a migration adapter.

## Boundaries

- `parsers/`: decode source formats and preserve source values. Source column names,
  accession parsing and exact source crosswalks belong here or in the resource module.
- Resource module: choose native Biolink classes, predicate slots, annotation slots
  and enum values. A source crosswalk must explain the meaning of each mapping;
  never infer biological meaning from substrings or similar labels.
- `tabular_builder` and `omnipath_core`: validate and serialize primitives. No source
  repair rules, biological alias dictionaries or fallback predicates belong here.
- Consolidation: merge identical qualified statements and retain evidence. Qualifiers
  are part of statement identity, not incidental metadata.
- API and UI: filter stored canonical values and use schema-derived vocabulary labels.
  Display projections must not become another biological vocabulary.

## Choosing a representation

1. Pick the most specific Biolink class justified by the source. For explicit source vocabulary concepts (currently GO terms and UniProt
   keywords), use `OntologyClass` as the primary type so users can filter them
   together in the existing entity browser. This is an intentional application
   convention: Biolink defines `OntologyClass` as a mixin and normally recommends
   a biological primary category. Preserve source identifiers and hierarchy;
   do not infer ontology membership merely from the presence of an identifier.
   SIGNOR and ChEBI both use `ChemicalEntity` for chemicals, including individual
   small molecules, so the browser exposes one chemical entity type. Do not
   split out `SmallMolecule` in these rewritten inputs.
2. Choose a predicate by its meaning, domain, range, symmetry and deprecation status
   in the pinned schema. Use `subclass_of` for ontology `is_a`; do not turn hierarchy
   into generic association. Explicit broad non-causal associations may use
   `associated_with`, but this does not assert causation, localization or function.
3. Use predicate plus qualifiers when the source supports them. SIGNOR uses `affects`
   with `DirectionQualifierEnum` and, where specified, `GeneOrGeneProductOrChemicalEntityAspectEnum`.
   Omit an unspecified qualifier; never invent a default aspect or direction.
4. Use slots for attributes: `name`, `synonym`, `description`, `in_taxon`,
   `publications`, `has_biological_sequence`, `has_chemical_formula`, etc.
   Narrative descriptions remain text; they are not entity identifiers or inferred edges.
5. Use `Namespace` for identifier namespaces. It is an application identifier registry,
   not an extension of Biolink's biological types. Add a member only for an actual
   identifier system. Translation datasets use the same members as graph inputs.
6. If no faithful Biolink representation exists, first adapt the representation
   (for example, preserve a narrative as `description`). Otherwise retain a published
   external ontology attribute CURIE with its original value and document it below.
   Do not mint an OmniPath term merely because an upstream label looks convenient.

Examples:

```python
from biolink_model.datamodel.model import Protein, slots, DirectionQualifierEnum
from omnipath_core.naming import Namespace
from pypath.internals.tabular_builder import CV, EntityBuilder, IdentifiersBuilder, AnnotationsBuilder

protein = EntityBuilder(
    entity_type=Protein,
    identifiers=IdentifiersBuilder(CV(term=Namespace.UNIPROT, value=lambda row: row['accession'])),
    annotations=AnnotationsBuilder(
        CV(term=slots.in_taxon, value=lambda row: f"NCBITaxon:{row['taxon']}"),
    ),
)
# In a relation builder:
direction = CV(term=slots.object_direction_qualifier, value=DirectionQualifierEnum.increased)
```

Ontology CURIE values retain their prefixes (`NCBITaxon:9606`, `PMID:123`, `GO:0005739`).
Namespace-specific accession parsing belongs to the input. Do not create redundant
namespace aliases. The `CV` helper name does not require a legacy CV enum.

## Decisions and remaining limits

| Source / field | Representation | Reason and path to a richer model |
| --- | --- | --- |
| ChEBI chemical branch | `ChemicalEntity` | Includes chemical classes, ions and mixtures; asserting `SmallMolecule` for every term would overstate the source. Inclusion follows ChEBI's chemical entity ancestry. Roles are not chemical entity nodes. |
| ChEBI `is_a`, `BFO:0000051` | `subclass_of`, `has_part` | Explicit source semantics; `BFO:0000051` is an exact mapping in the pinned schema. |
| Other ChEBI relationships | Original RO CURIE as an external attribute, target ChEBI CURIE as value | Preserve conjugate, stereochemical and role assertions without collapsing them into broad graph edges. Targets outside the chemical branch are also retained. Add typed role nodes and `has_chemical_role` if role traversal is needed; review specialized chemical relations individually before adding graph predicates. |
| ChEBI mass, monoisotopic mass, charge | `chemrof:mass`, `chemrof:monoisotopic_mass`, `chemrof:charge` | Published source ontology attributes, not new OmniPath terms. `chemrof` expands to `https://w3id.org/chemrof/`. Values retain the source's numeric strings and source-defined units. Formula uses Biolink `has_chemical_formula`. |
| UniProt GO / keyword links | `associated_with` to `OntologyClass` concepts | Explicit non-causal annotations. All source vocabulary concepts share the ontology class entity type for browsing. This deliberately uses Biolink’s mixin as a primary type. The slim TSV lacks GO branch metadata; more specific biological relationships would require an ontology join and structured source annotations. |
| UniProt keyword hierarchy | `subclass_of` between keyword concepts | Preserves source `is_a`. Keyword category tagging is `has_topic`, not a biological hierarchy assertion. Keyword identifiers use `Namespace.UNIPROT_KEYWORD`. |
| UniProt EC references | `has_topic` with EC CURIE | A classification annotation; do not invent a catalytic reaction from the EC accession. |
| UniProt narrative fields, length and mass | `description` values retaining source field headings | Avoid turning disease/localization prose into edges. Length is residue count and mass is the source's predicted sequence mass in Da. UniProt RDF places these properties on a Sequence node, so do not misapply them directly to Protein. A structured sequence/attribute representation can replace these descriptions if numeric queries are required. Sequence itself uses `has_biological_sequence`. |
| Source configuration | Existing resource/license/update enums | Operational metadata still uses `ResourceConfig`'s existing contract. These enums must not leak into biological records. Retire the package only after remaining inputs and configuration contracts are rewritten. |
| Serving `sign`, `is_directed` | Existing derived projections | Qualifier annotations already carry the authoritative aspect and effect direction. Sign projects effect direction; directedness concerns predicate symmetry. Removing the stored projections requires an explicit edge schema update and consumer changes. Do not interpret effect direction as endpoint orientation. |

Biolink's `original_predicate` is intended for mapping/test data, not production
knowledge graph assertions. Retain source predicate accessions in raw evidence instead.
External attributes currently use the annotation term/value storage representation;
this is not a claim that the complete graph instantiates every LinkML class directly.

## Verification before considering an input rewritten

- Exercise representative source rows through `SilverExtractor`, not only the builder.
- Check taxonomy/publication CURIEs, identifier domains, null fields and obsolete terms.
- Test what should *not* be asserted: unknown aspects, inferred locations, or chemical
  equivalence inferred from a broad mapping.
- Audit a full available source snapshot for unknown values and dropped records.
  Document exclusions and unresolved identifiers; do not report those as zero fallbacks.
- Rebuild immutable resource versions and publish a local release pinning those versions.

References: [Biolink modeling](https://biolink.github.io/biolink-model/),
[has_topic](https://biolink.github.io/biolink-model/has_topic/),
[UniProt vocabulary](https://www.uniprot.org/help/controlled_vocabulary),
[UniProt RDF schema](https://purl.uniprot.org/html/index-en.html).
