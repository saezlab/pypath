# Modeling inputs with Biolink

Use the installed, pinned `biolink-model` package (currently 4.4.4) as the biological
schema. `signor.py`, `chebi.py`, and `uniprot.py` are the rewritten examples.
Inputs are modeled individually; do not add a migration adapter.

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
4. Source names, synonyms and abbreviations belong in `IdentifiersBuilder`, using
   `Namespace.NAME` and `Namespace.SYNONYM`. This is the application convention
   requested for inputs_v2; keep real accession namespaces ahead of name aliases.
   Do not duplicate these values in name/synonym annotations. Display labels are
   selected from name identifiers.
   Use slots for attributes: `description`, `in_taxon`,
   `publications`, `has_biological_sequence`, `has_chemical_formula`, etc.
   Narrative descriptions remain text; they are not entity identifiers or inferred edges.
5. Use `Namespace` for identifier namespaces. It is an application identifier registry,
   not an extension of Biolink's biological types. Add a member only for an actual
   identifier system. Translation datasets use the same members as graph inputs.
6. Use `description` only for source prose, such as definitions and comments.
   Structured measurements, roles, categories and flags need an appropriate native
   slot or a reviewed published external ontology attribute CURIE. If neither fits,
   first consider a Biolink quantitative attribute with units and source-field
   context (see README); otherwise preserve them in parsed evidence payloads and document that they are not exposed
   as serving annotations. Do not turn them into labeled description strings or mint
   an OmniPath term merely because an upstream label looks convenient.
7. Ontology loaders (`ontology_term_to_entity`) emit only current terms: an
   obsolete term produces no entity and no hierarchy edges. Earlier CV-based
   inputs emitted obsolete terms with an `is_obsolete` flag; that is intentionally
   not carried over. Consequence: an annotation in another resource that still
   cites an obsolete GO/HPO/MONDO/EC/ChemOnt term keeps its CURIE but has no
   ontology name, parents or replacement link. `replaced_by`/`consider` hints are
   not used to redirect such references; doing so would need a reviewed mapping.

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
| UniProt narrative fields | `description` values retaining source field headings | Disease/localization/function prose remains prose. Numeric length/mass and protein-family classification stay in parsed payloads until a structured sequence/attribute model is provided; UniProt RDF places length/mass on Sequence rather than Protein. Sequence uses `has_biological_sequence`. |
| Source configuration | Existing resource/license/update enums | Operational metadata still uses `ResourceConfig`'s existing contract. These enums must not leak into biological records. Retire the package only after remaining inputs and configuration contracts are rewritten. |
| Serving `sign`, `is_directed` | Existing derived projections | Qualifier annotations already carry the authoritative aspect and effect direction. Sign projects effect direction; directedness concerns predicate symmetry. Removing the stored projections requires an explicit edge schema update and consumer changes. Do not interpret effect direction as endpoint orientation. |

Biolink's `original_predicate` is intended for mapping/test data, not production
knowledge graph assertions. Retain source predicate accessions in raw evidence instead.
External attributes currently use the annotation term/value storage representation;
this is not a claim that the complete graph instantiates every LinkML class directly.

## Annotation resources

Annotation inputs (classifications, localizations, markers, gene sets, disease and drug
annotations of genes and proteins) are modelled so that a term can be searched as an
entity and the molecules annotated with it are shown. The pattern is that of `uniprot.py`,
`hpo.py`, `mondo.py` and `macdb.py`. Examples: `cspa.py`, `tfcensus.py`, `kinasedotcom.py`,
`matrisome.py`.

**Statements.** The subject is a `Protein` when the source identifies proteins (UniProt,
PDB, entry name) and a `Gene` when it identifies genes (symbol, Entrez, Ensembl gene,
HGNC, MGI); when a row has both, `Protein` if there is a UniProt accession. Source
identifiers are kept; nothing is translated and no species is guessed. The object is an
`OntologyClass` entity. A statement with nothing attached is an `Association` on the
molecule's entity; one with a class, a score, a publication or a context is a `Relation`.
One statement per source statement.

**Term identity.** A source accession of a real vocabulary is written in that vocabulary's
`Namespace` exactly as its ontology module writes it (`GO:0005886`, `CL:0000312`,
`EC:3.4.21.4`, `IPR000001`), so the term is the same entity as the ontology's own. A term
the source gives only by name is not looked up in an ontology. A source-owned vocabulary
gets the namespace `<resource>_<field>` (one per field), the source's own value as the
identifier (label, code, number or path) and `Namespace.NAME` for the label. Levels of
one hierarchy share a namespace only where the ids are unique across levels
(`almen2009_class`, `gpcrdb_family`); otherwise there is one per level
(`kinasedotcom_group`, `_family`, `_subfamily`). A term is never identified by
`Namespace.NAME` alone.

**Terms datasets.** A source-owned vocabulary is also emitted by a dataset of its own,
called `terms`, with the name, the definition (`description`) where the source gives one,
and `subclass_of` for the levels the source states. The object of a statement and the
terms entity must have the same type, namespace and identifier; without the terms dataset
the term would be a code with no name. The dataset name has no meaning for the build; it
labels the evidence and selects with `--datasets`.

**Plain lists.** A resource that is a set of molecules with one property (a list of
transcription factors) is one term in the namespace `<resource>_list`, with a slug as the
identifier, and one statement per listed molecule. The list's own columns go on the
statement.

**Predicates.** `located_in` for the subcellular location or membrane of a molecule;
`expressed_in` for a marker or detection in a cell type or tissue; `associated_with` for
a class, role, list, gene set or disease. `associated_with` asserts nothing causal.
Domain and range are checked in the pinned model.

**Families.** A source organised as families with member lists (HGNC groups, kinase.com,
GPCRdb, the phosphatome tree, DNA-binding domain families, integrins) is emitted as
`GeneFamily` or `ProteinFamily` with `has_member`, nested as the source nests it. A source
that is a table with a class column per molecule (Almen, OPM, matrisome, LOCATE) is terms
with `subclass_of` and `associated_with`.

**Attributes.** Source fields on a statement are `<resource>:<field>`, lowercase, with
text as written. Every annotation has a value: the build leaves out a term without one.
Numbers the source states as measurements are a `Measurement` with the scale in a
comment. Yes/no flags are `'true'`/`'false'`. A confidence class is
`has_confidence_level`; a probability is `has_confidence_score`. A context that is an
ontology accession is an attribute value; a context that is a text label is kept as text.

**Species and provenance.** `in_taxon` is on the molecule, stated by the file where it
states one and a module constant where the resource covers one species (human), with a
comment. PubMed identifiers go through `pmid_curie`. Source databases named by the file go
to `supporting_data_source` under the source's own labels.

**Not ported.** OmniPath's own derived classifications (`intercell_annot`) and consensus
fields; rows the source marks as not members (TFcensus class x, Lambert "No", kinase
pseudogenes); collections the licence excludes (MSigDB BioCarta, KEGG legacy); and any
term-to-ontology matching by label.

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
