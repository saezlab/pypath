# Pathway, reaction and related inputs: native Biolink migration

The download and source parsing boundaries remain intact. Resource modules choose
native classes, predicates and annotation slots; shared builders perform validation.
Biological CV enums have been removed from these modules. Source configuration enums
remain operational metadata. The `kegg_metabolic` import remains an alias of KEGG.

Per the user’s explicit identifier convention, source names use `Namespace.NAME`
and synonyms/abbreviations use `Namespace.SYNONYM`, never annotation slots. Actual
accessions remain first; names and aliases supplement them. Source-qualified identifiers remain identifiers;
explicit unresolved names are retained for Recon3D inferred complexes and otherwise
unidentified WikiPathways participants. A classification accession is not automatically
an identifier of the protein or reaction it classifies. Taxonomy and publication
annotations contain `NCBITaxon:`, `PMID:` and `PMC:` CURIEs.

| Input | Representation and reviewed boundaries |
| --- | --- |
| BRENDA | UniProt proteins associated with explicit EC ontology concepts; EC hierarchy uses `subclass_of`, obsolete terms are excluded. No reaction is inferred from an EC classification. The existing parser's ligand and kinetic data remain in parsed source evidence; this module still exposes protein/EC annotations rather than a new ligand assay graph. |
| CORUM | `MacromolecularComplex`, protein members, GO concept associations, name identifiers, publications and explicit organism taxonomy. FunCat descriptions remain source data rather than being relabeled as an invented biological property. |
| IntAct | Native participant classes from exact PSI-MI interactor types. MITAB interaction observations use `interacts_with`; original MI interaction/method/role accessions are `has_topic` annotations, not arbitrary predicates. MI-score uses `intact_confidence_value`; stoichiometry remains typed. Unknown identifier systems are not relabeled as a known one. Negative/unavailable taxonomy values are omitted. |
| KEGG / alias | Molecular activities have `has_input`, `has_output`, `enabled_by` edges and typed stoichiometry. EC, KO and RCLASS are classification topics; KO/RCLASS are no longer reaction identifiers. Gene identifiers and orthology identifiers use distinct namespaces. KEGG's PubChem crosswalk contains **substance IDs**, not compound IDs. Pathways have reaction parts, molecular participants and broad links to linked pathways/orthology groups, rather than asserting every diagram item is a physical part. |
| MACdb | Direct noncausal chemical–source-trait `associated_with` relations. Trait IDs use the actual MACdb trait namespace; the EFO crosswalk is a topic, not assumed equivalence or subclassing. Contrast p-values use `p_value`; source study conclusions/conditions remain prose. Scalar case/control concentrations, bounds, confidence-interval columns, delta, log2FC and cohort sizes use native quantity values with the exact original source column retained as `source_field`. Comparators are retained; missing units are not invented. Non-scalar intervals and ambiguous text remain source evidence. In particular, confidence intervals are not mislabeled standard deviations. |
| PFOCR | Figure `InformationContentEntity` **mentions** `Gene` or `ChemicalEntity`. A figure is not an ontology term and co-mention is not a molecular interaction. Figure publications and gene taxonomy remain annotations. |
| Reactome | Source BioPAX types map explicitly to native classes. Reaction participants use explicit input/output predicates; unrecognized roles do not silently imply partonomy. Controls use `regulates` or explicit `catalyzes`; exact ACTIVATION/INHIBITION gives an effect direction, unknown controls do not. A degradation control retains its controlled process instead of synthesizing an effect on substrate quantity. Pathway hierarchy remains `part_of`. Parser cache version is 9 because the intermediate type tokens changed. |
| Recon3D | Native molecular activities with input/output/enzyme edges, stoichiometry, chemical formula and `chemrof:charge`. EC is a topic rather than reaction identity. Standalone model genes are `Gene`; GPR-derived enzyme members preserve the existing protein interpretation and inferred complex composition. Complexes use an explicitly unresolved name, not a fabricated accession. Transport and metabolic views use the same activity type for the same reaction ID. |
| Rhea | Master reactions are molecular activities with input/output/enzyme predicates; transport and metabolic views keep consistent entity identity. Equations are name identifiers, EC/cross-references are topics, GO and Reactome links remain broad associations. Conversion direction and compartment labels remain source evidence; they are not gene-effect direction qualifiers. Unexported legacy catalysis schemas with erroneous self-membership were removed. |
| SLC Tables | Taxon-scoped protein observations with gene identifiers and name/synonym identifiers. Substrates and tissue-expression prose do not become guessed chemical or anatomical entities. Family/transport descriptors remain source evidence; no structured values are dumped into `description`. |
| TCDB | UniProt proteins carry TC classification topics. Transporter–substrate observations are explicit broad interactions with ChEBI chemicals. A TC classification is neither a protein accession nor a unique transport-event identifier. Rows without a transporter accession cannot supply a protein endpoint. |
| WikiPathways | Explicit RDF participant types and exact interaction classes. Stimulation/inhibition supply direction only when unambiguous; conversion/catalysis diagram edges remain broad interactions until a reaction-node representation is introduced. Missing taxonomy is genuinely absent. Version identifiers remain distinct from names. Label-only participants remain explicitly unresolved name identifiers. |

`description` is reserved for real source prose. Structured measurements and source
roles are not flattened into labeled descriptions. Typed stoichiometry, chemical
formula/charge, classification topics, taxonomy, publications, IntAct confidence and
MACdb p-values and source-defined scalar quantities remain queryable attributes. Quantities use the shared optional typed `quantity` projection, including source field and comparator; legacy scalar annotations remain compatible. More elaborate assay and compartment
representations are a remaining boundary, not silently claimed complete support.

## Verification

Twenty-seven representative source-shaped and negative tests exercise each resource through `SilverExtractor`,
including negative tests for unknown roles/effect directions, missing taxonomy, obsolete
EC terms, identifier domains, figure co-mentions, and structured-versus-narrative fields.
A minimal BioPAX graph exercises the changed Reactome parser as well. These tests do
not use source downloads.

Read-only audits on the locally available source files produced the following counts
before the follow-up change classifying all names/synonyms as identifiers. Name-only
rows can now be retained as explicitly unresolved identifiers, so the historical empty-row
counts below are not claimed to be current build totals.
The audit constructs fresh extractors per row: counts measure mapping/semantic
validation, **not** globally consolidated entity/edge counts or a rebuilt release.

| Source dataset | Parsed rows | Mapped rows | Empty rows | Exceptions |
| --- | ---: | ---: | ---: | ---: |
| CORUM human complexes | 2,916 | 2,916 | 0 | 0 |
| SLC Tables | 464 | 464 | 0 | 0 |
| TCDB transporters | 24,627 | 24,616 | 11 | 0 |
| TCDB substrates | 20,050 | 20,039 | 11 | 0 |
| MACdb associations | 40,710 | 40,710 | 0 | 0 |
| PFOCR genes | 477,779 | 477,779 | 0 | 0 |
| PFOCR chemicals | 141,181 | 141,181 | 0 | 0 |
| Recon3D reactions | 10,600 | 10,600 | 0 | 0 |
| Recon3D genes | 1,883 | 1,883 | 0 | 0 |
| Recon3D inferred complexes | 182 | 182 | 0 | 0 |
| Rhea reactions | 18,343 | 18,343 | 0 | 0 |
| KEGG reactions | 12,389 | 12,389 | 0 | 0 |
| IntAct first 10,000 rows (capped) | 10,000 | 9,979 | 21 | 0 |
| WikiPathways pathways | 2,079 | 2,058 | 21 | 0 |
| WikiPathways interactions | 40,350 | 39,839 | 511 | 0 |
| Reactome human reactions | 14,566 | 14,566 | 0 | 0 |
| Reactome human pathways | 2,825 | 2,825 | 0 | 0 |
| Reactome human controls | 9,758 | 9,758 | 0 | 0 |
| Reactome human controller groups | 1,112 | 1,112 | 0 | 0 |

TCDB's 11 excluded records have blank UniProt accessions. WikiPathways has missing
pathway IDs or missing both identifiers and endpoint labels in the excluded rows.
IntAct is capped because its cached archive is approximately 1 GB; its excluded rows
lack a usable mapped endpoint. These are explicit exclusions, not successful fallbacks.
The full human Reactome RDF was reparsed (3,163,818 triples), bypassing old prepared
rows, before validating all four datasets. KEGG pathway/KGML outputs have fixture
coverage; the full KEGG audit above covers its reaction export.

The standalone Recon3D metabolite parser yields zero rows for this snapshot because
all metabolites are already emitted through reaction participants; its chemical fields
are exercised there. BRENDA's available `brenda.txt.tar.gz` is truncated (EOF before
the gzip end marker), preventing a complete source audit. BRENDA has fixture coverage,
including obsolete EC terms, but is not represented as a successful full-source run.

No source downloads, prepared caches, profiling artifacts or new resource versions were
written by these audits. Publishing rebuilt immutable resources remains a separate step.
