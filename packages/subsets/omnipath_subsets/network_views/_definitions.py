"""The datasets onboarded into the framework: MetalinksDB + LIANA.

Both are now **presets**: metadata over the interaction fact table — the
contributing sources, the interaction-class scope, the evidence scope, the
default and mandatory attributes, the labels, the curation thresholds, the
attribute sources, the mode the per-resource record folds under, the license a
resource must meet, and — where the dataset is not one query — the recipe that
assembles it. Registering the row is the whole build step, and it materialises
nothing at all.

The framework still carries the older **matview network** shape: a schema, a
combined relation and the curated SQL that fills them. No definition here uses
it any more, and the ``schema`` / ``combined_relation`` fields retire with the
registry columns behind them once the standing views are dropped.

**A database built before the conversion still holds both datasets' old
views, unmanaged.** Converting a dataset stops the framework *managing* its
views and does not drop them: ``apply_network`` returns early for a preset. So
they survive, frozen at their last refresh, serving data the build has since
corrected — they are not a live contract, and a fresh ``init-db
--drop-existing`` takes them. The SQL that built the MetalinksDB views is
deleted: it read annotation terms of the CV-based inputs that no build writes.
The copy under ``legacy/postgres`` stays as the parity reference.

Adding a dataset stays a declarative change: a definition here, served by the
same uniform API.
"""

from __future__ import annotations

from omnipath_subsets.network_views._framework import NetworkDefinition

# MetalinksDB: compound↔protein relations, a preset over the interaction
# fact table. It was fifteen materialized views and 1,571 lines of SQL; the
# recipe below is the whole of what that SQL said, expressed as parameters.
#
# **The resource names are the loaded ones, not the published ones.** Human-GEM
# loads under `metatlas` and the retiring views only read `humangem` because
# they label their own rows so by hand. A preset naming `humangem` would
# resolve it to nothing, and an empty contribution reads exactly like a
# resource that failed to load. `labels` carries the published name for output.
#
# **One of the thirteen contributes nothing, and it is worth stating so it does
# not read as a defect.** BindingDB is excluded by the recipe: its assertions
# stay in the record for any other query to find, and no row of this dataset
# counts them. ChEMBL's mechanism-of-action pairs do contribute, and no
# retiring view delivered them — the gate below says why.
#
# **The reaction-grain component delivers the transport half and not the rest.**
# Rhea, Recon3D and Human-GEM publish reactions as the source states them: a
# molecular activity with its participants and catalysts. The build projects the
# transport each states (``omnipath_build.projections``): a cargo a reaction
# takes in one compartment and gives in another, with each catalyst the same
# source row names, becomes an ordinary transport relation of that resource —
# catalyst ``affects`` cargo, aspect transport — so it reaches this component
# like any directly stated transporter-substrate pair. A conversion that moves
# nothing still has no pair reading: its enzyme reaches its substrates only
# through the reaction.
_CURATED_SOURCES = (
    "cellinker",
    "guidetopharma",
    "mrclinksdb",
    "stitch",
    "tcdb",
    "cellphonedb",
    "neuronchat",
    # Not in the original MetalinksDB: requested 2026-10-25. Only the pairs
    # MetaChatDB itself reports load, so it adds no copy of MRClinksDB's.
    "metachatdb",
)

# The classes the transport component scopes to. Named rather than left open
# because the reaction-derived contribution is a transport statement, and
# widening it would pull in the signaling rows of the same resources.
_TRANSPORT_CLASSES = ("transport", "ligand_receptor")

# The gate every component carries. It is on the entity, so no resource can
# route around it: one end of every row is a small molecule.
#
# **It gated on chemical class until this cycle, and now gates on entity type.**
# The class comes from `classify/chemical_class.yaml`, which reads a compound's
# class off the resources that contributed it and calls anything unmapped a
# metabolite. That is provenance standing in for chemistry, and it decides
# both ways. Gadolinium and olanzapine entered as metabolites because no rule
# spoke for their resources. BindingDB's catalog entered because one rule did.
# A dataset gated on that inherits both. The entity type asserts only what the
# record holds — this end is a small molecule — and leaves metabolite-or-not
# to the caller, who can still pass `chemical_classes` on any request.
#
# **The other end is a molecule, and the statement is an interaction.** One
# chemical end alone also admits a resource's annotations of its chemicals:
# CellPhoneDB's "secreted" role terms, MRClinksDB's ClassyFire classes (ChemOnt
# terms), Guide to Pharmacology's antibacterial role (``has_chemical_role`` to a
# ChEBI role, which the ChEBI load types as a chemical). So the partner is one
# of the molecular types — untyped targets (``named_thing``) and chemicals
# included, ontology terms and organisms not — and the predicate is one of the
# three the dataset's pairs are stated with.
_PARTNER_TYPES = [
    "protein",
    "gene",
    "polypeptide",
    "macromolecular_complex",
    "protein_family",
    "gene_family",
    "nucleic_acid_entity",
    "rna_product",
    "molecular_entity",
    "molecular_activity",
    "named_thing",
    "chemical_entity",
    "small_molecule",
]
_GATE = {
    "entity_types": ["chemical_entity", "small_molecule"],
    "partner_types": _PARTNER_TYPES,
    "predicates": ["interacts_with", "affects", "associated_with"],
}

# **Three filters the retiring views applied and this recipe does not**, the
# first two measured on the dev4 build 2026-08-25. None loses a row the views
# delivered: their output stays a subset of this preset's, pair for pair. But a
# consumer moving across meets a result several times larger, so this states
# the differences rather than leaving them to be met as a surprise.
#
# **Organism.** Every per-source view was human-only, filtering the protein
# mention on taxon 9606. This preset has none, for the reason the ligand-receptor
# drop has none: organism is a dimension a caller queries on, not a property of
# a dataset, and a human-scoped preset cannot be widened later without
# registering a second one. The effect is large and uneven — of STITCH's 202,669
# metabolite pairs, 78,495 are human and 93,636 are mouse — so **a consumer
# expecting the delivered human-only contract should pass an organism filter**.
#
# **Resolution status.** The combined view kept a pair only where both ends
# resolved cleanly, and the query surface has no parameter for that, so this
# recipe cannot express it. It costs little on most resources and a great deal
# on two: of Cellinker's 4,653 metabolite pairs only 41 have both ends resolved,
# and of MRC-LinkDB's 1,447 only 33. Those are pairs whose endpoints never
# reached a canonical identifier, and serving them under a curated dataset's
# name is the weaker half of this conversion. Tracked as a follow-up.
#
# **Chemical class.** Every per-source view kept a pair only where the compound
# classified as a metabolite, and so did this recipe. The gate above replaced
# it. The change admits ChEMBL, whose mechanism-of-action compounds classify
# as drugs. It also admits any compound of the other resources that a rule or
# the default classified away from metabolite.

METALINKSDB = NetworkDefinition(
    name="metalinksdb",
    kind="compound_protein",
    included_sources=(
        "chembl",
        "bindingdb",
        *_CURATED_SOURCES,
        "recon3d",
        "rhea",
        "metatlas",
    ),
    interaction_class_scope=(),
    default_attributes=("endpoints", "label", "evidence"),
    # The node classification the delivered contract carries inline. Mandatory
    # rather than default: a consumer of this dataset reads the compound and
    # protein annotations, and a request naming another attribute must not
    # take them away.
    #
    # **The name is `annotation`, not `intercell`.** The layer serves GO
    # cellular components, UniProt cellular-component keywords, UBERON tissues
    # and five chemical class vocabularies; `intercell` names a role
    # classification none of them states, and FR-046 forbids a changed meaning
    # behind an unchanged name. A build whose registry still says `intercell`
    # is still served — the API accepts the old spelling as a deprecated alias
    # and logs it — but it emits only the new column names, so re-registering
    # this product is what makes the registry agree with the response.
    mandatory_attributes=("annotation",),
    labels={
        "preset": "MetaLinksDB",
        "resources": {"metatlas": "humangem"},
    },
    curation={
        # Every threshold the retiring SQL held inline, as configuration.
        "entity_type_gate": "chemical_entity",
        "chembl_curation": "mechanism_of_action",
        "excluded_from_combined": ["bindingdb"],
    },
    attribute_sources={
        # Provenance for the node classification, and it says which stage
        # answered. The interim stage is the annotation content the resources
        # already publish; the Intercell rebuild replaces the vocabulary behind
        # the same output names, and `stage` is how a caller tells them apart.
        # It stays `interim` so the Phase 9 swap is a change to this one field.
        #
        # The source is named per layer rather than once for the lot. The old
        # value, "loaded resource annotations", described the
        # `entity_annotation_relation` store the API no longer reads and whose
        # five terms returned no rows at all — so it named the wrong store as
        # well as being too coarse to act on.
        "annotation": {
            "stage": "interim",
            "source": "omnipath.relation, per layer",
            "layers": {
                "location": (
                    "located_in facts over GO cellular components from LOCATE, "
                    "HMDB and iMM1415, over UBERON tissues and CL cell types, "
                    "beside the 71 UniProt cellular-component and topology "
                    "keywords as associated_with facts"
                ),
                "chemical_class": (
                    "associated_with facts over LIPID MAPS, the three RefMet "
                    "levels and ChemOnt, ranked in that order and rolled up "
                    "through the subclass_of tree in entity_ontology_relation"
                ),
                "metabolic_domain": (
                    "entity.metabolic_domain_id, a derived column of the build "
                    "rather than an asserted fact, so no resource is named"
                ),
            },
            "note": (
                "location and class terms as the contributing resources "
                "publish them, not a cross-resource consensus; the role "
                "classification the Intercell rebuild owns is not here, "
                "because the model states a role per assertion and an "
                "entity-level flag needs a consensus this stage has none of"
            ),
        },
    },
    composition={
        "operation": "union",
        "components": [
            {
                "parameters": {
                    "filters": {
                        "resources": ["chembl"],
                        # Presence of the mechanism annotation, folded onto the
                        # record as a flag. An affinity threshold is the
                        # alternative and it floods the set: pChEMBL above 6 is
                        # 1.6 M pairs, and mechanism and pChEMBL sit in
                        # different ChEMBL tables, so their conjunction is empty.
                        "curation_flags": ["mechanism_of_action"],
                        **_GATE,
                    },
                },
            },
            {
                "parameters": {
                    "filters": {
                        "resources": list(_CURATED_SOURCES),
                        **_GATE,
                    },
                },
            },
            {
                "parameters": {
                    "filters": {
                        "resources": ["recon3d", "rhea", "metatlas"],
                        "interaction_classes": list(_TRANSPORT_CLASSES),
                        **_GATE,
                    },
                },
            },
        ],
        "steps": [
            # Before the fold, so the dropped resource contributes no row and
            # no count. After it, its assertions would stay inside
            # source_count, the references and the sign flags.
            {"operation": "exclude", "resources": ["bindingdb"]},
            {"operation": "collapse"},
            {"operation": "annotate", "layer": "annotation"},
        ],
    },
)

# LIANA: ligand↔receptor pairs, a preset over the interaction fact table. It was
# a matview network over 5 cell-cell-communication resources; the reduced drop
# scopes to ConnectomeDB2025 alone, so the other four are out of scope here.
# Nothing materialises: registering this row is the whole build step. The SQL
# that built the old matview is deleted, so no build recreates it, and a
# database still holding one holds a relation nothing owns — frozen at its last
# refresh, its `sources` column reading a resource name the build has since
# corrected. Dropping that leftover is a one-line statement against the
# database and not a code change, which is why it is not here.
#
# Organism scope: all 14 taxa the resource loads, not human alone. That is a
# decision, not an oversight. Organism is a dimension a caller queries on, not a
# property of a preset: `subject_organism` and `object_organism` are populated on
# every one of the 44,455 rows and never disagree, so asking for organism 9606
# yields the human drop at no cost, whereas a human-scoped preset could not be
# widened later without registering a second one. Scoping to human here would
# throw away 92 per cent of a resource the build already holds. The consequence
# is worth stating plainly: **this preset's default result is all taxa**, so a
# consumer expecting the legacy human-only LIANA meets a result roughly twelve
# times larger, and should pass an organism filter to get the old shape.
#
# **`references` has left the default attributes, and so has it in the two
# definitions above.** FR-055 moves the four rebuildable columns off the
# default row and keeps them selectable by name, and a default attribute is
# exactly what would put `references` back on it: the API echoes a preset's
# declared attributes into the row, so declaring the column is the same thing
# as serving it unasked. A caller who wants it writes `attributes=references`
# and gets it, on this preset as on any other.
#
# For this preset it returns empty either way, and that is worth keeping on the
# record rather than in the defaults. ConnectomeDB2025 does publish PubMed ids,
# but the loader that produced this build
# (`pypath/inputs_v2/connectomedb.py`) never reads the `AI summary` column they
# sit in, so there is no reference annotation for the projection to collect and
# nothing to serve. That is an ingest gap, tracked as a follow-up in this
# cycle's task list, not a property of the dataset, and it fills in once the
# loader is fixed and the resource re-ingested. Nothing here claims otherwise:
# `attribute_sources` names no source for references.
LIANA = NetworkDefinition(
    name="liana",
    kind="ligand_receptor",
    included_sources=("connectomedb2025",),
    interaction_class_scope=("ligand_receptor",),
    default_attributes=("endpoints", "label", "evidence"),
)

# Reactions: metabolic reactions as native hyperedges, a preset over the same
# interaction fact table but at participant grain. One group is one reaction
# header whatever its arity, and its members come back as a participant list
# rather than as a first and a second endpoint. Nothing materialises: as with
# the two above, registering the row is the whole build step.
#
# **The grain identifies this dataset, not a class or a resource set.** Every
# other preset here says what it is by naming an interaction class or the
# resources that feed it. This one says it by how its groups are keyed, which
# is why `interaction_class_scope` is empty. That emptiness is a statement and
# not an omission: `vocab_interaction_class` holds no reaction slug, the
# `has_participant` predicate maps to `other`, and scoping to `other` would
# admit every unclassified pair in the build rather than narrowing anything. If
# the classification later grows a term that fits, the scope narrows to it
# without changing what the dataset means. `collapse_mode` records `none` for
# the reason the field's own documentation gives — at participant grain nothing
# reads it, because each mode folds an ordered endpoint pair and a reaction is
# not one, so `none` states no fold rather than the wrong one.
#
# **It ships four of the six resources the specification names.** Metabolic
# Atlas Human-GEM — loaded as `metatlas`, see the note above on loaded versus
# published names — RECON3D, Rhea and KEGG are here. Mouse-GEM and iMM1415 are
# not, and they are deliberately not carried as unloaded names either: the
# build has no source slug for either one, Human-GEM loading under `metatlas`
# and RECON3D under its own, so naming them would invent a loader rather than
# record a resource that failed to load. The dataset is short of whatever those
# two would have contributed, and that shortfall is an ingest gap tracked in
# this cycle's task list, not a property of the dataset.
#
# **Reactome is loaded, publishes reactions, and is deliberately out of scope.**
# This is a separate judgement from the COSMOS projection's, which excludes it
# too, and it rests on its own measurement against dev3. The reason is what its
# participants are. Across the four resources above, every one of the 722,435
# participant edges resolves to a `Chemical`. Reactome's resolve to 34.5 per
# cent `Chemical` and 65.5 per cent `Complex`, `Gene`, `Protein Family`,
# `Physical Entity`, `DNA` and `RNA`; 8,238 of its 14,778 reaction events —
# 55.8 per cent — contain no chemical participant at all. Those are not thinner
# descriptions of the same object, they are a different one, and merging them
# in would make the dataset's name false. It would also be an addition rather
# than a merge: Reactome shares 532 reaction entities with Rhea and 7 with
# Human-GEM, none with KEGG or RECON3D, so 96 per cent of it arrives as new
# events.
#
# **The thin annotation is not the reason, and saying so matters.** Reactome
# carries stoichiometry on 8.1 per cent of its participant edges and
# subcellular location on none, which looks like the obvious disqualification
# and cannot be one: Rhea publishes no stoichiometry on any of its 172,640
# edges, and KEGG publishes no compartment on any of its 365,481. A rule that
# turned either mandatory attribute into an admission bar would drop two of the
# four resources this dataset is built on. Coverage is a per-resource fact that
# `attribute_sources` below records; it decides nothing about scope.
#
# **Mandatory here means present, not known.** A mandatory attribute forces a
# name into the output projection and reaches no `WHERE` clause, so declaring
# role, side, stoichiometry and compartment mandatory guarantees the response
# shape and takes no row away — a participant whose stoichiometry nobody
# published is served as null. A caller who wants only the fully annotated
# reactions asks for them through `attribute_filters`, which is a request
# parameter and a different field.
REACTIONS = NetworkDefinition(
    name="reactions",
    kind="reaction",
    included_sources=("kegg", "rhea", "metatlas", "recon3d"),
    interaction_class_scope=(),
    default_attributes=("endpoints", "label", "evidence"),
    mandatory_attributes=("role", "side", "stoichiometry", "compartment"),
    grain="participant",
    collapse_mode="none",
    labels={
        "preset": "Reactions",
        "resources": {"metatlas": "humangem"},
    },
    curation={
        # The grain restated where a consumer reading configuration rather than
        # prose will find it, and the scope shortfall in machine-readable form.
        "participant_grain": "hyperedge",
        "excluded_from_scope": ["reactome"],
        "never_loaded": ["mouse-gem", "imm1415"],
    },
    # Every entry is a resource publishing the field itself, so no stage is
    # named: there is no interim vocabulary standing in for a later rebuild the
    # way the node annotation layer has one.
    attribute_sources={
        "role": {
            "sources": ["kegg", "rhea", "metatlas", "recon3d"],
            "note": "every participant of every reaction carries one",
        },
        "side": {
            "sources": ["kegg", "rhea", "metatlas", "recon3d"],
            "note": "derived from the reactant or product role",
        },
        "stoichiometry": {
            "sources": ["kegg", "metatlas", "recon3d"],
            "note": (
                "complete on every reaction these three publish; rhea states "
                "none on any of its reactions, and serves null"
            ),
        },
        "compartment": {
            "sources": ["metatlas", "recon3d"],
            "note": (
                "per-participant subcellular location, complete on both; kegg "
                "states none, and rhea states membrane sides that sit on "
                "transport records rather than on reactions"
            ),
        },
    },
)

NETWORKS: list[NetworkDefinition] = [METALINKSDB, LIANA, REACTIONS]
