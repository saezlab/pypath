"""A small hand-built canonical graph for the interaction-projection tests.

The interaction derive step projects ``relation`` and its
``relation_evidence`` provenance into ``interaction``, ``interaction_party`` and
``interaction_fact_resource``. Asserting its semantics against the full build
would make the assertions depend on whatever the resources happen to say this
week, so the projection tests build this fixture instead: a dozen relations in
a throwaway schema, each one carrying exactly the situation a requirement talks
about.

What the graph is built to exercise:

* ``ligand_receptor`` from **participant-role** evidence (tier 1),
* ``allosteric`` and ``orthosteric`` from **interaction-level** annotation
  (tier 2), ``signaling``/``transport`` from the **predicate** (tier 3) and
  ``other`` as the fallback,
* a pair whose resources **disagree** on sign (both flags true, cross-resource),
* a pair where **one** resource asserts both signs (single-resource conflict),
* a contributor asserting **neither** sign nor direction, so
  ``sign_source_count <= cardinality(sources)`` is a real inequality,
* a pair with **no sign at all**, so the sign columns stay NULL, and
* an **opposite-direction pair**, which must stay two rows.

Beside the ordered endpoint pairs the graph carries three **stars**, the shape a
metabolic resource publishes a reaction in: a parent entity and one
``has_participant`` relation per member, with the member's role, stoichiometry
and compartment riding on the relation evidence at **object** scope, and the
catalyst arriving as a sibling ``controls`` edge that points *at* the parent.

* two reaction parents with the **same** member multiset, reported by two
  different resources, so the hyperedge projection is asserted to merge them
  into one header rather than to mint one per resource,
* a transport parent whose members carry a membrane side rather than a
  subcellular location, and
* a **pathway** parent, which uses the same ``has_participant`` verb and must
  **not** become a hyperedge — the reaction projection keys on the parent's
  entity type, not on the predicate, and a pathway with two hundred members
  would otherwise mint a two-hundred-way interaction, and
* a reaction parent whose **whole** participant set is the two entities an
  ordinary pair already joins, which is the one shape where the two readings
  hash to the same header id.

Three further reaction stars are there for the **binary metabolic projection**,
which reads a reaction party row by party row and turns it into substrate and
product edges around a per-reaction enzyme node:

* an ordinary catalysed conversion, two substrates and one product, one of the
  substrates stated **without** a compartment,
* a conversion **no resource names an enzyme for**, which has to keep its edges
  rather than be dropped for want of a catalyst, and
* a conversion carrying a **cofactor** and a **regulator** beside its
  substrates, neither of which stands on a side of the arrow.

Five more carry a **direction**, which a resource states on the reaction event
as an entity annotation rather than on any one membership: two reversible
conversions, two one-way ones and a reversible one nobody names an enzyme for.
Each direction appears under both of the spellings the build holds, so a
projection that folded only one of them is caught rather than assumed correct.
Three of the reactions above carry **no** direction at all, which is what a
resource publishing only the chemistry looks like.

Two more stars hold the same metabolite on **both** sides of the arrow, which
is the shape a binary transport is read out of. One carries it across a
membrane — a reactant in the cytosol and a product outside — and the other
consumes and regenerates it in one compartment. Both statements ride on one
membership, because ``relation`` is unique on its endpoint triple and the role
lives on the evidence row rather than on the relation, so the two are told
apart by the evidence they sit on. What separates them is the compartment and
nothing else, so a projection reading the roles alone publishes a transporter
for the reaction that transports nothing.

One more star is **wide**: seventy members and their enzyme, with entity ids
that are hashes rather than the indexed ones the rest of the graph uses. The
hyperedge projection groups on a text multiset of one element per participant,
and a value that long is one a btree index row cannot hold — which is a
`CREATE INDEX` that fails rather than an assertion that misses. The ids have to
be hashes for it: a btree tuple holds the value compressed, and seventy ids
differing in a digit shrink back under the ceiling.

Beside those single-situation rows the graph carries a **coverage pair per
interaction class**: one ordered endpoint pair for every class the graph can
evidence, each reported by two resources that both publish a reference. Those
are what let the dedup and provenance rule be asserted class by class instead
of on ``signaling`` alone, and they are the reason a class reached only by a
verb no rule maps cannot pass unnoticed. One endpoint pair carries **two**
classes at once, which turns the class into a claim about rows rather than
about a column.

Every id is fixed, so a test can name the row it means.
"""

from __future__ import annotations

import uuid

from psycopg2 import sql
import psycopg2.extensions

# Sources. The ids are explicit so the partitioned evidence tables land in their
# default partition and the test can name a source by number.
SOURCE_A = 9001
SOURCE_B = 9002
SOURCE_C = 9003
SOURCE_LR = 9004
SOURCE_NAMES = {
    SOURCE_A: 'fixture_res_a',
    SOURCE_B: 'fixture_res_b',
    SOURCE_C: 'fixture_res_c',
    SOURCE_LR: 'fixture_res_lr',
}


# The predicate miRBase publishes its pre-miRNA to mature-miRNA processing
# under. The build stores the bare accession, because the controlled-vocabulary
# term is an (accession, label) pair and the loader writes the accession.
MATURATION_PREDICATE = 'OM:1257'

# One ordered endpoint pair per interaction class the graph can evidence:
# the class slug, the predicate the pair is reported under, and the
# interaction-level annotation that names the class where the verb cannot.
# `tf_target` is absent on purpose — no resource in the graph asserts it, and
# a fixture that invented one would hide the gap instead of surfacing it.
CLASS_PAIRS = (
    ('ligand_receptor', 'interacts_with', ()),
    ('signaling', 'controls', ()),
    ('transport', 'transports', ()),
    ('orthosteric', 'interacts_with', ('Agonist:OM:1001',)),
    ('allosteric', 'interacts_with', ('Allosteric Modulator:OM:1005',)),
    ('maturation', MATURATION_PREDICATE, ()),
    ('other', 'has_member', ()),
)

# The classes the coverage pairs above evidence, in the order they are built.
COVERED_CLASSES = tuple(slug for slug, _predicate, _terms in CLASS_PAIRS)

# The participants. The single-situation rows name theirs by the letter the
# docstrings above use. The coverage pairs name theirs after the class they
# carry, so a failure names the class that broke. Only hex digits are
# legal in a uuid, so the name indexes the entity and an ordinal carries it.
# The entity types the star rows need. Everything else in the graph is a
# protein, which is what the pair rows have always been.
PROTEIN_TYPE = 'protein:MI:0326'
CHEMICAL_TYPE = 'small molecule:MI:0328'
REACTION_TYPE = 'Reaction:OM:0015'
TRANSPORT_TYPE = 'Transport:OM:0035'
PATHWAY_TYPE = 'Pathway:OM:0014'

# The star participants. `rxn_a` and `rxn_b` are two **distinct** reaction
# entities over the same members — what two resources produce whenever the
# load side's reaction hash cannot merge them, which it cannot for a reaction
# that names only one side. The projection has to merge them anyway.
STAR_NAMES = (
    'rxn_a',
    'rxn_b',
    'trn_a',
    'pathway_a',
    'met_x',
    'met_y',
    'met_z',
    'met_in',
    'met_out',
    'enz_p',
    'enz_t',
    'pw_m1',
    'pw_m2',
    'pw_m3',
    'rxn_c',
)

# The reactions the binary metabolic projection is asserted on. They are
# separate from the stars above because that projection reads a reaction one
# party row at a time and each of these carries one situation it has to answer
# for: the ordinary catalysed conversion, a conversion nobody names an enzyme
# for, and a conversion whose party list holds a cofactor and a regulator
# beside its substrates. Every metabolite here is its own entity, so a test can
# scope an assertion to one reaction by naming the row it means.
COSMOS_NAMES = (
    # The ordinary case: two substrates, one product, one enzyme. `cos_met_b`
    # carries no compartment, which is the half of the label rule the other two
    # metabolites cannot exercise.
    'cos_rxn',
    'cos_met_a',
    'cos_met_b',
    'cos_met_c',
    'cos_enz',
    # A conversion with no catalyst. The resource states the chemistry and says
    # nothing about who runs it, which is the majority shape in a genome-scale
    # model and the one a projection keyed on the enzyme would drop.
    'cos_rxn_orphan',
    'cos_orphan_in',
    'cos_orphan_out',
    # A conversion whose party list holds roles that name no side of the arrow.
    # A cofactor is consumed and regenerated and a regulator is neither
    # consumed nor produced, so projecting either as a substrate would state
    # something the resource did not.
    'cos_rxn_helper',
    'cos_helper_in',
    'cos_helper_out',
    'cos_cofactor',
    'cos_regulator',
    'cos_enz_helper',
)

# The reactions that carry a **direction**. A resource states it on the
# reaction event itself, not on any one membership, so it arrives as an entity
# annotation and the fixture writes it as one. Four of these exist because each
# direction reaches the build under two spellings, and a projection that folded
# only the one it happened to meet first would read half the build as silent.
COSMOS_DIRECTION_NAMES = (
    # Reversible, stated in upper case. Two substrates and one product, so the
    # mirrored side has something to mirror.
    'cos_rxn_back',
    'cos_back_in_a',
    'cos_back_in_b',
    'cos_back_out',
    'cos_enz_back',
    # Reversible, stated in lower case.
    'cos_rxn_back_lower',
    'cos_back_lower_in',
    'cos_back_lower_out',
    'cos_enz_back_lower',
    # One way, stated with a hyphen and in upper case.
    'cos_rxn_fwd',
    'cos_fwd_in',
    'cos_fwd_out',
    'cos_enz_fwd',
    # One way, stated with an underscore and in lower case.
    'cos_rxn_fwd_lower',
    'cos_fwd_lower_in',
    'cos_fwd_lower_out',
    'cos_enz_fwd_lower',
    # Reversible and catalysed by nobody, which is where the two situations
    # meet: the placeholder node has to carry the reverse suffix as well.
    'cos_rxn_orphan_back',
    'cos_orphan_back_in',
    'cos_orphan_back_out',
)

# The reactions the **label translation** is asserted on. COSMOS reads UniProt
# on the gene side and ChEBI on the metabolite side, and the build hands the
# projection whatever identifier its resources agreed on, so a stage in between
# asks a mapping database and falls back where it gets no answer. Each of these
# reactions carries one population of that rule.
#
# The identifiers here are not `FIXTURE_*` strings, and they cannot be: the
# rule reads the shape of an identifier as well as its declared type, and a
# string no namespace would ever hold exercises none of it.
COSMOS_TRANSLATION_NAMES = (
    # Everything translates. Both metabolites are typed ChEBI and the build
    # stores them as bare numerals, which the label has to prefix; the catalyst
    # is an Entrez gene the mapping answers for.
    'cos_rxn_named',
    'cos_chebi_in',
    'cos_chebi_out',
    'cos_entrez_enz',
    # Nothing translates, and nothing is dropped. The substrate is an untyped
    # bare numeral no mapping database knows, the product is an untyped bare
    # numeral one **does** know as a ChEBI, and the catalyst is an Entrez gene
    # with no mapping at all. Three quarters of the untyped numerals on a real
    # build are the first kind, which is why shape alone never decides.
    'cos_rxn_unmapped',
    'cos_numeral_in',
    'cos_numeral_out',
    'cos_entrez_missing_enz',
    # The catalyst the build failed to type is a UniProt accession all the
    # same, which is the largest population on the gene side. Its label does
    # not change and its namespace is still `uniprot`, because the string is
    # one and a column claiming otherwise would hide coverage the network has.
    'cos_rxn_shaped',
    'cos_shaped_in',
    'cos_shaped_out',
    'cos_shaped_enz',
    # Two ambiguities, resolved in opposite directions. An InChIKey several
    # ChEBI entries claim keeps the InChIKey, because that disagreement is
    # about stereochemistry and an InChIKey is an identifier a consumer can
    # use. An Entrez gene several accessions answer for takes the lowest,
    # because a COSMOS gene node is one enzyme in one reaction and minting two
    # would state that the reaction runs twice.
    'cos_rxn_ambiguous',
    'cos_ambiguous_in',
    'cos_ambiguous_out',
    'cos_entrez_many_enz',
)

# The transport whose cargo crosses the membrane. A transport **is** a
# compartment change, and the resource states it by publishing the same
# metabolite twice on the same membership: once as the reactant, in the
# compartment it leaves, and once as the product, in the one it arrives in,
# with a stoichiometry on each side. The fuel is consumed and never
# regenerated, so it holds one role and stands beside the cargo as the
# unchanged half of the same rule.
TRANSPORT_SPLIT_NAMES = (
    'trn_b',
    'met_cargo',
    'met_fuel',
    'enz_tb',
)

# A metabolite a reaction consumes and regenerates without moving it. It holds
# both roles, exactly as a transported cargo does, and it ends up where it
# started — a cofactor turned over in the cytosol is the everyday case. What
# separates a transport from this is the compartment change and not the pair of
# roles, so a projection that read the roles alone would publish a transporter
# that transports nothing.
UNMOVED_MEMBER_NAMES = (
    'rxn_still',
    'met_still',
    'enz_still',
)

# A reaction wide enough that its participant set does not fit in a btree
# index row. The hyperedge projection groups on a text multiset of one element
# per participant — an entity uuid, a colon and a role word, some fifty bytes
# once the array header is counted — and a btree tuple cannot exceed 2704
# bytes, so a reaction past roughly fifty members makes `CREATE INDEX` fail
# outright on the staging tables that carry the signature. The build holds
# stars this wide and the binary reading never met one, so nothing in the
# graph above reproduces it. Seventy members puts the signature comfortably
# past the ceiling while still being a reaction a resource could publish.
WIDE_REACTION_MEMBERS = 70
WIDE_MEMBER_NAMES = tuple(
    f'wide_met_{index:03d}'
    for index in range(1, WIDE_REACTION_MEMBERS + 1)
)
WIDE_REACTION_NAMES = (
    'wide_rxn',
    'wide_enz',
    *WIDE_MEMBER_NAMES,
)

ENTITY_NAMES = (
    *'abcdefghijklmnopqr',
    *(
        f'{slug}_{side}'
        for slug in COVERED_CLASSES
        for side in ('s', 'o')
    ),
    'dual_s',
    'dual_o',
    # Appended, never inserted into: every id above is an index into this
    # tuple, and a test elsewhere names the row it means.
    *STAR_NAMES,
    *COSMOS_NAMES,
    *COSMOS_DIRECTION_NAMES,
    *TRANSPORT_SPLIT_NAMES,
    *COSMOS_TRANSLATION_NAMES,
    *UNMOVED_MEMBER_NAMES,
    *WIDE_REACTION_NAMES,
)

ENTITY_TYPES = {
    'rxn_a': REACTION_TYPE,
    'rxn_b': REACTION_TYPE,
    'rxn_c': REACTION_TYPE,
    'trn_a': TRANSPORT_TYPE,
    'pathway_a': PATHWAY_TYPE,
    'met_x': CHEMICAL_TYPE,
    'met_y': CHEMICAL_TYPE,
    'met_z': CHEMICAL_TYPE,
    'met_in': CHEMICAL_TYPE,
    'met_out': CHEMICAL_TYPE,
    'cos_rxn': REACTION_TYPE,
    'cos_rxn_orphan': REACTION_TYPE,
    'cos_rxn_helper': REACTION_TYPE,
    'cos_met_a': CHEMICAL_TYPE,
    'cos_met_b': CHEMICAL_TYPE,
    'cos_met_c': CHEMICAL_TYPE,
    'cos_orphan_in': CHEMICAL_TYPE,
    'cos_orphan_out': CHEMICAL_TYPE,
    'cos_helper_in': CHEMICAL_TYPE,
    'cos_helper_out': CHEMICAL_TYPE,
    'cos_cofactor': CHEMICAL_TYPE,
    # The regulator stays a protein, so a projection that turned every
    # non-enzyme protein party into a gene node would be caught here.
    'cos_rxn_back': REACTION_TYPE,
    'cos_rxn_back_lower': REACTION_TYPE,
    'cos_rxn_fwd': REACTION_TYPE,
    'cos_rxn_fwd_lower': REACTION_TYPE,
    'cos_rxn_orphan_back': REACTION_TYPE,
    'cos_back_in_a': CHEMICAL_TYPE,
    'cos_back_in_b': CHEMICAL_TYPE,
    'cos_back_out': CHEMICAL_TYPE,
    'cos_back_lower_in': CHEMICAL_TYPE,
    'cos_back_lower_out': CHEMICAL_TYPE,
    'cos_fwd_in': CHEMICAL_TYPE,
    'cos_fwd_out': CHEMICAL_TYPE,
    'cos_fwd_lower_in': CHEMICAL_TYPE,
    'cos_fwd_lower_out': CHEMICAL_TYPE,
    'cos_orphan_back_in': CHEMICAL_TYPE,
    'cos_orphan_back_out': CHEMICAL_TYPE,
    'trn_b': TRANSPORT_TYPE,
    'met_cargo': CHEMICAL_TYPE,
    'met_fuel': CHEMICAL_TYPE,
    'rxn_still': REACTION_TYPE,
    'met_still': CHEMICAL_TYPE,
    'cos_rxn_named': REACTION_TYPE,
    'cos_rxn_unmapped': REACTION_TYPE,
    'cos_rxn_shaped': REACTION_TYPE,
    'cos_rxn_ambiguous': REACTION_TYPE,
    'cos_chebi_in': CHEMICAL_TYPE,
    'cos_chebi_out': CHEMICAL_TYPE,
    'cos_numeral_in': CHEMICAL_TYPE,
    'cos_numeral_out': CHEMICAL_TYPE,
    'cos_shaped_in': CHEMICAL_TYPE,
    'cos_shaped_out': CHEMICAL_TYPE,
    'cos_ambiguous_in': CHEMICAL_TYPE,
    'cos_ambiguous_out': CHEMICAL_TYPE,
    'wide_rxn': REACTION_TYPE,
    **{name: CHEMICAL_TYPE for name in WIDE_MEMBER_NAMES},
}

# The identifier every other entity carries is `FIXTURE_<name>`, which is a
# string no namespace would claim and exactly what the rest of the graph wants.
# The translation rows need real ones, because the rule reads the shape of an
# identifier as well as its declared type: ChEBI numerals that have to gain a
# prefix, Entrez genes a mapping answers or does not, an accession the build
# failed to type. Written out rather than generated, so a test can name the
# string it expects to see in a node label.
ENTITY_IDENTIFIERS = {
    'cos_chebi_in': '15422',
    'cos_chebi_out': '16810',
    'cos_entrez_enz': '7157',
    'cos_numeral_in': '90000001',
    'cos_numeral_out': '17234',
    'cos_entrez_missing_enz': '90000002',
    'cos_shaped_out': '100',
    'cos_shaped_enz': 'P12345',
    'cos_ambiguous_in': 'FIXTUREINCHIKA-FIXTUREKEY-N',
    'cos_ambiguous_out': '200',
    'cos_entrez_many_enz': '7158',
}

# What the build canonicalised each of those to, as the name of a row in
# `vocab_identifier_type`. Every other entity in the graph carries no declared
# type at all, which is its own situation and the one the rest of the suite
# runs under.
#
# `cos_numeral_in` and `cos_numeral_out` are both untyped on purpose. They are
# the population the build could resolve to nothing, and the projection has to
# tell them apart by asking rather than by looking at the digits.
ENTITY_ID_TYPES = {
    'cos_chebi_in': 'Chebi:MI:0474',
    'cos_chebi_out': 'Chebi:MI:0474',
    'cos_entrez_enz': 'Entrez:MI:0477',
    'cos_numeral_in': 'omnipath:unresolved_entity_key',
    'cos_numeral_out': 'omnipath:unresolved_entity_key',
    'cos_entrez_missing_enz': 'Entrez:MI:0477',
    'cos_shaped_in': 'omnipath:unresolved_entity_key',
    'cos_shaped_out': 'Chebi:MI:0474',
    'cos_shaped_enz': 'omnipath:unresolved_entity_key',
    'cos_ambiguous_in': 'Standard Inchi Key:MI:1101',
    'cos_ambiguous_out': 'Chebi:MI:0474',
    'cos_entrez_many_enz': 'Entrez:MI:0477',
}

ENTITY = {
    name: f'e0000000-0000-4000-8000-{index:012d}'
    for index, name in enumerate(ENTITY_NAMES, start=1)
}

# The wide reaction's members are the one group that does **not** take the id
# above. What overruns a btree index row is the width the signature has after
# Postgres has tried to compress it, and an id differing from its neighbour in
# a single digit compresses by an order of magnitude: seventy of those make a
# signature that reports four thousand bytes and indexes without complaint,
# which is the failure passing unnoticed rather than being reproduced. The
# build's entity ids are content hashes — random-looking hex that pglz gives
# up on — so these are hashes too, and a hash of the name keeps them as fixed
# and as nameable as every other id here.
_WIDE_MEMBER_NAMESPACE = uuid.UUID('e0000000-0000-4000-8000-000000000000')
ENTITY.update(
    {
        name: str(uuid.uuid5(_WIDE_MEMBER_NAMESPACE, name))
        for name in WIDE_MEMBER_NAMES
    }
)

# Relations: (key, subject letter, predicate, object letter).
RELATIONS = (
    ('lr', 'a', 'interacts_with', 'b'),
    ('sig_pos', 'c', 'positively_regulates', 'd'),
    ('sig_neg', 'c', 'negatively_regulates', 'd'),
    ('sig_plain', 'c', 'controls', 'd'),
    ('unsigned', 'e', 'interacts_with', 'f'),
    ('forward', 'g', 'controls', 'h'),
    ('reverse', 'h', 'controls', 'g'),
    ('self_pos', 'i', 'positively_regulates', 'j'),
    ('self_neg', 'i', 'negatively_regulates', 'j'),
    ('orthosteric', 'k', 'interacts_with', 'l'),
    ('allosteric', 'm', 'interacts_with', 'n'),
    ('transport', 'o', 'transports', 'p'),
    ('other', 'q', 'has_member', 'r'),
    # One coverage pair per class, the maturation pair reversed so the ordered
    # key holds outside `signaling` too, and an endpoint pair carrying two
    # classes at once.
    *(
        (f'pair_{slug}', f'{slug}_s', predicate, f'{slug}_o')
        for slug, predicate, _terms in CLASS_PAIRS
    ),
    ('pair_maturation_reverse', 'maturation_o', MATURATION_PREDICATE,
     'maturation_s'),
    ('dual_signaling', 'dual_s', 'positively_regulates', 'dual_o'),
    ('dual_orthosteric', 'dual_s', 'interacts_with', 'dual_o'),
    # The stars. The parent is the **subject** of `has_participant` and the
    # member its object, which is the direction the membership projection
    # writes. The catalyst is the subject of a `controls` edge pointing at
    # the parent, which is the direction `predicate_for_membership` flips a
    # catalytic membership into.
    ('rxn_a_x', 'rxn_a', 'has_participant', 'met_x'),
    ('rxn_a_y', 'rxn_a', 'has_participant', 'met_y'),
    ('rxn_a_z', 'rxn_a', 'has_participant', 'met_z'),
    ('rxn_a_enzyme', 'enz_p', 'controls', 'rxn_a'),
    ('rxn_b_x', 'rxn_b', 'has_participant', 'met_x'),
    ('rxn_b_y', 'rxn_b', 'has_participant', 'met_y'),
    ('rxn_b_z', 'rxn_b', 'has_participant', 'met_z'),
    ('rxn_b_enzyme', 'enz_p', 'controls', 'rxn_b'),
    ('trn_a_in', 'trn_a', 'has_participant', 'met_in'),
    ('trn_a_out', 'trn_a', 'has_participant', 'met_out'),
    ('trn_a_enzyme', 'enz_t', 'controls', 'trn_a'),
    ('pathway_a_1', 'pathway_a', 'has_participant', 'pw_m1'),
    ('pathway_a_2', 'pathway_a', 'has_participant', 'pw_m2'),
    ('pathway_a_3', 'pathway_a', 'has_participant', 'pw_m3'),
    # A reaction whose whole participant set is the two entities the `other`
    # coverage pair already joins, so the shared identity scheme hands both
    # readings the same header id. The projection has to resolve that rather
    # than write four party rows under an `arity` of two.
    ('rxn_c_q', 'rxn_c', 'has_participant', 'q'),
    ('rxn_c_r', 'rxn_c', 'has_participant', 'r'),
    # The reactions the binary metabolic projection reads. Same shape as the
    # stars above: the parent is the subject of every membership, and the
    # catalyst points at the parent under `controls`.
    ('cos_reactant_a', 'cos_rxn', 'has_participant', 'cos_met_a'),
    ('cos_reactant_b', 'cos_rxn', 'has_participant', 'cos_met_b'),
    ('cos_product', 'cos_rxn', 'has_participant', 'cos_met_c'),
    ('cos_enzyme', 'cos_enz', 'controls', 'cos_rxn'),
    # No `controls` edge points at this one, so it reaches the header with no
    # enzyme party at all.
    ('cos_orphan_reactant', 'cos_rxn_orphan', 'has_participant',
     'cos_orphan_in'),
    ('cos_orphan_product', 'cos_rxn_orphan', 'has_participant',
     'cos_orphan_out'),
    ('cos_helper_reactant', 'cos_rxn_helper', 'has_participant',
     'cos_helper_in'),
    ('cos_helper_product', 'cos_rxn_helper', 'has_participant',
     'cos_helper_out'),
    ('cos_helper_cofactor', 'cos_rxn_helper', 'has_participant',
     'cos_cofactor'),
    ('cos_helper_regulator', 'cos_rxn_helper', 'has_participant',
     'cos_regulator'),
    ('cos_helper_enzyme', 'cos_enz_helper', 'controls', 'cos_rxn_helper'),
    # The reactions whose event entity carries a direction. The stars
    # themselves say nothing about it — the annotation rides on the event, and
    # `ENTITY_ANNOTATION` below is where it lives.
    ('cos_back_reactant_a', 'cos_rxn_back', 'has_participant', 'cos_back_in_a'),
    ('cos_back_reactant_b', 'cos_rxn_back', 'has_participant', 'cos_back_in_b'),
    ('cos_back_product', 'cos_rxn_back', 'has_participant', 'cos_back_out'),
    ('cos_back_enzyme', 'cos_enz_back', 'controls', 'cos_rxn_back'),
    ('cos_back_lower_reactant', 'cos_rxn_back_lower', 'has_participant',
     'cos_back_lower_in'),
    ('cos_back_lower_product', 'cos_rxn_back_lower', 'has_participant',
     'cos_back_lower_out'),
    ('cos_back_lower_enzyme', 'cos_enz_back_lower', 'controls',
     'cos_rxn_back_lower'),
    ('cos_fwd_reactant', 'cos_rxn_fwd', 'has_participant', 'cos_fwd_in'),
    ('cos_fwd_product', 'cos_rxn_fwd', 'has_participant', 'cos_fwd_out'),
    ('cos_fwd_enzyme', 'cos_enz_fwd', 'controls', 'cos_rxn_fwd'),
    ('cos_fwd_lower_reactant', 'cos_rxn_fwd_lower', 'has_participant',
     'cos_fwd_lower_in'),
    ('cos_fwd_lower_product', 'cos_rxn_fwd_lower', 'has_participant',
     'cos_fwd_lower_out'),
    ('cos_fwd_lower_enzyme', 'cos_enz_fwd_lower', 'controls',
     'cos_rxn_fwd_lower'),
    ('cos_orphan_back_reactant', 'cos_rxn_orphan_back', 'has_participant',
     'cos_orphan_back_in'),
    ('cos_orphan_back_product', 'cos_rxn_orphan_back', 'has_participant',
     'cos_orphan_back_out'),
    # The reactions the label translation reads. Same star shape as the rest:
    # the parent is the subject of every membership and the catalyst points at
    # the parent under `controls`.
    ('cos_named_reactant', 'cos_rxn_named', 'has_participant', 'cos_chebi_in'),
    ('cos_named_product', 'cos_rxn_named', 'has_participant', 'cos_chebi_out'),
    ('cos_named_enzyme', 'cos_entrez_enz', 'controls', 'cos_rxn_named'),
    ('cos_unmapped_reactant', 'cos_rxn_unmapped', 'has_participant',
     'cos_numeral_in'),
    ('cos_unmapped_product', 'cos_rxn_unmapped', 'has_participant',
     'cos_numeral_out'),
    ('cos_unmapped_enzyme', 'cos_entrez_missing_enz', 'controls',
     'cos_rxn_unmapped'),
    ('cos_shaped_reactant', 'cos_rxn_shaped', 'has_participant',
     'cos_shaped_in'),
    ('cos_shaped_product', 'cos_rxn_shaped', 'has_participant',
     'cos_shaped_out'),
    ('cos_shaped_enzyme', 'cos_shaped_enz', 'controls', 'cos_rxn_shaped'),
    ('cos_ambiguous_reactant', 'cos_rxn_ambiguous', 'has_participant',
     'cos_ambiguous_in'),
    ('cos_ambiguous_product', 'cos_rxn_ambiguous', 'has_participant',
     'cos_ambiguous_out'),
    ('cos_ambiguous_enzyme', 'cos_entrez_many_enz', 'controls',
     'cos_rxn_ambiguous'),
    # The transport whose cargo holds both roles. `trn_b_cargo` is **one**
    # relation — `relation` is unique on (subject, predicate, object), so the
    # two statements about the cargo cannot be two of them — and the role it
    # carries is whichever of the two evidence rows under
    # `SPLIT_ROLE_EVIDENCE` you are reading. The fuel and the transporter
    # take the ordinary path.
    ('trn_b_cargo', 'trn_b', 'has_participant', 'met_cargo'),
    ('trn_b_fuel', 'trn_b', 'has_participant', 'met_fuel'),
    ('trn_b_enzyme', 'enz_tb', 'controls', 'trn_b'),
    # The reaction that turns a cofactor over without moving it. Same shape as
    # the transport above, down to the two evidence rows on one membership,
    # and the compartment is the same on both of them.
    ('rxn_still_member', 'rxn_still', 'has_participant', 'met_still'),
    ('rxn_still_enzyme', 'enz_still', 'controls', 'rxn_still'),
    # The wide reaction. Ordinary star shape, and the only thing unusual about
    # it is how many spokes it has.
    *(
        (f'wide_member_{index:03d}', 'wide_rxn', 'has_participant', name)
        for index, name in enumerate(WIDE_MEMBER_NAMES, start=1)
    ),
    ('wide_enzyme', 'wide_enz', 'controls', 'wide_rxn'),
)

# `has_participant` is a membership, and the build files memberships under the
# association category rather than the interaction one. Every other relation in
# the graph is an interaction, which is what the default says.
RELATION_CATEGORY = {
    key: 1
    for key, _subject, predicate, _object in RELATIONS
    if predicate == 'has_participant'
}

# Evidence: (relation key, source id, relation-level annotation terms).
EVIDENCE = (
    ('lr', SOURCE_LR, ()),
    ('sig_pos', SOURCE_A, ('Activation:OM:0930',)),
    ('sig_neg', SOURCE_B, ('Inhibition:OM:0931',)),
    # The third resource on the same endpoint pair asserts neither sign nor
    # anything else — only a reference. It stays in `sources` all the same.
    ('sig_plain', SOURCE_C, ()),
    ('unsigned', SOURCE_C, ()),
    ('forward', SOURCE_A, ()),
    ('reverse', SOURCE_A, ()),
    # One resource asserting both signs, under two predicates that share the
    # `signaling` class. 93 per cent of both-flags rows are this rather than
    # a disagreement between resources.
    ('self_pos', SOURCE_A, ('Activation:OM:0930',)),
    ('self_neg', SOURCE_A, ('Inhibition:OM:0931',)),
    ('orthosteric', SOURCE_A, ('Agonist:OM:1001',)),
    ('allosteric', SOURCE_A, ('Allosteric Modulator:OM:1005',)),
    ('transport', SOURCE_A, ()),
    ('other', SOURCE_A, ()),
    # Two resources on every coverage pair, so folding one is folding a group
    # rather than copying a row. The annotation that names the class rides on
    # the first resource alone — it resolves the relation, not the record.
    *(
        (f'pair_{slug}', source, terms if source == SOURCE_A else ())
        for slug, _predicate, terms in CLASS_PAIRS
        for source in (SOURCE_A, SOURCE_B)
    ),
    ('pair_maturation_reverse', SOURCE_A, ()),
    ('dual_signaling', SOURCE_A, ()),
    ('dual_orthosteric', SOURCE_B, ('Agonist:OM:1001',)),
    # One resource per reaction parent, so the merge the projection performs is
    # a merge across resources rather than a fold inside one.
    ('rxn_a_x', SOURCE_A, ()),
    ('rxn_a_y', SOURCE_A, ()),
    ('rxn_a_z', SOURCE_A, ()),
    ('rxn_a_enzyme', SOURCE_A, ()),
    ('rxn_b_x', SOURCE_B, ()),
    ('rxn_b_y', SOURCE_B, ()),
    ('rxn_b_z', SOURCE_B, ()),
    ('rxn_b_enzyme', SOURCE_B, ()),
    ('trn_a_in', SOURCE_A, ()),
    ('trn_a_out', SOURCE_A, ()),
    ('trn_a_enzyme', SOURCE_A, ()),
    ('pathway_a_1', SOURCE_A, ()),
    ('pathway_a_2', SOURCE_A, ()),
    ('pathway_a_3', SOURCE_A, ()),
    ('rxn_c_q', SOURCE_C, ()),
    ('rxn_c_r', SOURCE_C, ()),
    # One resource per binary-projection reaction. The projection carries the
    # contributing resource onto every edge it emits, so a single publisher
    # makes that column a claim a test can name.
    ('cos_reactant_a', SOURCE_A, ()),
    ('cos_reactant_b', SOURCE_A, ()),
    ('cos_product', SOURCE_A, ()),
    ('cos_enzyme', SOURCE_A, ()),
    ('cos_orphan_reactant', SOURCE_A, ()),
    ('cos_orphan_product', SOURCE_A, ()),
    ('cos_helper_reactant', SOURCE_A, ()),
    ('cos_helper_product', SOURCE_A, ()),
    ('cos_helper_cofactor', SOURCE_A, ()),
    ('cos_helper_regulator', SOURCE_A, ()),
    ('cos_helper_enzyme', SOURCE_A, ()),
    ('cos_back_reactant_a', SOURCE_A, ()),
    ('cos_back_reactant_b', SOURCE_A, ()),
    ('cos_back_product', SOURCE_A, ()),
    ('cos_back_enzyme', SOURCE_A, ()),
    ('cos_back_lower_reactant', SOURCE_A, ()),
    ('cos_back_lower_product', SOURCE_A, ()),
    ('cos_back_lower_enzyme', SOURCE_A, ()),
    ('cos_fwd_reactant', SOURCE_A, ()),
    ('cos_fwd_product', SOURCE_A, ()),
    ('cos_fwd_enzyme', SOURCE_A, ()),
    ('cos_fwd_lower_reactant', SOURCE_A, ()),
    ('cos_fwd_lower_product', SOURCE_A, ()),
    ('cos_fwd_lower_enzyme', SOURCE_A, ()),
    ('cos_orphan_back_reactant', SOURCE_A, ()),
    ('cos_orphan_back_product', SOURCE_A, ()),
    ('cos_named_reactant', SOURCE_A, ()),
    ('cos_named_product', SOURCE_A, ()),
    ('cos_named_enzyme', SOURCE_A, ()),
    ('cos_unmapped_reactant', SOURCE_A, ()),
    ('cos_unmapped_product', SOURCE_A, ()),
    ('cos_unmapped_enzyme', SOURCE_A, ()),
    ('cos_shaped_reactant', SOURCE_A, ()),
    ('cos_shaped_product', SOURCE_A, ()),
    ('cos_shaped_enzyme', SOURCE_A, ()),
    ('cos_ambiguous_reactant', SOURCE_A, ()),
    ('cos_ambiguous_product', SOURCE_A, ()),
    ('cos_ambiguous_enzyme', SOURCE_A, ()),
    # The transport with the two-role cargo. `trn_b_cargo` is deliberately
    # absent here: one row per (relation, resource) is exactly what it cannot
    # be expressed as, and `SPLIT_ROLE_EVIDENCE` below writes its two.
    ('trn_b_fuel', SOURCE_A, ()),
    ('trn_b_enzyme', SOURCE_A, ()),
    # `rxn_still_member` is absent for the same reason `trn_b_cargo` is.
    ('rxn_still_enzyme', SOURCE_A, ()),
    *(
        (f'wide_member_{index:03d}', SOURCE_A, ())
        for index in range(1, WIDE_REACTION_MEMBERS + 1)
    ),
    ('wide_enzyme', SOURCE_A, ()),
)

# Participant-descriptive annotations: (relation key, source, term, value).
# These ride on the **relation** evidence at `object` scope, because the thing
# they describe is the object of `has_participant` — the member — and not the
# relation as a whole. The catalyst's role annotation is written at the same
# scope even though the party it describes is the *subject* of the `controls`
# edge: the membership projection stamps the scope before the catalytic flip
# swaps the endpoints, so the column says `object` on every one of the 810,240
# catalyst annotations in the build. The fixture reproduces that rather than
# correcting it, so a projection that trusted the scope column would fail here
# exactly as it fails on real data.
PARTICIPANT_EVIDENCE = (
    ('rxn_a_x', SOURCE_A, 'Reactant:OM:0310', None),
    ('rxn_a_x', SOURCE_A, 'Stoichiometry:OM:1226', '2'),
    ('rxn_a_x', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('rxn_a_y', SOURCE_A, 'Reactant:OM:0310', None),
    ('rxn_a_y', SOURCE_A, 'Stoichiometry:OM:1226', '1'),
    ('rxn_a_y', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('rxn_a_z', SOURCE_A, 'Product:OM:0311', None),
    ('rxn_a_z', SOURCE_A, 'Stoichiometry:OM:1226', '1'),
    ('rxn_a_z', SOURCE_A, 'Subcellular Location:OM:0604', 'm'),
    ('rxn_a_enzyme', SOURCE_A, 'Enzyme:MI:0501', None),
    ('rxn_b_x', SOURCE_B, 'Reactant:OM:0310', None),
    ('rxn_b_x', SOURCE_B, 'Stoichiometry:OM:1226', '2'),
    ('rxn_b_x', SOURCE_B, 'Subcellular Location:OM:0604', 'c'),
    # The second resource states the role and nothing else for this member, so
    # the merged party can only carry a stoichiometry and a compartment if the
    # merge reaches across the two resources.
    ('rxn_b_y', SOURCE_B, 'Reactant:OM:0310', None),
    ('rxn_b_z', SOURCE_B, 'Product:OM:0311', None),
    ('rxn_b_z', SOURCE_B, 'Stoichiometry:OM:1226', '1'),
    ('rxn_b_z', SOURCE_B, 'Subcellular Location:OM:0604', 'm'),
    ('rxn_b_enzyme', SOURCE_B, 'Enzyme:MI:0501', None),
    # A transport states where each side of the membrane is, not which
    # organelle the member sits in.
    ('trn_a_in', SOURCE_A, 'Reactant:OM:0310', None),
    ('trn_a_in', SOURCE_A, 'Membrane Side:OM:1231', 'in'),
    ('trn_a_out', SOURCE_A, 'Product:OM:0311', None),
    ('trn_a_out', SOURCE_A, 'Membrane Side:OM:1231', 'out'),
    ('trn_a_enzyme', SOURCE_A, 'Enzyme:MI:0501', None),
    # The pathway's members carry the pathway-component role, which is not a
    # reaction role and names no side of an arrow.
    ('pathway_a_1', SOURCE_A, 'Pathway Component:OM:0315', None),
    ('pathway_a_2', SOURCE_A, 'Pathway Component:OM:0315', None),
    ('pathway_a_3', SOURCE_A, 'Pathway Component:OM:0315', None),
    ('rxn_c_q', SOURCE_C, 'Reactant:OM:0310', None),
    ('rxn_c_q', SOURCE_C, 'Stoichiometry:OM:1226', '3'),
    ('rxn_c_r', SOURCE_C, 'Product:OM:0311', None),
    # The ordinary conversion. `cos_met_b` is stated without a compartment, so
    # the node string for it can carry no location and must not carry an empty
    # one either.
    ('cos_reactant_a', SOURCE_A, 'Reactant:OM:0310', None),
    ('cos_reactant_a', SOURCE_A, 'Stoichiometry:OM:1226', '2'),
    ('cos_reactant_a', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_reactant_b', SOURCE_A, 'Reactant:OM:0310', None),
    ('cos_product', SOURCE_A, 'Product:OM:0311', None),
    ('cos_product', SOURCE_A, 'Subcellular Location:OM:0604', 'm'),
    ('cos_enzyme', SOURCE_A, 'Enzyme:MI:0501', None),
    # The catalyst-free conversion. Both sides are stated, the enzyme is not.
    ('cos_orphan_reactant', SOURCE_A, 'Reactant:OM:0310', None),
    ('cos_orphan_reactant', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_orphan_product', SOURCE_A, 'Product:OM:0311', None),
    ('cos_orphan_product', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    # The conversion with helpers. The cofactor and the regulator hold roles
    # that name no side of the arrow, and the projection has to leave both
    # where they are rather than read them as substrates.
    ('cos_helper_reactant', SOURCE_A, 'Reactant:OM:0310', None),
    ('cos_helper_reactant', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_helper_product', SOURCE_A, 'Product:OM:0311', None),
    ('cos_helper_product', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_helper_cofactor', SOURCE_A, 'Cofactor:OM:0317', None),
    ('cos_helper_cofactor', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_helper_regulator', SOURCE_A, 'Regulator:MI:2274', None),
    ('cos_helper_enzyme', SOURCE_A, 'Enzyme:MI:0501', None),
    # The directional reactions. Each member sits in the cytosol, so a test can
    # name the node string it expects without a second compartment to track.
    ('cos_back_reactant_a', SOURCE_A, 'Reactant:OM:0310', None),
    ('cos_back_reactant_a', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_back_reactant_b', SOURCE_A, 'Reactant:OM:0310', None),
    ('cos_back_reactant_b', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_back_product', SOURCE_A, 'Product:OM:0311', None),
    ('cos_back_product', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_back_enzyme', SOURCE_A, 'Enzyme:MI:0501', None),
    ('cos_back_lower_reactant', SOURCE_A, 'Reactant:OM:0310', None),
    ('cos_back_lower_reactant', SOURCE_A,
     'Subcellular Location:OM:0604', 'c'),
    ('cos_back_lower_product', SOURCE_A, 'Product:OM:0311', None),
    ('cos_back_lower_product', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_back_lower_enzyme', SOURCE_A, 'Enzyme:MI:0501', None),
    ('cos_fwd_reactant', SOURCE_A, 'Reactant:OM:0310', None),
    ('cos_fwd_reactant', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_fwd_product', SOURCE_A, 'Product:OM:0311', None),
    ('cos_fwd_product', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_fwd_enzyme', SOURCE_A, 'Enzyme:MI:0501', None),
    ('cos_fwd_lower_reactant', SOURCE_A, 'Reactant:OM:0310', None),
    ('cos_fwd_lower_reactant', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_fwd_lower_product', SOURCE_A, 'Product:OM:0311', None),
    ('cos_fwd_lower_product', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_fwd_lower_enzyme', SOURCE_A, 'Enzyme:MI:0501', None),
    ('cos_orphan_back_reactant', SOURCE_A, 'Reactant:OM:0310', None),
    ('cos_orphan_back_reactant', SOURCE_A,
     'Subcellular Location:OM:0604', 'c'),
    ('cos_orphan_back_product', SOURCE_A, 'Product:OM:0311', None),
    ('cos_orphan_back_product', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    # The translation reactions. Every substrate sits in the cytosol, so a test
    # naming a node string has one compartment to track and the prefix rule is
    # asserted on a label that carries a location as well as an identifier.
    ('cos_named_reactant', SOURCE_A, 'Reactant:OM:0310', None),
    ('cos_named_reactant', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_named_product', SOURCE_A, 'Product:OM:0311', None),
    ('cos_named_product', SOURCE_A, 'Subcellular Location:OM:0604', 'm'),
    ('cos_named_enzyme', SOURCE_A, 'Enzyme:MI:0501', None),
    ('cos_unmapped_reactant', SOURCE_A, 'Reactant:OM:0310', None),
    ('cos_unmapped_reactant', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_unmapped_product', SOURCE_A, 'Product:OM:0311', None),
    ('cos_unmapped_product', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_unmapped_enzyme', SOURCE_A, 'Enzyme:MI:0501', None),
    ('cos_shaped_reactant', SOURCE_A, 'Reactant:OM:0310', None),
    ('cos_shaped_reactant', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_shaped_product', SOURCE_A, 'Product:OM:0311', None),
    ('cos_shaped_product', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_shaped_enzyme', SOURCE_A, 'Enzyme:MI:0501', None),
    ('cos_ambiguous_reactant', SOURCE_A, 'Reactant:OM:0310', None),
    ('cos_ambiguous_reactant', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_ambiguous_product', SOURCE_A, 'Product:OM:0311', None),
    ('cos_ambiguous_product', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('cos_ambiguous_enzyme', SOURCE_A, 'Enzyme:MI:0501', None),
    # The fuel of the transport below: consumed in the cytosol and never
    # produced, so it holds one role and one compartment.
    ('trn_b_fuel', SOURCE_A, 'Reactant:OM:0310', None),
    ('trn_b_fuel', SOURCE_A, 'Stoichiometry:OM:1226', '1'),
    ('trn_b_fuel', SOURCE_A, 'Subcellular Location:OM:0604', 'c'),
    ('trn_b_enzyme', SOURCE_A, 'Enzyme:MI:0501', None),
    ('rxn_still_enzyme', SOURCE_A, 'Enzyme:MI:0501', None),
    # The wide reaction's members, half on each side of the arrow. Each one
    # carries a role, because the role is half of what the signature is built
    # from and a member with none would shorten it.
    *(
        (
            f'wide_member_{index:03d}',
            SOURCE_A,
            'Reactant:OM:0310' if index % 2 else 'Product:OM:0311',
            None,
        )
        for index in range(1, WIDE_REACTION_MEMBERS + 1)
    ),
    *(
        (
            f'wide_member_{index:03d}',
            SOURCE_A,
            'Subcellular Location:OM:0604',
            'c',
        )
        for index in range(1, WIDE_REACTION_MEMBERS + 1)
    ),
    ('wide_enzyme', SOURCE_A, 'Enzyme:MI:0501', None),
)

# The two statements a resource makes about a metabolite it transports. Both
# ride on the **same** canonical relation and are told apart by the evidence
# row each sits on, which is the grain the resource publishes them at: the
# role and the compartment that qualifies it are on one row together, and
# reading either without the other says nothing. A fixture that keyed one
# evidence row per (relation, resource) — which is what `EVIDENCE` above does
# — cannot express this, so these rows carry an occurrence number and are
# written separately.
#
# The cargo leaves the cytosol two at a time and arrives outside one at a
# time. The stoichiometries differ so that a projection collapsing them to one
# value is caught on the number as well as on the compartment.
#
# (relation key, source, occurrence, ((term, value), ...))
SPLIT_ROLE_EVIDENCE = (
    (
        'trn_b_cargo',
        SOURCE_A,
        1,
        (
            ('Reactant:OM:0310', None),
            ('Subcellular Location:OM:0604', 'c'),
            ('Stoichiometry:OM:1226', '2'),
        ),
    ),
    (
        'trn_b_cargo',
        SOURCE_A,
        2,
        (
            ('Product:OM:0311', None),
            ('Subcellular Location:OM:0604', 'e'),
            ('Stoichiometry:OM:1226', '1'),
        ),
    ),
    # Consumed and regenerated in the cytosol. Both roles, one compartment,
    # and the reaction has a catalyst — everything a transport has except the
    # movement.
    (
        'rxn_still_member',
        SOURCE_A,
        1,
        (
            ('Reactant:OM:0310', None),
            ('Subcellular Location:OM:0604', 'c'),
        ),
    ),
    (
        'rxn_still_member',
        SOURCE_A,
        2,
        (
            ('Product:OM:0311', None),
            ('Subcellular Location:OM:0604', 'c'),
        ),
    ),
)

# The term a resource publishes the direction of a conversion under, and the
# four spellings it reaches the build in. Both halves of each pair are live on
# the build — the upper-case forms carry the KEGG and Reactome side, the
# lower-case forms the Human-GEM one — so a projection has to fold case and
# read `-` and `_` alike before it compares.
CONVERSION_DIRECTION_TERM = 'Conversion Direction:OM:1211'

# Annotations on an **entity**: (entity name, source, term, value). A
# conversion's direction describes the event, not any one of its memberships,
# so the resource states it here rather than on a spoke. It reaches the build
# through `entity_evidence_annotation`, and a consumer finds it by resolving
# the evidence back to the entity.
#
# `cos_rxn`, `cos_rxn_orphan` and `cos_rxn_helper` are deliberately absent. A
# resource that publishes no direction is not the same as one that publishes a
# one-way reaction, and something has to hold the first case.
ENTITY_ANNOTATION = (
    ('cos_rxn_back', SOURCE_A, CONVERSION_DIRECTION_TERM, 'REVERSIBLE'),
    ('cos_rxn_back_lower', SOURCE_A, CONVERSION_DIRECTION_TERM, 'reversible'),
    ('cos_rxn_fwd', SOURCE_A, CONVERSION_DIRECTION_TERM, 'LEFT-TO-RIGHT'),
    ('cos_rxn_fwd_lower', SOURCE_A, CONVERSION_DIRECTION_TERM,
     'left_to_right'),
    ('cos_rxn_orphan_back', SOURCE_A, CONVERSION_DIRECTION_TERM, 'REVERSIBLE'),
)

# Which evidence rows carry a PubMed reference, and which id.
REFERENCES = {
    ('sig_pos', SOURCE_A): '11111111',
    ('sig_plain', SOURCE_C): '33333333',
    ('lr', SOURCE_LR): '44444444',
    # Both contributors to a coverage pair publish a reference, and they
    # publish different ones, so a collapsed row reporting two of them is
    # reporting a collection rather than one resource's list.
    **{
        (f'pair_{slug}', source): f'{prefix}{index:07d}'
        for index, (slug, _predicate, _terms) in enumerate(CLASS_PAIRS, start=1)
        for source, prefix in ((SOURCE_A, '7'), (SOURCE_B, '8'))
    },
}

PREDICATES = (
    'interacts_with',
    'positively_regulates',
    'negatively_regulates',
    'controls',
    'transports',
    'has_member',
    'has_participant',
    MATURATION_PREDICATE,
)

# Participant-role annotations, and the endpoint evidence rows they hang off:
# (relation key, source, side, term). Role evidence is the one tier that lives
# at the entity-evidence grain, so a relation reaching a class through it needs
# an evidence row per endpoint rather than a bare entity id.
ROLE_EVIDENCE = (
    ('lr', SOURCE_LR, 'subject', 'Ligand:OM:7777'),
    ('lr', SOURCE_LR, 'object', 'Receptor:OM:7778'),
    ('pair_ligand_receptor', SOURCE_A, 'subject', 'Ligand:OM:7777'),
    ('pair_ligand_receptor', SOURCE_A, 'object', 'Receptor:OM:7778'),
)

TAXONOMY_ID = 9606

# Where the entity-annotation evidence rows start numbering, clear of the row
# ids the relation evidence above uses.
_ENTITY_ANNOTATION_ROW_BASE = 900000

# And where the split-role evidence rows start, clear of both.
_SPLIT_ROLE_ROW_BASE = 800000


def _uuid5(prefix: str, index: int) -> str:
    """A stable uuid: a hex prefix in the first block, the index in the last."""
    return f'{prefix}-0000-4000-8000-{index:012d}'


def build_interaction_fixture(
    conn: psycopg2.extensions.connection,
    schema: str,
) -> dict[str, object]:
    """Populate ``schema`` with the fixture graph; return its named ids."""
    schema_id = sql.Identifier(schema)

    def q(statement: str) -> sql.Composed:
        return sql.SQL(statement).format(schema_id)

    relation_ids = {
        key: _uuid5('11111111', index)
        for index, (key, *_rest) in enumerate(RELATIONS, start=1)
    }
    evidence_ids = {
        (key, source): _uuid5('22222222', index)
        for index, (key, source, _terms) in enumerate(EVIDENCE, start=1)
    }
    # Endpoint evidence ids exist only for the ligand/receptor relations: those
    # are the ones whose class comes from participant-role annotations, which
    # hang off the entity-evidence grain.
    endpoint_ids = {
        (key, source, side): _uuid5('33333333', index)
        for index, (key, source, side, _term) in enumerate(
            ROLE_EVIDENCE,
            start=1,
        )
    }

    with conn.cursor() as cur:
        type_names = sorted({PROTEIN_TYPE, *ENTITY_TYPES.values()})
        cur.executemany(
            q('INSERT INTO {}.vocab_entity_type (name) VALUES (%s) '
              'ON CONFLICT (name) DO NOTHING').as_string(cur.connection),
            [(name,) for name in type_names],
        )
        cur.execute(q('SELECT name, entity_type_id FROM {}.vocab_entity_type'))
        entity_type_ids = dict(cur.fetchall())
        entity_type_id = entity_type_ids[PROTEIN_TYPE]

        # The identifier types the translation rows declare. They are seeded
        # rows of the schema rather than fixture inventions, so this reads the
        # ids back rather than minting them; a type the schema does not carry
        # would be a typo in the fixture and the lookup says so at once.
        cur.execute(
            q('SELECT name, identifier_type_id FROM {}.vocab_identifier_type')
        )
        identifier_type_ids = dict(cur.fetchall())
        cur.executemany(
            q(
                'INSERT INTO {}.entity (entity_id, entity_type_id, '
                'taxonomy_id, canonical_identifier, '
                'canonical_identifier_type_id, resolution_status_id) '
                'VALUES (%s, %s, %s, %s, %s, 1) ON CONFLICT DO NOTHING'
            ).as_string(cur.connection),
            [
                (
                    entity_id,
                    entity_type_ids[ENTITY_TYPES.get(letter, PROTEIN_TYPE)],
                    TAXONOMY_ID,
                    ENTITY_IDENTIFIERS.get(letter, f'FIXTURE_{letter}'),
                    (
                        identifier_type_ids[ENTITY_ID_TYPES[letter]]
                        if letter in ENTITY_ID_TYPES
                        else None
                    ),
                )
                for letter, entity_id in ENTITY.items()
            ],
        )
        cur.executemany(
            q(
                'INSERT INTO {}.data_source (source_id, name) VALUES (%s, %s) '
                'ON CONFLICT (source_id) DO NOTHING'
            ).as_string(cur.connection),
            list(SOURCE_NAMES.items()),
        )
        cur.executemany(
            q(
                'INSERT INTO {}.dataset (dataset_id, source_id, name) '
                'VALUES (%s, %s, %s) ON CONFLICT (dataset_id) DO NOTHING'
            ).as_string(cur.connection),
            [
                (source_id, source_id, f'{name}_default')
                for source_id, name in SOURCE_NAMES.items()
            ],
        )
        cur.executemany(
            q(
                'INSERT INTO {}.vocab_relation_category '
                '(relation_category_id, name) VALUES (%s, %s) '
                'ON CONFLICT (relation_category_id) DO NOTHING'
            ).as_string(cur.connection),
            [(1, 'association'), (2, 'interaction')],
        )
        cur.executemany(
            q(
                'INSERT INTO {}.vocab_relation_predicate (name) VALUES (%s) '
                'ON CONFLICT (name) DO NOTHING'
            ).as_string(cur.connection),
            [(name,) for name in PREDICATES],
        )
        cur.execute(
            q('SELECT name, relation_predicate_id '
              'FROM {}.vocab_relation_predicate')
        )
        predicate_ids = dict(cur.fetchall())

        cur.executemany(
            q(
                'INSERT INTO {}.relation (relation_id, subject_entity_id, '
                'predicate_id, object_entity_id, relation_category_id) '
                'VALUES (%s, %s, %s, %s, %s) ON CONFLICT DO NOTHING'
            ).as_string(cur.connection),
            [
                (
                    relation_ids[key],
                    ENTITY[subject],
                    predicate_ids[predicate],
                    ENTITY[obj],
                    RELATION_CATEGORY.get(key, 2),
                )
                for key, subject, predicate, obj in RELATIONS
            ],
        )

        predicate_of = {key: predicate for key, _s, predicate, _o in RELATIONS}
        evidence_rows = []
        for row_id, (key, source, _terms) in enumerate(EVIDENCE, start=1):
            subject_evidence = endpoint_ids.get((key, source, 'subject'))
            object_evidence = endpoint_ids.get((key, source, 'object'))
            subject_letter = next(s for k, s, _p, _o in RELATIONS if k == key)
            object_letter = next(o for k, _s, _p, o in RELATIONS if k == key)
            evidence_rows.append(
                (
                    source,
                    evidence_ids[(key, source)],
                    source,
                    row_id,
                    subject_evidence,
                    None if subject_evidence else ENTITY[subject_letter],
                    predicate_ids[predicate_of[key]],
                    object_evidence,
                    None if object_evidence else ENTITY[object_letter],
                    RELATION_CATEGORY.get(key, 2),
                )
            )
        # The endpoint evidence rows the ligand/receptor annotations hang off.
        # `relation_evidence` carries a foreign key onto them.
        cur.executemany(
            q(
                'INSERT INTO {}.entity_evidence (source_id, '
                'entity_evidence_id, dataset_id, row_id, entity_role_id, '
                'entity_type_id, taxonomy_id) '
                'VALUES (%s, %s, %s, %s, 1, %s, %s) ON CONFLICT DO NOTHING'
            ).as_string(cur.connection),
            [
                (
                    source,
                    endpoint_id,
                    source,
                    ordinal,
                    entity_type_id,
                    TAXONOMY_ID,
                )
                for ordinal, ((_key, source, _side), endpoint_id) in enumerate(
                    endpoint_ids.items(),
                    start=1,
                )
            ],
        )
        cur.executemany(
            q(
                'INSERT INTO {}.relation_evidence (source_id, '
                'relation_evidence_id, dataset_id, row_id, '
                'subject_entity_evidence_id, subject_entity_id, predicate_id, '
                'object_entity_evidence_id, object_entity_id, '
                'relation_category_id) '
                'VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s) '
                'ON CONFLICT DO NOTHING'
            ).as_string(cur.connection),
            evidence_rows,
        )
        cur.executemany(
            q(
                'INSERT INTO {}.relation_evidence_relation '
                '(source_id, relation_id, relation_evidence_id) '
                'VALUES (%s, %s, %s) ON CONFLICT DO NOTHING'
            ).as_string(cur.connection),
            [
                (source, relation_ids[key], evidence_ids[(key, source)])
                for key, source, _terms in EVIDENCE
            ],
        )

        # Annotations. `annotation` is content-addressed by `annotation_key`, so
        # the fixture mints one key per (term, value) pair it needs.
        annotation_keys: dict[tuple[str, str | None], str] = {}

        def annotation_key(term: str, value: str | None = None) -> str:
            existing = annotation_keys.get((term, value))
            if existing is not None:
                return existing
            key = _uuid5('44444444', len(annotation_keys) + 1)
            annotation_keys[(term, value)] = key
            cur.execute(
                q(
                    'INSERT INTO {}.annotation (annotation_key, term, value) '
                    'VALUES (%s, %s, %s) ON CONFLICT DO NOTHING'
                ),
                [key, term, value],
            )
            return key

        relation_annotations = []
        for key, source, terms in EVIDENCE:
            for term in terms:
                relation_annotations.append(
                    (source, evidence_ids[(key, source)], annotation_key(term), 1)
                )
        for (key, source), pubmed in REFERENCES.items():
            relation_annotations.append(
                (
                    source,
                    evidence_ids[(key, source)],
                    annotation_key('Pubmed:MI:0446', pubmed),
                    1,
                )
            )
        # Scope 3 is `object` in `vocab_annotation_scope`. Everything above is
        # scope 1, `relation`.
        for key, source, term, value in PARTICIPANT_EVIDENCE:
            relation_annotations.append(
                (
                    source,
                    evidence_ids[(key, source)],
                    annotation_key(term, value),
                    3,
                )
            )
        cur.executemany(
            q(
                'INSERT INTO {}.relation_evidence_annotation (source_id, '
                'relation_evidence_id, annotation_key, annotation_scope_id) '
                'VALUES (%s, %s, %s, %s) ON CONFLICT DO NOTHING'
            ).as_string(cur.connection),
            relation_annotations,
        )

        # The second and further evidence rows a resource publishes about one
        # membership. `relation` holds a single row for the cargo of a
        # transport — it is unique on (subject, predicate, object) — so the
        # reactant statement and the product statement about that cargo can
        # only be told apart by the evidence row each of them rides on. These
        # rows are what a projection has to read at if it wants to keep the
        # role and the compartment together, and the fixture writes them the
        # way the loader does: one `relation_evidence` row per statement, all
        # of them resolving to the same relation.
        split_evidence_ids = {
            (key, source, occurrence): _uuid5('66666666', index)
            for index, (key, source, occurrence, _terms) in enumerate(
                SPLIT_ROLE_EVIDENCE,
                start=1,
            )
        }
        cur.executemany(
            q(
                'INSERT INTO {}.relation_evidence (source_id, '
                'relation_evidence_id, dataset_id, row_id, '
                'subject_entity_id, predicate_id, object_entity_id, '
                'relation_category_id) '
                'VALUES (%s, %s, %s, %s, %s, %s, %s, %s) '
                'ON CONFLICT DO NOTHING'
            ).as_string(cur.connection),
            [
                (
                    source,
                    split_evidence_ids[(key, source, occurrence)],
                    source,
                    _SPLIT_ROLE_ROW_BASE + ordinal,
                    ENTITY[
                        next(s for k, s, _p, _o in RELATIONS if k == key)
                    ],
                    predicate_ids[predicate_of[key]],
                    ENTITY[
                        next(o for k, _s, _p, o in RELATIONS if k == key)
                    ],
                    RELATION_CATEGORY.get(key, 2),
                )
                for ordinal, (key, source, occurrence, _terms) in enumerate(
                    SPLIT_ROLE_EVIDENCE,
                    start=1,
                )
            ],
        )
        cur.executemany(
            q(
                'INSERT INTO {}.relation_evidence_relation '
                '(source_id, relation_id, relation_evidence_id) '
                'VALUES (%s, %s, %s) ON CONFLICT DO NOTHING'
            ).as_string(cur.connection),
            [
                (
                    source,
                    relation_ids[key],
                    split_evidence_ids[(key, source, occurrence)],
                )
                for key, source, occurrence, _terms in SPLIT_ROLE_EVIDENCE
            ],
        )
        # Scope 3, `object`, exactly as the single-row participant evidence
        # above: what these annotations describe is the member.
        cur.executemany(
            q(
                'INSERT INTO {}.relation_evidence_annotation (source_id, '
                'relation_evidence_id, annotation_key, annotation_scope_id) '
                'VALUES (%s, %s, %s, %s) ON CONFLICT DO NOTHING'
            ).as_string(cur.connection),
            [
                (
                    source,
                    split_evidence_ids[(key, source, occurrence)],
                    annotation_key(term, value),
                    3,
                )
                for key, source, occurrence, terms in SPLIT_ROLE_EVIDENCE
                for term, value in terms
            ],
        )
        cur.executemany(
            q(
                'INSERT INTO {}.entity_evidence_annotation (source_id, '
                'entity_evidence_id, annotation_key) VALUES (%s, %s, %s) '
                'ON CONFLICT DO NOTHING'
            ).as_string(cur.connection),
            [
                (
                    source,
                    endpoint_ids[(key, source, side)],
                    annotation_key(term),
                )
                for key, source, side, term in ROLE_EVIDENCE
            ],
        )

        # The annotations a resource states about an entity rather than about
        # a relation. A conversion's direction is one: it describes the event,
        # and it reaches a consumer by resolving the evidence back to the
        # entity, which is the path the reaction event itself resolves along.
        entity_annotation_ids = {
            (name, source): _uuid5('55555555', index)
            for index, (name, source, _term, _value) in enumerate(
                ENTITY_ANNOTATION,
                start=1,
            )
        }
        cur.executemany(
            q(
                'INSERT INTO {}.entity_evidence (source_id, '
                'entity_evidence_id, dataset_id, row_id, entity_role_id, '
                'entity_type_id, taxonomy_id) '
                'VALUES (%s, %s, %s, %s, 1, %s, %s) ON CONFLICT DO NOTHING'
            ).as_string(cur.connection),
            [
                (
                    source,
                    entity_annotation_ids[(name, source)],
                    source,
                    _ENTITY_ANNOTATION_ROW_BASE + ordinal,
                    entity_type_ids[ENTITY_TYPES.get(name, PROTEIN_TYPE)],
                    TAXONOMY_ID,
                )
                for ordinal, (name, source, _term, _value) in enumerate(
                    ENTITY_ANNOTATION,
                    start=1,
                )
            ],
        )
        # The resolution is what ties the evidence to the canonical entity.
        # Without it the annotation describes a row nobody can reach from the
        # entity it is about.
        cur.executemany(
            q(
                'INSERT INTO {}.entity_evidence_resolution (source_id, '
                'entity_evidence_id, status_id, entity_id) '
                'VALUES (%s, %s, 1, %s) ON CONFLICT DO NOTHING'
            ).as_string(cur.connection),
            [
                (source, entity_annotation_ids[(name, source)], ENTITY[name])
                for name, source, _term, _value in ENTITY_ANNOTATION
            ],
        )
        cur.executemany(
            q(
                'INSERT INTO {}.entity_evidence_annotation (source_id, '
                'entity_evidence_id, annotation_key) VALUES (%s, %s, %s) '
                'ON CONFLICT DO NOTHING'
            ).as_string(cur.connection),
            [
                (
                    source,
                    entity_annotation_ids[(name, source)],
                    annotation_key(term, value),
                )
                for name, source, term, value in ENTITY_ANNOTATION
            ],
        )
    conn.commit()

    return {
        'entities': ENTITY,
        'relations': relation_ids,
        'evidence': evidence_ids,
        'sources': SOURCE_NAMES,
    }
