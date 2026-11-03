"""Per-resource MetSigDB extraction rules (cycle 010).

One rule per v1 resource. The rule says which build-database source it reads,
which extraction shape that source needs, which semantic it contributes, and
how its organism is derived. `contracts/mapping-rules.md` freezes all four.

The rules hold no identifier translation of their own. The core build already
resolved these evidences, and the extraction reads that result, so the
metabolite side stays identical to every other consumer of the canonical entity
layer.
"""

from __future__ import annotations

from dataclasses import dataclass

# The canonical metabolite entity type. Every published membership resolves to
# a `Chemical:OM:0037` entity, and to nothing else.
CHEMICAL_ENTITY_TYPE = 'Chemical:OM:0037'

# The set side is a pathway for three resources and a controlled-vocabulary
# term for two. Named here rather than by numeric id, because the build resolves
# the id from `vocab_entity_type` at run time.
PATHWAY_ENTITY_TYPE = 'Pathway:OM:0014'
CV_TERM_ENTITY_TYPE = 'Cv Term:OM:0012'

# The default projection, in the priority order the row contract publishes. The
# key is the column, the value the `vocab_identifier_type` name behind it.
#
# `inchi` is absent by design: no InChI identifier type exists in the schema.
PROJECTION_IDENTIFIERS: tuple[tuple[str, str], ...] = (
    ('inchikey', 'Standard Inchi Key:MI:1101'),
    ('smiles', 'Smiles:MI:0239'),
    ('hmdb', 'Hmdb:OM:0004'),
    ('pubchem', 'Pubchem Compound:OM:0002'),
    ('chebi', 'Chebi:MI:0474'),
    ('kegg', 'Kegg Compound:MI:2012'),
)


@dataclass(frozen=True)
class ResourceRule:
    """How one v1 resource becomes MetSigDB memberships.

    The organism is not a field here. It is read from the taxonomy the source
    recorded on its set evidence, in ``publish_membership.sql``, which is the
    same lookup for every resource. An earlier version derived it from the
    identifier — 9606 for anything matching ``R-HSA-`` — and that guessed where
    the data already answers, and published null for the 38 species
    WikiPathways covers.
    """

    name: str
    source_name: str
    set_type: str
    set_entity_type: str
    extraction: str
    # ClassyFire alone reads a second source: HMDB assigns the class, and
    # ChemOnt supplies the hierarchy the assignment expands over.
    hierarchy_source_name: str | None = None


REACTOME = ResourceRule(
    name='Reactome',
    source_name='reactome',
    set_type='pathway',
    set_entity_type=PATHWAY_ENTITY_TYPE,
    extraction='extract_onehop.sql',
)

WIKIPATHWAYS = ResourceRule(
    name='WikiPathways',
    source_name='wikipathways',
    set_type='pathway',
    set_entity_type=PATHWAY_ENTITY_TYPE,
    extraction='extract_onehop.sql',
)

# KEGG's global and overview maps: whole-metabolism diagrams rather than
# pathways. Eleven of them carry 2,332 of KEGG's 4,969 memberships, and a hit
# in one says a compound participates in metabolism. The list is KEGG's own
# BRITE grouping, declared here so no row's sub-type comes from a label.
KEGG_OVERVIEW_MAPS: tuple[str, ...] = (
    'rn01100',  # Metabolic pathways
    'rn01110',  # Biosynthesis of secondary metabolites
    'rn01120',  # Microbial metabolism in diverse environments
    'rn01200',  # Carbon metabolism
    'rn01210',  # 2-Oxocarboxylic acid metabolism
    'rn01212',  # Fatty acid metabolism
    'rn01220',  # Degradation of aromatic compounds
    'rn01230',  # Biosynthesis of amino acids
    'rn01232',  # Nucleotide metabolism
    'rn01240',  # Biosynthesis of cofactors
    'rn01250',  # Biosynthesis of nucleotide sugars
)


KEGG = ResourceRule(
    name='KEGG',
    source_name='kegg',
    set_type='pathway',
    set_entity_type=PATHWAY_ENTITY_TYPE,
    # KEGG publishes no compound-to-pathway edge. A compound reaches a pathway
    # through a reaction assigned to it, so this resource alone is two hops.
    extraction='extract_kegg.sql',
)

MACDB = ResourceRule(
    name='MACdb',
    source_name='macdb',
    set_type='disease',
    set_entity_type=CV_TERM_ENTITY_TYPE,
    # MACdb resolves its trait on the evidence side, so the subject needs the
    # same resolution join as the object.
    extraction='extract_macdb.sql',
)

CLASSYFIRE = ResourceRule(
    name='ClassyFire',
    source_name='hmdb',
    set_type='chemical_class',
    set_entity_type=CV_TERM_ENTITY_TYPE,
    extraction='extract_classyfire.sql',
    hierarchy_source_name='chemont',
)

RESOURCES: tuple[ResourceRule, ...] = (
    REACTOME,
    WIKIPATHWAYS,
    KEGG,
    MACDB,
    CLASSYFIRE,
)

_BY_NAME = {rule.name: rule for rule in RESOURCES}


def rule_for(name: str) -> ResourceRule:
    """The rule for one contract resource name.

    Raises ``KeyError`` for anything outside the v1 scope, because the substrate
    constrains `resource` to the same five names.
    """
    return _BY_NAME[name]
