"""Put the COSMOS labels into the namespaces the formalism reads.

COSMOS wants UniProt on the gene side and ChEBI on the metabolite side. The
build canonicalises an entity to whichever identifier its resources agreed on,
which for a metabolic graph is ChEBI a third of the time, an InChIKey or a
model-local accession the rest, and on the enzyme side an Entrez gene as often
as a UniProt accession. Translating the gap is this module's whole job.

Two things shape how it does that.

**The mappings live in another database.** They are in the utils Postgres, on
another server from the build, so no query joins the two. The stage collects
the identifiers the projection is about to label, pushes that list down through
a second connection in batches, brings the answers back and stages them beside
the build, which is the same shape `emit_build_manifest` uses to read the utils
capabilities. Pulling whole mapping tables is not an option: `id_mapping` holds
618 million rows, and the slice this step needs is a quarter of a million.

**A missing mapping is an answer, not a failure.** Entrez on the gene side and
an InChIKey on the chemical side are identifiers a consumer can work with, so
an entity nobody can translate keeps the identifier the build gave it and the
edge records which namespace that was. The alternative — dropping the node, or
labelling it as though it were in the wanted namespace — would either delete
chemistry the resources published or lie about it. The same rule covers the
utils database being unset or unreachable: every label falls back, the
namespace columns say so, and the projection is complete either way.

The two namespace columns on an edge are what makes the fallback readable. A
consumer that needs UniProt can take the gene endpoints whose namespace says
`uniprot` and leave the rest, without parsing a single label.
"""

from __future__ import annotations

import logging
import re
import time
from dataclasses import dataclass

import psycopg2
import psycopg2.extensions
from psycopg2 import sql
from psycopg2.extras import execute_values

logger = logging.getLogger(__name__)

# Where the mapping tables live in the utils database. A caller reading a copy
# under another name — a staging utils build, a harness that stages a mapping
# of its own — says so rather than editing this.
UTILS_SCHEMA = 'omnipath_utils'

# The staging tables, in the projection's own schema and dropped when the step
# ends, following the convention the interaction derive's `_if_*` tables set:
# unlogged rather than temporary, because an unlogged table uses the shared
# buffer pool and a temporary one exhausts a backend's local buffers.
LABEL_TABLE = '_cos_label'
NAMESPACE_TABLE = '_cos_namespace'

# How many identifiers go down the wire in one statement. Each batch is an
# index scan per identifier on the utils side, so the size trades round trips
# against the parameter array Postgres has to materialise; ten thousand keeps
# both small.
BATCH_SIZE = 10_000

# The namespaces COSMOS asks for, one per side of an edge.
GENE_NAMESPACE = 'uniprot'
CHEMICAL_NAMESPACE = 'chebi'

# What the build calls an entity whose identifier no resolver could place. The
# string is still an identifier — usually the one the resource published — so
# the stage reads its shape rather than trusting the type.
UNRESOLVED_TYPES = ('omnipath:unresolved_entity_key', 'Fallback')

# The build stores a ChEBI identifier as a bare numeral, `15422`. Both the
# utils database and the COSMOS network spell it `CHEBI:15422`, so the prefix
# is added in one place and both the lookup and the label read it from there.
# Getting this wrong is silent: three hundred build identifiers looked up
# unprefixed match nothing at all, and the same three hundred prefixed match a
# third of the time.
CHEBI_PREFIX = 'CHEBI:'

# A UniProt accession, in the pattern the UniProt knowledge base publishes,
# with an isoform suffix allowed. The unresolved population on the gene side is
# mostly accessions that merely failed to be typed as such, and this is what
# tells those apart from the rest of that population rather than assuming.
UNIPROT_ACCESSION = re.compile(
    r'^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]'
    r'|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?$'
)

# A bare numeral, which is what a ChEBI identifier looks like in the build and
# also what plenty of model-local accessions look like. Shape is not evidence,
# so a numeral is prefixed only where the utils database confirms it.
BARE_NUMERAL = re.compile(r'^[0-9]+$')

# The build's vocabulary names the identifier types after the controlled-
# vocabulary term they came from; the edge columns carry a plain namespace
# name. Most translate by stripping the accession and lowering the rest, and
# these are the ones that do not.
NAMESPACE_NAMES = {
    'Uniprot:MI:1097': 'uniprot',
    'Uniprot Trembl:MI:1099': 'uniprot',
    'Entrez:MI:0477': 'entrez',
    'Chebi:MI:0474': 'chebi',
    'Standard Inchi Key:MI:1101': 'inchikey',
    'Gene Name Primary:OM:0200': 'genesymbol',
    'Gene Name Synonym:OM:0201': 'genesymbol_synonym',
    'omnipath:unresolved_entity_key': 'unresolved',
    'Fallback': 'unresolved',
}

# What each of the build's identifier types is called in the utils database. A
# type absent here is one the utils mappings cannot be asked about, which is a
# fallback rather than an error: Human-GEM, MetaNetX and Reactome accessions
# reach ChEBI through no mapping table the utils build holds.
UTILS_ID_TYPES = {
    # The gene side.
    'Uniprot:MI:1097': 'uniprot',
    'Uniprot Trembl:MI:1099': 'uniprot',
    'Entrez:MI:0477': 'entrez',
    'Hgnc:MI:1095': 'hgnc',
    'Gene Name Primary:OM:0200': 'genesymbol',
    'Gene Name Synonym:OM:0201': 'genesymbol-syn',
    'Ensembl:MI:0476': 'ensg',
    'Refseq Protein:OM:0001': 'refseqp',
    # The chemical side.
    'Chebi:MI:0474': 'chebi',
    'Standard Inchi Key:MI:1101': 'inchikey',
    'Pubchem Compound:OM:0002': 'pubchem',
    'Pubchem:MI:0730': 'pubchem',
    'Hmdb:OM:0004': 'hmdb',
    'Kegg Compound:MI:2012': 'kegg',
    'Lipidmaps:OM:0003': 'lipidmaps',
    'Swisslipids:OM:0009': 'swisslipids',
    'Bigg Metabolite:OM:0233': 'bigg',
    'Cas:MI:2011': 'cas',
    'Drugbank:MI:2002': 'drugbank',
    'Chembl Compound:MI:0967': 'chembl',
}

# Everything after the last two colon-separated fields of a vocabulary name is
# the term accession — `Chebi:MI:0474` is the ChEBI type under MI:0474 — and
# the namespace is what is left once it goes.
_TERM_ACCESSION = re.compile(r':[A-Za-z][A-Za-z0-9]*:[0-9]+$')
_NOT_A_NAME = re.compile(r'[^a-z0-9]+')


@dataclass(frozen=True)
class LabelEntity:
    """One entity the projection is about to put into a node label."""

    entity_id: str
    identifier: str
    id_type: str | None
    taxonomy_id: int | None
    # `gene` for a catalyst, `chemical` for a reactant or a product. The role
    # decides it rather than the entity's own type, because the role is what
    # decides which side of a COSMOS edge the label lands on, and it is the
    # side that decides which namespace the label is wanted in.
    side: str


@dataclass(frozen=True)
class TranslatedLabel:
    """The identifier a label carries, and the namespace it belongs to."""

    identifier: str
    namespace: str
    # The utils database supplied the identifier. False where the build's own
    # identifier stands, whether it was already in the wanted namespace or
    # nothing could be found for it.
    mapped: bool


@dataclass(frozen=True)
class TranslationStats:
    """What one translation reached, per side of an edge.

    ``in_namespace`` counts the labels that carry the namespace COSMOS asks for
    — UniProt on the gene side, ChEBI on the chemical one — however they got
    there. ``mapped`` is the subset a utils lookup moved, and the difference
    between the two is the population that was already in the right namespace
    and needed nothing. ``fallback`` is the rest: labels left in the build's
    own namespace, each one recorded as such on the edge.
    """

    gene_entities: int
    gene_in_namespace: int
    gene_mapped: int
    gene_fallback: int
    chemical_entities: int
    chemical_in_namespace: int
    chemical_mapped: int
    chemical_fallback: int
    # The utils database answered. False where no URL was given or the
    # connection failed, which degrades every label to a fallback and is not an
    # error.
    utils_available: bool
    seconds: float


def namespace_name(id_type: str | None) -> str:
    """The plain namespace name behind one of the build's vocabulary names."""
    if not id_type:
        return 'unknown'
    known = NAMESPACE_NAMES.get(id_type)
    if known:
        return known
    stem = _TERM_ACCESSION.sub('', id_type)
    return _NOT_A_NAME.sub('_', stem.lower()).strip('_') or 'unknown'


def chebi_prefixed(identifier: str) -> str:
    """A ChEBI identifier as the utils database and COSMOS both spell it.

    The one place the two spellings meet, and it is used by the lookup and by
    the label alike rather than applied at each site: a lookup that forgot it
    matches nothing and reports a database that answers as one that does not,
    and a label that forgot it cannot be joined to the delivered network.

    Idempotent, because what comes back from a mapping is already prefixed and
    what comes out of the build is not.
    """
    bare = identifier.strip()
    if bare.upper().startswith(CHEBI_PREFIX):
        return CHEBI_PREFIX + bare[len(CHEBI_PREFIX):]
    return CHEBI_PREFIX + bare


class IdentifierMappings:
    """The utils database's mappings, asked one batch of identifiers at a time.

    Every method takes a set of identifiers and pushes it down as an array
    parameter, which is what keeps this a lookup rather than a download: the
    mapping table is indexed on exactly the key each query supplies, so a batch
    costs one index scan per identifier instead of a scan of 618 million rows.
    """

    def __init__(
        self,
        conn: psycopg2.extensions.connection,
        *,
        schema: str = UTILS_SCHEMA,
        batch_size: int = BATCH_SIZE,
    ) -> None:
        self._conn = conn
        self._schema = sql.Identifier(schema)
        self._batch_size = batch_size
        self._type_ids: dict[str, int] | None = None
        # Set when a query against this database raises. The caller reports it
        # rather than the bare fact that a URL was given, so a build whose
        # mapping database answered with an error is not recorded as one that
        # translated and happened to find nothing.
        self.failed = False

    def type_id(self, name: str) -> int | None:
        """The utils database's numeric id for one identifier type.

        Read by name rather than hard coded, so a utils build that renumbers
        its vocabulary does not silently translate into the wrong namespace.
        """
        if self._type_ids is None:
            with self._conn.cursor() as cur:
                cur.execute(
                    sql.SQL('SELECT name, id FROM {}.id_type').format(
                        self._schema
                    )
                )
                self._type_ids = {row[0]: int(row[1]) for row in cur.fetchall()}
        return self._type_ids.get(name)

    def _batches(self, identifiers: list[str]):
        for start in range(0, len(identifiers), self._batch_size):
            yield identifiers[start:start + self._batch_size]

    def map_to(
        self,
        source_type: str,
        target_type: str,
        taxonomy_id: int,
        identifiers: set[str],
    ) -> dict[str, tuple[str, int]]:
        """Translate identifiers of one type, under one taxonomy.

        Answers ``identifier -> (target, how many targets answered)``. The
        count is what lets the caller apply a rule to an ambiguous mapping
        rather than picking silently, and the target is the lowest by sort
        order so that a rebuild reading the same mappings mints the same label.

        Both directions of the mapping table are read. A pair is stored once,
        in whichever direction the resource published it — ChEBI to InChIKey
        exists and InChIKey to ChEBI does not, and the same holds for PubChem —
        so asking forward only would report half the coverage as missing. The
        reverse index covers the reverse question, so neither is a scan.
        """
        source_id = self.type_id(source_type)
        target_id = self.type_id(target_type)
        if source_id is None or target_id is None or not identifiers:
            return {}

        found: dict[str, tuple[str, int]] = {}
        forward = sql.SQL(
            """
            SELECT source_id, min(target_id), count(DISTINCT target_id)
            FROM {}.id_mapping
            WHERE source_type_id = %s
              AND target_type_id = %s
              AND ncbi_tax_id = %s
              AND source_id = ANY(%s)
            GROUP BY source_id
            """
        ).format(self._schema)
        reverse = sql.SQL(
            """
            SELECT target_id, min(source_id), count(DISTINCT source_id)
            FROM {}.id_mapping
            WHERE source_type_id = %s
              AND target_type_id = %s
              AND ncbi_tax_id = %s
              AND target_id = ANY(%s)
            GROUP BY target_id
            """
        ).format(self._schema)

        wanted = sorted(identifiers)
        with self._conn.cursor() as cur:
            for batch in self._batches(wanted):
                cur.execute(
                    forward, (source_id, target_id, taxonomy_id, batch)
                )
                for key, target, count in cur.fetchall():
                    found[key] = (target, int(count))
            missing = sorted(key for key in wanted if key not in found)
            for batch in self._batches(missing):
                cur.execute(
                    reverse, (target_id, source_id, taxonomy_id, batch)
                )
                for key, target, count in cur.fetchall():
                    found[key] = (target, int(count))
        return found

    def resolve_proteins(
        self,
        source_type: str,
        taxonomy_id: int,
        identifiers: set[str],
    ) -> dict[str, tuple[str, int]]:
        """Translate gene-side identifiers through the protein resolver.

        This is the gene side's first question, and the mapping table is asked
        only for what it leaves. Two measurements put it there, both taken over
        the 102,664 Entrez catalysts of the current build.

        **Coverage.** The mapping table's Entrez-to-UniProt pairs cover two
        organisms, mouse and human; the resolver covers 32,838. The catalysts
        of a metabolic build are not two organisms — Arabidopsis, zebrafish,
        rat, fly, worm and yeast together outnumber the human ones — and the
        split is 9,174 the mapping table answers against 99,163 the resolver
        does, 93,489 of them the resolver alone.

        **Agreement.** Where both answer they disagree on 3,478 of 9,174, and
        the disagreement has a direction: the mapping table lists every
        accession an Entrez gene ever reached, so picking one of them by sort
        order returns an unreviewed isoform — Entrez 14810 answers `A2AI14`
        there and `P35438` here. A published node label naming an isoform
        nobody works with is a worse answer than a narrower table.

        Answers in the same shape the mapping table answers in, so the caller
        merges the two without knowing which replied.
        """
        if not identifiers:
            return {}
        statement = sql.SQL(
            """
            SELECT source_id, min(uniprot), count(DISTINCT uniprot)
            FROM {}.resolver_protein
            WHERE ncbi_tax_id = %s
              AND source_type = %s
              AND source_id = ANY(%s)
              AND uniprot IS NOT NULL
            GROUP BY source_id
            """
        ).format(self._schema)
        found: dict[str, tuple[str, int]] = {}
        with self._conn.cursor() as cur:
            for batch in self._batches(sorted(identifiers)):
                cur.execute(statement, (taxonomy_id, source_type, batch))
                for key, target, count in cur.fetchall():
                    found[key] = (target, int(count))
        return found

    def known_chemicals(
        self,
        source_type: str,
        identifiers: set[str],
    ) -> set[str]:
        """Which of these identifiers the utils database knows as ``type``.

        The chemical resolver is the cheapest oracle for the question "is this
        string an identifier of that namespace at all", which is what an
        untyped bare numeral has to be asked before it can be called a ChEBI.
        """
        if not identifiers:
            return set()
        statement = sql.SQL(
            """
            SELECT DISTINCT source_id
            FROM {}.resolver_chemical
            WHERE source_type = %s AND source_id = ANY(%s)
            """
        ).format(self._schema)
        answered: set[str] = set()
        with self._conn.cursor() as cur:
            for batch in self._batches(sorted(identifiers)):
                cur.execute(statement, (source_type, batch))
                answered.update(row[0] for row in cur.fetchall())
        return answered


def collect_label_entities(
    cur: psycopg2.extensions.cursor,
) -> list[LabelEntity]:
    """The entities the projection is about to put into a node label.

    Exactly the party rows the projection reads — a reactant, a product or an
    enzyme — and one row per entity rather than per participation, because an
    enzyme catalysing five hundred reactions carries the same label identifier
    in all five hundred of them.

    An entity that is a catalyst anywhere takes the gene side even where it is
    also a reactant somewhere else. That case is vanishingly rare and the
    choice is arbitrary either way: the rules are driven by the declared type,
    and the side only decides which namespace is worth asking for.
    """
    cur.execute(
        """
        WITH labelled AS (
          SELECT
            party.entity_id,
            bool_or(role.name = 'enzyme') AS catalyst
          FROM interaction_party party
          JOIN vocab_relation_role role
            ON role.relation_role_id = party.role_id
          WHERE role.name IN ('reactant', 'product', 'enzyme')
          GROUP BY party.entity_id
        )
        SELECT
          labelled.entity_id,
          participant.canonical_identifier,
          id_type.name,
          participant.taxonomy_id,
          CASE WHEN labelled.catalyst THEN 'gene' ELSE 'chemical' END
        FROM labelled
        JOIN entity participant
          ON participant.entity_id = labelled.entity_id
        LEFT JOIN vocab_identifier_type id_type
          ON id_type.identifier_type_id
             = participant.canonical_identifier_type_id
        WHERE participant.canonical_identifier IS NOT NULL
        """
    )
    return [
        LabelEntity(
            entity_id=str(entity_id),
            identifier=identifier,
            id_type=id_type,
            taxonomy_id=None if taxonomy_id is None else int(taxonomy_id),
            side=side,
        )
        for entity_id, identifier, id_type, taxonomy_id, side in cur.fetchall()
    ]


def _fallback(entity: LabelEntity) -> TranslatedLabel:
    """The build's own identifier, under the namespace it actually came from.

    The one thing this does not do is take the declared type at its word on the
    gene side. A large part of the unresolved population is UniProt accessions
    that merely failed to be typed as such, and calling those `unresolved` when
    the string is an accession would hide coverage the network already has.
    """
    namespace = namespace_name(entity.id_type)
    if (
        entity.side == 'gene'
        and namespace == 'unresolved'
        and UNIPROT_ACCESSION.match(entity.identifier)
    ):
        namespace = GENE_NAMESPACE
    return TranslatedLabel(entity.identifier, namespace, mapped=False)


# The namespace each side of an edge wants its labels in, which is also the
# target every lookup from that side translates into.
TARGET_NAMESPACE = {'gene': GENE_NAMESPACE, 'chemical': CHEMICAL_NAMESPACE}


def _lookup_key(entity: LabelEntity) -> tuple[str, str, int] | None:
    """Which lookup answers for one entity, and ``None`` where none can.

    The key is the side, the utils database's name for the entity's identifier
    type, and the taxonomy to ask under. Chemical mappings are published
    without an organism and the utils build stores them under taxonomy zero, so
    the chemical side always asks under zero whatever the build's entity
    happens to say. The gene side asks under the entity's own taxonomy, because
    an Entrez gene names a different protein in a different organism and a
    mapping keyed on the wrong one would be wrong rather than merely absent.
    """
    utils_type = UTILS_ID_TYPES.get(entity.id_type or '')
    if utils_type is None:
        return None
    if entity.side == 'gene':
        if entity.taxonomy_id is None:
            return None
        return entity.side, utils_type, entity.taxonomy_id
    return entity.side, utils_type, 0


def _lookup_identifier(entity: LabelEntity) -> str:
    """The identifier as the utils database spells it.

    ChEBI is the only namespace the two databases disagree about, and they
    disagree on the prefix alone.
    """
    if UTILS_ID_TYPES.get(entity.id_type or '') == CHEMICAL_NAMESPACE:
        return chebi_prefixed(entity.identifier)
    return entity.identifier


def _grouped(
    entities: list[LabelEntity],
) -> dict[tuple[str, str, int], set[str]]:
    """The identifiers to ask for, one set per lookup."""
    grouped: dict[tuple[str, str, int], set[str]] = {}
    for entity in entities:
        key = _lookup_key(entity)
        if key is None:
            continue
        grouped.setdefault(key, set()).add(_lookup_identifier(entity))
    return grouped


def _translate_gene(
    entity: LabelEntity,
    mappings: dict[tuple[str, str, int], dict[str, tuple[str, int]]],
) -> TranslatedLabel:
    """One catalyst's label identifier.

    A UniProt accession is already the target and passes through. An Entrez
    gene is looked up, and where several accessions answer the lowest is taken
    rather than two nodes minted: a COSMOS gene node is one enzyme in one
    reaction, and splitting it would state that the reaction runs twice.
    Anything that does not answer keeps what the build holds.
    """
    if namespace_name(entity.id_type) == GENE_NAMESPACE:
        return TranslatedLabel(entity.identifier, GENE_NAMESPACE, mapped=False)

    key = _lookup_key(entity)
    if key is not None:
        hit = mappings.get(key, {}).get(_lookup_identifier(entity))
        if hit is not None:
            return TranslatedLabel(hit[0], GENE_NAMESPACE, mapped=True)
    return _fallback(entity)


def _translate_chemical(
    entity: LabelEntity,
    mappings: dict[tuple[str, str, int], dict[str, tuple[str, int]]],
    known_chebi: set[str],
) -> TranslatedLabel:
    """One reactant's or product's label identifier.

    A ChEBI identifier is already the target and only gains the prefix. Every
    other namespace is looked up, and an ambiguous answer keeps the build's
    identifier rather than picking one: an InChIKey that several ChEBI entries
    claim is a stereochemistry the two databases disagree about, and an
    InChIKey is an identifier a consumer can use, so the honest answer is to
    leave it alone and say `inchikey`.

    An untyped bare numeral is the one case where shape invites a guess. It is
    checked against the ChEBI identifiers the utils database knows and prefixed
    only where one answers; roughly three quarters of that population answers
    nothing, and calling those ChEBI because they are numbers would put a wrong
    identifier into a published label.
    """
    declared = namespace_name(entity.id_type)
    if declared == CHEMICAL_NAMESPACE:
        return TranslatedLabel(
            chebi_prefixed(entity.identifier), CHEMICAL_NAMESPACE, mapped=False
        )

    key = _lookup_key(entity)
    if key is not None:
        hit = mappings.get(key, {}).get(_lookup_identifier(entity))
        if hit is not None and hit[1] == 1:
            return TranslatedLabel(
                chebi_prefixed(hit[0]), CHEMICAL_NAMESPACE, mapped=True
            )
        return _fallback(entity)

    if (
        entity.id_type in UNRESOLVED_TYPES
        and BARE_NUMERAL.match(entity.identifier)
        and chebi_prefixed(entity.identifier) in known_chebi
    ):
        return TranslatedLabel(
            chebi_prefixed(entity.identifier), CHEMICAL_NAMESPACE, mapped=True
        )
    return _fallback(entity)


def _fallback_labels(
    entities: list[LabelEntity],
) -> dict[str, TranslatedLabel]:
    """Every label, with nothing asked of any mapping database.

    A complete answer rather than a failure: each label keeps the identifier
    the build canonicalised it to and each namespace column says which one that
    is. The ChEBI prefix is still applied, because it is a spelling of an
    identifier the build already holds rather than a translation of one it does
    not, and the same goes for reading an untyped accession as the accession it
    plainly is. Neither needs a database.
    """
    return {
        entity.entity_id: (
            TranslatedLabel(
                chebi_prefixed(entity.identifier),
                CHEMICAL_NAMESPACE,
                mapped=False,
            )
            if namespace_name(entity.id_type) == CHEMICAL_NAMESPACE
            else _fallback(entity)
        )
        for entity in entities
    }


def translate_identifiers(
    entities: list[LabelEntity],
    mappings: IdentifierMappings | None,
) -> dict[str, TranslatedLabel]:
    """Every entity's label identifier and the namespace it belongs to.

    ``mappings`` of ``None`` is the degraded path. So is a mapping database
    that answers with an error — an older utils build missing one of the tables
    this reads, a permission the connection does not hold — because the step
    cannot tell those from an absent database and must not treat either as a
    reason to fail a projection of content the build already wrote. Both land
    on the same fallback, and the edge columns record it.
    """
    if mappings is None:
        return _fallback_labels(entities)

    try:
        return _translate_against(entities, mappings)
    except psycopg2.Error as error:
        mappings.failed = True
        logger.warning(
            'cosmos: the identifier mappings could not be read (%s); '
            'every label keeps the identifier the build canonicalised it to',
            error,
        )
        return _fallback_labels(entities)


def _translate_against(
    entities: list[LabelEntity],
    mappings: IdentifierMappings,
) -> dict[str, TranslatedLabel]:
    """The translation proper, with a mapping database that answers."""
    answers: dict[tuple[str, str, int], dict[str, tuple[str, int]]] = {}
    for key, identifiers in _grouped(entities).items():
        side, utils_type, taxonomy = key
        target = TARGET_NAMESPACE[side]
        if utils_type == target:
            # Already in the namespace this side wants. Asking the mapping
            # table to translate UniProt into UniProt is a scan for nothing.
            continue
        if side == 'gene':
            # The protein resolver first, the mapping table for what it leaves.
            # Both the coverage and the quality of the answer put them in that
            # order, and `resolve_proteins` records the measurements.
            found = mappings.resolve_proteins(
                utils_type, taxonomy, identifiers
            )
            unanswered = {
                identifier
                for identifier in identifiers
                if identifier not in found
            }
            found.update(
                mappings.map_to(utils_type, target, taxonomy, unanswered)
            )
        else:
            found = mappings.map_to(utils_type, target, taxonomy, identifiers)
        answers[key] = found

    numerals = {
        chebi_prefixed(entity.identifier)
        for entity in entities
        if entity.side == 'chemical'
        and entity.id_type in UNRESOLVED_TYPES
        and BARE_NUMERAL.match(entity.identifier)
    }
    known_chebi = mappings.known_chemicals(CHEMICAL_NAMESPACE, numerals)

    return {
        entity.entity_id: (
            _translate_gene(entity, answers)
            if entity.side == 'gene'
            else _translate_chemical(entity, answers, known_chebi)
        )
        for entity in entities
    }


def open_mappings(
    utils_db_url: str | None,
    *,
    utils_schema: str = UTILS_SCHEMA,
) -> tuple[psycopg2.extensions.connection | None, IdentifierMappings | None]:
    """Connect to the utils database, or report that there is none.

    Both an unset URL and a refused connection answer ``(None, None)``. The
    projection degrades to the build's own identifiers on that answer, which is
    what lets the step run in an environment that has no utils database at all
    — a test harness, a build machine off the lab network — without the caller
    having to know whether one is reachable before it asks.
    """
    if not utils_db_url:
        return None, None
    try:
        conn = psycopg2.connect(utils_db_url)
    except psycopg2.Error as error:
        logger.warning(
            'cosmos: the identifier mappings are unreachable (%s); '
            'every label keeps the identifier the build canonicalised it to',
            error,
        )
        return None, None
    # Read only, and stated as such: nothing in this module writes to the utils
    # database, and a transaction that says so cannot start to.
    conn.set_session(readonly=True, autocommit=True)
    return conn, IdentifierMappings(conn, schema=utils_schema)


def _stage_namespaces(cur: psycopg2.extensions.cursor) -> None:
    """The plain namespace name behind every identifier type in the build.

    Small enough to stage whole. It exists so that the endpoints the
    translation does not touch — the pseudo-enzyme of a reaction nobody named a
    catalyst for, whose label is the reaction's own identifier — report their
    namespace in the same vocabulary as every other endpoint, rather than one
    column holding `chebi` and the next `Rhea Id:OM:0239`.
    """
    cur.execute(f'DROP TABLE IF EXISTS {NAMESPACE_TABLE}')
    cur.execute(
        f"""
        CREATE UNLOGGED TABLE {NAMESPACE_TABLE} (
          identifier_type_id integer PRIMARY KEY,
          namespace text NOT NULL
        )
        """
    )
    cur.execute('SELECT identifier_type_id, name FROM vocab_identifier_type')
    rows = [
        (int(type_id), namespace_name(name)) for type_id, name in cur.fetchall()
    ]
    if rows:
        execute_values(
            cur,
            f'INSERT INTO {NAMESPACE_TABLE} '
            '(identifier_type_id, namespace) VALUES %s',
            rows,
        )
    cur.execute(f'ANALYZE {NAMESPACE_TABLE}')


def stage_label_translation(
    cur: psycopg2.extensions.cursor,
    *,
    utils_db_url: str | None = None,
    utils_schema: str = UTILS_SCHEMA,
    progress: bool = False,
) -> TranslationStats:
    """Stage the label identifier and namespace of every entity to be labelled.

    The caller has already pointed the search path at the projection's schema,
    so the staging tables land beside the edge table and the projection SQL
    joins them unqualified. They are dropped when the step ends.

    Returns what the translation reached, per side, which the projection folds
    into its own statistics.
    """
    started = time.monotonic()
    _stage_namespaces(cur)
    entities = collect_label_entities(cur)

    utils_conn, mappings = open_mappings(
        utils_db_url, utils_schema=utils_schema
    )
    try:
        labels = translate_identifiers(entities, mappings)
    finally:
        if utils_conn is not None:
            utils_conn.close()

    cur.execute(f'DROP TABLE IF EXISTS {LABEL_TABLE}')
    cur.execute(
        f"""
        CREATE UNLOGGED TABLE {LABEL_TABLE} (
          entity_id uuid NOT NULL,
          identifier text NOT NULL,
          id_namespace text NOT NULL
        )
        """
    )
    if labels:
        execute_values(
            cur,
            f'INSERT INTO {LABEL_TABLE} '
            '(entity_id, identifier, id_namespace) VALUES %s',
            [
                (entity_id, label.identifier, label.namespace)
                for entity_id, label in labels.items()
            ],
        )
    cur.execute(
        f'CREATE INDEX {LABEL_TABLE}_idx ON {LABEL_TABLE} (entity_id)'
    )
    cur.execute(f'ANALYZE {LABEL_TABLE}')

    counted = {'gene': [0, 0, 0], 'chemical': [0, 0, 0]}
    for entity in entities:
        label = labels[entity.entity_id]
        tally = counted[entity.side]
        tally[0] += 1
        if label.namespace == TARGET_NAMESPACE[entity.side]:
            tally[1] += 1
        if label.mapped:
            tally[2] += 1

    stats = TranslationStats(
        gene_entities=counted['gene'][0],
        gene_in_namespace=counted['gene'][1],
        gene_mapped=counted['gene'][2],
        gene_fallback=counted['gene'][0] - counted['gene'][1],
        chemical_entities=counted['chemical'][0],
        chemical_in_namespace=counted['chemical'][1],
        chemical_mapped=counted['chemical'][2],
        chemical_fallback=counted['chemical'][0] - counted['chemical'][1],
        utils_available=mappings is not None and not mappings.failed,
        seconds=round(time.monotonic() - started, 1),
    )
    if progress:
        logger.info(
            '[cosmos] event=cosmos_translation utils=%s '
            'gene_entities=%s gene_uniprot=%s gene_mapped=%s '
            'chemical_entities=%s chemical_chebi=%s chemical_mapped=%s '
            'seconds=%s',
            stats.utils_available,
            stats.gene_entities, stats.gene_in_namespace, stats.gene_mapped,
            stats.chemical_entities, stats.chemical_in_namespace,
            stats.chemical_mapped, stats.seconds,
        )
    return stats


def drop_label_translation(cur: psycopg2.extensions.cursor) -> None:
    """Drop the staging tables, as the interaction derive drops its own."""
    cur.execute(f'DROP TABLE IF EXISTS {LABEL_TABLE}')
    cur.execute(f'DROP TABLE IF EXISTS {NAMESPACE_TABLE}')
