"""Project the reaction stars into the COSMOS prior-knowledge network.

The step reads what the interactions phase already wrote. The N-ary headers in
``interaction`` and their ``interaction_party`` rows hold the reactants, the
products and the catalysts, each in the role its resources published, so the
projection needs no upstream parse and no download: it turns one reaction into
the binary metabolite/enzyme edges COSMOS expects, and nothing else.

Its output is a table of its own. A COSMOS gene node is one enzyme in one
reaction, so an enzyme catalysing five hundred reactions is five hundred nodes;
those are artifacts of the formalism rather than entities of the graph, and
writing them into the canonical layer would move every count the rest of the
build reports.
"""

from __future__ import annotations

import logging
import time
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import psycopg2.extensions
from psycopg2 import sql

from omnipath_build.cosmos.translate import (
    UTILS_SCHEMA,
    drop_label_translation,
    stage_label_translation,
)

_SQL_DIR = Path(__file__).with_name('sql')

logger = logging.getLogger(__name__)

TABLE = 'cosmos_edge'

# The entity types an event node carries, mirroring what the interaction
# projection treats as the parent of a reaction star. Restated here rather than
# imported, because the two answer different questions of the same vocabulary
# and a shared private constant would tie this step to the derive's internals.
REACTION_ENTITY_TYPES = ('Reaction:OM:0015', 'Transport:OM:0035')

# A transport event moves a molecule rather than transforming it, and COSMOS
# distinguishes the two.
TRANSPORT_ENTITY_TYPE = 'Transport:OM:0035'

# The annotation that says which way a reaction runs. It sits on the event
# entity rather than on any of its participant relations, and both the
# genome-scale-model and the Recon inputs declare it.
DIRECTION_TERM = 'Conversion Direction:OM:1211'

# The resources whose reactions this network is made of. A positive list
# rather than an exclusion, so a reaction resource loaded later has to be
# named here before it reaches the output instead of arriving unannounced.
#
# **Reactome is deliberately absent.** It is a general reaction database rather
# than a metabolic reconstruction: its events are complex assembly, binding and
# modification, so 16,743 of its 18,903 participants are complexes, genes,
# protein families, DNA or RNA rather than small molecules, and it states no
# catalysis at all — it holds no `controls` edge, so every one of its reactions
# would reach the output through the unknown-enzyme path. The four named here
# contribute chemicals and nothing else, on every participant of every reaction
# they hold.
#
# A reaction two resources report keeps its place as long as one of them is
# named here, because the gate is an overlap rather than a containment: the 532
# reactions Reactome and Rhea both describe stay, on Rhea's word.
REACTION_SOURCES = ('kegg', 'rhea', 'metatlas', 'recon3d')

# The participant roles the projection reads. The rest of the vocabulary —
# cofactor, regulator, member — is counted and left out, for the reason the
# projection query records.
SKIPPED_ROLES = ('cofactor', 'regulator', 'member')


@dataclass(frozen=True)
class CosmosBuildStats:
    """What one projection published, and what it declined to publish."""

    build_id: str
    edges: int
    reactions: int
    orphan_reactions: int
    connectors: int
    reverse_edges: int
    reversible_reactions: int
    metabolite_nodes: int
    gene_nodes: int
    skipped_parties: int
    # Reaction headers that carry side roles but no reachable event entity.
    # Nothing labels them, so they are not projected; a non-zero figure means
    # the record lost a star its header kept, and is worth looking at.
    headers_without_event: int
    # What the label translation reached, counted over entities rather than
    # over edges: one enzyme in five hundred reactions is five hundred nodes
    # carrying one identifier, and the question these answer is how much of the
    # build reaches the namespaces COSMOS asks for.
    #
    # `translated` is a label that carries the wanted namespace -- UniProt on
    # the gene side, ChEBI on the chemical one -- whether a mapping put it
    # there or it was already in it. `mapped` is the subset a mapping moved.
    # `fallback` is the rest, each one recorded on its edges under the
    # namespace it actually came from.
    gene_labels_translated: int
    gene_labels_mapped: int
    gene_labels_fallback: int
    chemical_labels_translated: int
    chemical_labels_mapped: int
    chemical_labels_fallback: int
    # A mapping database answered. False where none was given or the
    # connection failed, which falls every label back and is not an error.
    identifier_translation: bool
    seconds: float


def _sql_text(name: str) -> str:
    return (_SQL_DIR / name).read_text(encoding='utf-8')


def _log(progress: bool, event: str, **fields: object) -> None:
    """One structured progress line, in the shape the derive's output uses."""
    if not progress:
        return
    details = ' '.join(f'{key}={value}' for key, value in fields.items())
    logger.info('[cosmos] event=%s%s', event, f' {details}' if details else '')


def _scalar(cur, query: str, params=None):
    cur.execute(query, params)
    row = cur.fetchone()
    return None if row is None else row[0]


def ensure_cosmos_edge_table(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
) -> None:
    """Create the COSMOS edge table and its indexes if they are absent.

    The DDL is idempotent. Applying it to a populated table leaves both the
    schema and the rows untouched, so a rebuild can call it first without a
    guard.
    """
    with conn.cursor() as cur:
        # Unqualified names in the DDL land in the target schema. public.*
        # still resolves, so the table can reference the canonical layer.
        cur.execute(
            sql.SQL('SET search_path = {}, public').format(sql.Identifier(schema))
        )
        cur.execute(_sql_text('cosmos_edge_table.sql'))
        cur.execute('RESET search_path')
    conn.commit()
    logger.info('cosmos: %s.%s is present', schema, TABLE)


def build_id(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
) -> str:
    """The build the edges are stamped with, read from the manifest.

    An edge records which build produced it, so a consumer can tell a stale row
    from a current one. The manifest holds exactly one row.

    Read through the target schema's search path rather than qualified,
    deliberately: only the derive's manifest step creates the table, so a
    namespace that holds a projection but no manifest of its own — a throwaway
    one built from the canonical schema — takes the stamp of the database it
    lives in instead of failing for want of a row it was never going to have.
    """
    with conn.cursor() as cur:
        cur.execute(
            sql.SQL('SET search_path = {}, public').format(sql.Identifier(schema))
        )
        stamp = _scalar(cur, 'SELECT build_id FROM build_manifest')
        cur.execute('RESET search_path')
    if stamp is None:
        raise LookupError('cosmos: build_manifest is empty; run the derive first')
    return stamp


def _node_counts(cur) -> tuple[int, int]:
    """Distinct metabolite and gene nodes across the reaction edges.

    Counted over both endpoints of every edge, because a metabolite that is
    only ever produced never appears on the source side and a gene node of a
    reaction with no products never appears on the target side.
    """
    cur.execute(
        f"""
        SELECT
          count(DISTINCT node.label)
            FILTER (WHERE node.node_type = 'metabolite'),
          count(DISTINCT node.label)
            FILTER (WHERE node.node_type <> 'metabolite')
        FROM {TABLE} edge
        CROSS JOIN LATERAL (
          VALUES
            (edge.source_label, edge.source_type),
            (edge.target_label, edge.target_type)
        ) AS node (label, node_type)
        WHERE edge.interaction_type <> 'connector'
        """
    )
    metabolites, genes = cur.fetchone()
    return int(metabolites), int(genes)


def _skipped_parties(cur) -> int:
    """Participants of a projected reaction that the projection left out.

    Restricted to the reactions that reached the table, and that restriction is
    load bearing: `member` is also the role a binary header gives an entity
    that appears on both of its ends, and counting those would report thirteen
    million skipped participants of reactions that do not exist.
    """
    return int(
        _scalar(
            cur,
            f"""
            SELECT count(*)
            FROM interaction_party party
            JOIN vocab_relation_role role
              ON role.relation_role_id = party.role_id
            WHERE role.name = ANY(%(skipped_roles)s)
              AND EXISTS (
                SELECT 1 FROM {TABLE} edge
                WHERE edge.interaction_id = party.interaction_id
              )
            """,
            {'skipped_roles': list(SKIPPED_ROLES)},
        )
    )


def _headers_without_event(cur) -> int:
    """Reaction headers the projection could not label.

    A header reaches the table only through its event entity, which is what
    carries the reaction's own identifier and what the reaction index ranks on.
    The question is asked from the participants inward, because the alternative
    — scanning every header for a side role — walks fourteen million rows to
    find a few tens of thousands.
    """
    return int(
        _scalar(
            cur,
            f"""
            WITH reaction_header AS (
              SELECT DISTINCT party.interaction_id
              FROM interaction_party party
              JOIN vocab_relation_role role
                ON role.relation_role_id = party.role_id
              WHERE role.name IN ('reactant', 'product')
            )
            SELECT count(*)
            FROM reaction_header header
            WHERE NOT EXISTS (
              SELECT 1 FROM {TABLE} edge
              WHERE edge.interaction_id = header.interaction_id
            )
            """,
        )
    )


def build_cosmos_projection(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
    utils_db_url: str | None = None,
    utils_schema: str = UTILS_SCHEMA,
    reaction_sources: Sequence[str] = REACTION_SOURCES,
    progress: bool = False,
) -> CosmosBuildStats:
    """Rebuild the whole COSMOS edge table from the reaction stars.

    The table is a projection, so it is replaced whole rather than merged: a
    reaction that lost a participant upstream has to lose the edges that
    participant carried, and the reaction index is ranked over the reactions
    this build sees, so a surviving row from an earlier one would carry a label
    from a different numbering.

    ``utils_db_url`` points at the identifier mappings, which live in another
    database on another server and are read through a second connection. It is
    optional, and so is the database answering: a step that cannot reach it
    labels every node with the identifier the build canonicalised it to and
    records that namespace on the edge. The output is then narrower, never
    absent, which is what lets this run on a machine that has no mapping
    database at all.
    """
    started = time.monotonic()
    ensure_cosmos_edge_table(conn, schema=schema)
    stamp = build_id(conn, schema=schema)
    _log(progress, 'cosmos_start', build=stamp, schema=schema)

    with conn.cursor() as cur:
        cur.execute(
            sql.SQL('SET search_path = {}, public').format(sql.Identifier(schema))
        )
        # RESTART IDENTITY so the surrogate does not climb by the size of the
        # projection on every rebuild.
        cur.execute(
            sql.SQL('TRUNCATE {}.{} RESTART IDENTITY').format(
                sql.Identifier(schema), sql.Identifier(TABLE)
            )
        )

        # The labels come out of the staging tables this leaves behind, so the
        # translation runs before either projection statement and the drop at
        # the end of the step is what takes them away again.
        translation = stage_label_translation(
            cur,
            utils_db_url=utils_db_url,
            utils_schema=utils_schema,
            progress=progress,
        )
        _log(
            progress,
            'cosmos_translation',
            utils=translation.utils_available,
            gene_uniprot=translation.gene_in_namespace,
            gene_fallback=translation.gene_fallback,
            chemical_chebi=translation.chemical_in_namespace,
            chemical_fallback=translation.chemical_fallback,
        )

        cur.execute(
            _sql_text('project_edges.sql'),
            {
                'build_id': stamp,
                'reaction_entity_types': list(REACTION_ENTITY_TYPES),
                'reaction_sources': list(reaction_sources),
                'transport_entity_type': TRANSPORT_ENTITY_TYPE,
                'direction_term': DIRECTION_TERM,
            },
        )
        edges = int(cur.rowcount)
        _log(progress, 'cosmos_edges', edges=edges)

        # The connectors read the edges back, so they are a second statement
        # rather than a branch of the first: the gene nodes they attach to are
        # exactly the distinct gene endpoints the projection just wrote, and
        # asking the table is both shorter and impossible to get out of step
        # with.
        cur.execute(_sql_text('project_connectors.sql'), {'build_id': stamp})
        connectors = int(cur.rowcount)
        _log(progress, 'cosmos_connectors', connectors=connectors)

        cur.execute(
            f"""
            SELECT
              count(DISTINCT reaction_entity_id),
              count(DISTINCT reaction_entity_id) FILTER (WHERE orphan),
              count(*) FILTER (WHERE reverse),
              count(DISTINCT reaction_entity_id)
                FILTER (WHERE direction = 'reversible')
            FROM {TABLE}
            WHERE interaction_type <> 'connector'
            """
        )
        (
            reactions,
            orphan_reactions,
            reverse_edges,
            reversible_reactions,
        ) = cur.fetchone()
        metabolite_nodes, gene_nodes = _node_counts(cur)
        skipped_parties = _skipped_parties(cur)
        headers_without_event = _headers_without_event(cur)

        drop_label_translation(cur)
        cur.execute('RESET search_path')
    conn.commit()

    stats = CosmosBuildStats(
        build_id=stamp,
        edges=edges,
        reactions=int(reactions),
        orphan_reactions=int(orphan_reactions),
        connectors=connectors,
        reverse_edges=int(reverse_edges),
        reversible_reactions=int(reversible_reactions),
        metabolite_nodes=metabolite_nodes,
        gene_nodes=gene_nodes,
        skipped_parties=skipped_parties,
        headers_without_event=headers_without_event,
        gene_labels_translated=translation.gene_in_namespace,
        gene_labels_mapped=translation.gene_mapped,
        gene_labels_fallback=translation.gene_fallback,
        chemical_labels_translated=translation.chemical_in_namespace,
        chemical_labels_mapped=translation.chemical_mapped,
        chemical_labels_fallback=translation.chemical_fallback,
        identifier_translation=translation.utils_available,
        seconds=round(time.monotonic() - started, 1),
    )
    _log(
        progress,
        'cosmos_done',
        edges=stats.edges,
        reactions=stats.reactions,
        orphan_reactions=stats.orphan_reactions,
        connectors=stats.connectors,
        reverse_edges=stats.reverse_edges,
        reversible_reactions=stats.reversible_reactions,
        metabolite_nodes=stats.metabolite_nodes,
        gene_nodes=stats.gene_nodes,
        skipped_parties=stats.skipped_parties,
        headers_without_event=stats.headers_without_event,
        gene_labels_translated=stats.gene_labels_translated,
        gene_labels_mapped=stats.gene_labels_mapped,
        gene_labels_fallback=stats.gene_labels_fallback,
        chemical_labels_translated=stats.chemical_labels_translated,
        chemical_labels_mapped=stats.chemical_labels_mapped,
        chemical_labels_fallback=stats.chemical_labels_fallback,
        identifier_translation=stats.identifier_translation,
        seconds=stats.seconds,
    )
    logger.info(
        'cosmos: build=%s edges=%s reactions=%s orphans=%s reversible=%s '
        'reverse_edges=%s connectors=%s metabolites=%s genes=%s skipped=%s '
        'unlabelled=%s translation=%s gene_uniprot=%s gene_fallback=%s '
        'chemical_chebi=%s chemical_fallback=%s seconds=%s',
        stats.build_id, stats.edges, stats.reactions, stats.orphan_reactions,
        stats.reversible_reactions, stats.reverse_edges, stats.connectors,
        stats.metabolite_nodes, stats.gene_nodes, stats.skipped_parties,
        stats.headers_without_event, stats.identifier_translation,
        stats.gene_labels_translated, stats.gene_labels_fallback,
        stats.chemical_labels_translated, stats.chemical_labels_fallback,
        stats.seconds,
    )
    return stats
