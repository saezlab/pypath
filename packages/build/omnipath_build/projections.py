"""Binary statements a resource's reactions make, projected as ordinary relations.

A reaction is published as the source states it: a ``molecular_activity`` with
its participants (``has_input`` / ``has_output``) and catalysts
(``enabled_by``, or ``catalyzes`` stated apart from it). Some of what it says
is a statement between two of its members, and a consumer asking for pairs
should not have to reassemble the reaction to find it. This module reads those
pairs off the reactions of one resource, in the resource's working tables,
before they are written; they then reach every product as relations of that
resource, with no logic of their own downstream.

**Transport** is the first. A participant that a reaction takes in one
compartment and gives in another moved, and each catalyst of that reaction
moved it: the catalyst ``affects`` the cargo, with the transport aspect and the
relocalization mechanism, which is how TCDB, Cellinker, Mebocost and MRClinksDB
state a transporter and its substrate. The compartments are the GO cellular
components of the participant links (Rhea: its in/out membrane sides) and ride
as ``transport_from`` / ``transport_to``.

**Everything is read per source row**, not per reaction entity: resolution
merges reactions (Reactome's PI3K step in fifty signalling contexts is one
Rhea reaction), and only the row says which participant sat where and which
catalyst it was. A catalyst stated in another row (Reactome's controls) pairs
with its reaction by the source's reaction id that both carry.

Only chemical cargo, and only reactions that name a catalyst: a movement
without one stays a reaction, since inventing a transporter would be inventing
the interaction.
"""

from __future__ import annotations

from omnipath_core import reactions as _r
from omnipath_core.biolink import (
    direction_sign_from_qualifiers,
    entity_type,
    is_symmetric,
    presentation_category,
    schema,
    statement_identity,
)
from omnipath_core.interaction_profiles import TRANSPORT_QUALIFIERS
from omnipath_core.source_attributes import (
    PROJECTED_FROM,
    SOURCE_REACTION,
    TRANSPORT_FROM,
    TRANSPORT_SIDE,
    TRANSPORT_TO,
)

TRANSPORT_DATASET = "projection:transport"
_PREDICATE = "affects"
_EVENT_TYPES = ("molecular_activity", "biological_process")
# Not catalysts: the events themselves, and what describes rather than does.
_NOT_CATALYST_TYPES = (*_EVENT_TYPES, "ontology_class", "pathway")


def _chemical_types() -> tuple[str, ...]:
    """Cargo the projection keeps: chemical entities, not nucleic acids."""
    return tuple(
        sorted(
            entity_type(name)
            for name in schema().class_descendants("chemical entity", reflexive=True)
            if entity_type(name) != "nucleic_acid_entity"
        )
    )


def _literal(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def _list(values) -> str:
    return "(" + ", ".join(_literal(v) for v in values) + ")"


def project_transport(db, first_event_id: int) -> int:
    """Append the transport pairs of the resource's reactions; returns the events added.

    Reads and extends the writer's working tables ``entities``, ``relations``,
    ``relation_annotations`` and ``entity_evidence``.
    """
    annotations = [{"term": term, "value": value} for term, value in TRANSPORT_QUALIFIERS]
    _, predicate, _, qualified = statement_identity("", _PREDICATE, "", annotations)
    sign = direction_sign_from_qualifiers(qualified)
    directed = not is_symmetric(predicate)
    interaction_class = (
        "undirected"
        if not directed
        else {1: "directed_stimulatory", -1: "directed_inhibitory", 0: "directed"}[sign]
    )
    chemical = _list(_chemical_types())
    events = _list(_EVENT_TYPES)
    not_catalyst = _list((*_NOT_CATALYST_TYPES, *_chemical_types()))

    db.execute(
        "CREATE OR REPLACE TEMP TABLE _p_types AS "
        "SELECT entity_key, min(entity_type) AS entity_type, "
        "min(reference_entity_key) AS reference_entity_key FROM entities GROUP BY entity_key"
    )
    # Each participant link of a reaction row, with where the participant is.
    db.execute(
        f"""CREATE OR REPLACE TEMP TABLE _p_member AS
        SELECT r.source, r.dataset, r.row_id, r.subject_entity_key AS reaction,
            r.object_entity_key AS participant, r.predicate,
            coalesce(
                max(a.value) FILTER (WHERE a.term = '{_r.COMPARTMENT}'),
                max(a.value) FILTER (WHERE a.term = '{TRANSPORT_SIDE}')
            ) AS place
        FROM relations r
        JOIN _p_types ev ON ev.entity_key = r.subject_entity_key AND ev.entity_type IN {events}
        JOIN _p_types p ON p.entity_key = r.object_entity_key AND p.entity_type IN {chemical}
        LEFT JOIN relation_annotations a ON a.event_id = r.event_id AND a.scope = 'relation'
        WHERE r.statement_kind = 'relation' AND r.predicate IN ('has_input', 'has_output')
        GROUP BY ALL"""
    )
    db.execute(
        """CREATE OR REPLACE TEMP TABLE _p_moved AS
        SELECT DISTINCT i.source, i.dataset, i.row_id, i.reaction, i.participant AS cargo,
            i.place AS place_from, o.place AS place_to
        FROM _p_member i JOIN _p_member o
          ON o.source = i.source AND o.dataset = i.dataset AND o.row_id = i.row_id
         AND o.reaction = i.reaction AND o.participant = i.participant
        WHERE i.predicate = 'has_input' AND o.predicate = 'has_output'
          AND i.place IS NOT NULL AND o.place IS NOT NULL AND i.place <> o.place"""
    )
    if not db.execute("SELECT count(*) FROM _p_moved").fetchone()[0]:
        _drop(db)
        return 0
    # The source's id of each moving reaction row, where it states one.
    db.execute(
        f"""CREATE OR REPLACE TEMP TABLE _p_row_reaction AS
        SELECT DISTINCT e.item.source AS source, e.item.dataset AS dataset,
            e.item.row_id AS row_id, e.entity_key AS reaction, a.value AS source_reaction
        FROM entity_evidence e, unnest(e.item.annotations) AS t(a)
        WHERE a.term = '{SOURCE_REACTION}'
          AND e.entity_key IN (SELECT reaction FROM _p_moved)"""
    )
    # Catalysts: members of the reaction row itself (Rhea's enzymes, a model's
    # gene rules), or catalysis the resource states apart from it.
    db.execute(
        f"""CREATE OR REPLACE TEMP TABLE _p_catalyst AS
        SELECT r.source, r.dataset, r.row_id, m.reaction,
            CASE WHEN r.subject_entity_key = m.reaction THEN r.object_entity_key
                 ELSE r.subject_entity_key END AS catalyst
        FROM (SELECT DISTINCT source, dataset, row_id, reaction FROM _p_moved) m
        JOIN relations r ON r.source = m.source AND r.dataset = m.dataset AND r.row_id = m.row_id
         AND starts_with(r.upstream_id, r.row_id || ':member:')
         AND (
            (r.predicate = 'enabled_by' AND r.subject_entity_key = m.reaction)
            OR (r.predicate = 'associated_with'
                AND m.reaction IN (r.subject_entity_key, r.object_entity_key))
         )
        UNION
        SELECT m.source, m.dataset, m.row_id, m.reaction, r.subject_entity_key
        FROM (SELECT DISTINCT source, dataset, row_id, reaction FROM _p_moved) m
        JOIN relations r ON r.source = m.source AND r.predicate = 'catalyzes'
         AND r.object_entity_key = m.reaction
        LEFT JOIN (
            SELECT event_id, max(value) AS source_reaction FROM relation_annotations
            WHERE term = '{SOURCE_REACTION}' GROUP BY event_id
        ) stated ON stated.event_id = r.event_id
        LEFT JOIN _p_row_reaction rr ON rr.source = m.source AND rr.dataset = m.dataset
         AND rr.row_id = m.row_id AND rr.reaction = m.reaction
        -- Where both say which source reaction they are about, they must agree.
        WHERE stated.source_reaction IS NULL OR rr.source_reaction IS NULL
           OR stated.source_reaction = rr.source_reaction"""
    )
    db.execute(
        f"""CREATE OR REPLACE TEMP TABLE _p_pair AS
        SELECT m.source, m.dataset, m.row_id, m.reaction, c.catalyst, m.cargo,
            list(DISTINCT m.place_from ORDER BY m.place_from) AS places_from,
            list(DISTINCT m.place_to ORDER BY m.place_to) AS places_to,
            {first_event_id} + row_number() OVER (
                ORDER BY m.source, m.dataset, m.row_id, m.reaction, c.catalyst, m.cargo
            ) - 1 AS event_id
        FROM _p_moved m
        JOIN _p_catalyst c USING (source, dataset, row_id, reaction)
        JOIN _p_types ct ON ct.entity_key = c.catalyst AND ct.entity_type NOT IN {not_catalyst}
        WHERE c.catalyst <> m.cargo
        GROUP BY m.source, m.dataset, m.row_id, m.reaction, c.catalyst, m.cargo"""
    )
    added = db.execute("SELECT count(*) FROM _p_pair").fetchone()[0]
    if not added:
        _drop(db)
        return 0
    qualified_text = str(qualified)
    db.execute(
        f"""INSERT INTO relations BY NAME
        SELECT sha256(p.catalyst || chr(0) || '{predicate}' || chr(0) || p.cargo || chr(0)
                || {_literal(qualified_text)} || chr(0)) AS relation_key,
            'relation' AS statement_kind,
            p.catalyst AS subject_entity_key, p.catalyst AS subject_label, 'protein' AS subject_type,
            '{predicate}' AS predicate,
            p.cargo AS object_entity_key, p.cargo AS object_label, 'protein' AS object_type,
            '' AS taxon, {str(directed).lower()} AS is_directed, {sign} AS sign,
            {_literal(presentation_category(predicate))} AS category,
            '{interaction_class}' AS interaction_class,
            p.event_id, p.source, '{TRANSPORT_DATASET}' AS dataset, p.row_id,
            p.row_id || ':transport' AS upstream_id,
            {_literal(qualified_text)} AS qualified,
            s.reference_entity_key AS subject_reference_entity_key,
            o.reference_entity_key AS object_reference_entity_key
        FROM _p_pair p
        JOIN _p_types s ON s.entity_key = p.catalyst
        JOIN _p_types o ON o.entity_key = p.cargo"""
    )
    qualifiers = " UNION ALL ".join(
        f"SELECT {_literal(term)} AS term, {_literal(value)} AS value, {i} AS position"
        for i, (term, value) in enumerate(TRANSPORT_QUALIFIERS)
    )
    db.execute(
        f"""INSERT INTO relation_annotations BY NAME
        WITH terms AS (
            SELECT p.event_id, p.source, q.term, q.value, q.position
            FROM _p_pair p, ({qualifiers}) q
            UNION ALL SELECT event_id, source, '{TRANSPORT_FROM}', unnest(places_from), 10 FROM _p_pair
            UNION ALL SELECT event_id, source, '{TRANSPORT_TO}', unnest(places_to), 20 FROM _p_pair
            UNION ALL SELECT event_id, source, '{PROJECTED_FROM}', reaction, 30 FROM _p_pair
        )
        SELECT r.relation_key, t.event_id,
            (row_number() OVER (PARTITION BY t.event_id ORDER BY t.position, t.value) - 1) AS ordinal,
            t.term, t.value, NULL AS quantity, t.source, '{TRANSPORT_DATASET}' AS dataset,
            'relation' AS scope
        FROM terms t JOIN (
            SELECT event_id, relation_key FROM relations WHERE dataset = '{TRANSPORT_DATASET}'
        ) r ON r.event_id = t.event_id"""
    )
    _drop(db)
    return added


def _drop(db) -> None:
    for table in ("_p_types", "_p_member", "_p_moved", "_p_row_reaction", "_p_catalyst", "_p_pair"):
        db.execute(f"DROP TABLE IF EXISTS {table}")
