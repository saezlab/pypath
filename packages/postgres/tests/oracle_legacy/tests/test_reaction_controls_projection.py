import duckdb

from omnipath_build.duckdb_load import (
    DuckDBEvidenceProjector,
    _create_duckdb_evidence_tables,
)
from pypath.internals.cv_terms import (
    BiologicalRoleCv,
    ControlEffectCv,
    EntityTypeCv,
    IdentifierNamespaceCv,
    InteractionMetadataCv,
    cv_term_label_accession,
)
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    Identifier,
    Membership,
)


def test_reaction_catalyst_projects_to_joinable_controls_relation():
    con = duckdb.connect(":memory:")
    _create_duckdb_evidence_tables(con)

    reaction = Entity(
        type=EntityTypeCv.REACTION,
        identifiers=[Identifier(IdentifierNamespaceCv.NAME, "A to B")],
        membership=[
            Membership(
                member=Entity(
                    type=EntityTypeCv.CHEMICAL,
                    identifiers=[Identifier(IdentifierNamespaceCv.CHEBI, "CHEBI:1")],
                ),
                annotations=[Annotation(BiologicalRoleCv.REACTANT)],
            ),
            Membership(
                member=Entity(
                    type=EntityTypeCv.CHEMICAL,
                    identifiers=[Identifier(IdentifierNamespaceCv.CHEBI, "CHEBI:2")],
                ),
                annotations=[Annotation(BiologicalRoleCv.PRODUCT)],
            ),
            Membership(
                member=Entity(
                    type=EntityTypeCv.PROTEIN,
                    identifiers=[Identifier(IdentifierNamespaceCv.UNIPROT, "P12345")],
                ),
                annotations=[Annotation(BiologicalRoleCv.CATALYST)],
            ),
        ],
    )

    stats = DuckDBEvidenceProjector(con).project_records(
        [reaction],
        source="test_source",
        dataset="reactions",
    )

    assert stats.entity_evidence == 4
    assert stats.relation_evidence == 3

    missing_endpoints = con.execute(
        """
        SELECT count(*)
        FROM relation_evidence_raw re
        LEFT JOIN entity_evidence_raw se
          ON se.source = re.source
         AND se.entity_evidence_id = re.subject_entity_evidence_id
        LEFT JOIN entity_evidence_raw oe
          ON oe.source = re.source
         AND oe.entity_evidence_id = re.object_entity_evidence_id
        WHERE se.entity_evidence_id IS NULL
           OR oe.entity_evidence_id IS NULL
        """
    ).fetchone()[0]
    assert missing_endpoints == 0

    controls = con.execute(
        """
        SELECT se.entity_type, re.predicate, oe.entity_type
        FROM relation_evidence_raw re
        JOIN entity_evidence_raw se
          ON se.source = re.source
         AND se.entity_evidence_id = re.subject_entity_evidence_id
        JOIN entity_evidence_raw oe
          ON oe.source = re.source
         AND oe.entity_evidence_id = re.object_entity_evidence_id
        WHERE re.predicate = 'controls'
        """
    ).fetchall()
    assert controls == [
        (
            cv_term_label_accession(EntityTypeCv.PROTEIN),
            "controls",
            cv_term_label_accession(EntityTypeCv.REACTION),
        )
    ]

    participant_count = con.execute(
        """
        SELECT count(*)
        FROM relation_evidence_raw
        WHERE predicate = 'has_participant'
        """
    ).fetchone()[0]
    assert participant_count == 2

    joined_pairs = con.execute(
        """
        WITH controls AS (
          SELECT subject_entity_evidence_id AS protein_id,
                 object_entity_evidence_id AS reaction_id
          FROM relation_evidence_raw
          WHERE predicate = 'controls'
        ),
        participants AS (
          SELECT subject_entity_evidence_id AS reaction_id,
                 object_entity_evidence_id AS metabolite_id
          FROM relation_evidence_raw
          WHERE predicate = 'has_participant'
        )
        SELECT count(*)
        FROM controls
        JOIN participants USING (reaction_id)
        """
    ).fetchone()[0]
    assert joined_pairs == 2

    relation_annotations = con.execute(
        """
        SELECT annotation_scope, term, value
        FROM relation_annotation_raw
        ORDER BY annotation_scope, term, value
        """
    ).fetchall()
    assert (
        "relation",
        cv_term_label_accession(InteractionMetadataCv.CONTROL_EFFECT),
        cv_term_label_accession(ControlEffectCv.CATALYSIS),
    ) in relation_annotations
    assert (
        "object",
        cv_term_label_accession(BiologicalRoleCv.CATALYST),
        None,
    ) in relation_annotations
