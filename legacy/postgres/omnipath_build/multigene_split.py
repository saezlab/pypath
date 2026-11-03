"""Multi-gene protein split.

The rare *identical gene copies* case: one UniProt accession maps to several
NCBI (Entrez) genes that produce the identical protein. Per the gene-anchored
principle the record must be **split into one gene-anchored record per gene** —
no anchor chosen, nothing dropped; the protein-centric view re-collapses them.

We realise the split the way protein-centric OmniPath did: keep every record
**1:1** (one evidence → one entity) and *duplicate* the multi-gene mention's
sub-graph, one copy per gene, each copy carrying its own provenance. This avoids
a 1→N evidence→entity cardinality change (and the PK/relation redesign that
would ripple through the whole pipeline); the existing 1:1 resolver / relation /
state machinery then handles each copy unchanged.

The explosion runs in DuckDB **before** ``entity_resolution_base`` (it needs
``needed_resolver_lookup`` to know the candidate genes). It:

1. finds protein mentions whose **UniProt** id resolves, per taxon, to >1 Entrez
   gene (``multigene_split``: one target gene per copy, fresh suffixed ids);
2. duplicates every raw row that references such a mention — the mention itself,
   its identifiers/annotations, and the relations/ontology-relations it
   participates in (both endpoints → cross-product, with regenerated
   ``relation_evidence_id``);
3. emits ``multigene_resolution`` so ``entity_resolution_base`` resolves each
   copy **directly** to its assigned gene (bypassing the ``candidate_count > 1``
   → *unresolved* branch). The retained UniProt then yields a per-gene protein
   ``state`` via the existing evidence-state logic, so the same AC ends up
   under each gene.

Restricted to the UniProt key type: the build collapses Ensembl gene/protein
into one ``Ensembl`` id-type, so an Ensembl multi-map cannot be told apart from
gene-level ambiguity; gene-symbol / Ensembl multi-maps are id *ambiguity*, not
identical copies, and must NOT be split (it would fabricate edges).
"""

from __future__ import annotations

from pypath.internals.cv_terms import (
    IdentifierNamespaceCv,
    cv_term_label_accession,
)

from omnipath_build.cv_terms import GENE_ENTITY_TYPE

UNIPROT_TYPE = cv_term_label_accession(IdentifierNamespaceCv.UNIPROT)
ENTREZ_TYPE = cv_term_label_accession(IdentifierNamespaceCv.ENTREZ)

def _lit(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def _content_uuid(expr: str) -> str:
    """SQL for a deterministic UUID from a string ``expr`` — mirrors the build's
    ``content_uuid`` macro (md5 → 8-4-4-4-12). Every evidence id is cast
    ``::UUID`` on the Postgres copy, so a split id MUST be a valid UUID, not a
    suffixed string. Kept inline (not the macro) so the module is self-contained
    and unit-testable on a bare DuckDB."""
    # Returns TEXT (raw evidence-id columns are VARCHAR; the ::UUID cast happens
    # only at the Postgres copy — the md5-hyphenated text is a valid UUID there).
    m = f'md5({expr})'
    h = chr(39) + '-' + chr(39)  # SQL literal '-'
    return (
        f'(substr({m},1,8)||{h}||substr({m},9,4)||{h}'
        f'||substr({m},13,4)||{h}||substr({m},17,4)||{h}'
        f'||substr({m},21,12))'
    )


def explode_multi_gene_protein_mentions(con, *, log=lambda *_: None) -> int:
    """Duplicate multi-gene UniProt mentions 1:1 per gene; return #copies added.

    Operates on the raw tables (``entity_evidence_raw`` etc.) + the prebuilt
    ``needed_resolver_lookup``; leaves ``multigene_split`` / ``multigene_resolution``
    behind for ``entity_resolution_base``. A no-op (returns 0) when no mention
    maps to multiple genes.
    """
    uniprot = _lit(UNIPROT_TYPE)

    # 1) Detect multi-gene UniProt mentions; one target-gene row per copy.
    con.execute(
        f"""
        CREATE OR REPLACE TABLE multigene_split AS
        WITH uniprot_type AS (
          SELECT identifier_type_id FROM identifier_type_all WHERE name = {uniprot}
        ),
        mention_gene AS (
          SELECT DISTINCT
            ee.source,
            ee.entity_evidence_id,
            coalesce(rl.taxonomy_id, ee.taxonomy_id) AS taxonomy_id,
            rl.canonical_identifier_type_id AS entrez_type_id,
            rl.canonical_identifier AS entrez
          FROM entity_evidence_raw ee
          JOIN entity_identifier_raw ei
            ON ei.source = ee.source
           AND ei.entity_evidence_id = ee.entity_evidence_id
          JOIN uniprot_type ut ON TRUE
          JOIN identifier_type_all kit
            ON kit.identifier_type_id = ut.identifier_type_id
           AND kit.name = ei.identifier_type
          JOIN needed_resolver_lookup rl
            ON rl.key_identifier_type_id = ut.identifier_type_id
           AND rl.key_value = ei.identifier
           AND rl.evidence_entity_type = ee.entity_type
           AND (
             rl.taxonomy_id = ee.taxonomy_id
             OR rl.taxonomy_id IS NULL
             OR rl.taxonomy_optional_match
           )
        ),
        multi AS (
          SELECT source, entity_evidence_id
          FROM mention_gene
          GROUP BY source, entity_evidence_id
          HAVING count(DISTINCT entrez) > 1
        )
        SELECT
          mg.source,
          mg.entity_evidence_id AS orig_entity_evidence_id,
          mg.taxonomy_id,
          mg.entrez_type_id,
          mg.entrez,
          {_content_uuid("mg.entity_evidence_id || '#mg=' || mg.entrez")}
            AS new_entity_evidence_id
        FROM mention_gene mg
        JOIN multi USING (source, entity_evidence_id)
        """
    )
    copies = con.execute('SELECT count(*) FROM multigene_split').fetchone()[0]
    if not copies:
        # Still expose the (empty) resolution table — same shape as the
        # populated one below — so entity_resolution_base can join it.
        con.execute(
            f"""
            CREATE OR REPLACE TABLE multigene_resolution AS
            SELECT
              source,
              new_entity_evidence_id AS entity_evidence_id,
              {_lit(GENE_ENTITY_TYPE)} AS entity_type,
              taxonomy_id,
              entrez_type_id AS canonical_identifier_type_id,
              entrez AS canonical_identifier
            FROM multigene_split
            """
        )
        return 0

    log(f'multigene split: {copies} per-gene copies')

    # 2) Explode the entity-side raw tables (mention, identifiers, annotations).
    #    Original multi-gene rows are dropped; replaced by one copy per gene.
    _explode_one(
        con, 'entity_evidence_raw', 'entity_evidence_id',
        cols=[
            'source', 'dataset', 'row_id',
            ('entity_evidence_id', 'new_entity_evidence_id'),
            'parent_entity_evidence_id', 'entity_role', 'entity_type',
            'taxonomy_id',
        ],
    )
    _explode_one(
        con, 'entity_identifier_raw', 'entity_evidence_id',
        cols=[
            'source', ('entity_evidence_id', 'new_entity_evidence_id'),
            'identifier_id', 'identifier_type', 'identifier',
            'identifier_normalized',
        ],
    )
    _explode_one(
        con, 'entity_annotation_raw', 'evidence_id',
        cols=[
            'source', ('evidence_id', 'new_entity_evidence_id'),
            'annotation_key', 'term', 'value', 'unit',
        ],
    )

    # 3) Explode relations. Both endpoints may be multi-gene → cross-product;
    #    regenerate relation_evidence_id from the (possibly) new endpoints.
    con.execute(
        f"""
        CREATE OR REPLACE TABLE relation_evidence_raw_mg AS
        SELECT
          r.source, r.dataset, r.row_id,
          CASE WHEN sm.new_entity_evidence_id IS NULL
                AND om.new_entity_evidence_id IS NULL
               THEN r.relation_evidence_id
               ELSE {_content_uuid("r.relation_evidence_id || '#mgs=' || coalesce(sm.entrez,'') || '#mgo=' || coalesce(om.entrez,'')")}
          END AS relation_evidence_id,
          r.relation_evidence_id AS orig_relation_evidence_id,
          coalesce(sm.new_entity_evidence_id, r.subject_entity_evidence_id)
            AS subject_entity_evidence_id,
          r.predicate,
          coalesce(om.new_entity_evidence_id, r.object_entity_evidence_id)
            AS object_entity_evidence_id,
          r.relation_category
        FROM relation_evidence_raw r
        LEFT JOIN multigene_split sm
          ON sm.source = r.source
         AND sm.orig_entity_evidence_id = r.subject_entity_evidence_id
        LEFT JOIN multigene_split om
          ON om.source = r.source
         AND om.orig_entity_evidence_id = r.object_entity_evidence_id
        """
    )
    # Map old→new relation ids for the relation-keyed annotation tables.
    con.execute(
        """
        CREATE OR REPLACE TABLE relation_id_map AS
        SELECT DISTINCT source, orig_relation_evidence_id, relation_evidence_id
        FROM relation_evidence_raw_mg
        WHERE relation_evidence_id <> orig_relation_evidence_id
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TABLE relation_evidence_raw AS
        SELECT source, dataset, row_id, relation_evidence_id,
               subject_entity_evidence_id, predicate,
               object_entity_evidence_id, relation_category
        FROM relation_evidence_raw_mg
        """
    )
    con.execute('DROP TABLE relation_evidence_raw_mg')

    # relation_annotation_raw is keyed on the relation id → follow the remap.
    con.execute(
        """
        CREATE OR REPLACE TABLE relation_annotation_raw AS
        SELECT ra.source, ra.evidence_id, ra.annotation_key,
               ra.annotation_scope, ra.term, ra.value, ra.unit
        FROM relation_annotation_raw ra
        WHERE NOT EXISTS (
          SELECT 1 FROM relation_id_map m
          WHERE m.source = ra.source
            AND m.orig_relation_evidence_id = ra.evidence_id
        )
        UNION ALL
        SELECT ra.source, m.relation_evidence_id, ra.annotation_key,
               ra.annotation_scope, ra.term, ra.value, ra.unit
        FROM relation_annotation_raw ra
        JOIN relation_id_map m
          ON m.source = ra.source
         AND m.orig_relation_evidence_id = ra.evidence_id
        """
    )

    # annotation_relation_evidence_raw: subject may be multi-gene; object is by
    # id (not an evidence id). Regenerate its relation id from the subject only.
    con.execute(
        f"""
        CREATE OR REPLACE TABLE annotation_relation_evidence_raw AS
        SELECT
          CASE WHEN sm.new_entity_evidence_id IS NULL
               THEN ar.relation_evidence_id
               ELSE {_content_uuid("ar.relation_evidence_id || '#mgs=' || sm.entrez")}
          END AS relation_evidence_id,
          ar.source, ar.dataset, ar.row_id,
          coalesce(sm.new_entity_evidence_id, ar.subject_entity_evidence_id)
            AS subject_entity_evidence_id,
          ar.predicate, ar.object_entity_type, ar.object_id_type,
          ar.object_id, ar.relation_category
        FROM annotation_relation_evidence_raw ar
        LEFT JOIN multigene_split sm
          ON sm.source = ar.source
         AND sm.orig_entity_evidence_id = ar.subject_entity_evidence_id
        """
    )

    # ontology_relation_raw: subject by evidence id (object is by identifier).
    _explode_one(
        con, 'ontology_relation_raw', 'subject_entity_evidence_id',
        cols=[
            'source', 'dataset',
            ('subject_entity_evidence_id', 'new_entity_evidence_id'),
            'ontology_id', 'subject_entity_type', 'subject_identifier_type',
            'subject_identifier', 'predicate', 'object_entity_type',
            'object_identifier_type', 'object_identifier',
        ],
    )

    # 4) Direct gene resolution for each copy (consumed by entity_resolution_base).
    con.execute(
        f"""
        CREATE OR REPLACE TABLE multigene_resolution AS
        SELECT
          source,
          new_entity_evidence_id AS entity_evidence_id,
          {_lit(GENE_ENTITY_TYPE)} AS entity_type,
          taxonomy_id,
          entrez_type_id AS canonical_identifier_type_id,
          entrez AS canonical_identifier
        FROM multigene_split
        """
    )
    return int(copies)


def _explode_one(con, table: str, ref_col: str, *, cols: list) -> None:
    """Rewrite ``table`` fanning out rows whose ``ref_col`` is a multi-gene
    mention into one copy per gene (``ref_col`` ← the gene-specific new id)."""

    def _proj(use_new: bool) -> str:
        parts = []
        for col in cols:
            if isinstance(col, tuple):
                src_col, new_expr = col
                parts.append(
                    f'ms.{new_expr} AS {src_col}' if use_new else f't.{src_col}'
                )
            else:
                parts.append(f't.{col}')
        return ', '.join(parts)

    kept = _proj(use_new=False)
    exploded = _proj(use_new=True)
    con.execute(
        f"""
        CREATE OR REPLACE TABLE {table}_mg AS
        SELECT {kept}
        FROM {table} t
        WHERE NOT EXISTS (
          SELECT 1 FROM multigene_split ms
          WHERE ms.source = t.source
            AND ms.orig_entity_evidence_id = t.{ref_col}
        )
        UNION ALL
        SELECT {exploded}
        FROM {table} t
        JOIN multigene_split ms
          ON ms.source = t.source
         AND ms.orig_entity_evidence_id = t.{ref_col}
        """
    )
    con.execute(f'DROP TABLE {table}')
    con.execute(f'ALTER TABLE {table}_mg RENAME TO {table}')
