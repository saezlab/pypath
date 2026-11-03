"""Non-lipid chemical fallback resolution.

When a chemical mention has no usable **structure** (InChIKey), the legacy
resolver dropped it to the opaque ``unresolved_entity_key`` md5 hash — ~44% of
chemicals. This module replaces the distrusted transitive co-occurrence
clustering with a **per-record, priority-ordered pick** (no transitivity, no
false chain-merges): each structure-less chemical canonicalises to its single
best identifier by a fixed priority, recording the producing
``resolution_mechanism``.

Priority (this module — STAGE 1):

  ChEBI → ChEMBL → PubChem → SwissLipids → HMDB → LIPID MAPS
  → keep-original (any other real external id, e.g. KEGG/FooDB/MetaNetX/CAS…)
  → name (exact, basic collision handling)

A mention's primary id is therefore never a lower-priority id when a higher one
is present (a ChEBI+PubChem record anchors on ChEBI). Same id across resources
merges; different ids stay distinct. **InChIKey** is intentionally absent —
those mentions resolve through the existing direct InChIKey path and never reach
here. **SMILES** is also intentionally absent: a structure-less lipid like
``PC(18:1_16:0)`` yields a placeholder SMILES that would collapse distinct
species, so SMILES is never a canonical merge key — it stays an attached
identifier only.

Deferred to a follow-up (the ~3.4% non-priority bulk): the **within-resource
id merge** (FooDB membership→compound, lifting structure-less membership rows
onto the structured compound that shares their FooDB id) and the **ChEBI-xref
1:1 translation** (CAS/KEGG→ChEBI, abort on 1→many). UniChem is not ingested.
Until then those ids resolve via *keep-original* (a stable id, the
``original_id`` mechanism — still far better than the md5 hash).

The result table ``chemical_fallback_resolution`` is consumed by
``entity_resolution_base``'s unresolved branch: it is applied only to chemical
mentions the structure + resolver-candidate paths left unresolved.
"""

from __future__ import annotations

from pypath.internals.cv_terms import (
    IdentifierNamespaceCv,
    cv_term_label_accession,
)

from omnipath_build.cv_terms import CHEMICAL_ENTITY_TYPE

NAME_TYPE = cv_term_label_accession(IdentifierNamespaceCv.NAME)


def normalize_name(value: str) -> str:
    """The one normalized form of a chemical name or synonym (spec 011 R5).

    Case- and whitespace-folded so the same molecule under different
    capitalization or incidental spacing collapses to one canonical entity,
    instead of two content-addressed entity ids for "Taurine" and "taurine".
    Mirrors ``omnipath_utils.mapping._id_types.normalize_name`` exactly (the
    established cross-repo duplication precedent for a small, pure function
    needed identically in both repos) -- keep the two in sync if either
    changes.

    Not called directly by ``build_chemical_fallback_resolution`` below (the
    actual per-row computation runs as DuckDB SQL, ``lower(trim(...))``, for
    bulk performance) -- this is the single documented, tested source of
    truth for what that SQL must compute, and the function any other Python
    caller should use.
    """

    return str(value).strip().lower()

# (identifier_type label, tier [lower=preferred], resolution_mechanism).
# Tiers are distinct so the per-mention pick is fully deterministic.
_TIERS: tuple[tuple[str, int, str], ...] = (
    # SMILES intentionally omitted — never a canonical merge key.
    ('Chebi:MI:0474', 2, 'chebi'),
    ('Chembl Compound:MI:0967', 3, 'chembl'),
    ('Pubchem Compound:OM:0002', 4, 'pubchem'),
    ('Pubchem:MI:0730', 5, 'pubchem'),
    ('Swisslipids:OM:0009', 6, 'swisslipids'),
    ('Hmdb:OM:0004', 7, 'hmdb'),
    ('Lipidmaps:OM:0003', 8, 'lipidmaps'),
    # keep-original: any other real external id → its own id is the canonical
    # identity (stable, mergeable across resources), mechanism 'original_id'.
    ('Kegg Compound:MI:2012', 20, 'original_id'),
    ('Metanetx:OM:0005', 21, 'original_id'),
    ('Bigg Metabolite:OM:0233', 22, 'original_id'),
    ('Human Gem Metabolite:OM:0243', 23, 'original_id'),
    ('Refmet:OM:0137', 24, 'original_id'),
    ('Cas:MI:2011', 25, 'original_id'),
    ('Drugbank:MI:2002', 26, 'original_id'),
    ('Drugcentral:OM:0242', 27, 'original_id'),
    ('Guidetopharma:OM:0008', 28, 'original_id'),
    ('Bindingdb:OM:0006', 29, 'original_id'),
    ('Reactome Stable Id:OM:0130', 30, 'original_id'),
    ('Ramp Id:OM:0132', 31, 'original_id'),
    ('Phenol Explorer:OM:0247', 32, 'original_id'),
    ('Foodb:OM:0213', 33, 'original_id'),
    ('Ptfi:OM:0214', 34, 'original_id'),
    ('Pubchem Substance:OM:0028', 35, 'original_id'),
    # name (no id at all): merge by exact name. Basic guard only — drop
    # empty/InChIKey-shaped/pure-numeric values; the unambiguous-name guard is
    # part of the deferred step.
    ('Name:OM:0202', 40, 'name'),
    ('Synonym:OM:0203', 41, 'name'),
)


def _values(rows: tuple) -> str:
    parts = []
    for name, tier, mech in rows:
        n = name.replace("'", "''")
        m = mech.replace("'", "''")
        parts.append(f"('{n}', {tier}, '{m}')")
    return ', '.join(parts)


def chemical_fallback_fires_sql(
    rcs_alias: str = 'rcs',
    cf_alias: str = 'cf',
) -> str:
    """SQL predicate: may the per-record chemical fallback supply the identity?

    The resolver wins at ``candidate_count = 1``; the fallback (``cf``) fires
    **only when the resolver produced no candidates** (``candidate_count`` 0 or
    NULL). When the resolver is genuinely ambiguous (``candidate_count > 1``)
    the entity stays **unresolved** — the fallback must not pick one of several
    distinct structures. This is the single source of truth for the gate,
    consumed by ``entity_resolution_base`` and unit-tested.
    """

    return (
        f'{cf_alias}.canonical_identifier IS NOT NULL '
        f'AND coalesce({rcs_alias}.candidate_count, 0) = 0'
    )


def build_chemical_anchor_map(con, *, log=lambda *_: None) -> int:
    """Build ``chemical_anchor_map`` — 1:1 translations of a non-structural id to
    a structure/ChEBI **anchor** (stage 2 of the fallback).

    For every chemical mention that carries an InChIKey or ChEBI, map its *other*
    ids (FooDB/KEGG/PubChem/CAS/…) to that mention's anchor (InChIKey preferred,
    else ChEBI). This is the resource's own per-record assertion (e.g. FooDB's
    compound row links its `public_id` to the compound's InChIKey; a ChEBI row
    links ChEBI to its KEGG/CAS xref) — **one hop, not transitive**. Keep only
    ids that map to exactly **one** distinct anchor (abort on 1→many).

    Consumed by ``build_chemical_fallback_resolution``: a structure-less mention
    whose id translates this way resolves to the structure/ChEBI (merging it onto
    the structured entity) instead of keeping its own id.
    """
    chem = CHEMICAL_ENTITY_TYPE.replace("'", "''")
    inchikey_re = "'^[A-Z]{14}-[A-Z]{10}-[A-Z]$'"
    # id types that are NOT translation *sources* (structures, names, formula, hash).
    non_src = (
        "'Standard Inchi Key:MI:1101','Smiles:MI:0239','Name:OM:0202',"
        "'Synonym:OM:0203','Iupac Name:OM:0210','Iupac Traditional Name:OM:0211',"
        "'Abbreviated Name:OM:0208','Inn:OM:0120','Molecular Formula:OM:0212',"
        "'omnipath:unresolved_entity_key'"
    )
    con.execute(
        f"""
        CREATE OR REPLACE TABLE chemical_anchor_map AS
        WITH ids AS (
          SELECT ee.source, ee.entity_evidence_id,
                 ei.identifier_type AS type_name, trim(ei.identifier) AS val
          FROM entity_evidence_raw ee
          JOIN entity_identifier_raw ei
            ON ei.source = ee.source
           AND ei.entity_evidence_id = ee.entity_evidence_id
          WHERE ee.entity_type = '{chem}'
            AND ei.identifier IS NOT NULL AND trim(ei.identifier) <> ''
        ),
        mention_anchor AS (   -- per mention: best anchor (InChIKey > ChEBI)
          SELECT source, entity_evidence_id,
            max(val) FILTER (
              WHERE type_name = 'Standard Inchi Key:MI:1101'
                AND val ~ {inchikey_re}
            ) AS inchikey,
            max(val) FILTER (WHERE type_name = 'Chebi:MI:0474') AS chebi
          FROM ids GROUP BY source, entity_evidence_id
        ),
        pair AS (             -- (non-anchor id) -> (anchor) on the same mention
          SELECT i.type_name AS src_type, i.val AS src_value,
                 ma.inchikey, ma.chebi
          FROM ids i
          JOIN mention_anchor ma
            ON ma.source = i.source AND ma.entity_evidence_id = i.entity_evidence_id
          WHERE (ma.inchikey IS NOT NULL OR ma.chebi IS NOT NULL)
            AND i.type_name NOT IN ({non_src})
        ),
        agg AS (
          SELECT src_type, src_value,
            count(DISTINCT inchikey) AS n_ik, max(inchikey) AS ik,
            count(DISTINCT chebi)    AS n_chebi, max(chebi) AS chebi
          FROM pair GROUP BY src_type, src_value
        )
        SELECT src_type, src_value,
          CASE WHEN n_ik = 1 THEN 'Standard Inchi Key:MI:1101' ELSE 'Chebi:MI:0474' END
            AS anchor_type,
          CASE WHEN n_ik = 1 THEN ik ELSE chebi END AS anchor_value,
          CASE WHEN n_ik = 1 THEN 0 ELSE 2 END AS anchor_tier,
          CASE WHEN n_ik = 1 THEN 'anchored_structure' ELSE 'anchored_chebi' END
            AS mechanism
        FROM agg
        WHERE n_ik = 1 OR (n_ik = 0 AND n_chebi = 1)   -- 1:1 only
        """
    )
    # Collision guard: a chemical name is ambiguous if, on
    # structure-bearing mentions, it appears with >1 distinct InChIKey (e.g. a
    # trivial name shared by L-/D-/racemic forms). Such names must NOT be used as
    # a canonical identity (they would false-merge distinct molecules).
    #
    # Grouped by the SAME normalize_name() fold build_chemical_fallback_resolution
    # applies to canonical_identifier below (lower(trim(...))) -- otherwise this
    # guard could disagree on a name's ambiguity status across case variants that
    # fold now merges (e.g. "Taurine" flagged ambiguous, "taurine" not, for what
    # becomes one canonical identity either way).
    con.execute(
        f"""
        CREATE OR REPLACE TABLE chemical_ambiguous_name AS
        SELECT name_val FROM (
          SELECT lower(trim(nm.identifier)) AS name_val,
                 count(DISTINCT trim(ik.identifier)) AS n_ik
          FROM entity_evidence_raw ee
          JOIN entity_identifier_raw nm
            ON nm.source = ee.source
           AND nm.entity_evidence_id = ee.entity_evidence_id
           AND nm.identifier_type IN ('Name:OM:0202', 'Synonym:OM:0203')
          JOIN entity_identifier_raw ik
            ON ik.source = ee.source
           AND ik.entity_evidence_id = ee.entity_evidence_id
           AND ik.identifier_type = 'Standard Inchi Key:MI:1101'
           AND trim(ik.identifier) ~ {inchikey_re}
          WHERE ee.entity_type = '{chem}'
            AND nm.identifier IS NOT NULL AND trim(nm.identifier) <> ''
          GROUP BY lower(trim(nm.identifier))
        ) WHERE n_ik > 1
        """
    )
    rows = con.execute('SELECT count(*) FROM chemical_anchor_map').fetchone()[0]
    amb = con.execute(
        'SELECT count(*) FROM chemical_ambiguous_name'
    ).fetchone()[0]
    log(f'chemical anchor map: {rows} 1:1 translations; {amb} ambiguous names')
    return int(rows)


def build_chemical_fallback_resolution(con, *, log=lambda *_: None) -> int:
    """Build ``chemical_fallback_resolution`` (one best-id row per chemical
    mention with a usable non-structure identifier). Returns the row count.

    Requires ``chemical_anchor_map`` (``build_chemical_anchor_map``) to exist —
    its translations enter the candidate set so a structure-less mention prefers
    a structure/ChEBI anchor over keeping its own id.
    """

    chem = CHEMICAL_ENTITY_TYPE.replace("'", "''")
    name_re = "'^[A-Z]{14}-[A-Z]{10}-[A-Z]$'"  # InChIKey-shaped guard for names
    con.execute(
        f"""
        CREATE OR REPLACE TABLE chemical_fallback_resolution AS
        WITH tier(type_name, tier, mechanism) AS (VALUES {_values(_TIERS)}),
        candidate AS (
          SELECT
            ee.source,
            ee.entity_evidence_id,
            t.tier,
            0 AS is_anchored,
            t.mechanism,
            it.identifier_type_id AS canonical_identifier_type_id,
            -- name tiers: case/whitespace-folded (normalize_name, spec 011
            -- R5) so "Taurine" and "taurine" from two different mentions
            -- merge into one canonical entity instead of two. Every other
            -- tier's identifier is a structured id (chebi/chembl/pubchem/
            -- kegg/...), never folded -- those must stay exactly as-is.
            CASE
              WHEN t.mechanism = 'name' THEN lower(trim(ei.identifier))
              ELSE trim(ei.identifier)
            END AS canonical_identifier
          FROM entity_evidence_raw ee
          JOIN entity_identifier_raw ei
            ON ei.source = ee.source
           AND ei.entity_evidence_id = ee.entity_evidence_id
          JOIN tier t ON t.type_name = ei.identifier_type
          JOIN identifier_type_all it ON it.name = ei.identifier_type
          WHERE ee.entity_type = '{chem}'
            AND ei.identifier IS NOT NULL
            AND trim(ei.identifier) <> ''
            -- name tiers: drop junk (InChIKey-shaped / pure-numeric)
            AND NOT (
              t.mechanism = 'name'
              AND (
                trim(ei.identifier) ~ {name_re}
                OR trim(ei.identifier) ~ '^[0-9]+$'
                OR length(trim(ei.identifier)) < 2
                -- collision guard: drop names that map to >1
                -- distinct structure (ambiguous → never a canonical identity).
                -- chemical_ambiguous_name is grouped by the same fold.
                OR lower(trim(ei.identifier)) IN (SELECT name_val FROM chemical_ambiguous_name)
              )
            )
          UNION ALL
          -- stage 2: translate a mention's id to its 1:1
          -- structure/ChEBI anchor; enters at the anchor's tier (InChIKey 0,
          -- ChEBI 2) so it beats keep-original and merges the mention onto the
          -- structured / ChEBI-keyed entity.
          SELECT
            ee.source, ee.entity_evidence_id, am.anchor_tier AS tier,
            1 AS is_anchored,
            am.mechanism,
            ait.identifier_type_id AS canonical_identifier_type_id,
            am.anchor_value AS canonical_identifier
          FROM entity_evidence_raw ee
          JOIN entity_identifier_raw ei
            ON ei.source = ee.source
           AND ei.entity_evidence_id = ee.entity_evidence_id
          JOIN chemical_anchor_map am
            ON am.src_type = ei.identifier_type
           AND am.src_value = trim(ei.identifier)
          JOIN identifier_type_all ait ON ait.name = am.anchor_type
          WHERE ee.entity_type = '{chem}'
        ),
        ranked AS (
          SELECT
            source, entity_evidence_id, mechanism,
            canonical_identifier_type_id, canonical_identifier,
            row_number() OVER (
              PARTITION BY source, entity_evidence_id
              ORDER BY tier, is_anchored, canonical_identifier
            ) AS rk
          FROM candidate
        )
        SELECT
          source, entity_evidence_id,
          canonical_identifier_type_id, canonical_identifier, mechanism
        FROM ranked
        WHERE rk = 1
        """
    )
    rows = con.execute(
        'SELECT count(*) FROM chemical_fallback_resolution'
    ).fetchone()[0]
    log(f'chemical fallback: {rows} chemical mentions with a fallback id')
    return int(rows)
