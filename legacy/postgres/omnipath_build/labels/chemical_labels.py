"""Brevity-first chemical label cascade.

Runs during ``derive`` **after** :func:`populate_entity_labels` (which has
already labelled genes by symbol and given every entity a universal
identifier fallback). This step overwrites the labels of **chemical** entities
(``Chemical:OM:0037``) with a human-readable cascade.

Selection follows three passes:

1. **Name cascade** — over the name-bearing identifier types collected from
   ChEBI / HMDB / ChEMBL / RaMP / RefMet / LIPID MAPS / SwissLipids (all
   consolidated into the generic ``Name`` / ``Synonym`` / IUPAC identifier
   types in the MAIN layer). Candidates are tiered by *recognisability* —
   curated primary ``Name`` first, INN / abbreviation / traditional-IUPAC next,
   the synonym soup below, and the systematic IUPAC name only as a last-resort
   name. **Within a tier the most broadly-attested name wins (recognisability),
   then the shortest (brevity), then alphabetical** for determinism. This is the
   ordering the user asked for: recognisable first (e.g. ``alanine``), short
   second, systematic IUPAC only when nothing better exists.

2. **The shortest name the authority supports** (spec 011 data-model.md
   section 4) — for whatever pass 1 still left unlabeled, the shortest name
   ``entity_name.py`` marked ``is_preferred`` (the minting resource's own
   recommended name, the authority relation applied to names, R1/FR-016).
   Additive only: it never overrides a pass-1 label.

3. **Identifier fallback** — chemicals with no usable name take their best real
   external identifier (ChEBI / KEGG / PubChem / HMDB / ChEMBL / … →
   molecular formula → SMILES → InChIKey as the true last resort). The opaque
   ``unresolved_entity_key`` hash that the universal fallback would otherwise
   leave is **never** used, and an InChIKey is used only when a chemical has
   literally no other identifier (none do in current builds).

Junk guards drop empty / single-character values, pure-numeric "names" (leaked
ids), InChIKey-shaped strings, anything equal to the canonical identifier, and
pathologically long (> 120 char) systematic names (those fall through to the id
fallback, which yields a clean ``CHEBI:..`` instead of a 200-char label).
"""

from __future__ import annotations

from dataclasses import dataclass

from psycopg2 import sql
import psycopg2.extensions


CHEMICAL_ENTITY_TYPE = 'Chemical:OM:0037'

CHEMICAL_NAME_RULE = 'chemical_name'
CHEMICAL_IUPAC_RULE = 'chemical_iupac_name'
CHEMICAL_PREFERRED_NAME_RULE = 'chemical_preferred_name'
CHEMICAL_IDENTIFIER_RULE = 'chemical_identifier'

# An InChIKey is 14-10-1 uppercase blocks — never a human-readable label.
INCHIKEY_RE = r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$'

# Name-bearing identifier types, tiered by recognisability (lower = preferred).
# Tier 6 (systematic IUPAC) is the last-resort *name* and is recorded under a
# distinct rule so QA can see how many chemicals fell that far.
_NAME_TIERS = (
    ('Name:OM:0202', 1),
    ('Inn:OM:0120', 2),
    ('Abbreviated Name:OM:0208', 3),
    ('Iupac Traditional Name:OM:0211', 4),
    ('Synonym:OM:0203', 5),
    ('Iupac Name:OM:0210', 6),
)
_IUPAC_LAST_RESORT_TIER = 6

# Real external identifiers for the no-name fallback (lower = preferred). Bare
# numeric namespaces get a readable prefix; self-identifying ids pass through.
# The ``unresolved_entity_key`` hash is deliberately absent — it is excluded.
_ID_TIERS = (
    ('Chebi:MI:0474', 1, 'CHEBI:'),
    ('Kegg Compound:MI:2012', 2, ''),
    ('Pubchem Compound:OM:0002', 3, 'CID:'),
    ('Pubchem:MI:0730', 4, 'CID:'),
    ('Hmdb:OM:0004', 5, ''),
    ('Chembl Compound:MI:0967', 6, ''),
    ('Lipidmaps:OM:0003', 7, ''),
    ('Swisslipids:OM:0009', 8, ''),
    ('Drugbank:MI:2002', 9, ''),
    ('Metanetx:OM:0005', 10, ''),
    ('Bigg Metabolite:OM:0233', 11, 'BiGG:'),
    ('Refmet:OM:0137', 12, ''),
    ('Ramp Id:OM:0132', 13, ''),
    ('Human Gem Metabolite:OM:0243', 14, ''),
    ('Cas:MI:2011', 15, 'CAS:'),
    ('Reactome Stable Id:OM:0130', 16, ''),
    ('Guidetopharma:OM:0008', 17, 'GtoPdb:'),
    ('Bindingdb:OM:0006', 18, 'BindingDB:'),
    ('Drugcentral:OM:0242', 19, 'DrugCentral:'),
    ('Foodb:OM:0213', 20, ''),
    ('Phenol Explorer:OM:0247', 21, 'PhenolExplorer:'),
    ('Ptfi:OM:0214', 22, ''),
    ('Molecular Formula:OM:0212', 55, ''),
    ('Smiles:MI:0239', 60, 'SMILES:'),
    ('Standard Inchi Key:MI:1101', 99, ''),
)


@dataclass(frozen=True)
class ChemicalLabelStats:
    chemical_name: int = 0
    chemical_iupac_name: int = 0
    chemical_preferred_name: int = 0
    chemical_identifier: int = 0
    chemical_without_real_label: int = 0


def _scalar(cur: psycopg2.extensions.cursor, query: sql.SQL, params: list) -> int:
    cur.execute(query, params)
    row = cur.fetchone()
    return int(row[0]) if row and row[0] is not None else 0


def _values_clause(rows: tuple) -> sql.SQL:
    """Render a literal ``VALUES`` list (trusted constants, not user input).

    Embeds the tier table as SQL literals so the queries carry only the named
    parameters (psycopg2 forbids mixing positional ``%s`` and named ``%(x)s``).
    """

    return sql.SQL(', ').join(
        sql.SQL('({})').format(
            sql.SQL(', ').join(sql.Literal(col) for col in row)
        )
        for row in rows
    )


def populate_chemical_labels(
    conn: psycopg2.extensions.connection,
    *,
    schema: str = 'public',
) -> ChemicalLabelStats:
    """Overwrite chemical labels via the recognisability-then-brevity cascade."""

    schema_id = sql.Identifier(schema)

    with conn.cursor() as cur:
        chem_type_id = _scalar(
            cur,
            sql.SQL(
                'SELECT entity_type_id FROM {}.vocab_entity_type WHERE name = %s'
            ).format(schema_id),
            [CHEMICAL_ENTITY_TYPE],
        )
        if not chem_type_id:
            return ChemicalLabelStats()

        # --- Pass 1: name cascade ------------------------------------------
        # Attestation is counted case-insensitively (so "Caffeine" + "caffeine"
        # reinforce each other); the displayed casing is the best-attested,
        # then shortest, then alphabetical exact form.
        cur.execute(
            sql.SQL(
                """
                WITH name_type(type_name, tier) AS (VALUES {name_values}),
                candidate AS (
                  -- One source name may pack several names with a ``|``
                  -- delimiter and a trailing ``(PTF#####)`` source-id suffix
                  -- (FooDB/PTFI); split and strip so each clean name competes
                  -- on its own merit.
                  SELECT
                    ei.entity_id,
                    nt.tier,
                    lower(c.cleaned) AS norm,
                    c.cleaned        AS val,
                    ei.source_id
                  FROM {schema}.entity e
                  JOIN {schema}.entity_identifier ei
                    ON ei.entity_id = e.entity_id
                  JOIN {schema}.identifier_evidence ie
                    ON ie.identifier_id = ei.identifier_id
                  JOIN {schema}.vocab_identifier_type it
                    ON it.identifier_type_id = ie.identifier_type_id
                  JOIN name_type nt ON nt.type_name = it.name
                  CROSS JOIN LATERAL (
                    SELECT btrim(
                      regexp_replace(part, '\\s*\\(PTF[0-9]+\\)\\s*$', '')
                    ) AS cleaned
                    FROM unnest(string_to_array(ie.value, '|')) AS part
                  ) c
                  WHERE e.entity_type_id = %(chem)s
                    AND c.cleaned <> ''
                    AND length(c.cleaned) BETWEEN 2 AND 120
                    AND c.cleaned !~ '^[0-9]+$'
                    AND c.cleaned !~ %(inchikey)s
                    AND c.cleaned IS DISTINCT FROM e.canonical_identifier
                ),
                exact AS (
                  SELECT entity_id, tier, norm, val,
                         count(DISTINCT source_id) AS srcs
                  FROM candidate
                  GROUP BY entity_id, tier, norm, val
                ),
                by_norm AS (
                  SELECT entity_id, tier, norm,
                         sum(srcs) AS attestations,
                         min(length(val)) AS len
                  FROM exact
                  GROUP BY entity_id, tier, norm
                ),
                display AS (
                  SELECT DISTINCT ON (entity_id, tier, norm)
                         entity_id, tier, norm, val
                  FROM exact
                  ORDER BY entity_id, tier, norm,
                           srcs DESC, length(val), val
                ),
                ranked AS (
                  SELECT
                    n.entity_id, n.tier, d.val,
                    row_number() OVER (
                      PARTITION BY n.entity_id
                      ORDER BY n.tier, n.attestations DESC, n.len, d.val
                    ) AS rk
                  FROM by_norm n
                  JOIN display d USING (entity_id, tier, norm)
                )
                UPDATE {schema}.entity e
                SET label = ranked.val,
                    label_rule = CASE
                      WHEN ranked.tier >= %(iupac_tier)s THEN %(iupac_rule)s
                      ELSE %(name_rule)s
                    END
                FROM ranked
                WHERE ranked.entity_id = e.entity_id
                  AND ranked.rk = 1
                  AND e.entity_type_id = %(chem)s
                """
            ).format(
                schema=schema_id,
                name_values=_values_clause(_NAME_TIERS),
            ),
            dict(
                chem=chem_type_id,
                inchikey=INCHIKEY_RE,
                iupac_tier=_IUPAC_LAST_RESORT_TIER,
                iupac_rule=CHEMICAL_IUPAC_RULE,
                name_rule=CHEMICAL_NAME_RULE,
            ),
        )

        # --- Pass 1b: the shortest name the authority supports -------------
        # entity_name (spec 011 data-model.md section 4) marks a name
        # is_preferred when the minting resource itself recommends it --
        # the same authority relation the resolver uses (R1), applied to
        # names. Pass 1's own tiered cascade already covers most of this
        # ground through entity_identifier directly; this pass only fires
        # for the entities it left unlabeled, catching a preferred name
        # attached through entity_name but not through the tiered
        # identifier types above (e.g. a resource entity_name reaches that
        # chemical_labels.py's own name-bearing type list does not).
        cur.execute(
            sql.SQL(
                """
                WITH candidate AS (
                  SELECT en.entity_id, en.name AS val, length(en.name) AS len
                  FROM {schema}.entity e
                  JOIN {schema}.entity_name en ON en.entity_id = e.entity_id
                  WHERE e.entity_type_id = %(chem)s
                    AND e.label_rule IS DISTINCT FROM %(name_rule)s
                    AND e.label_rule IS DISTINCT FROM %(iupac_rule)s
                    AND en.is_preferred
                ),
                ranked AS (
                  SELECT entity_id, val,
                    row_number() OVER (
                      PARTITION BY entity_id ORDER BY len, val
                    ) AS rk
                  FROM candidate
                )
                UPDATE {schema}.entity e
                SET label = ranked.val,
                    label_rule = %(preferred_name_rule)s
                FROM ranked
                WHERE ranked.entity_id = e.entity_id
                  AND ranked.rk = 1
                  AND e.entity_type_id = %(chem)s
                """
            ).format(schema=schema_id),
            dict(
                chem=chem_type_id,
                name_rule=CHEMICAL_NAME_RULE,
                iupac_rule=CHEMICAL_IUPAC_RULE,
                preferred_name_rule=CHEMICAL_PREFERRED_NAME_RULE,
            ),
        )

        # --- Pass 2: real-identifier fallback (no usable name) -------------
        cur.execute(
            sql.SQL(
                """
                WITH id_type(type_name, tier, prefix) AS (VALUES {id_values}),
                candidate AS (
                  SELECT
                    ei.entity_id,
                    idt.tier,
                    idt.prefix || btrim(ie.value) AS label_value,
                    length(btrim(ie.value)) AS len
                  FROM {schema}.entity e
                  JOIN {schema}.entity_identifier ei
                    ON ei.entity_id = e.entity_id
                  JOIN {schema}.identifier_evidence ie
                    ON ie.identifier_id = ei.identifier_id
                  JOIN {schema}.vocab_identifier_type it
                    ON it.identifier_type_id = ie.identifier_type_id
                  JOIN id_type idt ON idt.type_name = it.name
                  WHERE e.entity_type_id = %(chem)s
                    AND e.label_rule IS DISTINCT FROM %(name_rule)s
                    AND e.label_rule IS DISTINCT FROM %(iupac_rule)s
                    AND e.label_rule IS DISTINCT FROM %(preferred_name_rule)s
                    AND ie.value IS NOT NULL
                    AND btrim(ie.value) <> ''
                ),
                ranked AS (
                  SELECT entity_id, label_value,
                    row_number() OVER (
                      PARTITION BY entity_id
                      ORDER BY tier, len, label_value
                    ) AS rk
                  FROM candidate
                )
                UPDATE {schema}.entity e
                SET label = ranked.label_value,
                    label_rule = %(id_rule)s
                FROM ranked
                WHERE ranked.entity_id = e.entity_id
                  AND ranked.rk = 1
                  AND e.entity_type_id = %(chem)s
                """
            ).format(
                schema=schema_id,
                id_values=_values_clause(_ID_TIERS),
            ),
            dict(
                chem=chem_type_id,
                name_rule=CHEMICAL_NAME_RULE,
                iupac_rule=CHEMICAL_IUPAC_RULE,
                preferred_name_rule=CHEMICAL_PREFERRED_NAME_RULE,
                id_rule=CHEMICAL_IDENTIFIER_RULE,
            ),
        )

        # --- Stats ---------------------------------------------------------
        def _count(rule: str) -> int:
            return _scalar(
                cur,
                sql.SQL(
                    'SELECT count(*) FROM {}.entity '
                    'WHERE entity_type_id = %s AND label_rule = %s'
                ).format(schema_id),
                [chem_type_id, rule],
            )

        chemical_name = _count(CHEMICAL_NAME_RULE)
        chemical_iupac_name = _count(CHEMICAL_IUPAC_RULE)
        chemical_preferred_name = _count(CHEMICAL_PREFERRED_NAME_RULE)
        chemical_identifier = _count(CHEMICAL_IDENTIFIER_RULE)
        without_real_label = _scalar(
            cur,
            sql.SQL(
                "SELECT count(*) FROM {}.entity "
                "WHERE entity_type_id = %s "
                "AND (label IS NULL OR label = '' "
                "     OR label ~ %s "
                "     OR label_rule NOT IN (%s, %s, %s, %s))"
            ).format(schema_id),
            [
                chem_type_id,
                INCHIKEY_RE,
                CHEMICAL_NAME_RULE,
                CHEMICAL_IUPAC_RULE,
                CHEMICAL_PREFERRED_NAME_RULE,
                CHEMICAL_IDENTIFIER_RULE,
            ],
        )

    conn.commit()
    return ChemicalLabelStats(
        chemical_name=chemical_name,
        chemical_iupac_name=chemical_iupac_name,
        chemical_preferred_name=chemical_preferred_name,
        chemical_identifier=chemical_identifier,
        chemical_without_real_label=without_real_label,
    )
