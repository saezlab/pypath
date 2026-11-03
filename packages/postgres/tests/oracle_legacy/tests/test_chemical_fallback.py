"""Unit test for the non-lipid chemical fallback resolution.

Synthetic DuckDB raw tables → ``build_chemical_fallback_resolution`` → assert the
per-record priority pick (no transitive clustering): SMILES → ChEBI → ChEMBL →
PubChem → SwissLipids → HMDB → LIPID MAPS → keep-original → name, with the
producing ``resolution_mechanism`` recorded. Needs only ``duckdb``.
"""

from __future__ import annotations

import pytest

duckdb = pytest.importorskip('duckdb')

from omnipath_build.chemical_fallback import (  # noqa: E402
    build_chemical_anchor_map,
    build_chemical_fallback_resolution,
    normalize_name,
)


def test_normalize_name_folds_case_and_whitespace():
    assert normalize_name('Taurine') == 'taurine'
    assert normalize_name('  taurine  ') == 'taurine'
    assert normalize_name('TAURINE') == 'taurine'


def _resolve(con):
    """Stage-2 anchor map + the fallback pick (the real pipeline order)."""
    build_chemical_anchor_map(con)
    build_chemical_fallback_resolution(con)

CHEM = 'Chemical:OM:0037'
PROTEIN = 'Protein:MI:0326'
# (type label -> id) for identifier_type_all
TYPES = {
    'Smiles:MI:0239': 1,
    'Chebi:MI:0474': 2,
    'Chembl Compound:MI:0967': 3,
    'Pubchem Compound:OM:0002': 4,
    'Hmdb:OM:0004': 5,
    'Foodb:OM:0213': 6,
    'Name:OM:0202': 7,
    'Standard Inchi Key:MI:1101': 8,
    'Kegg Compound:MI:2012': 9,
}
IK1 = 'AAAAAAAAAAAAAA-BBBBBBBBBB-N'
IK2 = 'CCCCCCCCCCCCCC-DDDDDDDDDD-N'


def _con(rows):
    con = duckdb.connect(':memory:')
    con.execute('CREATE TABLE identifier_type_all (identifier_type_id BIGINT, name VARCHAR)')
    con.executemany('INSERT INTO identifier_type_all VALUES (?, ?)',
                    [(v, k) for k, v in TYPES.items()])
    con.execute("""CREATE TABLE entity_evidence_raw (
      source VARCHAR, dataset VARCHAR, row_id BIGINT, entity_evidence_id VARCHAR,
      parent_entity_evidence_id VARCHAR, entity_role VARCHAR, entity_type VARCHAR,
      taxonomy_id VARCHAR)""")
    con.execute("""CREATE TABLE entity_identifier_raw (
      source VARCHAR, entity_evidence_id VARCHAR, identifier_id VARCHAR,
      identifier_type VARCHAR, identifier VARCHAR)""")
    # rows: (ev_id, entity_type, [(id_type, id_value), ...])
    for ev, etype, ids in rows:
        con.execute('INSERT INTO entity_evidence_raw VALUES (?,?,?,?,?,?,?,?)',
                    ('s', 'd', 1, ev, None, 'r', etype, None))
        for it, val in ids:
            con.execute('INSERT INTO entity_identifier_raw VALUES (?,?,?,?,?)',
                        ('s', ev, ev + it, it, val))
    return con


def _pick(con, ev):
    r = con.execute(
        'SELECT canonical_identifier_type_id, canonical_identifier, mechanism '
        'FROM chemical_fallback_resolution WHERE entity_evidence_id = ?', [ev],
    ).fetchall()
    return r[0] if r else None


def test_priority_pick_and_mechanisms():
    con = _con([
        ('m_chebi_pubchem', CHEM, [('Chebi:MI:0474', '15365'), ('Pubchem Compound:OM:0002', '2244')]),
        ('m_pubchem', CHEM, [('Pubchem Compound:OM:0002', '999')]),
        ('m_foodb', CHEM, [('Foodb:OM:0213', 'FDB000123')]),
        ('m_smiles_chebi', CHEM, [('Smiles:MI:0239', 'CCO'), ('Chebi:MI:0474', '16236')]),
        ('m_name', CHEM, [('Name:OM:0202', 'aspirin')]),
    ])
    _resolve(con)
    # ChEBI beats PubChem (higher priority), not both → no clustering.
    assert _pick(con, 'm_chebi_pubchem') == (TYPES['Chebi:MI:0474'], '15365', 'chebi')
    assert _pick(con, 'm_pubchem') == (TYPES['Pubchem Compound:OM:0002'], '999', 'pubchem')
    # FooDB-only → keep-original (stable id, not a hash).
    assert _pick(con, 'm_foodb') == (TYPES['Foodb:OM:0213'], 'FDB000123', 'original_id')
    # SMILES is no longer a fallback tier → ChEBI wins.
    assert _pick(con, 'm_smiles_chebi') == (TYPES['Chebi:MI:0474'], '16236', 'chebi')
    # name-only → name.
    assert _pick(con, 'm_name') == (TYPES['Name:OM:0202'], 'aspirin', 'name')


def test_smiles_is_not_a_fallback_tier():
    # A placeholder SMILES of a structure-less lipid (e.g. unknown
    # sn-position) would false-merge distinct species, so SMILES is NOT a
    # canonical merge key — it stays an attached identifier only.
    from omnipath_build.chemical_fallback import _TIERS

    assert all(name != 'Smiles:MI:0239' for name, _, _ in _TIERS)


def test_smiles_only_chemical_has_no_fallback_pick():
    con = _con([
        # SMILES alone → no fallback row (SMILES is not a tier).
        ('m_smiles_only', CHEM, [('Smiles:MI:0239', 'CCO')]),
        # SMILES + ChEBI → ChEBI wins, SMILES ignored.
        ('m_smiles_chebi2', CHEM, [
            ('Smiles:MI:0239', 'CCO'), ('Chebi:MI:0474', '16236'),
        ]),
    ])
    _resolve(con)
    assert _pick(con, 'm_smiles_only') is None
    assert _pick(con, 'm_smiles_chebi2') == (
        TYPES['Chebi:MI:0474'], '16236', 'chebi',
    )


def test_chemical_fallback_gate_predicate():
    # The per-record fallback may supply the canonical identity ONLY
    # when the resolver produced NO candidates. An ambiguous chemical
    # (candidate_count > 1) must stay unresolved — never a fallback pick.
    from omnipath_build.chemical_fallback import chemical_fallback_fires_sql

    pred = chemical_fallback_fires_sql()
    con = duckdb.connect(':memory:')
    con.execute('CREATE TABLE rcs (candidate_count BIGINT)')
    con.execute('CREATE TABLE cf (canonical_identifier VARCHAR)')

    def fires(count, cf_val):
        con.execute('DELETE FROM rcs')
        con.execute('DELETE FROM cf')
        con.execute('INSERT INTO rcs VALUES (?)', [count])
        con.execute('INSERT INTO cf VALUES (?)', [cf_val])
        return con.execute(f'SELECT {pred} FROM rcs, cf').fetchone()[0]

    assert fires(0, 'CHEBI:1') is True       # no resolver candidates → fires
    assert fires(None, 'CHEBI:1') is True    # NULL (no candidates) → fires
    assert fires(2, 'CHEBI:1') is False      # ambiguous → stays unresolved
    assert fires(1, 'CHEBI:1') is False      # resolver handles =1 upstream
    assert fires(0, None) is False           # no fallback id → nothing to pick


def test_inchikey_and_junk_excluded():
    con = _con([
        # InChIKey present: not a fallback tier → no fallback row (handled by the
        # direct InChIKey path upstream).
        ('m_inchikey_only', CHEM, [('Standard Inchi Key:MI:1101', 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N')]),
        # junk name (InChIKey-shaped / numeric) → dropped, no row.
        ('m_junk_name', CHEM, [('Name:OM:0202', 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N')]),
        ('m_numeric_name', CHEM, [('Name:OM:0202', '12345')]),
        # non-chemical → ignored entirely.
        ('p_protein', PROTEIN, [('Chebi:MI:0474', '1')]),
    ])
    _resolve(con)
    for ev in ('m_inchikey_only', 'm_junk_name', 'm_numeric_name', 'p_protein'):
        assert _pick(con, ev) is None, f'{ev} should have no fallback row'


def test_name_falls_below_real_ids():
    con = _con([
        ('m_id_and_name', CHEM, [('Foodb:OM:0213', 'FDB1'), ('Name:OM:0202', 'foo')]),
    ])
    _resolve(con)
    # a real id (even keep-original FooDB) beats a name.
    assert _pick(con, 'm_id_and_name') == (TYPES['Foodb:OM:0213'], 'FDB1', 'original_id')


def test_anchor_translation_stage2():
    con = _con([
        # anchor mention: FooDB FDB9 co-occurs with a structure -> FDB9 -> IK1
        ('cmpd', CHEM, [('Foodb:OM:0213', 'FDB9'), ('Standard Inchi Key:MI:1101', IK1)]),
        # structure-less membership of the SAME FooDB id -> lifts onto the structure
        ('memb', CHEM, [('Foodb:OM:0213', 'FDB9')]),
        # ChEBI-anchor mention (no InChIKey): KEGG C1 co-occurs with ChEBI 100
        ('keggchebi', CHEM, [('Kegg Compound:MI:2012', 'C1'), ('Chebi:MI:0474', '100')]),
        # KEGG-only mention -> translates up to ChEBI 100 (anchored_chebi)
        ('keggonly', CHEM, [('Kegg Compound:MI:2012', 'C1')]),
        # ambiguous: FDBX maps to two different structures -> NOT translated
        ('ambA', CHEM, [('Foodb:OM:0213', 'FDBX'), ('Standard Inchi Key:MI:1101', IK1)]),
        ('ambB', CHEM, [('Foodb:OM:0213', 'FDBX'), ('Standard Inchi Key:MI:1101', IK2)]),
        ('ambmemb', CHEM, [('Foodb:OM:0213', 'FDBX')]),
    ])
    _resolve(con)
    # FooDB membership lifted onto the structure (merges with the compound).
    assert _pick(con, 'memb') == (TYPES['Standard Inchi Key:MI:1101'], IK1, 'anchored_structure')
    # KEGG-only lifted to ChEBI via the resource's own xref.
    assert _pick(con, 'keggonly') == (TYPES['Chebi:MI:0474'], '100', 'anchored_chebi')
    # ambiguous FooDB id (2 distinct structures) -> NOT translated, keep-original.
    assert _pick(con, 'ambmemb') == (TYPES['Foodb:OM:0213'], 'FDBX', 'original_id')


def test_ambiguous_name_guard():
    con = _con([
        # 'alanine' appears on two distinct structures -> ambiguous name
        ('al1', CHEM, [('Name:OM:0202', 'alanine'), ('Standard Inchi Key:MI:1101', IK1)]),
        ('al2', CHEM, [('Name:OM:0202', 'alanine'), ('Standard Inchi Key:MI:1101', IK2)]),
        # name-only 'alanine' -> guard drops it (no canonical-by-name) -> no row
        ('alname', CHEM, [('Name:OM:0202', 'alanine')]),
        # an unambiguous name (never on a conflicting structure) -> kept
        ('uniq', CHEM, [('Name:OM:0202', 'uniquechem')]),
    ])
    _resolve(con)
    assert _pick(con, 'alname') is None, 'ambiguous name must not be canonical'
    assert _pick(con, 'uniq') == (TYPES['Name:OM:0202'], 'uniquechem', 'name')


def test_case_variant_names_merge_into_one_canonical_identity():
    """spec 011 R5: 'Taurine' and 'taurine' from two different mentions must
    fold to the same canonical_identifier, not two separate content-addressed
    entities."""
    con = _con([
        ('m_upper', CHEM, [('Name:OM:0202', 'Taurine')]),
        ('m_lower', CHEM, [('Name:OM:0202', 'taurine')]),
        ('m_spaced', CHEM, [('Name:OM:0202', '  taurine  ')]),
    ])
    _resolve(con)
    upper = _pick(con, 'm_upper')
    lower = _pick(con, 'm_lower')
    spaced = _pick(con, 'm_spaced')
    assert upper == lower == spaced == (TYPES['Name:OM:0202'], 'taurine', 'name')


def test_ambiguous_name_guard_is_case_insensitive():
    """The collision guard must catch a case variant of an ambiguous name too
    -- otherwise 'Alanine' could slip through as canonical while 'alanine' is
    correctly dropped, for what becomes one merged identity either way."""
    con = _con([
        ('al1', CHEM, [('Name:OM:0202', 'alanine'), ('Standard Inchi Key:MI:1101', IK1)]),
        ('al2', CHEM, [('Name:OM:0202', 'alanine'), ('Standard Inchi Key:MI:1101', IK2)]),
        # case variant of the ambiguous name -- the guard must still drop it.
        ('alname_upper', CHEM, [('Name:OM:0202', 'Alanine')]),
    ])
    _resolve(con)
    assert _pick(con, 'alname_upper') is None, (
        'a case variant of an ambiguous name must not be canonical either'
    )


def test_non_name_tiers_are_never_case_folded():
    """Only the name/synonym tiers fold case -- a structured id (ChEBI here)
    must stay byte-exact."""
    con = _con([
        ('m1', CHEM, [('Chebi:MI:0474', 'CHEBI:15891')]),
    ])
    _resolve(con)
    assert _pick(con, 'm1') == (TYPES['Chebi:MI:0474'], 'CHEBI:15891', 'chebi')
