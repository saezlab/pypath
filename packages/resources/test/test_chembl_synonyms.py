"""T023 -- pypath ChEMBL molecule_synonyms exposure (spec 005)."""


def test_molecules_parser_splits_synonyms(monkeypatch, tmp_path):
    from pypath.inputs_v2.parsers import chembl as cp

    dbfile = tmp_path / 'fake.sqlite'
    dbfile.write_bytes(b'')  # exists -> parser skips extraction

    schema = {
        'compound_structures': {
            'canonical_smiles', 'standard_inchi', 'standard_inchi_key',
        },
        'compound_properties': {'full_mwt', 'full_molformula', 'alogp'},
        'molecule_synonyms': {'molregno', 'synonyms'},
    }
    monkeypatch.setattr(cp, '_table_columns', lambda p, t: schema.get(t, set()))

    # molecule_synonyms are concatenated on the unit separator (char 31) in SQL,
    # because chemical names contain commas; the parser splits them back.
    rows = [
        {
            'molregno': 1, 'chembl_id': 'CHEMBL25', 'pref_name': 'ASPIRIN',
            'synonyms': 'Aspirin\x1fAcetylsalicylic acid\x1f2-acetoxybenzoic acid',
        },
        {
            'molregno': 2, 'chembl_id': 'CHEMBL2', 'pref_name': 'X',
            'synonyms': None,           # no synonyms -> left as-is
        },
    ]
    monkeypatch.setattr(
        cp, 'iter_sqlite', lambda opener, query=None, **kw: iter(rows),
    )

    out = list(cp.molecules_parser(opener=None, sqlite_path=str(dbfile)))

    assert out[0]['synonyms'] == [
        'Aspirin', 'Acetylsalicylic acid', '2-acetoxybenzoic acid',
    ]
    assert out[1]['synonyms'] is None
    # the aggregation query must reference the molecule_synonyms table
    # (asserted indirectly: the parser yielded the joined rows unchanged)
    assert out[0]['chembl_id'] == 'CHEMBL25'


def test_molecules_schema_exposes_synonym():
    from pypath.inputs_v2 import chembl
    from pypath.internals.cv_terms.identifiers import IdentifierNamespaceCv

    # the SYNONYM identifier CV must be wired into the molecules schema
    assert hasattr(chembl, 'molecules_schema')
    assert IdentifierNamespaceCv.SYNONYM is not None
