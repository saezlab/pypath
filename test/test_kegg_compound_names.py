"""T021 -- pypath KEGG compound names input (spec 005)."""


def test_kegg_compound_names_parsing(monkeypatch):
    from pypath.inputs import kegg
    from pypath.inputs import kegg_api

    sample = [
        ['cpd:C00001', 'H2O; Water'],
        ['cpd:C00002', "ATP; Adenosine 5'-triphosphate"],
        ['cpd:C00031', 'D-Glucose; Grape sugar; Dextrose; Glucose'],
        ['cpd:C99999', ''],          # no names -> skipped
        ['cpd:C88888'],              # malformed row -> skipped
    ]
    monkeypatch.setattr(kegg_api, '_kegg_list', lambda db: sample)

    out = kegg.kegg_compound_names()

    assert out['C00001'] == ['H2O', 'Water']                 # id de-prefixed
    assert out['C00002'][0] == 'ATP'
    assert 'D-Glucose' in out['C00031'] and len(out['C00031']) == 4
    assert 'C99999' not in out                               # empty skipped
    assert 'C88888' not in out                               # malformed skipped
