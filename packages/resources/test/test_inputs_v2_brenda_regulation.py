"""BRENDA regulation: compounds affect enzyme activity, identified by BRENDA's own ligand table."""

from types import SimpleNamespace

from omnipath_build.silver import SilverExtractor
from omnipath_resources.inputs_v2 import brenda
from omnipath_resources.inputs_v2.parsers import brenda as parsers

NAMES = (
    '?id\t?kind\t?name\n'
    '1\t"label"\t"water"\n'
    '1\t"synonym"\t"H2O"\n'
    '2\t"label"\t"pyrazole"\n'
    '3\t"label"\t"2,3-dimethylpyrazine-1-sulfonamide"\n'
    '4\t"synonym"\t"same name"\n'
    '5\t"synonym"\t"same name"\n'
)
STRUCTURES = (
    '?id\t?inchikey\t?inchi\n'
    '2\t"InChIKey=WTKZEGDFNFYCGP-UHFFFAOYSA-N"\t"InChI=1S/C3H4N2/c1-2-4-5-3-1/h1-3H,(H,4,5)"\n'
)
FLAT_FILE = '\n'.join([
    'ID\t1.1.1.1',
    'PR\t#1# Mus musculus {Q9QYY9; source: UniProt} <1>',
    'PR\t#2# Homo sapiens <2>',
    'RF\t<1> Ref one {Pubmed:111}',
    'RF\t<2> Ref two {Pubmed:222}',
    'IN\t#1,2# pyrazole (#1# competitive <1>; #2# weak <2>) <1,2>',
    'IN\t#1# more <1>',
    'AC\t#2# 2,3-dimethylpyrazine-1- sulfonamide <2>',
    'CF\t#1# H2O <1>',
    'KI\t#1# 0.0051 {pyrazole} <1>',
    'KI\t#2# -999 {pyrazole} <2>',
    '///',
    '',
])


def _opener(text):
    return SimpleNamespace(result=text.splitlines(keepends=True))


def _records():
    return list(
        parsers.iter_regulation(
            SimpleNamespace(result={'brenda.txt': _opener(FLAT_FILE).result}),
            names_download=SimpleNamespace(open=lambda **_: _opener(NAMES)),
            structures_download=SimpleNamespace(open=lambda **_: _opener(STRUCTURES)),
        )
    )


def test_names_find_brenda_ligands_only_when_unique():
    compounds = parsers.BrendaCompounds(_opener(NAMES), _opener(STRUCTURES))
    assert compounds.ligand('water') == '1'
    assert compounds.ligand('H2O') == '1'
    # The flat file wraps long names and gains a space where it does.
    assert compounds.ligand('2,3-dimethylpyrazine-1- sulfonamide') == '3'
    assert compounds.ligand('same name') is None
    assert compounds.ligand('unknown') is None
    assert compounds.structures['2'][0] == 'WTKZEGDFNFYCGP-UHFFFAOYSA-N'


def test_each_protein_keeps_its_notes_references_and_ki():
    records = {(r['role'], r['compound_name'], r['protein_record_id']): r for r in _records()}
    assert set(records) == {
        ('inhibitor', 'pyrazole', '1'),
        ('inhibitor', 'pyrazole', '2'),
        ('activator', '2,3-dimethylpyrazine-1- sulfonamide', '2'),
        ('cofactor', 'H2O', '1'),
    }
    mouse = records['inhibitor', 'pyrazole', '1']
    assert (mouse['UniProt'], mouse['organism'], mouse['notes'], mouse['Refs']) == (
        'Q9QYY9', 'Mus musculus', ['competitive'], ['111'],
    )
    assert (mouse['brenda_ligand'], mouse['inchikey']) == ('2', 'WTKZEGDFNFYCGP-UHFFFAOYSA-N')
    assert mouse['ki'] == ['0.0051']
    human = records['inhibitor', 'pyrazole', '2']
    assert (human['UniProt'], human['notes'], human['Refs'], human['ki']) == (None, ['weak'], ['222'], [])
    assert records['activator', '2,3-dimethylpyrazine-1- sulfonamide', '2']['brenda_ligand'] == '3'


def _relation(row):
    ex = SilverExtractor(brenda.config.id.value, 'regulation')
    ex.process_record(brenda.regulation_schema.build(row), row, 'regulation:1', 1)
    (relation,) = ex.relations
    return ex, relation


def test_roles_are_qualified_effects_on_enzyme_activity():
    rows = {(r['role'], r['protein_record_id']): r for r in _records()}
    ex, inhibition = _relation(rows['inhibitor', '1'])
    annotations = {(a['term'], str(a['value'])) for a in inhibition.annotations}
    assert inhibition.predicate == 'affects'
    assert annotations >= {
        ('object_aspect_qualifier', 'activity'),
        ('object_direction_qualifier', 'decreased'),
        ('causal_mechanism_qualifier', 'inhibition'),
        ('brenda:role', 'inhibitor'),
        ('brenda:note', 'competitive'),
        ('publications', 'PMID:111'),
    }
    (ki,) = [a for a in inhibition.annotations if a['term'] == 'BAO:0000192']
    assert ki['quantity']['has_numeric_value'] == 0.0051 and ki['quantity']['has_unit'] == 'mM'
    compound = ex.entities[inhibition.subject_entity_key]
    assert [(str(i['ns']), i['id']) for i in compound.identifiers][:3] == [
        ('inchikey', 'WTKZEGDFNFYCGP-UHFFFAOYSA-N'),
        ('inchi', 'InChI=1S/C3H4N2/c1-2-4-5-3-1/h1-3H,(H,4,5)'),
        ('brenda_ligand', '2'),
    ]
    _, activation = _relation(rows['activator', '2'])
    assert {(a['term'], str(a['value'])) for a in activation.annotations} >= {
        ('object_direction_qualifier', 'increased'),
        ('causal_mechanism_qualifier', 'activation'),
    }
    _, cofactor = _relation(rows['cofactor', '1'])
    terms = {a['term']: str(a['value']) for a in cofactor.annotations}
    assert terms['causal_mechanism_qualifier'] == 'cofactor'
    assert 'object_direction_qualifier' not in terms


def test_trembl_accessions_are_uniprot_and_other_sources_are_not():
    assert parsers.extract_uniprot_ids('#54# Daucus carota {Q2VEY0; source: TREMBL} <44>') == ['Q2VEY0']
    assert parsers.extract_uniprot_ids('#3# Homo sapiens {P00330; source: SwissProt} <1>') == ['P00330']
    assert parsers.extract_uniprot_ids('#7# Bos taurus {GQ403010; source: GenBank} <2>') == []
