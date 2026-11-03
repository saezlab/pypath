"""PDB complexes from the SIFTS chain-to-UniProt table: reader, parsed rows and mapping."""

from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.parsers.pdb import COLUMNS, iter_pdb_entries
from pypath.inputs_v2.pdb import map_complex, resource

RELEASE = '# 2026/10/04 - 08:26 | PDB: 40.26 | UniProt: 2026.03'

# Rows of pdb_chain_uniprot.tsv.gz (release 2026/10/04), as written: CRLF line ends and
# empty cells where a residue number is missing.
TABLE = [
    # one protein, one chain
    '101m\tA\tP02185\t1\t154\t0\t153\t1\t154',
    # three proteins on three chains
    '1a02\tF\tP01100\t1\t56\t\t\t138\t193',
    '1a02\tJ\tP05412\t1\t56\t\t318\t253\t308',
    '1a02\tN\tQ13469\t15\t301\t\t678\t392\t678',
    # two proteins on two chains
    '10ed\tA\tQ88DX3\t1\t598\t1\t\t1\t598',
    '10ed\tX\tA0A140FWQ6\t1\t58\t\t\t1\t58',
    # a protein that is on a chimeric chain and a protein on another chain
    '1fav\tA\tP03069\t1\t29\t\t29\t252\t280',
    '1fav\tA\tP03377\t30\t79\t30\t79\t546\t595',
    '1fav\tC\tP04578\t4\t33\t125\t\t636\t665',
    # one chain made of two proteins: a fusion
    '10zf\tA\tP08515\t9\t225\t\t\t2\t218',
    '10zf\tA\tP22681\t237\t625\t47\t435\t47\t435',
    # two chains, each the same fusion
    '11af\tA\tP68431\t1\t12\t1\t\t2\t13',
    '11af\tA\tQ92613\t18\t190\t\t\t198\t370',
    '11af\tB\tP68431\t1\t12\t1\t\t2\t13',
    '11af\tB\tQ92613\t18\t190\t\t\t198\t370',
    # one protein on two chains
    '10ah\tA\tA7YT55\t6\t364\t\t\t288\t646',
    '10ah\tB\tA7YT55\t6\t364\t\t\t288\t646',
]


def text(rows=TABLE, comment=RELEASE, header='\t'.join(COLUMNS)):
    lines = ([comment] if comment else []) + ([header] if header else []) + rows
    return ''.join(line + '\r\n' for line in lines)


def opener(content=None):
    return SimpleNamespace(result={'pdb_chain_uniprot.tsv': content if content is not None else text()})


def entries(**kwargs):
    return {e['pdb_id']: e for e in iter_pdb_entries(opener(), **kwargs)}


def extract(row):
    record = map_complex(row)
    assert record is not None
    out = SilverExtractor('pdb', 'complexes')
    out.process_record(record, row, 'fixture:1', 1)
    return record, out


def test_module_imports_with_native_schema():
    assert set(resource.datasets()) == {'complexes'}
    assert resource.config.name == 'PDB'


def test_reader_yields_entries_with_distinct_proteins_on_different_chains():
    found = entries()
    assert list(found) == ['1a02', '10ed', '1fav']
    assert found['1a02']['accessions'] == ['P01100', 'P05412', 'Q13469']
    assert found['1a02']['chains'] == ['F', 'J', 'N']
    assert found['1a02']['release'] == RELEASE.lstrip('# ')
    # the chimeric chain is kept visible in the parsed row
    assert found['1fav']['accessions'] == ['P03069', 'P03377', 'P04578']
    assert found['1fav']['chimeric_chains'] == ['A']


def test_reader_leaves_out_fusions_homo_oligomers_and_single_proteins():
    everything = entries(complexes_only=False)
    assert len(everything) == 7
    for pdb_id in ('101m', '10zf', '11af', '10ah'):
        assert pdb_id in everything and pdb_id not in entries()
    assert everything['10ah']['accessions'] == ['A7YT55'] and everything['10ah']['chains'] == ['A', 'B']
    assert everything['10zf']['accessions'] == ['P08515', 'P22681']
    assert everything['11af']['chimeric_chains'] == ['A', 'B']


def test_reader_rejects_malformed_tables():
    with pytest.raises(ValueError, match='release line'):
        list(iter_pdb_entries(opener(text(comment=None))))
    with pytest.raises(ValueError, match='Unexpected SIFTS header'):
        list(iter_pdb_entries(opener(text(header='PDB\tCHAIN\tUNIPROT'))))
    with pytest.raises(ValueError, match='no entries'):
        list(iter_pdb_entries(opener(text([]))))
    with pytest.raises(ValueError, match='another number of fields'):
        list(iter_pdb_entries(opener(text(['1abc\tA\tP02185']))))
    with pytest.raises(ValueError, match='Malformed'):
        list(iter_pdb_entries(opener(text(['1abc\tA\tTP53\t1\t2\t1\t2\t1\t2']))))
    with pytest.raises(ValueError, match='not contiguous'):
        list(iter_pdb_entries(opener(text([TABLE[0], TABLE[1], TABLE[0]]))))


def test_complex_is_the_pdb_entry_with_its_proteins_as_given():
    _, out = extract(entries()['1a02'])
    (cplx,) = [e for e in out.entities.values() if e.entity_type == 'macromolecular_complex']
    assert (cplx.namespace, cplx.identifier) == (Namespace.PDB, '1a02')
    assert cplx.annotations == []
    members = [e for e in out.entities.values() if e.entity_type == 'protein']
    assert sorted(m.identifier for m in members) == ['P01100', 'P05412', 'Q13469']
    # the table does not state a species
    assert all(m.namespace == Namespace.UNIPROT and m.taxon is None for m in members)
    assert [r.predicate for r in out.relations] == ['has_member'] * 3


def test_no_stoichiometry_or_pairs_are_asserted():
    _, out = extract(entries()['1fav'])
    assert len(out.relations) == 3
    assert {r.predicate for r in out.relations} == {'has_member'}
    assert all(not r.annotations for r in out.relations)


def test_entries_with_the_same_proteins_stay_two_complexes():
    first = entries()['10ed']
    second = {**first, 'pdb_id': '9zzz'}
    out = SilverExtractor('pdb', 'complexes')
    for i, row in enumerate((first, second)):
        out.process_record(map_complex(row), row, f'fixture:{i}', i)
    complexes = [e for e in out.entities.values() if e.entity_type == 'macromolecular_complex']
    assert sorted(e.identifier for e in complexes) == ['10ed', '9zzz']
    assert len(out.relations) == 4


def test_entry_with_one_protein_gives_no_record():
    assert map_complex(entries(complexes_only=False)['10ah']) is None
