"""Source contract of the Integrins input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.integrins import FAMILY_NAME, FAMILY_NAMESPACE, resource
from pypath.inputs_v2.parsers.integrins import COLUMNS, iter_family, iter_subunits

# Rows of Table 1 of Takada et al. 2007 as Europe PMC serves it (JATS XML).
ITGA1 = ('ITGA1', 'α1', 'CD49a', 'NM_181501', '5q11.2', 'P56199')
ITGA4 = ('ITGA4', 'α4', 'CD49d, α4 subunit of VLA-4', 'NM_000885', '2q31.3', 'AAB25486')
ITGA7 = ('ITGA7', 'α7', '', 'NM_002206', '12q13', 'Q86W93')
ITGB1 = ('ITGB1', 'β1', 'Fibronectin receptor β, CD29, MDF2, MSK12', 'BC020057', '10p11.2', 'P05556')


def cell(text, *, italic=False, link=False):
    if not text:
        return '<td colspan="1" rowspan="1"/>'
    if italic:
        text = f'<italic>{text}</italic>'
    if link:
        text = f'<ext-link ext-link-type="gen" xlink:href="{text}">{text}</ext-link>'
    return f'<td align="left" colspan="1" rowspan="1">{text}</td>'


def article(rows, caption='Human integrin subunits', header=COLUMNS, table_id='T1'):
    head = ''.join(f'<td align="left">{name}</td>' for name in header)
    body = ''.join(
        '<tr>'
        + cell(r[0], italic=True)
        + ''.join(cell(c) for c in r[1:3])
        + cell(r[3])
        + cell(r[4])
        + cell(r[5])
        + '</tr>'
        for r in rows
    )
    return SimpleNamespace(
        result=io.BytesIO(
            (
                '<?xml version="1.0"?><article xmlns:xlink="http://www.w3.org/1999/xlink">'
                f'<body><table-wrap position="float" id="{table_id}"><label>Table 1</label>'
                f'<caption><p>{caption}</p></caption><table><thead><tr>{head}</tr></thead>'
                f'<tbody>{body}</tbody></table></table-wrap></body></article>'
            ).encode('utf-8')
        )
    )


def family(*rows, **kwargs):
    (record,) = iter_family(article(rows or [ITGA1, ITGA4, ITGA7, ITGB1], **kwargs))
    return record


def extract():
    record = family()
    out = SilverExtractor('integrins', 'families')
    out.process_record(resource.families.mapper(record), record, 'fixture:1', 1)
    return out


def test_the_dataset_of_families_is_called_families():
    assert list(resource.datasets()) == ['families']


def test_parser_reads_the_subunits_of_the_table():
    subunits = list(iter_subunits(article([ITGA1, ITGA7])))
    assert subunits == [
        {
            'symbol': 'ITGA1', 'name': 'α1', 'synonyms': 'CD49a',
            'gene_accession': 'NM_181501', 'location': '5q11.2', 'accession': 'P56199',
        },
        {
            'symbol': 'ITGA7', 'name': 'α7', 'synonyms': None,
            'gene_accession': 'NM_002206', 'location': '12q13', 'accession': 'Q86W93',
        },
    ]


def test_the_family_has_the_subunits_as_members():
    out = extract()
    (group,) = [e for e in out.entities.values() if e.entity_type == 'protein_family']
    # A namespace of the source, not the bare name: another resource's family
    # called "Integrins" must not merge with this one.
    assert {(i['ns'], i['id']) for i in group.identifiers} == {
        (FAMILY_NAMESPACE, FAMILY_NAME),
        (Namespace.NAME, FAMILY_NAME),
    }
    assert (group.namespace, group.identifier) == ('integrins_family', 'Integrins')
    assert group.taxon == '9606'
    members = [
        (out.entities[r.object_entity_key], r)
        for r in out.relations
        if r.predicate == 'has_member' and r.subject_entity_key == group.entity_key
    ]
    assert len(members) == 4
    by_symbol = {
        next(i['id'] for i in m.identifiers if i['ns'] == Namespace.GENESYMBOL): m
        for m, _ in members
    }
    assert {m.entity_type for m in by_symbol.values()} == {'protein'}
    assert {m.taxon for m in by_symbol.values()} == {'9606'}
    assert {(i['ns'], i['id']) for i in by_symbol['ITGA1'].identifiers} == {
        (Namespace.UNIPROT, 'P56199'),
        (Namespace.GENESYMBOL, 'ITGA1'),
    }


def test_a_subunit_without_a_uniprot_accession_keeps_its_genbank_accession():
    out = extract()
    itga4 = next(
        e for e in out.entities.values()
        if (Namespace.GENBANK, 'AAB25486') in {(i['ns'], i['id']) for i in e.identifiers}
    )
    assert (Namespace.GENESYMBOL, 'ITGA4') in {(i['ns'], i['id']) for i in itga4.identifiers}
    assert Namespace.UNIPROT not in {i['ns'] for i in itga4.identifiers}


def test_nothing_is_asserted_between_the_subunits_or_for_pairs():
    out = extract()
    assert {r.predicate for r in out.relations} == {'has_member'}
    # The alpha and beta subunits are not nested: the table does not group them.
    assert sum(1 for e in out.entities.values() if e.entity_type == 'protein_family') == 1


def test_parser_raises_on_another_table_header_or_caption_and_an_empty_table():
    with pytest.raises(ValueError, match='no table'):
        family(table_id='T2')
    with pytest.raises(ValueError, match='not'):
        family(caption='Mouse integrin subunits')
    with pytest.raises(ValueError, match='header'):
        family(header=COLUMNS[::-1])
    with pytest.raises(ValueError, match='occurs twice'):
        family(ITGA1, ITGA1)
    with pytest.raises(ValueError, match='no rows'):
        list(iter_subunits(article([])))
    with pytest.raises(ValueError, match='No integrin article'):
        list(iter_subunits(SimpleNamespace(result=None)))
