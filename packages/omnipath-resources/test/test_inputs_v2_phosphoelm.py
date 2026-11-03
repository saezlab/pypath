"""Phospho.ELM kinase-substrate sites: parser, mapper and extractor on source rows."""

import io
import tarfile
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import phosphoelm
from pypath.inputs_v2.parsers.phosphoelm import (
    COLUMNS,
    iter_phosphoelm,
    kinase_accessions,
)

HEADER = '\t'.join(COLUMNS) + '\n'
# Lines copied from phosphoELM_vertebrate_2015-04.dump. A sequence of more than
# 130 residues is cut a few residues after the site (the file has the whole
# protein). In file order of the kinds: a kinase with an accession, a site
# without a kinase and with no publication ("N.N."), a line with an empty
# publication, a high-throughput line with a kinase, a group, an isoform, an
# Ensembl protein of human and of mouse, a mouse kinase site, a group with a
# slash, a site without a kinase of another species, and a repeated line.
ROWS = (
    (
        'P61925',
        'MTDVETTYADFIASGRTGRRNAIHDILVSSASGNSNELALKLAGLDINKTEGEEDAQRSSTEQSGEAQGEAAKSES',
        '8', 'Y', '1956339', 'EGFR', 'LTP', 'Homo sapiens', '2004-12-31 00:00:00+01',
    ),
    (
        'P00740',
        'MQRVNMIMAESPGLITICLLGYLLSAECTVFLDHENANKILNRPKRYNSGKLEEFVQGNLERECMEEKCSFEEAREVFENTERTTEFWKQYVDGDQCESNPCLNGGSCKDDINSYEC',
        '114', 'S', 'N.N.', '', 'LTP', 'Homo sapiens', '2004-12-31 00:00:00+01',
    ),
    (
        'P04637',
        'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAM',
        '37', 'S', '', 'DNA-PK', 'LTP', 'Homo sapiens', '2005-11-21 09:44:58.784817+01',
    ),
    (
        'Q9UJX3',
        'MNVIDHVRDMAAAGLHSNVR',
        '17', 'S', '14657031', 'PLK1', 'HTP', 'Homo sapiens', '2005-11-21 09:38:36.993443+01',
    ),
    (
        'P62807',
        'MPEPAKSAPAPKKGSKKAVTKAQKKDGKKRKRSRKESYSVYVYKVLKQVHPDTGISSKAMGIMNSFVNDIFERIAGEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSK',
        '37', 'S', '192223', 'PKA_group', 'LTP', 'Homo sapiens', '2004-12-31 00:00:00+01',
    ),
    (
        'O43521-2',
        'MAKQPSDVSSECDREGRQLQPAERPPQLRPGAPTSLQTEPQDRSPAPMSCDKSTQTPSP',
        '56', 'T', '12591950', 'MAPK9', 'LTP', 'Homo sapiens', '2006-04-07 10:18:23.324613+02',
    ),
    (
        'ENSP00000267569',
        'MMPGQIPDPSVTTGSLPGLGPLTGLPSSALTVEELKYADIRNLGAMIAPLHFLEVKLGKRPQPVKSELDEEEERRKRRREKNKVAAARCRNKKKERTEFLQRESERLELMNAELKTQIEELKQERQQLILMLNRHRPTCIVRTDSVKTPES',
        '148', 'T', '12225289', 'MAPK9', 'LTP', 'Homo sapiens', '2006-10-17 12:06:37.35763+02',
    ),
    (
        'ENSMUSP00000043430',
        'MMPGQIPDPSVTAGSLPGLGPLTGLPSSALTTEELKYADIRNIGAMIAPLHFLEVKLGKRPQPVKSELDEEEERRKRRREKNKVAAARCRNKKKERTEFLQRESERLELMNAELKTQIEELKLERQQLILMLNRHRPTCIVRTDSVRTPES',
        '148', 'T', '11602244', 'MAPK8', 'LTP', 'Mus musculus', '2007-07-13 15:16:00.704002+02',
    ),
    (
        'P61014',
        'MEKVQYLTRSAIRRASTIEMPQQARQNLQNLFINFCLILICLLLICIIVMLL',
        '17', 'T', '2544595', 'CaM-KII_alpha', 'LTP', 'Mus musculus', '2007-07-13 15:13:51.834628+02',
    ),
    (
        'P05204',
        'MPKRKAEGDAKGDKAKVKDEPQRRSARLSAKPAPPKPEPKPKKAPAKKGEKVPKGKKGKADAGKEGNNPAENGDAKTDQAQKAEGAGDAK',
        '29', 'S', '6279643', 'PKG/cGK_group', 'LTP', 'Homo sapiens', '2004-12-31 00:00:00+01',
    ),
    (
        'P17306',
        'MSTSRKLKSHGMRRGKNRAPHKGVKRGGSKRKYRKGSLKSRKRCDDANRNYRSHL',
        '40', 'S', '8019440', '', 'LTP', 'Sus scrofa', '2004-12-31 00:00:00+01',
    ),
)
# The file repeats some lines; this one is the group line above.
ROWS = ROWS + (ROWS[4],)

# Rows of the table of http://phospho.elm.eu.org/kinases.html.
KINASE_PAGE = """<html><body>
<table id="statistics"><tr><td>Kinases</td><td>299</td></tr></table>
<table id="kinaselist" class="results sortable" width="500px;">
<thead>
<tr class="clickable" > <th>Name</th> <th>#Instances</th> <th>Uniprot</th> <th>HUGO&nbsp;Id</th> <th>Classification</th> <th style="text-align:left; width=100px;">Description</th></tr>
</thead>
<tbody>
<tr> <td style="text-align:left;"><a name="Abl"></a><a href="/byKinase/Abl.html">Abl</a></td> <td>48</td> <td><a href="http://www.uniprot.org/uniprot/P00519" target="_blank">P00519</a></td> <td><a href="http://www.genenames.org/data/hgnc_data.php?hgnc_id=76" target="_blank">76</a></td> <td>TK</td> <td style="text-align:left;">Proto-oncogene tyrosine-protein kinase ABL1</td> </tr>
<tr> <td style="text-align:left;"><a name="CaM-KII_alpha"></a><a href="/byKinase/CaM-KII_alpha.html">CaM-KII_alpha</a></td> <td>36</td> <td><a href="http://www.uniprot.org/uniprot/Q9UQM7" target="_blank">Q9UQM7</a></td> <td><a href="http://www.genenames.org/data/hgnc_data.php?hgnc_id=1460" target="_blank">1460</a></td> <td>CAMK</td> <td style="text-align:left;">calcium/calmodulin-dependent protein kinase (CaM kinase) II alpha</td> </tr>
<tr> <td style="text-align:left;"><a name="DNA-PK"></a><a href="/byKinase/DNA-PK.html">DNA-PK</a></td> <td>21</td> <td><a href="http://www.uniprot.org/uniprot/P78527" target="_blank">P78527</a></td> <td><a href="http://www.genenames.org/data/hgnc_data.php?hgnc_id=9413" target="_blank">9413</a></td> <td>Atypical</td> <td style="text-align:left;">DNA-dependent protein kinase catalytic subunit</td> </tr>
<tr> <td style="text-align:left;"><a name="EGFR"></a><a href="/byKinase/EGFR.html">EGFR</a></td> <td>55</td> <td><a href="http://www.uniprot.org/uniprot/P00533" target="_blank">P00533</a></td> <td><a href="http://www.genenames.org/data/hgnc_data.php?hgnc_id=3236" target="_blank">3236</a></td> <td>TK</td> <td style="text-align:left;">Epidermal growth factor receptor (Receptor tyrosine-protein kinase ErbB-1)</td> </tr>
<tr> <td style="text-align:left;"><a name="MAPK8"></a><a href="/byKinase/MAPK8.html">MAPK8</a></td> <td>33</td> <td><a href="http://www.uniprot.org/uniprot/P45983" target="_blank">P45983</a></td> <td><a href="http://www.genenames.org/data/hgnc_data.php?hgnc_id=6881" target="_blank">6881</a></td> <td>CMGC</td> <td style="text-align:left;">mitogen-activated protein kinase 8 (JNK1)</td> </tr>
<tr> <td style="text-align:left;"><a name="MAPK9"></a><a href="/byKinase/MAPK9.html">MAPK9</a></td> <td>15</td> <td><a href="http://www.uniprot.org/uniprot/P45984" target="_blank">P45984</a></td> <td><a href="http://www.genenames.org/data/hgnc_data.php?hgnc_id=6886" target="_blank">6886</a></td> <td>CMGC</td> <td style="text-align:left;">mitogen-activated protein kinase 9 (JNK3</td> </tr>
<tr> <td style="text-align:left;"><a name="PKA_group"></a><a href="/byKinase/PKA_group.html">PKA_group</a></td> <td>349</td> <td></td> <td></td> <td></td> <td style="text-align:left;"></td> </tr>
<tr> <td style="text-align:left;"><a name="PKG/cGK_group"></a><a href="/byKinase/PKG/cGK_group.html">PKG/cGK_group</a></td> <td>20</td> <td></td> <td></td> <td></td> <td style="text-align:left;"></td> </tr>
<tr> <td style="text-align:left;"><a name="PLK1"></a><a href="/byKinase/PLK1.html">PLK1</a></td> <td>42</td> <td><a href="http://www.uniprot.org/uniprot/P53350" target="_blank">P53350</a></td> <td><a href="http://www.genenames.org/data/hgnc_data.php?hgnc_id=9077" target="_blank">9077</a></td> <td>PLK:PLK1</td> <td style="text-align:left;">Polo like kinase 1</td> </tr>
</tbody></table></body></html>
"""
KINASES = SimpleNamespace(
    open=lambda force_refresh=False: SimpleNamespace(result=KINASE_PAGE)
)


def dump(rows=ROWS, header=HEADER, member='phosphoELM_vertebrate_2015-04.dump'):
    text = header + ''.join('\t'.join(row) + '\n' for row in rows)
    return SimpleNamespace(result={member: io.BytesIO(text.encode())})


def parse(rows=ROWS, kinases=KINASES, **kwargs):
    return list(
        iter_phosphoelm(
            dump(rows, **kwargs),
            organisms=frozenset(phosphoelm.ORGANISMS),
            kinases_download=kinases,
        )
    )


def extract(rows):
    out = SilverExtractor('phosphoelm', 'enzyme_substrate')
    for index, row in enumerate(rows):
        record = phosphoelm.resource.enzyme_substrate.mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def annotations(relation, term, scope='relation'):
    return [
        a['value']
        for a in relation.annotations
        if a['term'] == term and a['scope'] == scope
    ]


def site(out, relation):
    return out.entities[relation.object_entity_key].molecular_form['modifications'][0]


def test_parser_reads_the_lines_with_a_kinase_and_leaves_out_the_sequence():
    rows = parse()
    # Ten of the twelve lines name a kinase; the last one repeats the fifth.
    assert [(r['accession'], r['position'], r['enzyme']) for r in rows] == [
        ('P61925', 8, 'EGFR'),
        ('P04637', 37, 'DNA-PK'),
        ('Q9UJX3', 17, 'PLK1'),
        ('P62807', 37, 'PKA_group'),
        ('O43521-2', 56, 'MAPK9'),
        ('ENSP00000267569', 148, 'MAPK9'),
        ('ENSMUSP00000043430', 148, 'MAPK8'),
        ('P61014', 17, 'CaM-KII_alpha'),
        ('P05204', 29, 'PKG/cGK_group'),
        ('P62807', 37, 'PKA_group'),
    ]
    assert rows[0] == {
        'acc': 'P61925',
        'position': 8,
        'code': 'Y',
        'pmids': '1956339',
        'kinases': 'EGFR',
        'source': 'LTP',
        'species': 'Homo sapiens',
        'entry_date': '2004-12-31 00:00:00+01',
        'accession': 'P61925',
        'accession_database': 'uniprot',
        'pmid': '1956339',
        'enzyme': 'EGFR',
        'enzyme_is_group': False,
        'enzyme_uniprot': 'P00533',
    }
    assert [(r['accession_database']) for r in rows[4:6]] == ['uniprot', 'ensembl_protein']
    assert [r['enzyme_is_group'] for r in rows[2:4]] == [False, True]
    # A group has no accession in the kinase list.
    assert rows[3]['enzyme_uniprot'] is None
    # Without the kinase list no accession is added.
    assert parse(kinases=None)[0]['enzyme_uniprot'] is None


def test_parser_can_yield_the_sites_without_a_kinase():
    rows = list(
        iter_phosphoelm(
            dump(),
            organisms=frozenset(phosphoelm.ORGANISMS),
            kinase_required=False,
        )
    )
    assert len(rows) == 12
    unassigned = [r for r in rows if r['enzyme'] is None]
    assert [(r['accession'], r['pmid'], r['species']) for r in unassigned] == [
        ('P00740', None, 'Homo sapiens'),
        ('P17306', '8019440', 'Sus scrofa'),
    ]
    assert not any(r['enzyme_is_group'] or r['enzyme_uniprot'] for r in unassigned)


def test_publication_is_a_pubmed_id_or_nothing():
    rows = {r['accession'] + str(r['position']): r for r in parse()}
    assert rows['P0463737']['pmid'] is None
    assert rows['P0463737']['pmids'] == ''
    changed = tuple(
        cell if index != 4 else 'N.N.' for index, cell in enumerate(ROWS[0])
    )
    (row,) = parse([changed])
    assert (row['pmids'], row['pmid']) == ('N.N.', None)
    # A number with a leading zero or nine digits is no PubMed ID.
    for number in ('0123456', '123456789'):
        changed = tuple(
            cell if index != 4 else number for index, cell in enumerate(ROWS[0])
        )
        (row,) = parse([changed])
        assert (row['pmids'], row['pmid']) == (number, None)


def test_parser_rejects_a_changed_file_and_unknown_values():
    row = list(ROWS[0])

    def changed(index, value):
        return [tuple(value if i == index else c for i, c in enumerate(row))]

    with pytest.raises(ValueError, match='columns'):
        parse(header=HEADER.replace('sequence', 'seq'))
    with pytest.raises(ValueError, match='fields'):
        parse(header='acc\tsequence\tposition\n')
    with pytest.raises(ValueError, match='HTML'):
        parse(header='<html>\t<head><title>404 Not Found</title></head>\n')
    with pytest.raises(ValueError, match='no rows'):
        parse([])
    with pytest.raises(ValueError, match='not an archive'):
        list(
            iter_phosphoelm(
                SimpleNamespace(result=io.StringIO(HEADER)), organisms=frozenset()
            )
        )
    with pytest.raises(ValueError, match='dump files'):
        parse(member='readme.txt')
    with pytest.raises(ValueError, match='Unknown Phospho.ELM species'):
        parse(changed(7, 'Danio rerio'))
    with pytest.raises(ValueError, match='accession not understood'):
        parse(changed(0, 'NP_005219'))
    with pytest.raises(ValueError, match='is not at position'):
        parse(changed(3, 'S'))
    with pytest.raises(ValueError, match='is not at position'):
        parse(changed(2, '9999'))
    with pytest.raises(ValueError, match='position not understood'):
        parse(changed(2, '0'))
    with pytest.raises(ValueError, match='throughput'):
        parse(changed(6, 'MTP'))
    with pytest.raises(ValueError, match='reference not understood'):
        parse(changed(4, '1956339;2544595'))
    with pytest.raises(ValueError, match='date not understood'):
        parse(changed(8, '2004-12-31'))
    # A line with a field fewer is an error as well.
    with pytest.raises(ValueError):
        list(
            iter_phosphoelm(
                SimpleNamespace(
                    result={
                        'phosphoELM_vertebrate_2015-04.dump': io.BytesIO(
                            (HEADER + 'P61925\tMTDV\t1\tM\n').encode()
                        )
                    }
                ),
                organisms=frozenset(phosphoelm.ORGANISMS),
            )
        )


def test_kinase_list_gives_the_accessions_of_the_names_that_have_one():
    accessions = kinase_accessions(SimpleNamespace(result=KINASE_PAGE))
    assert accessions == {
        'Abl': 'P00519',
        'CaM-KII_alpha': 'Q9UQM7',
        'DNA-PK': 'P78527',
        'EGFR': 'P00533',
        'MAPK8': 'P45983',
        'MAPK9': 'P45984',
        'PLK1': 'P53350',
    }
    for page, message in (
        ('<html>404</html>', 'kinaselist'),
        (KINASE_PAGE.replace('<th>Uniprot</th>', '<th>Accession</th>'), 'other columns'),
        (KINASE_PAGE.replace('>MAPK9</a>', '>MAPK8</a>'), 'repeated'),
        (KINASE_PAGE.replace('>P78527<', '>NP_000'), 'accession'),
    ):
        with pytest.raises(ValueError, match=message):
            kinase_accessions(SimpleNamespace(result=page))


def test_kinase_affects_the_phosphorylation_of_the_substrate_at_the_site():
    record = phosphoelm.resource.enzyme_substrate.mapper(parse()[0])
    assert slot_name(record.predicate) == 'affects'
    # The kinase by its Phospho.ELM name, with the accession of the kinase list.
    assert record.subject.identifiers == [
        (Namespace.UNIPROT, 'P00533'),
        (Namespace.NAME, 'EGFR'),
    ]
    assert record.object.identifiers == [(Namespace.UNIPROT, 'P61925')]
    assert record.subject.molecular_form is None
    assert record.object.molecular_form['modifications'] == [
        {
            'term': 'phosphorylated residue',
            'residue': 'Y',
            'position': 8,
            'end_position': 8,
            'coordinate_reference': {
                'identifier': {'ns': 'uniprot', 'id': 'P61925'},
                'coordinate_system': 'protein',
                'position_base': 1,
            },
            'description': None,
        }
    ]

    out = extract(parse()[:1])
    (relation,) = out.relations
    assert relation.predicate == 'affects'
    assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']
    assert annotations(relation, 'object_direction_qualifier') == ['increased']
    assert annotations(relation, 'publications') == ['PMID:1956339']
    assert annotations(relation, 'phosphoelm:throughput') == ['Low Throughput']
    # The species is the substrate's; the dump does not give the kinase's.
    assert annotations(relation, 'in_taxon', scope='object') == ['NCBITaxon:9606']
    assert annotations(relation, 'in_taxon', scope='subject') == []
    kinase = out.entities[relation.subject_entity_key]
    assert (kinase.entity_type, kinase.namespace, kinase.identifier) == (
        'protein',
        Namespace.UNIPROT,
        'P00533',
    )
    assert kinase.label == 'EGFR'
    substrate = out.entities[relation.object_entity_key]
    assert substrate.molecular_form['sequence_identifiers'] is None


def test_kinase_without_an_accession_is_identified_by_its_name():
    rows = parse(kinases=None)
    out = extract(rows[:1])
    kinase = out.entities[out.relations[0].subject_entity_key]
    assert (kinase.entity_type, kinase.namespace, kinase.identifier) == (
        'protein',
        Namespace.NAME,
        'EGFR',
    )


def test_kinase_group_is_a_family_and_is_not_expanded():
    rows = parse()
    out = extract([rows[3], rows[8]])
    subjects = [out.entities[r.subject_entity_key] for r in out.relations]
    assert [(e.entity_type, e.namespace, e.identifier, e.taxon) for e in subjects] == [
        ('protein_family', Namespace.NAME, 'PKA_group', None),
        ('protein_family', Namespace.NAME, 'PKG/cGK_group', None),
    ]
    # Only the family and the substrates: no member of a group is invented.
    assert len(out.entities) == 4


def test_isoform_and_ensembl_substrates_keep_their_accession_and_sequence():
    rows = parse()
    out = extract(rows[4:7])
    substrates = [out.entities[r.object_entity_key] for r in out.relations]
    assert [(s.namespace, s.identifier, s.taxon) for s in substrates] == [
        (Namespace.UNIPROT, 'O43521-2', '9606'),
        (Namespace.ENSP, 'ENSP00000267569', '9606'),
        (Namespace.ENSP, 'ENSMUSP00000043430', '10090'),
    ]
    assert substrates[0].molecular_form['isoform_identifier'] == {
        'ns': 'uniprot',
        'id': 'O43521-2',
    }
    assert substrates[1].molecular_form['sequence_identifiers'] == [
        {'ns': 'ensp', 'id': 'ENSP00000267569'}
    ]
    assert [site(out, r)['coordinate_reference']['identifier'] for r in out.relations] == [
        {'ns': 'uniprot', 'id': 'O43521-2'},
        {'ns': 'ensp', 'id': 'ENSP00000267569'},
        {'ns': 'ensp', 'id': 'ENSMUSP00000043430'},
    ]
    # The same kinase on the same position of two sequences: two substrates.
    assert substrates[1].entity_key != substrates[2].entity_key


def test_high_throughput_lines_are_kept_with_their_throughput():
    out = extract(parse())
    values = [annotations(r, 'phosphoelm:throughput') for r in out.relations]
    assert values.count(['High Throughput']) == 1
    assert values.count(['Low Throughput']) == 9
    htp = out.relations[2]
    assert out.entities[htp.subject_entity_key].label == 'PLK1'
    assert annotations(htp, 'publications') == ['PMID:14657031']


def test_line_without_a_publication_has_none_and_nn_is_not_a_pmid():
    out = extract(parse()[1:2])
    (relation,) = out.relations
    assert annotations(relation, 'publications') == []
    assert annotations(relation, 'phosphoelm:throughput') == ['Low Throughput']


def test_repeated_lines_are_separate_observations_of_one_relation():
    out = extract(parse())
    assert len(out.relations) == 10
    assert len({r.relation_key for r in out.relations}) == 9
    assert out.relations[3].relation_key == out.relations[9].relation_key


def test_species_crosswalk_is_exact():
    assert phosphoelm.ORGANISMS == {
        'Homo sapiens': '9606',
        'Mus musculus': '10090',
        'Rattus norvegicus': '10116',
        'Bos taurus': '9913',
        'Gallus gallus': '9031',
        'Oryctolagus cuniculus': '9986',
        'Xenopus laevis': '8355',
        'Mesocricetus auratus': '10036',
        'Sus scrofa': '9823',
        'Ovis aries': '9940',
        'Canis familiaris': '9615',
        'Didelphis virginiana': '9267',
        'Cricetulus griseus': '10029',
    }
    out = extract(parse()[5:7])
    assert [e.taxon for e in out.entities.values() if e.identifier.startswith('ENS')] == [
        '9606',
        '10090',
    ]


def test_unknown_throughput_raises_and_nothing_is_invented():
    row = dict(parse()[0], source='MTP')
    with pytest.raises(ValueError, match='throughput'):
        phosphoelm.resource.enzyme_substrate.mapper(row)
    out = extract(parse())
    names = {e.identifier for e in out.entities.values()}
    # No site without a kinase reaches the output.
    assert 'P00740' not in names and 'P17306' not in names
    for relation in out.relations:
        # Phosphorylation is not an effect on activity.
        assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']
        assert annotations(relation, 'object_direction_qualifier') == ['increased']
        assert annotations(relation, 'in_taxon', scope='subject') == []
    # Kinase names are not turned into gene symbols.
    assert not [i for e in out.entities.values() for i in e.identifiers if i['ns'] == Namespace.GENESYMBOL]


def test_archive_member_is_found_by_its_name_in_a_real_archive():
    data = (HEADER + '\t'.join(ROWS[0]) + '\n').encode()
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode='w:gz') as archive:
        info = tarfile.TarInfo('phosphoELM_vertebrate_2015-04.dump')
        info.size = len(data)
        archive.addfile(info, io.BytesIO(data))
    buffer.seek(0)
    with tarfile.open(fileobj=buffer, mode='r:gz') as archive:
        handle = archive.extractfile('phosphoELM_vertebrate_2015-04.dump')
        (row,) = iter_phosphoelm(
            SimpleNamespace(result={'phosphoELM_vertebrate_2015-04.dump': handle}),
            organisms=frozenset(phosphoelm.ORGANISMS),
        )
    assert (row['accession'], row['enzyme']) == ('P61925', 'EGFR')
