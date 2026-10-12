"""MIMP kinase-substrate sites: reader, mapper and extractor on source rows."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import mimp
from pypath.inputs_v2.parsers import mimp as mimp_parser

HEADER = '\t'.join(mimp_parser.COLUMNS) + '\n'
# Rows copied from phosphorylation_data.tab: sites without a kinase, several
# kinases in one cell, a name of a source database, a kinase group, a
# replaced gene symbol, a fusion protein, a histidine and a padded window.
SOURCE = HEADER + (
    'C6orf57\tNM_145267\t\t47\tKSELVKQSLKKPKLP\t1\tPhosphoSitePlus\n'
    'YWHAB\tNM_003404\tPDKC;PRKCD\t60\tVVGARRSSWRVISSI\t8\tHPRD;PhosphoELM;PhosphoSitePlus\n'
    'YWHAB\tNM_003404\tPRKCZ\t132\tGDYFRYLSEVASGDN\t1\tHPRD;PhosphoELM;PhosphoSitePlus\n'
    'AKT1\tNM_005163\tAKT1\t72\tTERPRPNTFIIRCLQ\t1\tPhosphoSitePlus\n'
    'AANAT\tNM_001166579\tPKA_GROUP;PRKACA\t250\tHPFLRRNSGC-----\t1\tHPRD;PhosphoELM;PhosphoSitePlus\n'
    'ADRA2A\tNM_000681\tADRBK1\t311\tDALDLEESSSSDHAE\t1\tHPRD;PhosphoELM;PhosphoSitePlus\n'
    'BCR\tNM_004327\tBCR/ABL;FES\t283\tYQPYQSIYVGGMMEG\t91\tHPRD;PhosphoELM;PhosphoSitePlus\n'
    'NME2\tNM_002512\tNME2\t118\tQVGRNIIHGSDSVKS\t3\tHPRD;PhosphoSitePlus\n'
    'C6orf57\tNM_145267\t\t64\tRFDAPEDSHLEKEPL\t1\tPhosphoSitePlus\n'
)


def parse(text=SOURCE):
    return list(mimp_parser.iter_kinase_sites(SimpleNamespace(result=io.StringIO(text))))


def extract(rows):
    out = SilverExtractor('mimp', 'enzyme_substrate')
    for index, row in enumerate(rows):
        record = mimp.resource.enzyme_substrate.mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def annotations(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def site(out, relation):
    return out.entities[relation.object_entity_key].molecular_form['modifications'][0]


def test_reader_yields_one_row_per_kinase_and_skips_sites_without_one():
    rows = parse()
    # Seven of the nine sites name a kinase; three cells name two.
    assert [(r['gene'], r['enzyme'], r['site_position'], r['residue']) for r in rows] == [
        ('YWHAB', 'PDKC', 60, 'S'),
        ('YWHAB', 'PRKCD', 60, 'S'),
        ('YWHAB', 'PRKCZ', 132, 'S'),
        ('AKT1', 'AKT1', 72, 'T'),
        ('AANAT', 'PKA_GROUP', 250, 'S'),
        ('AANAT', 'PRKACA', 250, 'S'),
        ('ADRA2A', 'ADRBK1', 311, 'S'),
        ('BCR', 'BCR/ABL', 283, 'Y'),
        ('BCR', 'FES', 283, 'Y'),
        ('NME2', 'NME2', 118, 'H'),
    ]
    assert rows[0] == {
        'gene': 'YWHAB',
        'gene_acc': 'NM_003404',
        'kinase': 'PDKC;PRKCD',
        'position': '60',
        'flank': 'VVGARRSSWRVISSI',
        'npmid': '8',
        'source': 'HPRD;PhosphoELM;PhosphoSitePlus',
        'site_position': 60,
        'residue': 'S',
        'sources': ['HPRD', 'PhosphoELM', 'PhosphoSitePlus'],
        'enzyme': 'PDKC',
        'enzyme_genesymbol': None,
        'enzyme_name': 'PDKC',
        'enzyme_group': None,
    }
    kinds = {
        row['enzyme']: [
            kind for kind in ('genesymbol', 'name', 'group') if row[f'enzyme_{kind}']
        ]
        for row in rows
    }
    assert kinds == {
        'PDKC': ['name'],
        'PRKCD': ['genesymbol'],
        'PRKCZ': ['genesymbol'],
        'AKT1': ['genesymbol'],
        'PKA_GROUP': ['group'],
        'PRKACA': ['genesymbol'],
        # The symbol of 2015, now GRK2, stays as written.
        'ADRBK1': ['genesymbol'],
        'BCR/ABL': ['name'],
        'FES': ['genesymbol'],
        'NME2': ['genesymbol'],
    }


def test_reader_rejects_changed_layout_and_malformed_rows():
    with pytest.raises(ValueError, match='columns'):
        parse('gene\tkinase\tposition\nAKT1\tAKT1\t72\n')
    # A file that has the header and no rows is not an empty dataset.
    with pytest.raises(ValueError, match='no rows'):
        parse(HEADER)
    with pytest.raises(ValueError, match='HTML'):
        parse('<html>\n<head><title>404 Not Found</title></head>\n')
    row = '{}\t{}\t{}\t{}\t{}\t1\t{}\n'
    good = ('AKT1', 'NM_005163', 'AKT1', '72', 'TERPRPNTFIIRCLQ', 'PhosphoSitePlus')
    assert len(parse(HEADER + row.format(*good))) == 1
    for column, value in ((1, 'P31749'), (3, 'T72'), (4, 'TERPRPNTFII')):
        bad = list(good)
        bad[column] = value
        with pytest.raises(ValueError, match='Malformed'):
            parse(HEADER + row.format(*bad))
    with pytest.raises(ValueError, match='kinase label'):
        parse(HEADER + row.format(*good[:2], 'CaM-KII', *good[3:]))
    with pytest.raises(ValueError, match='source database'):
        parse(HEADER + row.format(*good[:5], 'NewDB'))


def test_kinase_increases_substrate_phosphorylation_at_the_site():
    row = parse()[2]
    record = mimp.resource.enzyme_substrate.mapper(row)
    assert slot_name(record.predicate) == 'affects'
    assert record.subject.identifiers == [(Namespace.GENESYMBOL, 'PRKCZ')]
    # The gene symbol and the RefSeq transcript, as the source gives them.
    assert record.object.identifiers == [
        (Namespace.GENESYMBOL, 'YWHAB'),
        (Namespace.REFSEQ, 'NM_003404'),
    ]
    assert record.subject.molecular_form is None
    assert record.object.molecular_form['modifications'] == [
        {
            'term': 'phosphorylated residue',
            'residue': 'S',
            'position': 132,
            'end_position': 132,
            'coordinate_reference': {
                'identifier': {'ns': 'refseq', 'id': 'NM_003404'},
                'coordinate_system': 'protein',
                'position_base': 1,
            },
            'description': 'GDYFRYLSEVASGDN',
        }
    ]

    out = extract([row])
    (relation,) = out.relations
    assert relation.predicate == 'affects'
    assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']
    assert annotations(relation, 'object_direction_qualifier') == ['increased']
    assert annotations(relation, 'supporting_data_source') == [
        'HPRD',
        'PhosphoELM',
        'PhosphoSitePlus',
    ]
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert {e.taxon for e in out.entities.values()} == {'9606'}
    substrate = out.entities[relation.object_entity_key]
    assert (substrate.namespace, substrate.identifier) == (Namespace.GENESYMBOL, 'YWHAB')
    # The transcript is the reference of the coordinate, not an isoform claim.
    assert substrate.molecular_form['isoform_identifier'] is None
    assert substrate.molecular_form['sequence_identifiers'] is None


def test_sites_and_kinases_stay_distinct():
    out = extract(parse()[:3])
    assert len(out.relations) == 3
    assert len({r.relation_key for r in out.relations}) == 3
    # Two kinases at S60 of YWHAB; PRKCZ at S132.
    assert out.relations[0].object_entity_key == out.relations[1].object_entity_key
    assert out.relations[1].object_entity_key != out.relations[2].object_entity_key
    assert [site(out, r)['position'] for r in out.relations] == [60, 60, 132]


def test_kinase_group_is_not_expanded_to_its_members():
    rows = parse()
    out = extract(rows[4:6])
    subjects = [out.entities[r.subject_entity_key] for r in out.relations]
    # One relation for the group label and one for the kinase the cell names.
    assert [(e.entity_type, e.namespace, e.identifier, e.taxon) for e in subjects] == [
        ('protein_family', Namespace.NAME, 'PKA_GROUP', None),
        ('protein', Namespace.GENESYMBOL, 'PRKACA', '9606'),
    ]
    assert {e.identifier for e in out.entities.values()} == {'PKA_GROUP', 'PRKACA', 'AANAT'}
    # A window padded at the C-terminus is kept as written.
    assert site(out, out.relations[0])['description'] == 'HPFLRRNSGC-----'


def test_other_kinase_names_are_not_gene_symbols():
    rows = parse()
    out = extract([rows[0], rows[7]])
    subjects = [out.entities[r.subject_entity_key] for r in out.relations]
    assert [(e.entity_type, e.namespace, e.identifier, e.taxon) for e in subjects] == [
        ('protein', Namespace.NAME, 'PDKC', None),
        ('protein', Namespace.NAME, 'BCR/ABL', None),
    ]


def test_source_residue_is_kept_as_reported():
    out = extract([parse()[9]])
    (relation,) = out.relations
    assert (site(out, relation)['residue'], site(out, relation)['position']) == ('H', 118)


def test_no_publication_from_a_count_and_no_invented_values():
    rows = parse()
    out = extract(rows)
    assert len(out.relations) == 10
    for relation in out.relations:
        # npmid is a number of publications, not a PubMed ID.
        assert annotations(relation, 'publications') == []
        assert not [a for a in relation.annotations if a['value'] in {'8', '91'}]
        # Phosphorylation is not an effect on activity.
        assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']
        assert annotations(relation, 'object_direction_qualifier') == ['increased']
    # No site without a kinase reaches the output.
    assert 'C6orf57' not in {e.identifier for e in out.entities.values()}
