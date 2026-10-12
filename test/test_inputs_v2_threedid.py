"""Source contract of the 3did domain-domain and domain-motif inputs."""

from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import entity_type
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import threedid
from pypath.inputs_v2.parsers import threedid as parser

# Records of 3did_flat: a homodimer pair, a heterodimer pair, and an instance
# without the topology field. The lines below an instance are residue contacts.
DDI = (
    '#=ID\t1-cysPrx_C\t1-cysPrx_C\t (PF10417.14@Pfam\tPF10417.14@Pfam)\n'
    '#=3D\t1n8j\tE:153-185\tO:153-185\t0.99\t1.35582\t2:2\n'
    'A\tA\t158 \t159 \tmm\n'
    'A\tH\t158 \t160 \tmm\n'
    '#=3D\t5jcg\tA:157-192\tH:157-192\t1.1\t1.53738\t2:2\n'
    'Q\tE\t159 \t162 \tsm\n'
    '//\n'
    '#=ID\t1-cysPrx_C\tAhpC-TSA\t (PF10417.14@Pfam\tPF00578.26@Pfam)\n'
    '#=3D\t1e2y\tB:162-176\tB:8-142\t4.07\t3.17739\t12:2\n'
    'T\tD\t163 \t10 \tss\n'
    '//\n'
    '#=ID\tProteasome\tProteasome\t (PF00227.31@Pfam\tPF00227.31@Pfam)\n'
    '#=3D\t1e94\tC:1-171\tD:1-171\t3.81\t2.92685\n'
    'V\tQ\t26 \t130 \tms\n'
    '//\n'
)
# Records of 3did_dmi_flat.
DMI = (
    '#=ID\tPCNA_C\t(Pfam)\tPCNA_C_LIG_0-0\t(PLoS_CB_2010)\n'
    '#=PT\tQ..[ILM]..F[FY] (Feb 2025)\n'
    '#=3D\t7m5l\tB:127-254\tE:2-11\tQCSMTCFY\t0\t0\n'
    '#=3D\t7m5l\tC:127-254\tF:2-10\tQCSMTCFY\t3\t0\n'
    '//\n'
    '#=ID\tPCNA_C\t(Pfam)\tPCNA_C_LIG_0-1\t(PLoS_CB_2010)\n'
    '#=PT\tQ..[ILM]..F[FY] (Feb 2025)\n'
    '#=3D\t5mlo\tA:127-254\tB:519-526\tQHDIRSFF\t2\t1\n'
    '//\n'
)


def opened(text):
    return SimpleNamespace(result=text)


def rows(reader, text):
    return list(reader(opened(text)))


def extract(dataset, parsed):
    mapper = threedid.resource.datasets()[dataset].mapper
    out = SilverExtractor('threedid', dataset)
    for index, row in enumerate(parsed):
        record = mapper(row)
        assert record is not None
        out.process_record(record, row, f'fixture:{index}', index)
    return out


def identifiers(entity):
    return {(i['ns'], i['id']) for i in entity.identifiers}


def test_domain_domain_instances_are_read_with_their_pairs_and_regions():
    parsed = rows(parser.iter_domain_domain, DDI)
    assert len(parsed) == 4
    assert parsed[0] == {
        'name_1': '1-cysPrx_C',
        'name_2': '1-cysPrx_C',
        'pfam_1': 'PF10417',
        'pfam_version_1': '14',
        'pfam_2': 'PF10417',
        'pfam_version_2': '14',
        'pdb_id': '1n8j',
        'region_1': {'chain': 'E', 'start': '153', 'end': '185'},
        'region_2': {'chain': 'O', 'start': '153', 'end': '185'},
        'score': 0.99,
        'z_score': 1.35582,
        'topology': '2:2',
    }
    assert (parsed[2]['name_2'], parsed[2]['pfam_2']) == ('AhpC-TSA', 'PF00578')
    # The topology field is missing on some lines of the file.
    assert parsed[3]['topology'] is None
    assert parsed[3]['pdb_id'] == '1e94'


def test_domain_motif_instances_carry_the_motif_and_its_pattern():
    parsed = rows(parser.iter_domain_motif, DMI)
    assert len(parsed) == 3
    assert parsed[0] == {
        'domain': 'PCNA_C',
        'motif': 'PCNA_C_LIG_0-0',
        'motif_source': 'PLoS_CB_2010',
        'pattern': 'Q..[ILM]..F[FY]',
        'pattern_date': 'Feb 2025',
        'pdb_id': '7m5l',
        'domain_region': {'chain': 'B', 'start': '127', 'end': '254'},
        'motif_region': {'chain': 'E', 'start': '2', 'end': '11'},
        'peptide': 'QCSMTCFY',
        'contacts': 0,
        'topology': 0,
    }
    assert parsed[2]['motif'] == 'PCNA_C_LIG_0-1'
    assert parsed[2]['topology'] == 1


@pytest.mark.parametrize(
    'reader,text,message',
    [
        (parser.iter_domain_domain, '', 'empty'),
        (parser.iter_domain_domain, '#=ID\tA\tB\t (PF1.1@Pfam\tPF00017.1@Pfam)\n', 'malformed domain pair'),
        (parser.iter_domain_domain, '#=3D\t1n8j\tE:1-2\tO:1-2\t0.99\t1.3\t2:2\n', 'without a domain pair'),
        (
            parser.iter_domain_domain,
            DDI.replace('E:153-185', 'E153-185'),
            'malformed chain range',
        ),
        (parser.iter_domain_domain, DDI.replace('\t0.99\t', '\tx\t'), 'not a number'),
        (parser.iter_domain_domain, DDI.replace('1n8j', '1N8J!'), 'malformed structural instance'),
        (parser.iter_domain_domain, '#=ID\tA\tB\t (PF00017.1@Pfam\tPF00017.1@Pfam)\n//\n', 'no structural instances'),
        (parser.iter_domain_motif, DMI.replace('(Pfam)', '(InterPro)', 1), 'malformed domain-motif'),
        (parser.iter_domain_motif, DMI.replace('\t0\t0\n', '\t0\n', 1), 'malformed structural instance'),
        (parser.iter_domain_motif, '#=3D\t7m5l\tB:1-2\tE:1-2\tQ\t0\t0\n', 'without a motif'),
    ],
)
def test_readers_reject_files_that_are_not_as_expected(reader, text, message):
    with pytest.raises(ValueError, match=message):
        rows(reader, text)


def test_domain_domain_is_a_contact_of_two_pfam_domains_not_of_proteins():
    parsed = rows(parser.iter_domain_domain, DDI)
    record = threedid.resource.domain_domain.mapper(parsed[2])
    assert entity_type(record.subject.type) == 'protein_domain'
    assert record.subject.identifiers == [
        ('pfam', 'PF10417'),
        ('pfam_id', '1-cysPrx_C'),
        (Namespace.NAME, '1-cysPrx_C'),
    ]
    assert record.object.identifiers[0] == ('pfam', 'PF00578')

    out = extract('domain_domain', parsed)
    assert {r.predicate for r in out.relations} == {'directly_physically_interacts_with'}
    assert {e.entity_type for e in out.entities.values()} == {'protein_domain'}
    # No protein, no gene and no species is given to a domain family.
    assert {e.taxon for e in out.entities.values()} == {None}
    assert {e.namespace for e in out.entities.values()} == {'pfam'}
    # Two instances of the homodimer pair, one of the heterodimer pair and one
    # of the proteasome pair give three statements and four observations.
    assert len({r.relation_key for r in out.relations}) == 3
    assert [r.upstream_id for r in out.relations] == ['1n8j', '5jcg', '1e2y', '1e94']


def test_instance_attributes_keep_the_chains_ranges_and_scores():
    parsed = rows(parser.iter_domain_domain, DDI)
    out = extract('domain_domain', parsed[2:3])
    (relation,) = out.relations
    assert relation.upstream_id == '1e2y'
    assert sorted((a['term'], a['value']) for a in relation.annotations) == [
        ('threedid:domain_1_region', 'B:162-176'),
        ('threedid:domain_2_region', 'B:8-142'),
        ('threedid:interprets_score', '4.07'),
        ('threedid:interprets_z_score', '3.17739'),
        ('threedid:topology', '12:2'),
    ]
    # An instance without a topology field has no topology attribute.
    (no_topology,) = extract('domain_domain', parsed[3:4]).relations
    assert 'threedid:topology' not in {a['term'] for a in no_topology.annotations}


def test_domain_motif_is_a_domain_and_a_motif_class_with_its_pattern():
    parsed = rows(parser.iter_domain_motif, DMI)
    out = extract('domain_motif', parsed)
    assert {r.predicate for r in out.relations} == {'directly_physically_interacts_with'}
    domains = [e for e in out.entities.values() if e.entity_type == 'protein_domain']
    motifs = [e for e in out.entities.values() if e.entity_type == 'ontology_class']
    # The file names the domain, not its accession; both motifs bind the same domain.
    assert len(domains) == 1
    assert identifiers(domains[0]) == {('pfam_id', 'PCNA_C'), (Namespace.NAME, 'PCNA_C')}
    # The motif is a concept of 3did's own vocabulary.
    assert {e.namespace for e in motifs} == {'threedid_motif'}
    assert {e.identifier for e in motifs} == {'PCNA_C_LIG_0-0', 'PCNA_C_LIG_0-1'}
    assert len({r.relation_key for r in out.relations}) == 2
    first = out.relations[0]
    assert first.upstream_id == '7m5l'
    assert {
        (a['term'], a['value']) for a in first.annotations if a['scope'] == 'relation'
    } == {
        ('threedid:domain_region', 'B:127-254'),
        ('threedid:motif_region', 'E:2-11'),
        ('threedid:motif_instance', 'QCSMTCFY'),
        ('threedid:contextual_contacts', '0'),
        ('threedid:topology', '0'),
    }
    assert {
        (a['term'], a['value']) for a in first.annotations if a['scope'] == 'object'
    } == {('threedid:motif_pattern', 'Q..[ILM]..F[FY]')}


def test_no_protein_is_asserted():
    # The file has chains and PDB entries, not proteins: nothing is made of a
    # chain, and no species, UniProt accession or gene is invented.
    out = extract('domain_domain', rows(parser.iter_domain_domain, DDI))
    out_motif = extract('domain_motif', rows(parser.iter_domain_motif, DMI))
    for result in (out, out_motif):
        assert not [e for e in result.entities.values() if e.entity_type == 'protein']
        assert all(
            i['ns'] not in (Namespace.UNIPROT, Namespace.GENESYMBOL)
            for e in result.entities.values()
            for i in e.identifiers
        )
