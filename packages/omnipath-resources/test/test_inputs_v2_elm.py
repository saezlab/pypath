"""Source contract of the ELM motif-domain interaction input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.elm import resource
from pypath.inputs_v2.parsers.elm import (
    join_evidence,
    parse_position,
    read_features,
    split_methods,
)
from pypath.inputs_v2.parsers.base import iter_tsv


# Header and rows of http://elm.eu.org/interactions/as_tsv ...
TABLE_HEADER = (
    'Elm\tDomain\tinteractorElm\tinteractorDomain\tStartElm\tStopElm\t'
    'StartDomain\tStopDomain\tAffinityMin\tAffinityMax\tPMID\ttaxonomyElm\t'
    'taxonomyDomain\t'
)
# ... and of the same records in http://elm.eu.org/interactions/as_mitab.
# The method and interaction type cells of the last record are shortened.
EVIDENCE_HEADER = (
    'Unique_identifier_A\tUnique_identifier_B\tAlternative_identifier_A\t'
    'Alternative_identifier_B\tAliases for A\tAliases for B\tMethods\t'
    'First_author\tpublication\tTaxon_A\tTaxon_B\tInteraction_types\t'
    'Source_databases\tInteraction_identifier\tConfidence score'
)
RECORDS = {
    'separase': (
        'CLV_Separin_Fungi\tPF03568\tQ12158\tQ03018\t175\t181\t1171\t1571\t'
        'None\tNone\t10403247,14585836\t"559292"(Saccharomyces cerevisiae S288c)\t'
        '"559292"(Saccharomyces cerevisiae S288c)',
        'uniprotkb:Q12158\tuniprotkb:Q03018\t\t\tuniprotkb:MCD1(gene name)\t'
        'uniprotkb:ESP1(gene name)\tpsi-mi:"MI:0363"(inferred by author);'
        'psi-mi:"MI:0019"(coimmunoprecipitation)\t'
        'Sullivan et al. (2004)|Uhlmann et al. (1999)\t'
        'pubmed:14585836|pubmed:10403247\t'
        'taxid:559292(Saccharomyces cerevisiae S288c)\t'
        'taxid:559292(Saccharomyces cerevisiae S288c)\t'
        'psi-mi:"MI:0407"(direct interaction)|psi-mi:"MI:0914"(association)\t'
        'ELM\tELMINT000001\t',
    ),
    # A measured affinity, and partners of two species.
    'affinity': (
        'LIG_PDZ_Class_1\tPF00595\tP35222\tQ9DBG9\t776\t781\t17\t110\t0.14\t'
        '0.24\t18835279\t"9606"(Homo sapiens)\t"10090"(Mus musculus)',
        'uniprotkb:P35222\tuniprotkb:Q9DBG9\t\t\tuniprotkb:CTNNB1(gene name)\t'
        'uniprotkb:Tax1bp3(gene name)\t'
        'psi-mi:"MI:0065"(isothermal titration calorimetry);'
        'psi-mi:"MI:0114"(x-ray crystallography)\tZhang et al. (2008)\t'
        'pubmed:18835279\ttaxid:9606(Homo sapiens)|taxid:9606(human)\t'
        'taxid:10090(Mus musculus)|taxid:10090(house mouse)\t'
        'psi-mi:"MI:0914"(association)\tELM\tELMINT000015\t',
    ),
    # The motif is in an isoform; the position of the domain is not given.
    'isoform': (
        'LIG_PDZ_Class_1\tPF00595\tP46937-2\tQ9UDY2\t483\t488\tNone\tNone\t'
        'None\tNone\t20868367\t"9606"(Homo sapiens)\t"9606"(Homo sapiens)',
        'uniprotkb:P46937-2\tuniprotkb:Q9UDY2\t\t\tuniprotkb:YAP1(gene name)\t'
        'uniprotkb:TJP2(gene name)\tpsi-mi:"MI:0363"(inferred by author);'
        'psi-mi:"MI:0019"(coimmunoprecipitation);'
        'psi-mi:"MI:0065"(isothermal titration calorimetry)\t'
        'Oka et al. (2010)\tpubmed:20868367\t'
        'taxid:9606(Homo sapiens)|taxid:9606(human)\t'
        'taxid:9606(Homo sapiens)|taxid:9606(human)\t'
        'psi-mi:"MI:0914"(association)\tELM\tELMINT000102\t',
    ),
    # The domain position is written 0, and the domain is an InterPro entry.
    'zero': (
        'LIG_WD40_WDR5_VDV_1\tIPR017986\tQ9H9L4\tP61964\t410\t413\t0\t0\t8.6\t'
        '8.6\t24788516\t"9606"(Homo sapiens)\t"9606"(Homo sapiens)',
        'uniprotkb:Q9H9L4\tuniprotkb:P61964\t\t\tuniprotkb:KANSL2(gene name)\t'
        'uniprotkb:WDR5(gene name)\tpsi-mi:"MI:0363"(inferred by author);'
        'psi-mi:"MI:0065"(isothermal titration calorimetry);'
        'psi-mi:"MI:0114"(x-ray crystallography)\tDias et al. (2014)\t'
        'pubmed:24788516\ttaxid:9606(Homo sapiens)|taxid:9606(human)\t'
        'taxid:9606(Homo sapiens)|taxid:9606(human)\t'
        'psi-mi:"MI:0914"(association)\tELM\tELMINT000980\t',
    ),
    # No publication, no method.
    'unreferenced': (
        'LIG_eIF4E_2\tPF01652\tQ9BY44\tP06730\t446\t452\t38\t200\tNone\tNone\t'
        '\t"9606"(Homo sapiens)\t"9606"(Homo sapiens)',
        'uniprotkb:Q9BY44\tuniprotkb:P06730\t\t\tuniprotkb:EIF2A(gene name)\t'
        'uniprotkb:EIF4E(gene name)\t\t\t\t'
        'taxid:9606(Homo sapiens)|taxid:9606(human)\t'
        'taxid:9606(Homo sapiens)|taxid:9606(human)\t'
        'psi-mi:"MI:0914"(association)\tELM\tELMINT000573\t',
    ),
    # An empty gene name.
    'unnamed': (
        'LIG_Actin_WH2_2\tPF00022\tP18281\tP68135\t2\t18\t4\t377\tNone\tNone\t'
        '1464307\t"5755"(Acanthamoeba castellanii)\t'
        '"9986"(Oryctolagus cuniculus)',
        'uniprotkb:P18281\tuniprotkb:P68135\t\t\tuniprotkb:(gene name)\t'
        'uniprotkb:ACTA1(gene name)\tpsi-mi:"MI:0363"(inferred by author);'
        'psi-mi:"MI:0030"(cross linking study)\t'
        'Vancompernolle et al. (1992)\tpubmed:1464307\t'
        'taxid:5755(Acanthamoeba castellanii)\t'
        'taxid:9986(Oryctolagus cuniculus)\t'
        'psi-mi:"MI:0407"(direct interaction)\tELM\tELMINT000156\t',
    ),
    # A class identifier with hyphens.
    'hyphens': (
        'LIG_14-3-3_CanoR_1\tPF00244\tB7UM88\tP63104\t22\t31\t9\t229\tNone\t'
        'None\t15339273\t"574521"(Escherichia coli O127:H6 str. E2348/69)\t'
        '"9606"(Homo sapiens)',
        'uniprotkb:B7UM88\tuniprotkb:P63104\t\t\tuniprotkb:espF(gene name)\t'
        'uniprotkb:YWHAZ(gene name)\tpsi-mi:"MI:0019"(coimmunoprecipitation)\t'
        'Viswanathan et al. (2004)\tpubmed:15339273\t'
        'taxid:574521(Escherichia coli O127:H6 str. E2348/69)\t'
        'taxid:9606(Homo sapiens)\tpsi-mi:"MI:0914"(association)\tELM\t'
        'ELMINT001682\t',
    ),
    'smart': (
        'DEG_APCC_TPR_1\tSM00028\tP51955\tP30260\t443\t445\tNone\tNone\tNone\t'
        'None\t23288039,16648845\t"9606"(Homo sapiens)\t"9606"(Homo sapiens)',
        'uniprotkb:P51955\tuniprotkb:P30260\t\t\tuniprotkb:NEK2(gene name)\t'
        'uniprotkb:CDC27(gene name)\t'
        'psi-mi:"MI:0997"(ubiquitinase assay);'
        'psi-mi:"MI:0019"(coimmunoprecipitation)\t'
        'Sedgwick et al. (2013)|Hayes et al. (2006)\t'
        'pubmed:23288039|pubmed:16648845\t'
        'taxid:9606(Homo sapiens)|taxid:9606(human)\t'
        'taxid:9606(Homo sapiens)|taxid:9606(human)\t'
        'psi-mi:"MI:0414"(enzymatic reaction)|psi-mi:"MI:0403"(colocalization)'
        '|psi-mi:"MI:0914"(association)\tELM\tELMINT000718\t',
    ),
}

# Two records of the same proteins, species and publication that differ only in
# the domain: the PSI-MI TAB export gives them the same cells but their numbers.
RECORDS['calmodulin_45'] = (
    'LIG_CaM_IQ_9\tPF13499\tP38377\tP62161\t13\t31\t85\t146\tNone\tNone\t21102557\t'
    '"9615"(Canis lupus familiaris)\t"10116"(Rattus norvegicus)',
    'uniprotkb:P38377\tuniprotkb:P62161\t\t\tuniprotkb:SEC61A1(gene name)\t'
    'uniprotkb:Calm1(gene name)\tpsi-mi:"MI:0363"(inferred by author);'
    'psi-mi:"MI:0051"(fluorescence technology)\tErdmann et al. (2011)\t'
    'pubmed:21102557\ttaxid:9615(Canis lupus familiaris)|taxid:9615(dog)\t'
    'taxid:10116(Rattus norvegicus)|taxid:10116(Norway rat)\t'
    'psi-mi:"MI:0407"(direct interaction)\tELM\tELMINT000045\t',
)
RECORDS['calmodulin_46'] = (
    RECORDS['calmodulin_45'][0].replace('\t85\t146\t', '\t8\t73\t'),
    RECORDS['calmodulin_45'][1].replace('ELMINT000045', 'ELMINT000046'),
)

PROTEIN = {'identifier': None, 'coordinate_system': 'protein', 'position_base': 1}


def opener(header, lines):
    return SimpleNamespace(result=io.StringIO('\n'.join([header, *lines]) + '\n'))


def feature_xml(table_lines, evidence_lines):
    """The PSI-MI XML export of the same records, in the layout of
    http://elm.eu.org/interactions/as_psimi: a list of interactors, and for each
    record its number and the motif and the domain as features of the participants."""
    interactors, interactions = [], []

    def feature(name, start, stop):
        positions = (
            f'<featureRangeList><featureRange><startStatus/><begin position="{start}"/>'
            f'<endStatus/><end position="{stop}"/></featureRange></featureRangeList>'
            if start not in ('None', '0')
            else '<featureRangeList><featureRange><startStatus/><endStatus/>'
            '</featureRange></featureRangeList>'
        )
        return f'<featureList><feature><names><shortLabel>{name}</shortLabel></names>{positions}</feature></featureList>'

    for number, (line, evidence_line) in enumerate(zip(table_lines, evidence_lines), 1):
        cells = line.split('\t')
        elm_class, domain, motif_protein, domain_protein = cells[:4]
        record = evidence_line.split('\t')[13]
        interactors += [
            f'<interactor id="m{number}"><xref><primaryRef id="{motif_protein}"/></xref></interactor>',
            f'<interactor id="d{number}"><xref><primaryRef id="{domain_protein}"/></xref></interactor>',
        ]
        interactions.append(
            f'<interaction id="{number}"><xref><primaryRef db="elm" id="{record}"/></xref>'
            f'<participantList><participant><interactorRef>m{number}</interactorRef>'
            f'{feature(elm_class, cells[4], cells[5])}</participant>'
            f'<participant><interactorRef>d{number}</interactorRef>'
            f'{feature(domain, cells[6], cells[7])}</participant></participantList></interaction>'
        )
    return SimpleNamespace(
        result=io.BytesIO(
            (
                '<entrySet xmlns="net:sf:psidev:mi" level="2"><entry><interactorList>'
                + ''.join(interactors)
                + '</interactorList><interactionList>'
                + ''.join(interactions)
                + '</interactionList></entry></entrySet>'
            ).encode()
        )
    )


def joined(*names, table=None, evidence=None, features=None):
    """Rows of the two exports, read and joined as the dataset does."""
    table = [RECORDS[name][0] for name in names] if table is None else table
    evidence = [RECORDS[name][1] for name in names] if evidence is None else evidence
    if features is None:
        features = read_features(feature_xml(table, evidence))
    return list(
        join_evidence(
            iter_tsv(opener(TABLE_HEADER, table)),
            iter_tsv(opener(EVIDENCE_HEADER, evidence)),
            features,
        )
    )


def extract(source_row):
    record = resource.datasets()['interactions'].mapper(source_row)
    assert record is not None
    out = SilverExtractor('elm', 'interactions')
    out.process_record(record, source_row, 'fixture:1', 1)
    (relation,) = out.relations
    return (
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
    )


def record(name):
    (source_row,) = joined(name)
    return extract(source_row)


def identifiers(entity):
    return {(i['ns'], i['id']) for i in entity.identifiers}


def annotations(relation, scope='relation'):
    return sorted(
        (a['term'], a['value'])
        for a in relation.annotations
        if a['scope'] == scope
    )


def region(entity):
    (item,) = entity.molecular_form['regions']
    return (
        item['type'],
        item['identifier'],
        item['position'],
        item['end_position'],
        item['coordinate_reference'],
    )


def test_interaction_is_motif_protein_to_domain_protein():
    relation, motif, domain = record('separase')
    assert relation.predicate == 'directly_physically_interacts_with'
    assert relation.upstream_id == 'ELMINT000001'
    assert identifiers(motif) == {
        (Namespace.UNIPROT, 'Q12158'),
        (Namespace.GENESYMBOL, 'MCD1'),
    }
    assert identifiers(domain) == {
        (Namespace.UNIPROT, 'Q03018'),
        (Namespace.GENESYMBOL, 'ESP1'),
    }
    assert (motif.namespace, motif.identifier) == ('uniprot', 'Q12158')
    assert {motif.entity_type, domain.entity_type} == {'protein'}
    assert annotations(relation) == [
        ('has_evidence_of_type', 'MI:0019'),
        ('has_evidence_of_type', 'MI:0363'),
        ('publications', 'PMID:10403247'),
        ('publications', 'PMID:14585836'),
    ]
    assert annotations(relation, 'subject') == [('in_taxon', 'NCBITaxon:559292')]
    assert annotations(relation, 'object') == [('in_taxon', 'NCBITaxon:559292')]


def test_motif_instance_and_domain_are_binding_regions_of_the_participants():
    _, motif, domain = record('separase')
    assert region(motif) == (
        'direct binding region',
        {'ns': 'elm', 'id': 'CLV_Separin_Fungi'},
        175,
        181,
        PROTEIN,
    )
    assert region(domain) == (
        'direct binding region',
        {'ns': 'pfam', 'id': 'PF03568'},
        1171,
        1571,
        PROTEIN,
    )


def test_each_participant_has_its_own_taxon():
    relation, motif, domain = record('affinity')
    assert (motif.taxon, domain.taxon) == ('9606', '10090')
    assert annotations(relation, 'object') == [('in_taxon', 'NCBITaxon:10090')]


def test_affinity_and_interaction_types_are_not_asserted():
    # The table does not say which constant the affinity is, and the
    # interaction types of the PSI-MI TAB export describe the experiments.
    (source_row,) = joined('affinity')
    assert (source_row['AffinityMin'], source_row['AffinityMax']) == ('0.14', '0.24')
    relation, _, _ = extract(source_row)
    assert {a['term'] for a in relation.annotations} == {
        'has_evidence_of_type',
        'in_taxon',
        'publications',
    }
    assert all(a['quantity'] is None for a in relation.annotations)
    assert all('MI:0914' != a['value'] for a in relation.annotations)


def test_isoform_is_the_reference_of_the_motif_position():
    _, motif, domain = record('isoform')
    isoform = {'ns': 'uniprot', 'id': 'P46937-2'}
    assert (Namespace.UNIPROT, 'P46937-2') in identifiers(motif)
    assert motif.molecular_form['isoform_identifier'] == isoform
    assert region(motif) == (
        'direct binding region',
        {'ns': 'elm', 'id': 'LIG_PDZ_Class_1'},
        483,
        488,
        {**PROTEIN, 'identifier': isoform},
    )
    # ``None`` for the domain position: the domain is named, not located.
    assert region(domain) == (
        'direct binding region',
        {'ns': 'pfam', 'id': 'PF00595'},
        None,
        None,
        None,
    )


def test_zero_is_not_a_domain_position():
    _, _, domain = record('zero')
    assert region(domain) == (
        'direct binding region',
        {'ns': 'interpro', 'id': 'IPR017986'},
        None,
        None,
        None,
    )


def test_domain_namespace_follows_the_accession_and_is_not_guessed():
    _, _, domain = record('smart')
    assert region(domain)[1] == {'ns': 'smart', 'id': 'SM00028'}
    (source_row,) = joined('separase')
    _, _, domain = extract(dict(source_row, Domain='SH3_domain'))
    assert region(domain) == ('direct binding region', None, 1171, 1571, PROTEIN)
    _, _, domain = extract(
        dict(source_row, Domain='SH3_domain', StartDomain='None', StopDomain='None')
    )
    assert domain.molecular_form is None


def test_class_identifier_with_hyphens_names_the_motif():
    _, motif, _ = record('hyphens')
    assert region(motif)[1] == {'ns': 'elm', 'id': 'LIG_14-3-3_CanoR_1'}


def test_record_without_publication_or_method_has_neither():
    relation, _, _ = record('unreferenced')
    assert annotations(relation) == []
    assert relation.upstream_id == 'ELMINT000573'


def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication():
    (source_row,) = joined('separase')
    row = dict(source_row, PMID='0123456,123456789,10403247')
    relation, _, _ = extract(row)
    assert [v for t, v in annotations(relation) if t == 'publications'] == [
        'PMID:10403247'
    ]


def test_empty_gene_name_is_not_an_identifier():
    _, motif, domain = record('unnamed')
    assert identifiers(motif) == {(Namespace.UNIPROT, 'P18281')}
    assert (Namespace.GENESYMBOL, 'ACTA1') in identifiers(domain)


def test_method_written_without_the_mi_prefix_is_read():
    (source_row,) = joined('hyphens')
    methods = (
        'psi-mi:"MI:0363"(inferred by author);psi-mi:"1260"(pyrene label );'
        'psi-mi:"MI:1088"(phenotype-based detection assay)'
    )
    relation, _, _ = extract(dict(source_row, Methods=methods))
    assert [v for t, v in annotations(relation) if t == 'has_evidence_of_type'] == [
        'MI:0363',
        'MI:1088',
        'MI:1260',
    ]


@pytest.mark.parametrize('column', ['interactorElm', 'interactorDomain'])
@pytest.mark.parametrize('value', ['', 'YAP1', 'P46937,Q9UDY2'])
def test_no_relation_without_an_accession_of_each_protein(column, value):
    (source_row,) = joined('separase')
    mapper = resource.datasets()['interactions'].mapper
    assert mapper(dict(source_row, **{column: value})) is None


def test_join_pairs_the_two_exports_by_position():
    rows = joined('separase', 'affinity')
    assert [row['Interaction_identifier'] for row in rows] == [
        'ELMINT000001',
        'ELMINT000015',
    ]
    # The table reader strips the quotes around the taxon number.
    assert rows[0]['taxonomyElm'] == '559292(Saccharomyces cerevisiae S288c)'
    assert rows[1]['Methods'].startswith('psi-mi:"MI:0065"')
    assert '' not in rows[0] and None not in rows[0]


def test_join_rejects_exports_in_a_different_order():
    with pytest.raises(ValueError, match='Row 1 .* not the same record'):
        joined(
            table=[RECORDS['separase'][0], RECORDS['affinity'][0]],
            evidence=[RECORDS['affinity'][1], RECORDS['separase'][1]],
        )


def test_join_rejects_a_record_whose_publications_differ():
    evidence = RECORDS['separase'][1].replace('pubmed:14585836|', '')
    with pytest.raises(ValueError, match='not the same record'):
        joined(table=[RECORDS['separase'][0]], evidence=[evidence])


def test_join_rejects_exports_of_different_length():
    with pytest.raises(ValueError, match='differ in length'):
        joined(
            table=[RECORDS['separase'][0], RECORDS['affinity'][0]],
            evidence=[RECORDS['separase'][1]],
        )


def test_join_rejects_changed_columns_and_empty_exports():
    table = opener(TABLE_HEADER.replace('StartElm', 'Start'), [RECORDS['separase'][0]])
    evidence = opener(EVIDENCE_HEADER, [RECORDS['separase'][1]])
    with pytest.raises(ValueError, match='StartElm'):
        list(join_evidence(iter_tsv(table), iter_tsv(evidence), {}))
    with pytest.raises(ValueError, match='no rows'):
        joined(features={})


@pytest.mark.parametrize(
    ('value', 'position'),
    [('175', 175), ('0', None), ('None', None), ('', None), (None, None)],
)
def test_parse_position_reads_the_null_conventions_of_the_table(value, position):
    assert parse_position(value) == position


def test_split_methods_splits_at_semicolons_not_pipes():
    cell = 'psi-mi:"MI:0363"(inferred by author);psi-mi:"MI:0019"(coimmunoprecipitation)'
    assert split_methods(cell) == [
        'psi-mi:"MI:0363"(inferred by author)',
        'psi-mi:"MI:0019"(coimmunoprecipitation)',
    ]
    assert split_methods(None) == ['']


def test_features_are_read_by_record_number_with_null_positions():
    features = read_features(
        feature_xml(
            [RECORDS['separase'][0], RECORDS['zero'][0]],
            [RECORDS['separase'][1], RECORDS['zero'][1]],
        )
    )
    assert features['ELMINT000001'] == {
        ('Q12158', 'CLV_Separin_Fungi', 175, 181): 1,
        ('Q03018', 'PF03568', 1171, 1571): 1,
    }
    # A range the export leaves out, or writes as 0, is no position.
    assert features['ELMINT000980'][('P61964', 'IPR017986', None, None)] == 1
    with pytest.raises(ValueError, match='not PSI-MI XML'):
        read_features(SimpleNamespace(result=io.BytesIO(b'<html/>')))
    with pytest.raises(ValueError, match='not PSI-MI XML'):
        read_features(None)
    twice = feature_xml([RECORDS['separase'][0]] * 2, [RECORDS['separase'][1]] * 2)
    with pytest.raises(ValueError, match='twice'):
        read_features(twice)


def test_join_rejects_records_of_another_species():
    # Real record ELMINT000015 with the second species of the PSI-MI TAB export
    # changed: the two exports no longer describe the same record.
    evidence = RECORDS['affinity'][1].replace('taxid:10090', 'taxid:9606')
    with pytest.raises(ValueError, match='not the same record'):
        joined(table=[RECORDS['affinity'][0]], evidence=[evidence])


def test_records_that_differ_only_in_the_domain_are_paired_by_their_features():
    # ELMINT000045 and ELMINT000046: the same proteins, species and publication;
    # the domain lies at 85-146 in the first and at 8-73 in the second.
    first, second = joined('calmodulin_45', 'calmodulin_46')
    assert (first['Interaction_identifier'], first['StartDomain']) == ('ELMINT000045', '85')
    assert (second['Interaction_identifier'], second['StartDomain']) == ('ELMINT000046', '8')


def test_join_rejects_a_pairing_that_swaps_two_records_of_one_group():
    # The table rows in the other order than the PSI-MI TAB rows: the key
    # (proteins, species, publications) is the same for both, only the domain
    # positions tell that the pairs are wrong.
    with pytest.raises(ValueError, match='ELMINT000045.*PSI-MI XML'):
        joined(
            table=[RECORDS['calmodulin_46'][0], RECORDS['calmodulin_45'][0]],
            evidence=[RECORDS['calmodulin_45'][1], RECORDS['calmodulin_46'][1]],
            features=read_features(
                feature_xml(
                    [RECORDS['calmodulin_45'][0], RECORDS['calmodulin_46'][0]],
                    [RECORDS['calmodulin_45'][1], RECORDS['calmodulin_46'][1]],
                )
            ),
        )


def test_join_rejects_a_record_the_xml_export_does_not_have():
    with pytest.raises(ValueError, match='ELMINT000015'):
        joined(
            'separase',
            'affinity',
            features=read_features(
                feature_xml([RECORDS['separase'][0]], [RECORDS['separase'][1]])
            ),
        )
