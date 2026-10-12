"""Source contract of the SPIKE input."""

import functools
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.parsers.spike import (
    QUALITY_LEVELS,
    iter_interactions,
    iter_spike_groups,
    iter_spike_interactions,
    iter_spike_regulations,
)
from pypath.inputs_v2.parsers.spike import iter_regulations as _iter_regulations
from pypath.inputs_v2.spike import (
    ASSAYS,
    EFFECTS,
    MECHANISMS,
    SPIKE_NAMESPACE,
    map_group,
    map_interaction,
    map_regulation,
    regulation_statements,
    resource,
)

iter_regulations = functools.partial(
    _iter_regulations, statements_of=regulation_statements
)

_COMMON = 'owner="public" status="Draft" isDeleted="0"'


def gene(internal_id, accession, name, entrez, long_name):
    return (
        f'<Gene spikeAccession="{accession}" name="{name}" id="{internal_id}" '
        f'dataSource="LocusLink" site="TAU" owner="public" integrity="1" '
        f'status="Approved" type="Gene" isDeleted="0" longName="{long_name}">'
        f'<XRef db="1" id="{entrez}"/><Display color="Purple" shape="Rectangle" '
        f'symbol="{name}"/><Description>Automatically obtained from Entrez Gene '
        f'on 23-05-11</Description></Gene>'
    )


def references(*pmids):
    return ''.join(f'<Reference pmid="{pmid}" firstAuthor="" journal=""/>' for pmid in pmids)


# Elements of LatestSpikeDB.xml (25 January 2012), without the aliases of the
# genes, the display settings of the groups and relations and the citation
# details of the references; long reference lists are shortened.
DOCUMENT = f"""<?xml version="1.0" encoding="UTF-8"?>
<SpikeDatabase version="2" creation-date="2012-01-25T14:49:15.248" author="spike_curator" origin="CEN">
<SiteBlock><Site name="TAU" description="Tel-Aviv University, Israel"/></SiteBlock>
<BuildingBlock>
{gene('487', 'TAU:P-000456', 'ATM', '472', 'ataxia telangiectasia mutated')}
{gene('882', 'TAU:P-000851', 'CCNB1', '891', 'cyclin B1')}
{gene('972', 'TAU:P-000941', 'CDK1', '983', 'cyclin-dependent kinase 1')}
{gene('1013', 'TAU:P-000982', 'CDKN1A', '1026', 'cyclin-dependent kinase inhibitor 1A (p21, Cip1)')}
{gene('3902', 'TAU:P-003871', 'LYN', '4067', 'v-yes-1 Yamaguchi sarcoma viral related oncogene homolog')}
{gene('3998', 'TAU:P-003967', 'MCM2', '4171', 'minichromosome maintenance complex component 2')}
{gene('4019', 'TAU:P-003988', 'MDM2', '4193', 'Mdm2 p53 binding protein homolog (mouse)')}
{gene('4419', 'TAU:P-004388', 'MYC', '4609', 'v-myc myelocytomatosis viral oncogene homolog (avian)')}
{gene('4892', 'TAU:P-004861', 'PCNA', '5111', 'proliferating cell nuclear antigen')}
{gene('5196', 'TAU:P-005165', 'POLB', '5423', 'polymerase (DNA directed), beta')}
{gene('5289', 'TAU:P-005258', 'PPP2R1A', '5518', 'protein phosphatase 2, regulatory subunit A, alpha')}
{gene('5368', 'TAU:P-005337', 'MAPK8', '5599', 'mitogen-activated protein kinase 8')}
{gene('5370', 'TAU:P-005339', 'MAPK9', '5601', 'mitogen-activated protein kinase 9')}
{gene('5371', 'TAU:P-005340', 'MAPK10', '5602', 'mitogen-activated protein kinase 10')}
{gene('5641', 'TAU:P-005610', 'RAC1', '5879', 'ras-related C3 botulinum toxin substrate 1 (rho family, small GTP binding protein Rac1)')}
{gene('5681', 'TAU:P-005650', 'RASA1', '5921', 'RAS p21 protein activator (GTPase activating protein) 1')}
{gene('5688', 'TAU:P-005657', 'RBBP4', '5928', 'retinoblastoma binding protein 4')}
{gene('6878', 'TAU:P-006847', 'TP53', '7157', 'tumor protein p53')}
{gene('7625', 'TAU:P-007594', 'CHAF1B', '8208', 'chromatin assembly factor 1, subunit B (p60)')}
{gene('7756', 'TAU:P-007725', 'BCAR3', '8412', 'breast cancer anti-estrogen resistance 3')}
{gene('8222', 'TAU:P-008191', 'HAP1', '9001', 'huntingtin-associated protein 1')}
{gene('8315', 'TAU:P-008284', 'CCNB2', '9133', 'cyclin B2')}
{gene('8824', 'TAU:P-008793', 'CHAF1A', '10036', 'chromatin assembly factor 1, subunit A (p150)')}
{gene('9600', 'TAU:P-009569', 'CHEK2', '11200', 'CHK2 checkpoint homolog (S. pombe)')}
{gene('100189', 'TAU:P-015181', 'ASF1A', '25842', 'ASF1 anti-silencing function 1 homolog A (S. cerevisiae)')}
{gene('1506225', 'TAU:P-023908', 'MIR24-1', '407012', 'microRNA 24-1')}
{gene('1506226', 'TAU:P-023909', 'MIR24-2', '407013', 'microRNA 24-2')}
{gene('1506254', 'TAU:P-023937', 'MIR34C', '407042', 'microRNA 34c')}
<SignallingMolecule spikeAccession="CEN:S-000033" name="S1P" id="1509499" dataSource="ChEBI" site="CEN" owner="public" integrity="1" status="Draft" type="Signaling Molecule" isDeleted="0" longName="Sphingosine 1-phosphate"><XRef db="4" id="37550"/></SignallingMolecule>
<Event spikeAccession="TAU:E-000011" name="Cell-cycle progression" id="1524067" dataSource="Published research" site="TAU" owner="public" integrity="1" status="Draft" type="Event" isDeleted="0"></Event>
<Group spikeAccession="TAU:C-000003" name="CycB-Cdc2" id="15079" dataSource="Published research" site="TAU" owner="public" integrity="1" status="Draft" type="Complex" isDeleted="0" longName="CyclinB-Cdc2 (CDK1)- M-phase-promoting factor (MPF)"><Description/><Member ref="972"/><Member ref="1509699"/></Group>
<Group spikeAccession="CEN:F-329953" name="CycB" id="1509699" dataSource="Published research" site="CEN" owner="public" integrity="1" status="Draft" type="Family" isDeleted="0" longName="Cyclin B"><Member ref="882"/><Member ref="8315"/></Group>
<Group spikeAccession="TAU:F-000008" name="JNK" id="101111" dataSource="Published research" site="TAU" owner="public" integrity="1" status="Draft" type="Family" isDeleted="0" longName="JNK1-3"><Description>JNK1=Mapk8, JNK2=Mapk9, JNK3=Mapk10</Description><Member ref="5368"/><Member ref="5370"/><Member ref="5371"/></Group>
<Group spikeAccession="TAU:F-329917" name="mir-24" id="1502630" dataSource="Published research" site="TAU" owner="public" integrity="1" status="Draft" type="Family" isDeleted="0"><Member ref="1506226"/><Member ref="1506225"/></Group>
<Group spikeAccession="ELK:C-000029" name="CAF-1" id="326506" dataSource="Published research" site="ELK" owner="public" integrity="1" status="Draft" type="Complex" isDeleted="0" longName="chromatin assembly factor 1"><Description/><Member ref="8824"/><Member ref="7625"/><Member ref="5688"/></Group>
<Group spikeAccession="TAU:C-000017" name="PIC" id="324770" dataSource="Published research" site="TAU" owner="public" integrity="1" status="Draft" type="Complex" isDeleted="0" longName="Transcription PreInitiationComplex"><Description>Initiation of transcription by RNA polymerase II requires the activities of more than 70 polypeptides. </Description></Group>
<Group spikeAccession="TAU:C-000039" name="CUL4-DDB1-2" id="1502120" dataSource="Published research" site="TAU" owner="public" integrity="1" status="Draft" type="Complex" isDeleted="1" longName="CUL4-DDB1-2 E3 ubiqitin ligase"></Group>
</BuildingBlock>
<RegulationBlock>
<Regulation spikeAccession="TAU:R-000001" id="15030" dataSource="Published research" site="TAU" integrity="1" type="Regulation" mechanism="Phosphorylation" effect="1" {_COMMON}>{references('9843217', '9733515')}<Description>Ser15</Description><Target ref="6878"/><PhysicalTarget ref="6878"/><Source ref="487"/></Regulation>
<Regulation spikeAccession="TAU:R-000253" id="324789" dataSource="Published research" site="TAU" integrity="2" type="Regulation" mechanism="Phosphorylation" effect="3" {_COMMON}>{references('17525332', '18309293', '15448142')}<Description>The 2nd and 3rd references are bazsed on Xenopus laevis</Description><Target ref="3998"/><PhysicalTarget ref="3998"/><Source ref="487"/></Regulation>
<Regulation spikeAccession="TAU:R-000009" id="15042" dataSource="Published research" site="TAU" integrity="1" type="Regulation" mechanism="Phosphorylation" effect="2" {_COMMON}>{references('10710310')}<Description>Phosphorylation on p53-Ser20</Description><Target ref="15041"/><PhysicalTarget ref="6878"/><Source ref="9600"/></Regulation>
<Regulation spikeAccession="TAU:R-000008" id="15041" dataSource="Published research" site="TAU" integrity="1" type="Regulation" mechanism="Physical Interaction" effect="2" {_COMMON}>{references('19579266', '20639885', '11331603')}<Target ref="6878"/><PhysicalTarget ref="6878"/><Source ref="4019"/></Regulation>
<Regulation spikeAccession="ELK:R-330365" id="330374" dataSource="Published research" site="ELK" integrity="1" type="Regulation" mechanism="Physical Interaction" effect="2" {_COMMON}>{references('15077146')}<Target ref="15079"/><PhysicalTarget ref="972"/><Source ref="1013"/></Regulation>
<Regulation spikeAccession="TAU:R-000004" id="15033" dataSource="Published research" site="TAU" integrity="1" type="Regulation" mechanism="Transcription Regulation" effect="1" {_COMMON}>{references('12748190', '19749794', '21163364')}<Target ref="1013"/><PhysicalTarget ref="1013"/><Source ref="6878"/></Regulation>
<Regulation spikeAccession="TAU:R-331071" id="1524807" dataSource="Shiloh's lab" site="TAU" owner="public" integrity="3" status="Draft" type="Regulation" isDeleted="1" mechanism="Phosphorylation" effect="3"><Description>Shelly</Description><Target ref="6878"/><PhysicalTarget ref="6878"/><Source ref="487"/></Regulation>
<Regulation spikeAccession="WEB:0000711" id="1526464" dataSource="Published research" site="WEB" integrity="1" type="Regulation" mechanism="mRNA stability" effect="2" {_COMMON}>{references('20212154')}<Target ref="4419"/><PhysicalTarget ref="4419"/><Source ref="1506254"/></Regulation>
<Regulation spikeAccession="WEB:0003765" id="1532576" dataSource="Published research" site="WEB" integrity="2" type="Regulation" mechanism="Other" effect="1" {_COMMON}>{references('21304987')}<Target ref="5641"/><PhysicalTarget ref="5641"/><Source ref="1509499"/></Regulation>
<Regulation spikeAccession="TAU:R-332275" id="1524263" dataSource="Published research" site="TAU" integrity="1" type="Regulation" mechanism="N/A" effect="2" {_COMMON}>{references('18590585')}<Target ref="1524067"/><PhysicalTarget ref="1524067"/><Source ref="1013"/></Regulation>
</RegulationBlock>
<InteractionBlock>
<Interaction spikeAccession="ELK:I-010444" id="1009080" dataSource="IntAct DB" site="ELK" integrity="4" type="Interaction" biologicalAssay="Y2H" {_COMMON}>{references('17145710')}<Description>intact:EBI-1074875,intact:EBI-1074875</Description><ProteinB ref="5289"/><ProteinA ref="7756"/></Interaction>
<Interaction spikeAccession="CEN:I-011523" id="1510365" dataSource="Published research" site="CEN" integrity="2" type="Interaction" biologicalAssay="Other" {_COMMON}>{references('8395016')}<ProteinB ref="5681"/><ProteinA ref="3902"/></Interaction>
<Interaction spikeAccession="CEN:I-011237" id="1507048" dataSource="Published research" site="CEN" integrity="1" type="Interaction" biologicalAssay="IP" {_COMMON}>{references('18378699', '11897662')}<ProteinB ref="326506"/><ProteinA ref="100189"/></Interaction>
<Interaction spikeAccession="CEN:I-012967" id="1511837" dataSource="Published research" site="CEN" owner="public" integrity="1" status="Draft" type="Interaction" isDeleted="1" biologicalAssay="Y2H">{references('9207062')}<ProteinB ref="5196"/><ProteinA ref="8222"/></Interaction>
</InteractionBlock>
</SpikeDatabase>
"""


def document(text=DOCUMENT):
    return SimpleNamespace(result={'LatestSpikeDB.xml': text.encode('utf-8')})


def regulations(accession):
    """The parsed rows (one per statement) of a regulation, at any level."""
    return [
        row
        for level in QUALITY_LEVELS
        for row in iter_regulations(document(), level=level)
        if row['accession'] == accession
    ]


def interaction(accession):
    (row,) = [
        row
        for level in QUALITY_LEVELS
        for row in iter_interactions(document(), level=level)
        if row['accession'] == accession
    ]
    return row


GROUPS = {row['name']: row for row in iter_spike_groups(document())}


def extract(record, row, dataset='regulations_quality_1'):
    assert record is not None
    out = SilverExtractor('spike', dataset)
    out.process_record(record, row, 'fixture:1', 1)
    return out


def relation_of(row, mapper=map_regulation):
    out = extract(mapper(row), row)
    (relation,) = out.relations
    return (
        relation,
        out.entities[relation.subject_entity_key],
        out.entities[relation.object_entity_key],
    )


def identifiers(entity):
    return {(i['ns'], i['id']) for i in entity.identifiers}


def annotations(relation, scope='relation'):
    return sorted(
        (a['term'], a['value']) for a in relation.annotations if a['scope'] == scope
    )


def qualifiers(relation):
    return {
        term: value for term, value in annotations(relation) if 'qualifier' in term
    }


def test_reader_resolves_the_references_of_a_regulation():
    (row,) = [
        r for r in iter_spike_regulations(document()) if r['accession'] == 'ELK:R-330365'
    ]
    assert (row['mechanism'], row['effect'], row['integrity']) == (
        'Physical Interaction', '2', '1',
    )
    assert row['pubmed'] == ['15077146']
    assert (row['source']['kind'], row['source']['name']) == ('Gene', 'CDKN1A')
    assert (row['source']['xref_db'], row['source']['xref_id']) == ('1', '1026')
    assert (row['target']['kind'], row['target']['group_type']) == ('Group', 'Complex')
    assert row['target']['accession'] == 'TAU:C-000003'
    assert row['physical_target']['name'] == 'CDK1'
    assert row['is_deleted'] is False


def test_reader_marks_a_regulation_that_targets_a_regulation():
    (row,) = [
        r for r in iter_spike_regulations(document()) if r['accession'] == 'TAU:R-000009'
    ]
    assert row['target']['kind'] == 'Regulation'
    assert row['target']['accession'] == 'TAU:R-000008'
    assert row['physical_target']['name'] == 'TP53'


def test_reader_reads_groups_and_interactions():
    complex_row = GROUPS['CycB-Cdc2']
    assert complex_row['group_type'] == 'Complex'
    assert [(m['kind'], m['name']) for m in complex_row['members']] == [
        ('Gene', 'CDK1'),
        ('Group', 'CycB'),
    ]
    assert GROUPS['JNK']['description'] == 'JNK1=Mapk8, JNK2=Mapk9, JNK3=Mapk10'
    assert GROUPS['CUL4-DDB1-2']['is_deleted'] is True
    rows = list(iter_spike_interactions(document()))
    assert [r['biological_assay'] for r in rows] == ['Y2H', 'Other', 'IP', 'Y2H']
    assert (rows[0]['protein_a']['name'], rows[0]['protein_b']['name']) == (
        'BCAR3', 'PPP2R1A',
    )


def test_reader_rejects_another_document():
    with pytest.raises(ValueError, match='SPIKE'):
        list(iter_spike_regulations(document('<sbml><model/></sbml>')))
    with pytest.raises(ValueError, match='InteractionBlock'):
        list(
            iter_spike_groups(
                document('<SpikeDatabase><BuildingBlock/><RegulationBlock/></SpikeDatabase>')
            )
        )


def test_activating_phosphorylation_gives_the_effect_and_the_modification():
    effect, modification = regulations('TAU:R-000001')
    assert (effect['STATEMENT'], modification['STATEMENT']) == ('effect', 'modification')
    relation, source, target = relation_of(effect)
    assert relation.predicate == 'affects'
    assert relation.upstream_id == 'TAU:R-000001'
    assert identifiers(source) == {
        (Namespace.ENTREZ, '472'),
        (Namespace.GENESYMBOL, 'ATM'),
        (Namespace.NAME, 'ataxia telangiectasia mutated'),
    }
    assert (Namespace.ENTREZ, '7157') in identifiers(target)
    assert {source.entity_type, target.entity_type} == {'protein'}
    assert {source.taxon, target.taxon} == {'9606'}
    # The effect is on TP53 as a whole: no aspect.
    assert annotations(relation) == [
        ('description', 'Ser15'),
        ('has_confidence_level', '1'),
        ('object_direction_qualifier', 'increased'),
        ('publications', 'PMID:9733515'),
        ('publications', 'PMID:9843217'),
    ]
    relation, _, target = relation_of(modification)
    assert qualifiers(relation) == {
        'object_direction_qualifier': 'increased',
        'object_aspect_qualifier': 'phosphorylation',
    }
    assert (Namespace.GENESYMBOL, 'TP53') in identifiers(target)


def test_phosphorylation_with_unknown_effect_gives_the_modification_only():
    (row,) = regulations('TAU:R-000253')
    assert row['STATEMENT'] == 'modification'
    relation, _, target = relation_of(row, map_regulation)
    assert (Namespace.GENESYMBOL, 'MCM2') in identifiers(target)
    assert qualifiers(relation) == {
        'object_direction_qualifier': 'increased',
        'object_aspect_qualifier': 'phosphorylation',
    }
    assert ('has_confidence_level', '2') in annotations(relation)


def test_inhibition_of_a_regulation_is_not_an_inhibition_of_the_physical_target():
    # CHEK2 inhibits the inhibition of TP53 by MDM2, by phosphorylating TP53.
    (row,) = regulations('TAU:R-000009')
    assert row['STATEMENT'] == 'modification'
    relation, source, target = relation_of(row)
    assert (Namespace.GENESYMBOL, 'CHEK2') in identifiers(source)
    assert (Namespace.GENESYMBOL, 'TP53') in identifiers(target)
    assert qualifiers(relation) == {
        'object_direction_qualifier': 'increased',
        'object_aspect_qualifier': 'phosphorylation',
    }
    # With a mechanism that is no modification, only "acts on" remains.
    other = dict(row, mechanism='Physical Interaction')
    assert regulation_statements(other) == ['action']
    relation, _, _ = relation_of(dict(other, STATEMENT='action'))
    assert qualifiers(relation) == {}


def test_inhibition_by_binding_has_a_direction_and_no_aspect():
    (row,) = regulations('TAU:R-000008')
    assert row['STATEMENT'] == 'effect'
    relation, _, _ = relation_of(row)
    assert qualifiers(relation) == {'object_direction_qualifier': 'decreased'}


def test_complex_target_stays_the_complex():
    (row,) = regulations('ELK:R-330365')
    out = extract(map_regulation(row), row)
    (relation,) = out.relations
    target = out.entities[relation.object_entity_key]
    assert target.entity_type == 'macromolecular_complex'
    assert identifiers(target) == {
        (SPIKE_NAMESPACE, 'TAU:C-000003'),
        (Namespace.NAME, 'CycB-Cdc2'),
        (Namespace.SYNONYM, 'CyclinB-Cdc2 (CDK1)- M-phase-promoting factor (MPF)'),
    }
    # Neither the physical target CDK1 nor any other member gets a relation.
    assert len(out.entities) == 2


def test_transcription_regulation_is_an_effect_on_expression():
    (row,) = regulations('TAU:R-000004')
    relation, _, target = relation_of(row)
    assert (Namespace.GENESYMBOL, 'CDKN1A') in identifiers(target)
    assert qualifiers(relation) == {
        'object_direction_qualifier': 'increased',
        'object_aspect_qualifier': 'expression',
    }


def test_microrna_small_molecule_and_event_participants():
    (row,) = regulations('WEB:0000711')
    relation, source, _ = relation_of(row)
    assert source.entity_type == 'microrna'
    assert (Namespace.GENESYMBOL, 'MIR34C') in identifiers(source)
    assert qualifiers(relation) == {
        'object_direction_qualifier': 'decreased',
        'object_aspect_qualifier': 'stability',
    }
    (row,) = regulations('WEB:0003765')
    relation, source, _ = relation_of(row)
    assert source.entity_type == 'chemical_entity'
    assert identifiers(source) == {
        (Namespace.CHEBI, '37550'),
        (Namespace.NAME, 'S1P'),
        (Namespace.SYNONYM, 'Sphingosine 1-phosphate'),
    }
    assert source.taxon is None
    (row,) = regulations('TAU:R-332275')
    relation, _, target = relation_of(row)
    assert target.entity_type == 'named_thing'
    assert identifiers(target) == {
        (SPIKE_NAMESPACE, 'TAU:E-000011'),
        (Namespace.NAME, 'Cell-cycle progression'),
    }


def test_deleted_regulation_gives_no_relation():
    (row,) = regulations('TAU:R-331071')
    assert row['STATEMENT'] is None
    assert map_regulation(row) is None
    assert map_regulation(dict(row, STATEMENT='modification')) is None


def test_datasets_are_the_quality_levels():
    names = set(resource.datasets())
    assert names == {
        'complexes',
        'protein_families',
        *(f'regulations_quality_{level}' for level in '1234'),
        *(f'interactions_quality_{level}' for level in '1234'),
    }
    level_2 = list(iter_regulations(document(), level='2'))
    assert {row['accession'] for row in level_2} == {'TAU:R-000253', 'WEB:0003765'}
    assert [row['accession'] for row in iter_interactions(document(), level='4')] == [
        'ELK:I-010444'
    ]
    with pytest.raises(ValueError, match='quality level'):
        list(iter_interactions(document(), level='5'))
    with pytest.raises(ValueError, match='quality level'):
        list(iter_interactions(document(DOCUMENT.replace('integrity="4"', 'integrity="7"')), level='1'))


def test_interaction_keeps_method_level_and_imported_database():
    row = interaction('ELK:I-010444')
    relation, first, second = relation_of(row, map_interaction)
    assert relation.predicate == 'interacts_with'
    assert (Namespace.GENESYMBOL, 'BCAR3') in identifiers(first)
    assert (Namespace.GENESYMBOL, 'PPP2R1A') in identifiers(second)
    assert annotations(relation) == [
        ('description', 'intact:EBI-1074875,intact:EBI-1074875'),
        ('has_confidence_level', '4'),
        ('has_evidence_of_type', 'MI:0018'),
        ('publications', 'PMID:17145710'),
        ('supporting_data_source', 'IntAct DB'),
    ]


def test_interaction_with_unspecific_assay_has_no_method():
    relation, _, _ = relation_of(interaction('CEN:I-011523'), map_interaction)
    assert annotations(relation) == [
        ('has_confidence_level', '2'),
        ('publications', 'PMID:8395016'),
    ]


def test_interaction_with_a_complex_is_with_the_complex():
    row = interaction('CEN:I-011237')
    out = extract(map_interaction(row), row)
    (relation,) = out.relations
    assert out.entities[relation.object_entity_key].entity_type == 'macromolecular_complex'
    assert (SPIKE_NAMESPACE, 'ELK:C-000029') in identifiers(
        out.entities[relation.object_entity_key]
    )
    assert ('has_evidence_of_type', 'MI:0019') in annotations(relation)
    assert len(out.entities) == 2


def test_deleted_interaction_gives_no_relation():
    assert map_interaction(interaction('CEN:I-012967')) is None


def test_complex_has_its_stated_members_and_a_family_stays_a_family():
    row = GROUPS['CycB-Cdc2']
    out = extract(map_group(row), row, 'complexes')
    assert [r.predicate for r in out.relations] == ['has_member', 'has_member']
    by_type = {
        out.entities[r.object_entity_key].entity_type: identifiers(
            out.entities[r.object_entity_key]
        )
        for r in out.relations
    }
    assert (Namespace.ENTREZ, '983') in by_type['protein']
    assert by_type['protein_family'] == {
        (SPIKE_NAMESPACE, 'CEN:F-329953'),
        (Namespace.NAME, 'CycB'),
        (Namespace.SYNONYM, 'Cyclin B'),
    }
    # The cyclins of the family are not members of the complex.
    assert len(out.entities) == 3


def test_family_has_its_members_and_its_description():
    row = GROUPS['JNK']
    out = extract(map_group(row), row, 'protein_families')
    (family,) = [e for e in out.entities.values() if e.entity_type == 'protein_family']
    assert identifiers(family) == {
        (SPIKE_NAMESPACE, 'TAU:F-000008'),
        (Namespace.NAME, 'JNK'),
        (Namespace.SYNONYM, 'JNK1-3'),
    }
    assert [(a['term'], a['value']) for a in family.annotations] == [
        ('description', 'JNK1=Mapk8, JNK2=Mapk9, JNK3=Mapk10')
    ]
    members = {
        next(v for ns, v in identifiers(out.entities[r.object_entity_key]) if ns == Namespace.GENESYMBOL)
        for r in out.relations
    }
    assert members == {'MAPK8', 'MAPK9', 'MAPK10'}


def test_family_of_microrna_genes_is_not_a_protein_family():
    row = GROUPS['mir-24']
    out = extract(map_group(row), row, 'protein_families')
    types = sorted(e.entity_type for e in out.entities.values())
    assert types == ['gene_family', 'microrna', 'microrna']


def test_group_without_members_or_deleted_gives_no_record():
    assert map_group(GROUPS['PIC']) is None
    assert map_group(GROUPS['CUL4-DDB1-2']) is None


def test_group_datasets_split_complexes_and_families():
    datasets = resource.datasets()
    assert {r['name'] for r in datasets['complexes']._raw_parser(document())} == {
        'CycB-Cdc2', 'CAF-1', 'PIC', 'CUL4-DDB1-2',
    }
    assert {r['name'] for r in datasets['protein_families']._raw_parser(document())} == {
        'CycB', 'JNK', 'mir-24',
    }


def test_unknown_codes_are_errors():
    (row,) = regulations('TAU:R-000008')
    with pytest.raises(ValueError, match='mechanism'):
        regulation_statements(dict(row, mechanism='Citrullination'))
    with pytest.raises(ValueError, match='effect'):
        regulation_statements(dict(row, effect='4'))
    with pytest.raises(ValueError, match='data source'):
        map_regulation(dict(row, data_source='STRING'))
    with pytest.raises(ValueError, match='assay'):
        map_interaction(dict(interaction('CEN:I-011523'), biological_assay='FRET'))
    assert set(EFFECTS) == {'1', '2', '3'} and EFFECTS['3'] is None
    assert ASSAYS['Other'] is None
    assert 'Neddylation' in MECHANISMS
