"""Pathway Commons extended SIF: parser, mapper and extractor on source lines."""

import gzip

import pytest
from cachedir._open import Opener
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import pathwaycommons as pc
from pypath.inputs_v2.parsers import pathwaycommons as parser

EDGE_HEADER = (
    'PARTICIPANT_A\tINTERACTION_TYPE\tPARTICIPANT_B\tINTERACTION_DATA_SOURCE\t'
    'INTERACTION_PUBMED_ID\tPATHWAY_NAMES\tMEDIATOR_IDS'
)
NODE_HEADER = (
    'PARTICIPANT\tPARTICIPANT_TYPE\tPARTICIPANT_NAME\tUNIFICATION_XREF\t'
    'RELATIONSHIP_XREF'
)
# Lines copied from pc-hgnc.txt.gz (Pathway Commons release 14): one edge of
# each relation type, edges with RNA participants and the participants'
# lines of the second section. ``EDGES`` are in the order of the first
# column of the table in ``test_each_relation_type_is_mapped_from_its_documentation``.
EDGES = [
    'A1BG\tcontrols-expression-of\tA2M\tpid\t12456685;7678052;9794795\tIL6-mediated signaling events\tpid:pid_74043;pid:pid_74038',
    'A1BG\tinteracts-with\tABCC6\tBioGRID\t21988832\t\tbiogrid:MolecularInteraction_bebab183-9762-4322-b91f-e9a69a0c9e97___null__162374_',
    'A1CF\tin-complex-with\tAPOBEC1\tIntAct;Reactome\t11134005;11352648;11727199;12881431;14681455;21496894;8999813\tFormation of the Editosome;mRNA Editing: C to U Conversion\treactome:Complex4742;reactome:Complex4741;intact.complex:Complex_d0de9be0-bce1-44a0-917e-a1ba2804f93d___null__19983_',
    'A1CF\tcontrols-production-of\tchebi:17775\tCTD\t23263486\t\tctdbase:ABU_4127922;ctdbase:ABU_MESH_D014527',
    'A4GALT\tcatalysis-precedes\tB3GNT3\tRecon X\t10747952;10748143\t\treconx:R_A4GALTg;reconx:control_-2029732522;reconx:R_B3GNT31g;reconx:control_-614896134',
    'AAAS\tcontrols-transport-of\tAGO2\tReactome\t15643423;18842624;19167051;22064859;23150874;23882114;24388755;24596094;9214382\tTranscriptional regulation by small RNAs\thttp://bioregistry.io/reactome:R-HSA-5578744;reactome:Catalysis4930',
    'AAAS\tcontrols-transport-of-chemical\tchebi:1604\tReactome\t10531003\tsnRNP Assembly\thttp://bioregistry.io/reactome:R-HSA-191830;http://bioregistry.io/reactome:R-HSA-191812',
    'CHEBI:10357\tchemical-affects\tACACA\tCTD\t27234712\t\tctdbase:PHO_7129021;ctdbase:PHO_7129034;ctdbase:PHO_GENE_31;ctdbase:PHO_7129030;ctdbase:PHO_7129009;ctdbase:PHO_7129016;ctdbase:PHO_7129039',
    'CHEBI:11684\tconsumption-controlled-by\tDIO1\tCTD\t21565810\t\tctdbase:ABU_MESH_C030629;ctdbase:MET_MESH_D014285;ctdbase:ABU_3478634;ctdbase:ABU_3478639;ctdbase:ACT_GENE_1733;ctdbase:MET_3478624',
    'CHEBI:136982\tused-to-produce\tchebi:16189\tRecon X\t11886493\t\treconx:R_RE3218L;reconx:R_RE3218R;reconx:R_RE3218C',
    'CHEBI:143180\treacts-with\tchebi:15824\tRecon X\t385176;8528071;890964\t\treconx:R_RE0926C',
    'AHR\tcontrols-expression-of\tCCDC177\tMSigDB\t\t\tmsigdb:control_7ff6df05-a25f-419d-a30d-18b181553071;msigdb:transcription_CCDC177_aa0d49ef-a919-45a3-b09e-3d9de9eec613',
    'MIR1228\tcontrols-state-change-of\tCASP3\tCTD\t26867688\t\tctdbase:CLV_GENE_836;ctdbase:RXN_6679439',
    'AGO2\tinteracts-with\tCASC2\tBioGRID\t33174006\t\tbiogrid:MolecularInteraction_fa83a73f-e48f-4c39-8885-79af4280a7f7___null__786531_',
    'CYP2A13\tcontrols-production-of\tCHEBI:PubChem:25245599\tPANTHER\t\tNicotine_degradation\tpanther:_CATALYSIS___CYP2A6_re7m1_re7;panther:STATE_TRANSITION_LEFT__nicotine-_NADPH-_H2O-RIGHT__NADP_PLUS_-_O2-_2_hydroxynicotine-_re7',
    'ABCF3\tinteracts-with\tMOCS2\tIntAct;HPRD\t16169070\t\tintact:MolecularInteraction_4e444cb4-67f1-4dc6-9e86-876958aa66a3__IM-16517-2581_230039_;hprd:MolecularInteraction_5d05a54e-49d8-4ecb-8b69-2dd5b503f54e___null__29318_',
]
NODES = [
    'A1BG\tDnaReference;ProteinReference;RnaReference\tA1BG;A1BG gene;A1BG mRNA;A1BG_HUMAN;FHOS78;NRON;ssDNA\trefseq:NM_130786;refseq:NP_570602;uniprot:P04217\thgnc.symbol:A1BG;ncbigene:1;omim:138670;uniprot:P04217',
    'A2M\tDnaReference;ProteinReference;RnaReference;RnaRegionReference\tA2M;A2M exon;A2M mRNA;A2MG_HUMAN;APRE\tensembl:Z11711;uniprot:P01023\thgnc.symbol:A2M;ncbigene:2',
    'ABCC6\tDnaReference;ProteinReference;RnaReference;RnaRegionReference\tABCC6;ABCC6 intron;ABCC6 mRNA;MRP6_HUMAN\trefseq:NM_001171;uniprot:O95255\thgnc.symbol:ABCC6;ncbigene:368;omim:603234',
    'A1CF\tDnaReference;DnaRegionReference;ProteinReference;RnaReference\tA1CF;A1CF gene SNP;A1CF mRNA;A1CF promoter;A1CF_HUMAN;ACF;Apobec 1 complementation factor\thprd:10627;refseq:NM_014576;refseq:XP_005269776;uniprot:Q9NQ94\thgnc.symbol:A1CF;ncbigene:29974;uniprot:Q9NQ94',
    'APOBEC1\tDnaRegionReference;ProteinReference;RnaReference\tABEC1_HUMAN;APOBEC1;APOBEC1 mRNA;APOBEC1 promoter\tuniprot:P41238\thgnc.symbol:APOBEC1;ncbigene:339',
    'chebi:17775\tSmallMoleculeReference\t7,9-dihydro-1H-purine-2,6,8(3H)-trione\tchebi:CHEBI:17775\t',
    'A4GALT\tDnaRegionReference;ProteinReference;RnaReference\tA4GALT;A4GALT mRNA;A4GALT promoter;A4GAT_HUMAN\trefseq:NP_059132;uniprot:Q9NPC4\thgnc.symbol:A4GALT;ncbigene:53947;uniprot:Q9NPC4',
    'B3GNT3\tProteinReference;RnaReference;RnaRegionReference\tB3GN3_HUMAN;B3GNT3 intron;B3GNT3 mRNA\tuniprot:Q9Y2A9\thgnc.symbol:B3GNT3;ncbigene:10331',
    'AAAS\tProteinReference;RnaReference\tAAAS;AAAS gene;AAAS mRNA;AAAS_HUMAN\tuniprot:Q9NRG9\thgnc.symbol:AAAS;ncbigene:8086',
    'AGO2\tProteinReference;RnaReference\tAGO2;AGO2 gene;AGO2 mRNA;AGO2_HUMAN\trefseq:XP_016868806;uniprot:Q9UKV8\thgnc.symbol:AGO2;ncbigene:27161;uniprot:Q9UKV8',
    'chebi:1604\tSmallMoleculeReference\t3-methylguanine\tchebi:CHEBI:1604\t',
    'CHEBI:10357\tSmallMoleculeReference\tcaryophyllene\tmesh:C024714\t',
    'ACACA\tProteinReference;RnaReference\tACACA;ACACA gene;ACACA mRNA;ACACA_HUMAN;Acetyl-CoA carboxylase alpha\thprd:01938;refseq:XP_006725383;uniprot:Q13085\thgnc.symbol:ACACA;ncbigene:31;uniprot:Q13085',
    'CHEBI:11684\tSmallMoleculeReference\tTriiodothyronine, Reverse\tmesh:D014285\t',
    'DIO1\tProteinReference;RnaReference\tDIO1;DIO1 gene SNP;DIO1 mRNA;IOD1_HUMAN\tuniprot:P49895\thgnc.symbol:DIO1;ncbigene:1733',
    'CHEBI:136982\tSmallMoleculeReference\t5alpha-Dihydrotestosterone sulfate\thmdb:HMDB06278\t',
    'chebi:16189\tSmallMoleculeReference\tsulfate\tchebi:CHEBI:16189\t',
    'CHEBI:143180\tSmallMoleculeReference\tMaltotetraose\thmdb:HMDB01296;kegg.compound:C02052;pubchem.substance:5141\t',
    'chebi:15824\tSmallMoleculeReference\tD-fructose\tchebi:CHEBI:15824\t',
    'AHR\tProteinReference;RnaReference\tAHR;AHR gene mutant form;AHR mRNA;AHR mutant form;AHR_HUMAN\trefseq:NP_001612;uniprot:P35869\thgnc.symbol:AHR;ncbigene:196;uniprot:P35869',
    'CCDC177\tRnaReference\tCCDC177;CCDC177 gene;CCDC177 mRNA\t\t',
    'MIR1228\tDnaRegionReference;RnaReference\tMIR1228 mRNA;MIR1228 promoter\t\t',
    'CASP3\tDnaReference;ProteinReference;RnaReference\tCASP3;CASP3 gene;CASP3_HUMAN;YAMA;casp3_human_gene\tensembl:ENSG00000164305;refseq:XP_011530603;uniprot:P42574\thgnc.symbol:CASP3;intact:EBI-16686578;ncbigene:836;uniprot:P42574',
    'CASC2\tRnaReference\tCASC2;CASC2 mRNA\t\t',
    'CYP2A13\tDnaRegionReference;ProteinReference;RnaReference\tCP2AD_HUMAN;CYP2A13 gene polymorphism;CYP2A13 mRNA;CYP2A13 promoter\tuniprot:Q16696\thgnc.symbol:CYP2A13;ncbigene:1553',
    'CHEBI:PubChem:25245599\tSmallMoleculeReference\t2-hydroxynicotine\tchebi:CHEBI:PubChem:25245599\t',
    'ABCF3\tDnaReference;ProteinReference;RnaReference\tABCF3;ABCF3 gene;ABCF3 mRNA;ABCF3_HUMAN;FLJ11198\trefseq:NM_018358;refseq:XP_011511261;uniprot:Q9NUQ8\thgnc.symbol:ABCF3;ncbigene:55324;uniprot:Q9NUQ8',
    'MOCS2\tDnaReference;ProteinReference;RnaReference\tGLUD1P3;MOC2A_HUMAN;MOC2B_HUMAN;MOCS2;MOCS2 gene;MOCS2 mRNA;MOCS2 promoter;Molybdenum cofactor synthesis 2\thprd:04750;ncbigene:4338;uniprot:O96007;uniprot:O96033\thgnc.symbol:MOCS2;ncbigene:4338;omim:603708;refseq:NP_004522;refseq:NP_789776;uniprot:O96007;uniprot:O96033',
]


def write(tmp_path, edges=EDGES, nodes=NODES, edge_header=EDGE_HEADER):
    path = tmp_path / 'pc-hgnc.txt.gz'
    text = '\n'.join([edge_header, *edges, '', NODE_HEADER, *nodes]) + '\n'
    path.write_bytes(gzip.compress(text.encode()))
    return Opener(path=str(path), ext='gz', large=True, default_mode='r', encoding='utf-8')


def read(tmp_path, **kwargs):
    return list(parser.iter_interactions(write(tmp_path, **kwargs)))


def edge(source, kind, target):
    for line in EDGES:
        fields = line.split('\t')
        if fields[:3] == [source, kind, target]:
            return line
    raise KeyError((source, kind, target))


def map_edge(tmp_path, source, kind, target):
    (row,) = read(tmp_path, edges=[edge(source, kind, target)])
    return row, pc.resource.interactions.mapper(row)


def extract(tmp_path, rows):
    out = SilverExtractor('pathwaycommons', 'interactions')
    for index, row in enumerate(rows):
        record = pc.resource.interactions.mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def values(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term]


def identifiers(entity):
    return [(i.type, i.value) for i in entity.identifiers]


def test_module_import_defines_the_one_dataset_of_the_release_14_file():
    assert list(pc.resource.datasets()) == ['interactions']
    assert pc.resource.interactions.download.url == (
        'https://download.baderlab.org/PathwayCommons/PC2/v14/pc-hgnc.txt.gz'
    )


def test_parser_joins_the_second_section_to_the_edges(tmp_path):
    rows = read(tmp_path)
    assert len(rows) == len(EDGES)
    first = rows[0]
    assert first['interaction_type'] == 'controls-expression-of'
    assert first['data_sources'] == ['pid']
    assert first['pubmed'] == ['12456685', '7678052', '9794795']
    assert first['pathway_names'] == 'IL6-mediated signaling events'
    assert first['participant_a'] == {
        'id': 'A1BG',
        'types': ['DnaReference', 'ProteinReference', 'RnaReference'],
        'names': ['A1BG', 'A1BG gene', 'A1BG mRNA', 'A1BG_HUMAN', 'FHOS78', 'NRON', 'ssDNA'],
        'unification_xrefs': [
            ('refseq', 'NM_130786'), ('refseq', 'NP_570602'), ('uniprot', 'P04217')
        ],
        'relationship_xrefs': [
            ('hgnc.symbol', 'A1BG'), ('ncbigene', '1'), ('omim', '138670'),
            ('uniprot', 'P04217'),
        ],
    }
    # A cross-reference accession can have colons itself.
    chebi = {r['participant_b']['id']: r for r in rows}['chebi:1604']
    assert chebi['participant_b']['unification_xrefs'] == [('chebi', 'CHEBI:1604')]


def test_parser_raises_on_a_changed_header_a_missing_participant_and_a_missing_section(
    tmp_path,
):
    with pytest.raises(ValueError, match='lacks expected columns'):
        read(tmp_path, edge_header=EDGE_HEADER.replace('MEDIATOR_IDS', 'MEDIATORS'))
    with pytest.raises(ValueError, match="'A2M' is not described"):
        read(tmp_path, nodes=[n for n in NODES if not n.startswith('A2M\t')])
    path = tmp_path / 'edges-only.txt.gz'
    path.write_bytes(gzip.compress((EDGE_HEADER + '\n' + EDGES[0] + '\n').encode()))
    opener = Opener(path=str(path), ext='gz', large=True, default_mode='r', encoding='utf-8')
    with pytest.raises(ValueError, match='no participant section'):
        list(parser.iter_interactions(opener))


def test_each_relation_type_is_mapped_from_its_documentation(tmp_path):
    expected = {
        ('A1BG', 'interacts-with', 'ABCC6'): ('interacts_with', None),
        ('A1CF', 'in-complex-with', 'APOBEC1'): ('in_complex_with', None),
        ('MIR1228', 'controls-state-change-of', 'CASP3'): ('affects', None),
        ('AAAS', 'controls-transport-of', 'AGO2'): ('affects', 'localization'),
        ('A1BG', 'controls-expression-of', 'A2M'): ('affects', 'expression'),
        ('CHEBI:10357', 'chemical-affects', 'ACACA'): ('affects', None),
        ('A1CF', 'controls-production-of', 'chebi:17775'): ('affects', None),
        ('AAAS', 'controls-transport-of-chemical', 'chebi:1604'): ('affects', 'localization'),
    }
    for key, (predicate, aspect) in expected.items():
        _, record = map_edge(tmp_path, *key)
        relation = extract(tmp_path, read(tmp_path, edges=[edge(*key)])).relations[0]
        assert relation.predicate == predicate, key
        assert values(relation, 'object_aspect_qualifier') == ([aspect] if aspect else []), key
        # The source's own type stays on the relation.
        assert values(relation, 'pathwaycommons:interaction_type') == [key[1]], key
        # No type states a sign or a physical contact.
        assert not values(relation, 'object_direction_qualifier'), key
        assert relation.predicate != 'directly_physically_interacts_with'


def test_types_between_reactions_or_without_a_faithful_predicate_are_left_out(tmp_path):
    for key in (
        ('A4GALT', 'catalysis-precedes', 'B3GNT3'),
        ('CHEBI:11684', 'consumption-controlled-by', 'DIO1'),
        ('CHEBI:136982', 'used-to-produce', 'chebi:16189'),
        ('CHEBI:143180', 'reacts-with', 'chebi:15824'),
    ):
        row, record = map_edge(tmp_path, *key)
        assert record is None, key
        assert pc.exclusion_reason(row) == 'relation type without a faithful predicate'
        assert not extract(tmp_path, [row]).relations


def test_an_unknown_relation_type_raises(tmp_path):
    line = edge('A1BG', 'interacts-with', 'ABCC6').replace('interacts-with', 'neighbor-of')
    (row,) = read(tmp_path, edges=[line])
    with pytest.raises(ValueError, match="Unknown Pathway Commons relation type: 'neighbor-of'"):
        pc.resource.interactions.mapper(row)


def test_data_sources_become_resources_and_publications_pubmed_ids(tmp_path):
    row, _ = map_edge(tmp_path, 'A1CF', 'in-complex-with', 'APOBEC1')
    (relation,) = extract(tmp_path, [row]).relations
    assert values(relation, 'supporting_data_source') == ['IntAct', 'Reactome']
    assert values(relation, 'publications') == [
        'PMID:11134005', 'PMID:11352648', 'PMID:11727199', 'PMID:12881431',
        'PMID:14681455', 'PMID:21496894', 'PMID:8999813',
    ]
    # Every source is carried as the file spells it, "pid" and "CTD" too.
    ctd, _ = map_edge(tmp_path, 'CHEBI:10357', 'chemical-affects', 'ACACA')
    assert values(extract(tmp_path, [ctd]).relations[0], 'supporting_data_source') == ['CTD']
    pid, _ = map_edge(tmp_path, 'A1BG', 'controls-expression-of', 'A2M')
    assert values(extract(tmp_path, [pid]).relations[0], 'supporting_data_source') == ['pid']
    # A value that is not a PubMed ID is not a publication.
    bad = {**row, 'pubmed': ['11134005', '888800001068', '', '0123']}
    assert [a.value for a in pc.map_interaction(bad).annotations
            if a.term == pc.slots.publications] == ['PMID:11134005']
    with pytest.raises(ValueError, match="Unknown Pathway Commons data source: 'NewDB'"):
        pc.map_interaction({**row, 'data_sources': ['NewDB']})


def test_participants_are_proteins_genes_or_chemicals_by_their_reference_types(tmp_path):
    protein = pc.map_interaction(read(tmp_path, edges=[EDGES[1]])[0]).subject
    assert protein.type.__name__ == 'Protein'
    assert identifiers(protein) == [
        (Namespace.UNIPROT, 'P04217'), (Namespace.GENESYMBOL, 'A1BG')
    ]
    rna = pc.map_interaction(read(tmp_path, edges=[edge('AGO2', 'interacts-with', 'CASC2')])[0])
    # Only RNA references: a gene known by its symbol, with no protein accession.
    assert rna.object.type.__name__ == 'Gene'
    assert identifiers(rna.object) == [(Namespace.GENESYMBOL, 'CASC2')]
    chemical = pc.map_interaction(
        read(tmp_path, edges=[edge('CHEBI:10357', 'chemical-affects', 'ACACA')])[0]
    )
    assert chemical.subject.type.__name__ == 'ChemicalEntity'
    # ``CHEBI:10357`` and ``chebi:10357`` are one ChEBI identifier; no taxon.
    assert identifiers(chemical.subject) == [
        (Namespace.CHEBI, 'CHEBI:10357'), (Namespace.NAME, 'caryophyllene')
    ]
    assert not chemical.subject.annotations
    out = extract(tmp_path, read(tmp_path, edges=[EDGES[1], EDGES[7]]))
    assert {e.taxon for e in out.entities.values() if e.entity_type != 'chemical_entity'} == {'9606'}


def test_participant_with_several_uniprot_accessions_keeps_only_its_symbol(tmp_path):
    record = pc.map_interaction(
        read(tmp_path, edges=[edge('ABCF3', 'interacts-with', 'MOCS2')])[0]
    )
    assert identifiers(record.object) == [(Namespace.GENESYMBOL, 'MOCS2')]
    assert identifiers(record.subject) == [
        (Namespace.UNIPROT, 'Q9NUQ8'), (Namespace.GENESYMBOL, 'ABCF3')
    ]


def test_a_chebi_form_that_is_not_a_number_gives_the_name_only(tmp_path):
    record = pc.map_interaction(
        read(tmp_path, edges=[edge('CYP2A13', 'controls-production-of', 'CHEBI:PubChem:25245599')])[0]
    )
    assert identifiers(record.object) == [(Namespace.NAME, '2-hydroxynicotine')]


def test_a_chemical_where_a_gene_is_expected_raises(tmp_path):
    (row,) = read(tmp_path, edges=[edge('A1BG', 'interacts-with', 'ABCC6')])
    row = {**row, 'interaction_type': 'chemical-affects'}
    with pytest.raises(ValueError, match='relation with a protein as participant_a'):
        pc.map_interaction(row)


def test_parser_reads_a_cell_longer_than_the_default_csv_limit(tmp_path):
    # The longest MEDIATOR_IDS cell of release 14 has 252648 characters.
    mediators = ';'.join(f'biogrid:MolecularInteraction_{n:08d}' for n in range(8000))
    assert len(mediators) > 131072
    line = EDGES[1].rsplit('\t', 1)[0] + '\t' + mediators
    (row,) = read(tmp_path, edges=[line])
    assert row['mediator_ids'] == mediators
