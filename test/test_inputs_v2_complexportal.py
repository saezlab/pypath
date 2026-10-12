"""Complex Portal ComplexTAB: parser, mapper and extractor on rows of the files."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import complexportal
from pypath.inputs_v2.parsers.complexportal import COLUMNS, iter_complextab

# Rows of 9606.tsv (human, curated) and of 9606_predicted.tsv (human,
# predicted) of the release of 14 January 2026, with the columns the parser
# reads, in file order.
ROWS = {
    # 9606.tsv
    'CPX-17': (
        'CPX-17',
        'Telomerase catalytic core complex',
        '-',
        '9606',
        'O14746(2)|URS00004A7003_9606(2)',
        'ECO:0000353(physical interaction evidence used in manual assertion)',
        'intact:EBI-8789396',
        'GO:0000722(telomere maintenance via recombination)|GO:0007004(telomere maintenance via telomerase)|GO:0000333(telomerase catalytic core complex)|GO:0003720(telomerase activity)',
        'qsproteome:QS241426(identity)|efo:MONDO:0004992(see-also)|efo:MONDO:0019499(see-also)|efo:EFO:0000768(see-also)|efo:MONDO:0015780(see-also)|pubmed:22057212(see-also)|intenz:2.7.7.49(identity)|pubmed:11432839(see-also)|complex portal:CPX-17(complex-primary)',
        "A ribonucleoprotein complex essential for the replication of chromosome termini in most eukaryotes. Catalytic component of the telomerase holoenzyme complex whose main activity is the elongation of telomeres. Acts as a reverse transcriptase that adds simple sequence repeats to chromosome ends by copying a template sequence within the RNA component of the enzyme. Catalyzes the RNA-dependent extension of 3'-chromosomal termini with the 6-nucleotide telomeric repeat unit, 5'-TTAGGG-3'. The catalytic cycle involves primer binding, primer extension and release of product once the template boundary has been reached or nascent product translocation followed by further extension. Overexpressed telomerase results in telomere lengthening which may lead to cell immortalization and cancer cell pathogenesis.",
        'psi-mi:"MI:0469"(IntAct)',
    ),
    # 9606.tsv
    'CPX-37': (
        'CPX-37',
        'Calprotectin heterotetramer',
        '-',
        '9606',
        'CHEBI:29035(1)|CHEBI:29105(0)|CHEBI:29108(0)|P05109(2)|P06702(2)',
        'ECO:0000353(physical interaction evidence used in manual assertion)',
        'intact:EBI-10097914',
        'GO:0034121(regulation of toll-like receptor signaling pathway)|GO:0035662(Toll-like receptor 4 binding)|GO:1990660(calprotectin complex)|GO:0050786(RAGE receptor binding)',
        'qsproteome:QS49793(identity)|signor:SIGNOR-C293(identity)|efo:MONDO:0009061(see-also)|efo:EFO:0000685(see-also)|efo:EFO:0006505(see-also)|pubmed:17050004(see-also)|wwpdb:4GGF(identity)|pubmed:25597447(see-also)|complex portal:CPX-37(complex-primary)|wwpdb:1XK4(identity)|wwpdb:8sjc(identity)|wwpdb:5w1f(identity)|wwpdb:4xjk(identity)|wwpdb:8sjb(identity)|wwpdb:6ds2(identity)|pubmed:29890074(see-also)',
        'A Ca(2+), Mn(2+) and Zn(2+)-binding complex used by the innate immune system in a metal-withholding strategy that limits Mn(2+) and Zn(2+) availability at sites of infection. Calprotectin is expressed and released by neutrophils and epithelial cells, and exhibits broad-spectrum antimicrobial activity attributed to its metal-binding properties. Upon neutrophil activation or endothelial adhesion of monocytes, Calprotectin becomes secreted via a microtubule-mediated pathway and can thus serve as a marker for the influx of mononuclear phagocytes into the site of inflammation. Calprotectin is an endogenous ligand of toll-like receptor 4 (TLR4) and of the receptor for advanced glycation end products (RAGE) initiating signal transduction through NF-kappa-B pathways. ',
        'psi-mi:"MI:0469"(IntAct)',
    ),
    # 9606.tsv
    'CPX-39': (
        'CPX-39',
        'Calprotectin heterodimer',
        '-',
        '9606',
        'CHEBI:29035(1)|CHEBI:29105(0)|CHEBI:29108(0)|P05109(1)|P06702(1)',
        'ECO:0005546(biological system reconstruction evidence based on paralogy evidence used in manual assertion)',
        '-',
        'GO:0034121(regulation of toll-like receptor signaling pathway)|GO:0035662(Toll-like receptor 4 binding)|GO:0050786(RAGE receptor binding)|GO:1990660(calprotectin complex)',
        'qsproteome:QS49792(identity)|efo:EFO:0000685(see-also)|efo:EFO:0003767(see-also)|efo:HP:0000347(see-also)|efo:MONDO:0009061(see-also)|pubmed:17050004(see-also)|wwpdb:1XK4(identity)|pubmed:25597447(see-also)|complex portal:CPX-37(inferred-from)|complex portal:CPX-39(complex-primary)',
        'A Ca(2+), Mn(2+) and Zn(2+)-binding complex used by the innate immune system in a metal-withholding strategy that limits Mn(2+) and Zn(2+) availability at sites of infection. Calprotectin is expressed and released by neutrophils and epithelial cells, and exhibits broad-spectrum antimicrobial activity attributed to its metal-binding properties. Upon neutrophil activation or endothelial adhesion of monocytes, Calprotectin becomes secreted via a microtubule-mediated pathway and can thus serve as a marker for the influx of mononuclear phagocytes into the site of inflammation. Calprotectin is an endogenous ligand of toll-like receptor 4 (TLR4) and of the receptor for advanced glycation end products (RAGE) initiating signal transduction through NF-kappa-B pathways. ',
        'psi-mi:"MI:0469"(IntAct)',
    ),
    # 9606.tsv
    'CPX-102': (
        'CPX-102',
        'DAPK1 - calmodulin complex',
        'Death-associated protein kinase 1 - calmodulin complex',
        '9606',
        'CHEBI:29108(1)|P53355(1)|[P0DP23,P0DP24,P0DP25](1)',
        'ECO:0000353(physical interaction evidence used in manual assertion)',
        'intact:EBI-9828233',
        'GO:0004683(calmodulin-dependent protein kinase activity)|GO:0004674(protein serine/threonine kinase activity)|GO:0005524(ATP binding)|GO:0002834(regulation of response to tumor cell)|GO:1904094(positive regulation of autophagic cell death)|GO:0043065(positive regulation of apoptotic process)|GO:1990722(DAPK1 - calmodulin complex)',
        'wwpdb:2x0g(identity)|wwpdb:2y4v(identity)|wwpdb:1yr5(identity)|intenz:2.7.11.1(identity)|complex portal:CPX-102(complex-primary)',
        'A serine/threonine protein kinase complex involved in cell survival, apoptosis and autophagic cell death pathways. DAPK1 is activated by dephosphorylation of Ser-308 and calcium-calmodulin binding. Complex activity is regulated via various phosphorylation sites, two of which (Ser-298 and Ser-308) are located on the autoregulatory domain (ARD) of DAPK1.',
        'psi-mi:"MI:2228"(ceitec)',
    ),
    # 9606.tsv
    'CPX-168': (
        'CPX-168',
        'Neuronal nicotinic acetylcholine receptor complex, 3xalpha4-2xbeta2',
        'neuronal acetylcholine receptor; alpha4-beta2|neuronal nicotinic acetylcholine-gated ion channel complex|neuronal nicotinic cholinergic receptor complex, alpha4-beta2|neuronal nicotinic acetylcholine-gated receptor-channel complex|neuronal nicotinic acetylcholine-gated cation channel complex|neuronal nicotinic acetylcholine-gated receptor complex|ionotropic acetylcholine receptor complex',
        '9606',
        'P17787-PRO_0000000379(2)|P43681-PRO_0000000351(3)',
        'ECO:0000353(physical interaction evidence used in manual assertion)',
        'intact:EBI-10690292',
        'GO:0005892(nicotinic acetylcholine-gated receptor-channel complex)|GO:0022848(acetylcholine-gated cation channel activity)|GO:0051899(membrane depolarization)|GO:0042166(acetylcholine binding)',
        'qsproteome:QS122603(identity)|efo:MONDO:0004975(see-also)|efo:MONDO:0005180(see-also)|efo:EFO:0003768(see-also)|efo:EFO:0000474(see-also)|pubmed:15470079(see-also)|pubmed:25869137(see-also)|pubmed:19126755(see-also)|pubmed:10353988(see-also)|pubmed:19481063(see-also)|pubmed:25460185(see-also)|reactome:R-HSA-532625(identity)|chembl target:CHEMBL1907589(identity)|complex portal:CPX-168(complex-primary)|rhea:RHEA:29463(identity)|rhea:RHEA:29671(identity)|rhea:RHEA:34963(identity)',
        'A ligand-gated ion channel receptor complex that is sensitive to endogenous achetylcholine and exogenous nicotinic agents. Ligand binding causes a conformational change of the receptor resulting in the formation of a pore that is permeable to Na+ and/or Ca2+ (inward) and K+ (outward) thereby depolarising the neuron. Mediates fast, short-lived, mainly pre-synaptic transmission of neurotransmitters. Major receptor in Central Nervous System and predominantly found in cerebellum, cortex, forebrain, hippocampus, mesencephalon, striatum, superior colliculus and thalamus. Up-regulated by pro-inflammatory cytokines, for example TNF-alpha.',
        'psi-mi:"MI:0469"(IntAct)',
    ),
    # 9606.tsv
    'CPX-578': (
        'CPX-578',
        'Cyclin-dependent protein kinase-activating kinase complex',
        'CDK-activating kinase complex|CAK complex',
        '9606',
        'P50613(1)|P51946(1)|P51948(1)',
        'ECO:0000353(physical interaction evidence used in manual assertion)',
        'intact:EBI-1245996',
        'GO:0097472(cyclin-dependent protein kinase activity)|GO:0005524(ATP binding)|GO:2000045(regulation of G1/S transition of mitotic cell cycle)|GO:0000307(cyclin-dependent protein kinase holoenzyme complex)',
        'wwpdb:6XBZ(identity)|emdb:EMD-22123(identity)|pubmed:32855301(see-also)|wwpdb:6xd3(identity)|emdb:EMD-22131(identity)|qsproteome:QS220938(identity)|reactome:R-HSA-69221(identity)|intenz:2.7.11.22(identity)|intenz:2.7.11.23(identity)|wwpdb:8p76(identity)|emdb:EMD-17513(identity)|wwpdb:8p70(identity)|emdb:EMD-17514(identity)|wwpdb:8p72(identity)|emdb:EMD-17512(identity)|emdb:EMD-17520(identity)|emdb:EMD-17518(identity)|emdb:8p75(identity)|emdb:EMD-17517(identity)|emdb:EMD-17754(identity)|wwpdb:8plz(identity)|emdb:EMD-12042(identity)|wwpdb:7b5q(identity)|emdb:EMD-17516(identity)|wwpdb:8p73(identity)|wwpdb:8p74(identity)|emdb:EMD-17519(identity)|wwpdb:7b5o(identity)|pubmed:8p6y(see-also)|emdb:EMD-17511(identity)|emdb:EMD-17522(identity)|wwpdb:8p79(identity)|pubmed:38480681(see-also)|emdb:EMD-17510(identity)|wwpdb:8p6x(identity)|wwpdb:8p77(identity)|emdb:EMD-17508(identity)|wwpdb:8p6v(identity)|emdb:EMD-17521(identity)|wwpdb:8p78(identity)|emdb:EMD-17509(identity)|wwpdb:8p6w(identity)|wwpdb:8p6z(identity)|emdb:EMD-17515(identity)|emdb:EMD-17129(identity)|wwpdb:8p71(identity)|wwpdb:8orm(identity)|rhea:Rhea 10216(identity)|rhea:Rhea 17989(identity)|rhea:Rhea 46608(identity)|complex portal:CPX-578(complex-primary)|signor:SIGNOR-C456(identity)|chembl target:CHEMBL3038473(identity)|chembl target:CHEMBL3430895(identity)',
        'Cyclin-dependent kinase (CDK) activation minimally depends on two events: binding to a cyclin and phosphorylation of a conserved Thr residue in the activation (T) loop. CAK activates the cyclin-associated kinases CDK1 (P06493), CDK2 (P24941), CDK4 (P11802) and CDK6 Q00534) by threonine phosphorylation, thus regulating cell cycle progression. CAK activity is itself regulated throughout the cell cycle by T-loop phosphorylation of CDK7 on Thr-170. Phosphorylation of Ser-164 during mitosis inactivates the enzyme. In the transcription cycle, CAK serine phosphorylates the carboxyl-terminal domain (CTD) of RNA polymerase II (POLR2A, P24928) and other proteins, as part of the general transcription factor TFIIH. Phosphorylation of POLR2A in complex with DNA promotes transcription initiation by triggering dissociation from DNA. CAK also phosphorylates CDK9 (P-TEFb, P50750), which releases POLR2A from the promoter and and enables elongation of the transcripts.',
        'psi-mi:"MI:0469"(IntAct)',
    ),
    # 9606.tsv
    'CPX-1265': (
        'CPX-1265',
        'Laminin111-nidogen complex',
        'Laminin-entactin complex|Laminin111-nidogen1 complex|Laminin-nidogen complex',
        '9606',
        'CPX-1770(1)|P14543(1)',
        'ECO:0005547(biological system reconstruction evidence based on inference from background scientific knowledge used in manual assertion)',
        '-',
        'GO:0110011(regulation of basement membrane organization)|GO:0005509(calcium ion binding)|GO:0098637(protein complex involved in cell-matrix adhesion)|GO:0005604(basement membrane)|GO:0051149(positive regulation of muscle cell differentiation)|GO:0045785(positive regulation of cell adhesion)|GO:2001046(positive regulation of integrin-mediated signaling pathway)|GO:0005518(collagen binding)|GO:0005178(integrin binding)|GO:0005201(extracellular matrix structural constituent)',
        'pubmed:25387007(see-also)|pubmed:7556571(see-also)|pubmed:24727304(see-also)|pubmed:1350783(see-also)|matrixdb:MULT_2_human(identity)|pubmed:3109910(see-also)|complex portal:CPX-1265(complex-primary)',
        'An extracellular matrix complex responsible for basement membrane stabilization and cell-matrix adhesion. Crosslinking of the nidogen subunit from one complex to a laminin arm of a second by tissue transglutaminases forms large assemblies. Facilitates cell adhesion by binding collagens (e.g. collagen type I, CPX-1650 or collagen type IV, CPX-1723) and integrins (e.g. alpha3beta1, CPX-1797 or alphavbeta3, CPX-1795). May modify type I collagen scaffolds and thereby enhance myotube formation. ',
        'psi-mi:"MI:0469"(IntAct)',
    ),
    # 9606.tsv
    'CPX-1282': (
        'CPX-1282',
        'Laminin211-nidogen complex',
        'Laminin-entactin complex|Laminin211-nidogen1 complex|Laminin-nidogen complex',
        '9606',
        'CPX-1771(0)|P14543(0)',
        'ECO:0005547(biological system reconstruction evidence based on inference from background scientific knowledge used in manual assertion)',
        '-',
        'GO:0005509(calcium ion binding)|GO:0110011(regulation of basement membrane organization)|GO:0098637(protein complex involved in cell-matrix adhesion)|GO:0005604(basement membrane)|GO:0051149(positive regulation of muscle cell differentiation)|GO:0045785(positive regulation of cell adhesion)|GO:2001046(positive regulation of integrin-mediated signaling pathway)|GO:0005518(collagen binding)|GO:0005178(integrin binding)|GO:0005201(extracellular matrix structural constituent)',
        'pubmed:25387007(see-also)|pubmed:7556571(see-also)|pubmed:24727304(see-also)|pubmed:1350783(see-also)|pubmed:3109910(see-also)|complex portal:CPX-1282(complex-primary)',
        'An extracellular matrix complex responsible for basement membrane stabilization and cell-matrix adhesion. Crosslinking of the nidogen subunit from one complex to a laminin arm of a second by tissue transglutaminases forms large assemblies. Facilitates cell adhesion by binding collagens (e.g. collagen type I, CPX-1650 or collagen type IV, CPX-1723) and integrins (e.g. alpha3beta1, CPX-1797 or alphavbeta3, CPX-1795). May modify type I collagen scaffolds and thereby enhance myotube formation. ',
        'psi-mi:"MI:0469"(IntAct)',
    ),
    # 9606.tsv
    'CPX-5983': (
        'CPX-5983',
        'Phosphatidylinositol 3-kinase complex class IA, p110delta/p85alpha',
        '-',
        '9606',
        'O00329(1)|P27986(1)',
        'ECO:0000353(physical interaction evidence used in manual assertion)',
        'intact:EBI-15939463',
        'GO:0030183(B cell differentiation)|GO:0006955(immune response)|GO:0030217(T cell differentiation)|GO:0005943(phosphatidylinositol 3-kinase complex, class IA)|GO:0035005(1-phosphatidylinositol-4-phosphate 3-kinase activity)',
        'qsproteome:QS119719(identity)|humap:huMAP3_05061.1(superset)|chembl:CHEMBL2111432(identity)|reactome:R-HSA-1045152(identity)|efo:MONDO:0018338(see-also)|pubmed:20379207(see-also)|pubmed:24587488(see-also)|intenz:2.7.1.153(identity)|rhea:RHEA:21292(identity)|complex portal:CPX-5983(complex-primary)|wwpdb:6pyr(identity)|wwpdb:6pyu(identity)|wwpdb:5ubt(identity)|wwpdb:5m6u(identity)|wwpdb:5vlr(identity)',
        'Uses 1-phosphatidyl-1D-myo-inositol 4,5-bisphosphate/PI(4,5)P2 (CHEBI:18348) as a substrate to generate the product 1-phosphatidyl-1D-myo-inositol 3,4,5-trisphosphate/PI(3,4,5)P3 (CHEBI:16618) which initiates major signaling pathways downstream. Performs a central role in cellular signaling, mediating responses to growth factors and extracellular stimuli. Expressed predominantly in leukocytes, where it mediates immune responses. Gain-of-function mutations in the phosphoinositide 3-kinase (PI3K) genes PIK3CD and PIK3R1 can cause a combined immunodeficiency syndrome, referred to as activated PI3Kdelta syndrome (APDS).',
        'psi-mi:"MI:0469"(IntAct)',
    ),
    # 9606.tsv
    'CPX-5984': (
        'CPX-5984',
        'Phosphatidylinositol 3-kinase complex class IA, p110delta/p55alpha',
        '-',
        '9606',
        'O00329(1)|P27986(1)',
        'ECO:0000353(physical interaction evidence used in manual assertion)',
        '-',
        'GO:0030183(B cell differentiation)|GO:0006955(immune response)|GO:0030217(T cell differentiation)|GO:0005943(phosphatidylinositol 3-kinase complex, class IA)|GO:0035005(1-phosphatidylinositol-4-phosphate 3-kinase activity)',
        'qsproteome:QS119719(identity)|humap:huMAP3_05061.1(superset)|chembl:CHEMBL2111432(identity)|reactome:R-HSA-1045152(identity)|efo:MONDO:0018338(see-also)|pubmed:20379207(see-also)|pubmed:24587488(see-also)|intenz:2.7.1.153(identity)|rhea:RHEA:21292(identity)|wwpdb:6pyr(identity)|wwpdb:6pyu(identity)|wwpdb:5ubt(identity)|wwpdb:5m6u(identity)|wwpdb:5vlr(identity)|complex portal:CPX-5984(complex-primary)',
        'Uses PI(4,5)P2 as a substrate to generate the product PI(3,4,5)P3 which initiates major signaling pathways downstream. Performs a central role in cellular signaling, mediating responses to growth factors and extracellular stimuli. Expressed predominantly in leukocytes, where it mediates immune responses. Gain-of-function mutations in the phosphoinositide 3-kinase (PI3K) genes PIK3CD and PIK3R1 can cause a combined immunodeficiency syndrome, referred to as activated PI3Kdelta syndrome (APDS).',
        'psi-mi:"MI:0469"(IntAct)',
    ),
    # 9606.tsv
    'CPX-4747': (
        'CPX-4747',
        'E3 ligase (RANBP2) complex',
        'RanBP2/RanGAP1/Ubc9 complex|RanBP2 E3 SUMO ligase',
        '9606',
        'P46060(1)|P49792(1)|P63279(1)',
        'ECO:0000353(physical interaction evidence used in manual assertion)',
        'wwpdb:3UIP',
        'GO:0005643(nuclear pore)|GO:0106068(SUMO ligase complex)|GO:0061665(SUMO ligase activity)|GO:0051168(nuclear export)|GO:0016925(protein sumoylation)',
        'qsproteome:QS141811(identity)|wwpdb:5FQ2(subset)|complex portal:CPX-4747(complex-primary)|intact:EBI-22028918(secondary-ac)|pubmed:22464730(see-also)|pubmed:27160050(see-also)|pubmed:22194619(see-also)|wwpdb:2GRO(subset)|wwpdb:2GRN(subset)|wwpdb:2GRQ(subset)|wwpdb:2GRP(subset)|wwpdb:2GRR(subset)|wwpdb:2PE6(subset)|intact:EBI-22181174(secondary-ac)|complex portal:CPX-4922(secondary-ac)',
        'SUMOylated E3 ligase complex that takes part in the termination of CRM1-mediated export by facilitating the hydrolysis of GTP by RAN (P62826), resulting in export complex disassembly. The subunit RANBP2/NUP358 links the E3 ligase complex to the Nuclear Pore Complex (CPX-873). SUMOylation, mediated by the enzymatic activity of RANBP2, might play a role in the directionality of nucleo-cytoplasmic transport through the nuclear pore.',
        'psi-mi:"MI:0469"(IntAct)',
    ),
    # 9606.tsv
    'CPX-2392': (
        'CPX-2392',
        'U1 small nuclear ribonucleoprotein complex',
        'U1 snRNP complex',
        '9606',
        'M14387(1)|P08621(1)|P09012(1)|P09234(1)|P14678(1)|P62304(1)|P62306(1)|P62308(1)|P62314(1)|P62316(1)|P62318(1)',
        'ECO:0005547(biological system reconstruction evidence based on inference from background scientific knowledge used in manual assertion)',
        '-',
        'GO:0005685(U1 snRNP)|GO:0005634(nucleus)|GO:0030619(U1 snRNA binding)|GO:0000398(mRNA splicing, via spliceosome)|GO:0005681(spliceosomal complex)',
        'qsproteome:QS121923(identity)|pubmed:30765414(see-also)|wwpdb:3cw1(subset)|complex portal:CPX-2392(complex-primary)|wwpdb:3pgw(subset)|pubmed:21113136(see-also)|pubmed:33677607(see-also)|reactome:R-HSA-71917(identity)|signor:SIGNOR-C480(identity)',
        "Recognizes the 5' splice site within precursor messenger RNAs and initiates the assembly of the spliceosome. U1 snRNP interacts with the pre-mRNA, via base-pairing between the 5' end of U1 snRNA and the 5' splice site in pre-mRNA. The spliceosome is a highly dynamic structure, assembled by sequential binding and release of the small nuclear RNAs and protein factors which removes intronic sequence from pre-mRNA.",
        'psi-mi:"MI:0469"(IntAct)',
    ),
    # 9606_predicted.tsv
    'CPX-10560': (
        'CPX-10560',
        '-',
        '-',
        '9606',
        'Q13349(0)|Q9BW72(0)',
        'ECO:0008004(machine learning method evidence used in automatic assertion)',
        '-',
        '-',
        'qsproteome:QS201849(identity)|complex portal:CPX-10560(complex-primary)|humap:huMAP3_05475.1(identity)',
        '-',
        'psi-mi:"MI:2424"(HuMap)',
    ),
    # 9606_predicted.tsv
    'CPX-26327': (
        'CPX-26327',
        'U2-type spliceosomal complex family',
        '-',
        '9606',
        'O60508(0)|O75494(0)|P09234(0)|P26368(0)|Q00839(0)|Q15287(0)|Q8IYB3(0)|Q9UQ35(0)',
        'ECO:0008004(machine learning method evidence used in automatic assertion)',
        '-',
        '-',
        'complex portal:CPX-26327(complex-primary)|music:750(identity)|cell line ontology:CLO:0037372(identity)|pubmed:34819669(see-also)',
        '-',
        'psi-mi:"-"(MuSIC)',
    ),
}


def table(*accessions):
    lines = ['\t'.join(COLUMNS)]
    for accession in accessions:
        lines.append('\t'.join(ROWS[accession]))
    return '\n'.join(lines) + '\n'


def parse(*accessions, **kwargs):
    text = table(*accessions, **kwargs)
    return list(iter_complextab(SimpleNamespace(result=io.StringIO(text))))


def extract(*accessions, dataset='human_complexes'):
    out = SilverExtractor('complexportal', dataset)
    for index, row in enumerate(parse(*accessions)):
        record = getattr(complexportal.resource, dataset).mapper(row)
        out.process_record(record, row, f'fixture:{index}', index)
    return out


def identifiers(entity):
    return [(i['ns'], i['id']) for i in entity.identifiers]


def terms(entity, term):
    return [a['value'] for a in entity.annotations if a['term'] == term]


def complex_of(out, accession):
    (entity,) = [
        e
        for e in out.entities.values()
        if e.entity_type == 'macromolecular_complex'
        and (Namespace.COMPLEXPORTAL, accession) in identifiers(e)
        and (e.evidence or e.annotations)
    ]
    return entity


def members(out, entity):
    """(type, first identifier, stoichiometry annotations) of the members of an entity."""
    return [
        (
            out.entities[r.object_entity_key].entity_type,
            identifiers(out.entities[r.object_entity_key])[0],
            [a['value'] for a in r.annotations if a['term'] == 'stoichiometry'],
        )
        for r in out.relations
        if r.predicate == 'has_member' and r.subject_entity_key == entity.entity_key
    ]


def test_parser_decodes_members_with_stoichiometry_and_sets():
    row = parse('CPX-102')[0]
    assert row['accession'] == 'CPX-102'
    assert row['name'] == 'DAPK1 - calmodulin complex'
    assert row['synonyms'] == ['Death-associated protein kinase 1 - calmodulin complex']
    assert row['taxon'] == '9606'
    assert row['members'] == [
        {'accession': 'CHEBI:29108', 'set': None, 'stoichiometry': 1},
        {'accession': 'P53355', 'set': None, 'stoichiometry': 1},
        {'accession': None, 'set': ['P0DP23', 'P0DP24', 'P0DP25'], 'stoichiometry': 1},
    ]
    assert row['evidence_code'] == 'ECO:0000353'
    assert row['experimental_evidence'] == [('intact', 'EBI-9828233')]
    assert row['go_terms'][0] == 'GO:0004683'
    assert row['source'] == {'name': 'ceitec', 'mi': 'MI:2228'}
    assert {'database': 'wwpdb', 'accession': '2x0g', 'qualifier': 'identity'} in row['xrefs']


def test_parser_keeps_zero_stoichiometry_and_missing_values():
    row = parse('CPX-1282', 'CPX-10560')[0]
    assert [m['stoichiometry'] for m in row['members']] == [0, 0]
    predicted = parse('CPX-10560')[0]
    assert predicted['name'] is None
    assert predicted['synonyms'] == []
    assert predicted['description'] is None
    assert predicted['go_terms'] == []
    assert predicted['experimental_evidence'] == []


def test_parser_reads_quote_characters_as_data():
    values = list(ROWS['CPX-17'])
    values[1] = 'Telomerase "core" complex'
    text = '\t'.join(COLUMNS) + '\n' + '\t'.join(values) + '\n'
    (row,) = iter_complextab(SimpleNamespace(result=io.StringIO(text)))
    assert row['name'] == 'Telomerase "core" complex'


def test_parser_reads_a_curation_source_without_psi_mi_term():
    (row,) = parse('CPX-26327')
    assert row['source'] == {'name': 'MuSIC', 'mi': None}


def test_parser_raises_on_missing_column_empty_file_and_malformed_member():
    short = '\t'.join(COLUMNS[:-1]) + '\n' + '\t'.join(ROWS['CPX-17'][:-1]) + '\n'
    with pytest.raises(ValueError, match='lacks columns'):
        list(iter_complextab(SimpleNamespace(result=io.StringIO(short))))
    with pytest.raises(ValueError, match='no complexes'):
        parse()
    values = list(ROWS['CPX-17'])
    values[4] = 'O14746'
    text = '\t'.join(COLUMNS) + '\n' + '\t'.join(values) + '\n'
    with pytest.raises(ValueError, match='malformed molecule'):
        list(iter_complextab(SimpleNamespace(result=io.StringIO(text))))


def test_complex_has_every_member_with_its_type_and_stoichiometry():
    out = extract('CPX-102')
    entity = complex_of(out, 'CPX-102')
    assert entity.entity_type == 'macromolecular_complex'
    assert identifiers(entity) == [
        (Namespace.COMPLEXPORTAL, 'CPX-102'),
        (Namespace.NAME, 'DAPK1 - calmodulin complex'),
        (Namespace.SYNONYM, 'Death-associated protein kinase 1 - calmodulin complex'),
    ]
    assert members(out, entity) == [
        ('chemical_entity', (Namespace.CHEBI, '29108'), ['1']),
        ('protein', (Namespace.UNIPROT, 'P53355'), ['1']),
        ('protein_family', (Namespace.NAME, '[P0DP23,P0DP24,P0DP25]'), ['1']),
    ]
    assert terms(entity, 'in_taxon') == ['NCBITaxon:9606']
    # The species of the molecules is not stated.
    assert all(
        e.taxon is None for e in out.entities.values() if e.entity_type != 'macromolecular_complex'
    )


def test_set_of_molecules_is_a_family_of_its_proteins():
    out = extract('CPX-102')
    (family,) = [e for e in out.entities.values() if e.entity_type == 'protein_family']
    assert members(out, family) == [
        ('protein', (Namespace.UNIPROT, accession), [])
        for accession in ('P0DP23', 'P0DP24', 'P0DP25')
    ]
    # The proteins are not members of the complex.
    assert len(members(out, complex_of(out, 'CPX-102'))) == 3


def test_stoichiometry_zero_is_unknown_and_not_emitted():
    out = extract('CPX-1282')
    entity = complex_of(out, 'CPX-1282')
    assert members(out, entity) == [
        ('macromolecular_complex', (Namespace.COMPLEXPORTAL, 'CPX-1771'), []),
        ('protein', (Namespace.UNIPROT, 'P14543'), []),
    ]
    assert not [
        a for r in out.relations for a in r.annotations if a['term'] == 'stoichiometry'
    ]


def test_sub_complex_stays_a_member_and_is_not_flattened():
    out = extract('CPX-17', 'CPX-2392')
    entity = complex_of(out, 'CPX-17')
    assert members(out, entity) == [
        ('protein', (Namespace.UNIPROT, 'O14746'), ['2']),
        ('rna_product', (Namespace.RNACENTRAL, 'URS00004A7003_9606'), ['2']),
    ]
    assert not [r for r in out.relations if r.predicate not in {'has_member', 'associated_with'}]


def test_rna_has_the_taxon_of_its_identifier_and_nucleotide_record_is_a_nucleic_acid():
    out = extract('CPX-17', 'CPX-2392')
    (rna,) = [e for e in out.entities.values() if e.entity_type == 'rna_product']
    assert rna.taxon == '9606'
    (record,) = [e for e in out.entities.values() if e.entity_type == 'nucleic_acid_entity']
    assert identifiers(record) == [(Namespace.GENBANK, 'M14387')]
    assert record.taxon is None


def test_isoform_and_chain_accessions_are_kept():
    out = extract('CPX-168')
    entity = complex_of(out, 'CPX-168')
    assert [(m[1], m[2]) for m in members(out, entity)] == [
        ((Namespace.UNIPROT, 'P17787-PRO_0000000379'), ['2']),
        ((Namespace.UNIPROT, 'P43681-PRO_0000000351'), ['3']),
    ]


def test_sub_complex_reference_of_nested_complex():
    out = extract('CPX-1265')
    entity = complex_of(out, 'CPX-1265')
    assert [m[:2] for m in members(out, entity)] == [
        ('macromolecular_complex', (Namespace.COMPLEXPORTAL, 'CPX-1770')),
        ('protein', (Namespace.UNIPROT, 'P14543')),
    ]


def test_complexes_with_the_same_members_stay_distinct():
    out = extract('CPX-5983', 'CPX-5984', 'CPX-37', 'CPX-39')
    first, second = complex_of(out, 'CPX-5983'), complex_of(out, 'CPX-5984')
    assert first.entity_key != second.entity_key
    assert [m[:2] for m in members(out, first)] == [m[:2] for m in members(out, second)]
    # The accession of each is the identifier of its complex.
    assert (Namespace.COMPLEXPORTAL, 'CPX-5983') not in identifiers(second)
    # The members of CPX-37 and CPX-39 are the same; their stoichiometry differs.
    heterotetramer, heterodimer = complex_of(out, 'CPX-37'), complex_of(out, 'CPX-39')
    assert heterotetramer.entity_key != heterodimer.entity_key
    assert [m[2] for m in members(out, heterotetramer)][-2:] == [['2'], ['2']]
    assert [m[2] for m in members(out, heterodimer)][-2:] == [['1'], ['1']]
    # The shared proteins and small molecules are single entities.
    assert len([e for e in out.entities.values() if e.entity_type == 'protein']) == 4


def test_evidence_go_publications_enzyme_and_identical_objects():
    out = extract('CPX-17')
    entity = complex_of(out, 'CPX-17')
    assert terms(entity, 'has_evidence_of_type') == ['ECO:0000353']
    assert terms(entity, 'supporting_data_source') == ['IntAct']
    assert terms(entity, 'publications') == ['PMID:22057212', 'PMID:11432839']
    assert terms(entity, 'has_topic') == ['EC:2.7.7.49']
    assert terms(entity, 'xref') == ['QSPROTEOME:QS241426']
    assert terms(entity, 'description').pop().startswith('A ribonucleoprotein complex')
    go = [
        out.entities[r.object_entity_key]
        for r in out.relations
        if r.predicate == 'associated_with'
    ]
    assert [identifiers(e)[0] for e in go] == [
        (Namespace.GO, 'GO:0000722'),
        (Namespace.GO, 'GO:0007004'),
        (Namespace.GO, 'GO:0000333'),
        (Namespace.GO, 'GO:0003720'),
    ]
    assert {e.entity_type for e in go} == {'ontology_class'}


def test_pubmed_number_with_a_leading_zero_or_nine_digits_is_not_a_publication(
    monkeypatch,
):
    row = tuple(
        cell.replace('pubmed:22057212', 'pubmed:0123456').replace(
            'pubmed:11432839', 'pubmed:123456789'
        )
        for cell in ROWS['CPX-17']
    )
    monkeypatch.setitem(ROWS, 'CPX-17', row)
    assert terms(complex_of(extract('CPX-17'), 'CPX-17'), 'publications') == []


def test_partial_matches_and_disease_links_are_not_cross_references():
    out = extract('CPX-2392', 'CPX-39')
    xrefs = terms(complex_of(out, 'CPX-2392'), 'xref')
    # Structures of part of the complex (qualifier subset) are not identical objects.
    assert xrefs == ['QSPROTEOME:QS121923', 'REACTOME:R-HSA-71917', 'SIGNOR:SIGNOR-C480']
    assert terms(complex_of(out, 'CPX-39'), 'xref') == ['QSPROTEOME:QS49792', 'PDB:1XK4']
    # The disease and phenotype links (efo) and inferred-from links are not emitted.
    assert not [
        a for e in out.entities.values() for a in e.annotations if 'MONDO' in str(a['value'])
    ]


def test_secondary_accessions_are_identifiers_of_the_complex():
    out = extract('CPX-4747')
    entity = complex_of(out, 'CPX-4747')
    ids = identifiers(entity)
    assert ids[0] == (Namespace.COMPLEXPORTAL, 'CPX-4747')
    # The accession of a merged complex and the two IntAct accessions.
    assert ids[1:4] == [
        (Namespace.INTACT, 'EBI-22028918'),
        (Namespace.INTACT, 'EBI-22181174'),
        (Namespace.COMPLEXPORTAL, 'CPX-4922'),
    ]


def test_predicted_complex_without_name_has_only_its_accession():
    out = extract('CPX-10560', 'CPX-26327', dataset='human_predicted_complexes')
    entity = complex_of(out, 'CPX-10560')
    assert identifiers(entity) == [(Namespace.COMPLEXPORTAL, 'CPX-10560')]
    assert terms(entity, 'has_evidence_of_type') == ['ECO:0008004']
    assert terms(entity, 'supporting_data_source') == ['HuMap']
    assert terms(entity, 'xref') == ['QSPROTEOME:QS201849', 'HUMAP:huMAP3_05475.1']
    assert terms(entity, 'description') == []
    # Stoichiometry 0 of the predicted complexes is not a stoichiometry.
    assert [m[2] for m in members(out, entity)] == [[], []]
    music = complex_of(out, 'CPX-26327')
    assert terms(music, 'supporting_data_source') == ['MuSIC']
    assert terms(music, 'xref') == ['MUSIC:750']
    assert terms(music, 'publications') == ['PMID:34819669']
    # The cell line the map was built from is not an identical object.
    assert not [v for v in terms(music, 'xref') if 'CLO' in v]


def test_unknown_source_and_unknown_molecule_raise():
    row = parse('CPX-17')[0]
    with pytest.raises(ValueError, match='source'):
        complexportal.map_complex({**row, 'source': {'name': 'Unknown', 'mi': None}})
    with pytest.raises(ValueError, match='molecule'):
        complexportal.map_complex(
            {**row, 'members': [{'accession': 'X', 'set': None, 'stoichiometry': 1}]}
        )


def test_datasets_are_the_two_human_files():
    assert set(complexportal.resource.datasets()) == {
        'human_complexes',
        'human_predicted_complexes',
    }
    download = complexportal.resource.human_complexes.download
    assert download.url.endswith('/complextab/9606.tsv')
