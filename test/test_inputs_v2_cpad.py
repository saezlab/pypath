"""CPAD entries: parser, mapper and extractor on source rows."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import cpad
from pypath.inputs_v2.parsers.cpad import iter_entries, iter_terms

# The header ends in three empty column names; the lines end in three tabs and CRLF.
HEADER = (
    'Regulator\tRegulator_Type\tPathway\tPathway_Category\tKEGG_ID\tRegulation_Type\tCancer\tID\tOutcome_Description'
    '\tDescription\tMethod\tSpecies\tTissue\tCell line\tPMID\tTitle\tYear\t\t\t'
)
# Lines of CPAD-data.txt (descriptions and titles shortened).
ENTRIES = {
    'regulator': 'AEG-1\tUpregulation\tWnt/beta-catenin signaling pathway\tWnt signaling pathway\tpath:04310\tActivating\tBreast cancer\tH00031\tpromote the proliferation and invasion\tAEG1 was significantly overexpressed in breast cancer cells it promoted cell proliferation and invasion via activating the Wnt/β-catenin signaling pathway.\tRT-qPCR; CCK-8 assay; Colony formation assay\tHomo sapiens\tBreast\tMCF-7;MDA-MB-231;SK-BR-3\t28454408\tAstrocyte elevated gene-1 promotes the proliferation and invasion\t2017\t\t\t',
    'no regulator': 'NULL\tNULL\t5-LO/LTA4 hydrolase pathway\tOthers\tNULL\tInhibiting\tGlioma\tH00042\tpartial suppression of tumor growth\tWe confirmed the expression of 5-LO in various human brain tumors.\tImmunohistochemistry; RT-PCR\tHomo sapiens\tBrain\tA-172;U-251MG;U-138MG\t19919819\t5-lipoxygenase pathway promotes cell proliferation\t2010\t\t\t',
    'two regulators': 'AXIN and APC\tNULL\tWnt signaling pathway\tWnt signaling pathway\tpath:04310\tInhibiting\tColorectal cancer\tH00020\tinhibit cell growth\tCo-expression of APC5 and Axin genes.\tMTT assay; Western blot\tHomo sapiens\tColorectum\tSW480\t28731177\tCo-expression of Axin and APC\t2017\t\t\t',
    # The expression pattern 0, and a quoted description.
    'defect': '8-Bromo-7-methoxychrysin\t0\tAkt/FoxO3a signaling pathway\tPI3K-Akt signaling pathway; FoxO signaling pathway\tpath:04151; path:04068\tActivating\tOvarian cancer\tH00027\tinduce apoptosis\t"BrMC-induced apoptotic cell death occurred mainly by the activation of Akt, which was accompanied by the overexpression of FOXO3a."\tELISA; MTT assay\tHomo sapiens\tOvary\tA2780\t26151347\t8-bromo-7-methoxychrysin induces apoptosis\t2015\t\t\t',
    # Values named by regulator.
    'named values': 'MEN1 and DNMT1\tmenin(Upregulation);Dnmt1(Downregulation)\tHedgehog signaling pathway\tHedgehog signaling pathway\tpath:04340\tInhibiting\tPancreatic cancer\tH00019\tinhibit growth\t"Menin is a negative effector of Hedgehog signaling."\treal-time qPCR; Western blot\tHomo sapiens\tPancreas\tPANC-1;AsPC-1\t26454216\tInterplay between menin and Dnmt1\t2016\t\t\t',
    # A genus, not a species.
    'genus': 'FMR\tUpregulation\tJNK signaling pathway\tMAPK signaling pathway\tpath:04010\tActivating\tMalignant pleural mesothelioma\tH00015\tresult in cell death\t"A novel link between FMR gene and the JNK pathway."\tBeta-galactosidase assay; Immunostaining\tDrosophila\tPleura\tNULL\t26425438\tA novel link between FMR gene and the JNK pathway\t2015\t\t\t',
    # A pathway state named by regulator, no publication.
    'named state': 'Cisplatin and Trastuzumab\tNULL\tmiR-125-HER2 pathway\tOthers\tNULL\tmicroRNA-125(Activating);HER2(Inhibiting)\tSmall cell lung cancer\tH00013\tinhibit growth\t"Treatment of the HER2-positive SCLC cells."\tWestern blot; CCK-8 assay\tHomo sapiens\tLung\tH69;NCI-H82\tNULL\tChemotherapy regulated microRNA-125/HER2 pathway\t2015\t\t\t',
}


def opener(names=None, header=HEADER):
    lines = [ENTRIES[name] for name in (ENTRIES if names is None else names)]
    return SimpleNamespace(result=io.StringIO('\r\n'.join([header, *lines]) + '\r\n', newline=''))


def entries(names=None, **kwargs):
    return list(iter_entries(opener(names, **kwargs)))


def extract(records):
    out = SilverExtractor('cpad', 'entries')
    for index, record in enumerate(records):
        out.process_record(cpad.resource.entries.mapper(record), record, f'fixture:{index}', index)
    return out


def values(relation, term):
    return [a['value'] for a in relation.annotations if a['term'] == term and a.get('scope') == 'relation']


def test_an_entry_gives_a_record_per_pair_it_relates():
    pairs = [r['pair'] for r in entries(['regulator'])]
    assert pairs == ['regulator_pathway', 'regulator_cancer', 'pathway_cancer']
    # Without a regulator the entry relates the pathway and the cancer.
    assert [r['pair'] for r in entries(['no regulator'])] == ['pathway_cancer']
    first = entries(['regulator'])[0]
    assert first['regulator'] == 'AEG-1' and first['cancer'] == 'Breast cancer'
    assert first['regulator_expression'] == 'Upregulation' and first['pathway_state'] == 'Activating'
    assert first['kegg_pathways'] == ['path:04310'] and first['pathway_categories'] == ['Wnt signaling pathway']
    assert first['disease_id'] == 'H00031' and first['outcome'] == 'promote the proliferation and invasion'
    assert first['methods'] == ['RT-qPCR', 'CCK-8 assay', 'Colony formation assay']
    assert first['cell_lines'] == ['MCF-7', 'MDA-MB-231', 'SK-BR-3'] and first['pmid'] == '28454408'
    # The Greek letters of the GBK file are read as such.
    assert 'β-catenin' in first['description']
    # A regulator of two molecules is one regulator, as the website shows it.
    assert entries(['two regulators'])[0]['regulator'] == 'AXIN and APC'


def test_parser_reads_quoted_cells_and_repairs_the_one_defect():
    defect = entries(['defect'])[0]
    assert defect['description'].startswith('BrMC-induced apoptotic cell death occurred mainly by the activation of Akt, which')
    assert defect['regulator_expression'] is None
    assert defect['kegg_pathways'] == ['path:04151', 'path:04068']
    named = entries(['named values'])[0]
    assert named['regulator_expression'] == 'menin(Upregulation);Dnmt1(Downregulation)'
    assert entries(['named state'])[0]['pathway_state'] == 'microRNA-125(Activating);HER2(Inhibiting)'
    assert entries(['named state'])[0]['pmid'] is None


def test_parser_rejects_unknown_values_a_changed_header_and_an_empty_table():
    def parse_line(old, new):
        text = HEADER + '\r\n' + ENTRIES['regulator'].replace(old, new) + '\r\n'
        return list(iter_entries(SimpleNamespace(result=io.StringIO(text, newline=''))))

    with pytest.raises(ValueError, match='unknown expression pattern'):
        parse_line('Upregulation', 'Up')
    with pytest.raises(ValueError, match='unknown pathway state'):
        parse_line('Activating', 'Promoting')
    with pytest.raises(ValueError, match='unexpected PubMed identifier'):
        parse_line('28454408', '028454408')
    with pytest.raises(ValueError, match='lacks expected columns'):
        entries(header=HEADER.replace('Pathway_Category', 'Category'))
    with pytest.raises(ValueError, match='has no rows'):
        entries([])


def test_pathways_and_cancers_are_listed_once():
    # The pathways, then the cancers, each once, in order of first appearance.
    assert list(iter_terms(opener(['regulator', 'two regulators', 'regulator']))) == [
        {'kind': 'pathway', 'name': 'Wnt/beta-catenin signaling pathway'},
        {'kind': 'pathway', 'name': 'Wnt signaling pathway'},
        {'kind': 'cancer', 'name': 'Breast cancer'},
        {'kind': 'cancer', 'name': 'Colorectal cancer'},
    ]
    assert [t['name'] for t in iter_terms(opener(['regulator', 'regulator'])) if t['kind'] == 'cancer'] == [
        'Breast cancer'
    ]
    assert set(cpad.resource.datasets()) == {'terms', 'entries'}
    pathway = cpad.resource.terms.mapper({'kind': 'pathway', 'name': 'Wnt signaling pathway'})
    cancer = cpad.resource.terms.mapper({'kind': 'cancer', 'name': 'Breast cancer'})
    assert pathway.type == 'pathway' and cancer.type == 'ontology_class'
    assert pathway.identifiers == [('cpad_pathway', 'Wnt signaling pathway'), (Namespace.NAME, 'Wnt signaling pathway')]
    assert cancer.identifiers == [('cpad_cancer', 'Breast cancer'), (Namespace.NAME, 'Breast cancer')]


def test_pairs_are_associations_without_a_direction():
    records = entries(['regulator'])
    out = extract(records)
    assert [r.predicate for r in out.relations] == ['associated_with'] * 3
    assert not any(a['term'].endswith('_qualifier') for r in out.relations for a in r.annotations)
    regulator_pathway, regulator_cancer, pathway_cancer = out.relations
    # The two observed states are attributes of the pairs they describe.
    assert values(regulator_pathway, 'cpad:regulator_expression') == ['Upregulation']
    assert values(regulator_pathway, 'cpad:pathway_state') == ['Activating']
    assert values(regulator_cancer, 'cpad:regulator_expression') == ['Upregulation']
    assert values(regulator_cancer, 'cpad:outcome') == ['promote the proliferation and invasion']
    assert values(regulator_cancer, 'cpad:pathway_state') == []
    assert values(pathway_cancer, 'cpad:pathway_state') == ['Activating']
    assert values(pathway_cancer, 'cpad:disease_id') == ['H00031']
    for relation in out.relations:
        assert values(relation, 'publications') == ['PMID:28454408']
        assert values(relation, 'cpad:method') == ['RT-qPCR', 'CCK-8 assay', 'Colony formation assay']
    types = {(e.entity_type, e.namespace) for e in out.entities.values()}
    assert types == {('named_thing', 'name'), ('pathway', 'cpad_pathway'), ('ontology_class', 'cpad_cancer')}


def test_objects_are_the_entities_of_the_term_datasets():
    out = extract(entries(['regulator']))
    pathway = cpad.resource.terms.mapper({'kind': 'pathway', 'name': 'Wnt/beta-catenin signaling pathway'})
    cancer = cpad.resource.terms.mapper({'kind': 'cancer', 'name': 'Breast cancer'})
    statement = cpad.resource.entries.mapper(entries(['regulator'])[2])
    assert statement.subject.type == pathway.type and statement.subject.identifiers[0] == pathway.identifiers[0]
    assert statement.object.type == cancer.type and statement.object.identifiers[0] == cancer.identifiers[0]
    assert len(out.relations) == 3


def test_unknown_term_kind_raises():
    with pytest.raises(ValueError, match='Unknown CPAD term kind'):
        cpad.resource.terms.mapper({'kind': 'regulator', 'name': 'AEG-1'})


def test_species_and_unknown_values():
    mapped = extract(entries(['regulator', 'genus']))
    taxa = [[a['value'] for a in r.annotations if a['term'] == 'in_taxon'] for r in mapped.relations]
    # Human for the first entry; no taxon for a genus.
    assert taxa[:3] == [['NCBITaxon:9606']] * 3 and all(t == [] for t in taxa[3:])
    with pytest.raises(ValueError, match='Unknown CPAD species'):
        cpad.resource.entries.mapper(entries(['regulator'])[0] | {'species': 'Sus scrofa'})
    with pytest.raises(ValueError, match='Unknown CPAD pair'):
        cpad.resource.entries.mapper(entries(['regulator'])[0] | {'pair': 'cancer_cancer'})
