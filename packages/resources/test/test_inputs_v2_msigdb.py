"""Source contract of the MSigDB gene set input."""

import io
import json
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor

from pypath.internals.cv_terms import LicenseCV, ResourceCv
from pypath.inputs_v2 import msigdb
from pypath.inputs_v2.msigdb import COLLECTION, COLLECTIONS, GENE_SET, resource, source_xref
from pypath.inputs_v2.parsers.msigdb import iter_genes, iter_terms

# Entries of the JSON files of release 2026.1.Hs (name, fields): the gene sets of
# the files read, copied from the data server. BIOCARTA is from a file that the
# module does not read.
H_NOTCH = (
    'HALLMARK_NOTCH_SIGNALING',
    {'collection': 'H',
     'systematicName': 'M5903',
     'pmid': '26771021',
     'exactSource': '',
     'externalDetailsURL': '',
     'msigdbURL': 'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_NOTCH_SIGNALING',
     'geneSymbols': ['APH1A', 'ARRB1', 'CCND1', 'CUL1', 'DLL1', 'DTX1', 'DTX2', 'DTX4', 'FBXW11',
                     'FZD1', 'FZD5', 'FZD7', 'HES1', 'HEYL', 'JAG1', 'KAT2A', 'LFNG', 'MAML2',
                     'NOTCH1', 'NOTCH2', 'NOTCH3', 'PPARD', 'PRKCA', 'PSEN2', 'PSENEN', 'RBX1',
                     'SAP30', 'SKP1', 'ST3GAL6', 'TCF7L2', 'WNT2', 'WNT5A'],
     'filteredBySimilarity': [],
     'externalNamesForSimilarTerms': []}
)

C2_CGP = (
    'SHARMA_ASTROCYTOMA_WITH_NF1_SYNDROM',
    {'collection': 'C2:CGP',
     'systematicName': 'M216',
     'pmid': '17283119',
     'exactSource': 'Fig. 2A',
     'externalDetailsURL': '',
     'msigdbURL': 'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/SHARMA_ASTROCYTOMA_WITH_NF1_SYNDROM',
     'geneSymbols': ['ADGRV1', 'CLEC3B', 'MAP3K7CL', 'SLC1A3', 'STN1'],
     'filteredBySimilarity': [],
     'externalNamesForSimilarTerms': []}
)

C2_CP_REACTOME = (
    'REACTOME_FBXW7_MUTANTS_AND_NOTCH1_IN_CANCER',
    {'collection': 'C2:CP:REACTOME',
     'systematicName': 'M27197',
     'pmid': '',
     'exactSource': 'R-HSA-2644605',
     'externalDetailsURL': 'https://www.reactome.org/content/detail/R-HSA-2644605',
     'msigdbURL': 'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/REACTOME_FBXW7_MUTANTS_AND_NOTCH1_IN_CANCER',
     'geneSymbols': ['CUL1', 'FBXW7', 'NOTCH1', 'RBX1', 'SKP1'],
     'filteredBySimilarity': [],
     'externalNamesForSimilarTerms': []}
)

C2_CP_WIKIPATHWAYS = (
    'WP_TCA_CYCLE_NUTRIENT_USE_AND_INVASIVENESS_OF_OVARIAN_CANCER',
    {'collection': 'C2:CP:WIKIPATHWAYS',
     'systematicName': 'M39621',
     'pmid': '',
     'exactSource': 'WP2868',
     'externalDetailsURL': 'https://www.wikipathways.org/pathways/WP2868',
     'msigdbURL': 'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/WP_TCA_CYCLE_NUTRIENT_USE_AND_INVASIVENESS_OF_OVARIAN_CANCER',
     'geneSymbols': ['EGFR', 'JAK1', 'MAPK1', 'MAPK3', 'STAT3'],
     'filteredBySimilarity': [],
     'externalNamesForSimilarTerms': []}
)

C2_CP_KEGG_MEDICUS = (
    'KEGG_MEDICUS_ENV_FACTOR_NNK_TO_DNA_ADDUCTS',
    {'collection': 'C2:CP:KEGG_MEDICUS',
     'systematicName': 'M47808',
     'pmid': '',
     'exactSource': 'N01381',
     'externalDetailsURL': 'https://www.kegg.jp/entry/N01381',
     'msigdbURL': 'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/KEGG_MEDICUS_ENV_FACTOR_NNK_TO_DNA_ADDUCTS',
     'geneSymbols': ['AKR1C2', 'CBR1', 'CYP2A13', 'CYP2A6', 'HSD11B1'],
     'filteredBySimilarity': ['N01382'],
     'externalNamesForSimilarTerms': ['NNK to DNA adducts']}
)

C3_TFT_GTRD = (
    'SOD1_TARGET_GENES',
    {'collection': 'C3:TFT:GTRD',
     'systematicName': 'M46044',
     'pmid': '30445619',
     'exactSource': 'P00441',
     'externalDetailsURL': 'http://gtrd.biouml.org/',
     'msigdbURL': 'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/SOD1_TARGET_GENES',
     'geneSymbols': ['GAS5', 'SNORA25', 'SNORA32', 'SNORD12', 'SNORD6'],
     'filteredBySimilarity': [],
     'externalNamesForSimilarTerms': []}
)

C5_GO_BP = (
    'GOBP_SINGLE_STRAND_BREAK_REPAIR',
    {'collection': 'C5:GO:BP',
     'systematicName': 'M22128',
     'pmid': '',
     'exactSource': 'GO:0000012',
     'externalDetailsURL': 'http://amigo.geneontology.org/amigo/term/GO:0000012',
     'msigdbURL': 'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_SINGLE_STRAND_BREAK_REPAIR',
     'geneSymbols': ['APLF', 'APTX', 'ERCC6', 'ERCC8', 'LIG4', 'NEIL3', 'SIRT1', 'TDP1', 'TERF2',
                     'TNP1', 'UNG', 'WHR1', 'XNDC1N', 'XRCC1'],
     'filteredBySimilarity': [],
     'externalNamesForSimilarTerms': []}
)

C5_GO_CC = (
    'GOCC_CUL7_RING_UBIQUITIN_LIGASE_COMPLEX',
    {'collection': 'C5:GO:CC',
     'systematicName': 'M50596',
     'pmid': '',
     'exactSource': 'GO:0031467',
     'externalDetailsURL': 'http://amigo.geneontology.org/amigo/term/GO:0031467',
     'msigdbURL': 'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOCC_CUL7_RING_UBIQUITIN_LIGASE_COMPLEX',
     'geneSymbols': ['CUL7', 'CUL9', 'FBXW8', 'RBX1', 'SKP1'],
     'filteredBySimilarity': [],
     'externalNamesForSimilarTerms': []}
)

C5_HPO = (
    'HP_2_4_FINGER_CUTANEOUS_SYNDACTYLY',
    {'collection': 'C5:HPO',
     'systematicName': 'M46181',
     'pmid': '',
     'exactSource': 'HP:0010709',
     'externalDetailsURL': 'https://hpo.jax.org/app/browse/term/HP:0010709',
     'msigdbURL': 'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HP_2_4_FINGER_CUTANEOUS_SYNDACTYLY',
     'geneSymbols': ['CDC45', 'CDH3', 'NECTIN1', 'NECTIN4', 'TWIST1'],
     'filteredBySimilarity': [],
     'externalNamesForSimilarTerms': []}
)

C1 = (
    'chr6q11',
    {'collection': 'C1',
     'systematicName': 'M22108',
     'pmid': '',
     'exactSource': 'Ensembl 115 human gene annotation on the GRCh38.p14 reference chromosomes '
                    'only.',
     'externalDetailsURL': 'https://ensembl.org/Homo_sapiens/Location/Chromosome?r=6',
     'msigdbURL': 'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/chr6q11',
     'geneSymbols': ['DHFRP5', 'GAPDHP41', 'KHDRBS2', 'PRIM2BP', 'RBBP4P3'],
     'filteredBySimilarity': [],
     'externalNamesForSimilarTerms': []}
)

BIOCARTA = (
    'BIOCARTA_BOTULIN_PATHWAY',
    {'collection': 'C2:CP:BIOCARTA',
     'systematicName': 'M22086',
     'pmid': '',
     'exactSource': '',
     'externalDetailsURL': 'https://data.broadinstitute.org/gsea-msigdb/msigdb/biocarta/human/h_botulinPathway.gif',
     'msigdbURL': 'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/BIOCARTA_BOTULIN_PATHWAY',
     'geneSymbols': ['CHRM1', 'CHRNA1', 'SNAP25', 'STX1A', 'VAMP2'],
     'filteredBySimilarity': [],
     'externalNamesForSimilarTerms': []}
)


class Json:
    """A stand-in for the Download of one JSON file: the sets of one collection."""

    def __init__(self, *entries, text=None):
        self.text = text if text is not None else json.dumps(dict(entries))

    def open(self, force_refresh=False, **_kwargs):
        return SimpleNamespace(result=io.StringIO(self.text))


def collections(**files):
    """The real collection table with the given JSON files in place of the downloads."""
    result = {code: dict(collection, download=None) for code, collection in COLLECTIONS.items()}
    for code, file in files.items():
        result[code.replace('__', ':')]['download'] = file
    return result


SAMPLE = collections(
    H=Json(H_NOTCH),
    C2__CGP=Json(C2_CGP),
    C2__CP__REACTOME=Json(C2_CP_REACTOME),
    C2__CP__WIKIPATHWAYS=Json(C2_CP_WIKIPATHWAYS),
    C2__CP__KEGG_MEDICUS=Json(C2_CP_KEGG_MEDICUS),
    C3__TFT__GTRD=Json(C3_TFT_GTRD),
    C5__GO__BP=Json(C5_GO_BP),
    C5__GO__CC=Json(C5_GO_CC),
    C5__HPO=Json(C5_HPO),
    C1=Json(C1),
)
ENTRIES = [
    H_NOTCH, C2_CGP, C2_CP_REACTOME, C2_CP_WIKIPATHWAYS, C2_CP_KEGG_MEDICUS, C3_TFT_GTRD,
    C5_GO_BP, C5_GO_CC, C5_HPO, C1,
]
NAMES = {name for name, _ in ENTRIES}


def terms(sample=SAMPLE):
    return list(iter_terms(collections=sample))


def term(name, sample=SAMPLE):
    return next(r for r in terms(sample) if r['id'] == name)


def genes(sample=SAMPLE):
    return {row['gene']: row['gene_sets'] for row in iter_genes(collections=sample)}


def extract(row, dataset):
    record = resource.datasets()[dataset].mapper(row)
    assert record is not None
    out = SilverExtractor('msigdb', dataset)
    out.process_record(record, row, 'fixture:1', 1)
    return out


def edges(out):
    return sorted(
        (
            out.entities[r.subject_entity_key].identifier,
            r.predicate,
            out.entities[r.object_entity_key].namespace,
            out.entities[r.object_entity_key].identifier,
        )
        for r in out.relations
    )


def gene_set(name):
    """The entity of a gene set and the extractor that produced it."""
    out = extract(term(name), 'terms')
    (entity,) = (e for e in out.entities.values() if e.label)
    return entity, out


def annotation_values(entity, name):
    return [a['value'] for a in entity.annotations if a['term'] == name]


def identifiers(entity):
    return [(i['ns'], i['id']) for i in entity.identifiers]


def test_collections_come_first_and_gene_sets_follow_with_their_collection():
    rows = terms()
    kinds = [r['kind'] for r in rows]
    assert kinds == ['collection'] * len(COLLECTIONS) + ['gene_set'] * 10
    assert term('HALLMARK_NOTCH_SIGNALING') == {
        'kind': 'gene_set',
        'id': 'HALLMARK_NOTCH_SIGNALING',
        'systematic_name': 'M5903',
        'url': 'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_NOTCH_SIGNALING',
        'pmid': '26771021',
        'exact_source': None,
        'collection': 'H',
    }
    assert {r['id']: r['collection'] for r in rows if r['kind'] == 'gene_set'} == {
        'HALLMARK_NOTCH_SIGNALING': 'H',
        'SHARMA_ASTROCYTOMA_WITH_NF1_SYNDROM': 'C2:CGP',
        'REACTOME_FBXW7_MUTANTS_AND_NOTCH1_IN_CANCER': 'C2:CP:REACTOME',
        'WP_TCA_CYCLE_NUTRIENT_USE_AND_INVASIVENESS_OF_OVARIAN_CANCER': 'C2:CP:WIKIPATHWAYS',
        'KEGG_MEDICUS_ENV_FACTOR_NNK_TO_DNA_ADDUCTS': 'C2:CP:KEGG_MEDICUS',
        'SOD1_TARGET_GENES': 'C3:TFT:GTRD',
        'GOBP_SINGLE_STRAND_BREAK_REPAIR': 'C5:GO:BP',
        'GOCC_CUL7_RING_UBIQUITIN_LIGASE_COMPLEX': 'C5:GO:CC',
        'HP_2_4_FINGER_CUTANEOUS_SYNDACTYLY': 'C5:HPO',
        'chr6q11': 'C1',
    }


def test_gene_set_is_a_source_term_below_its_collection():
    entity, out = gene_set('REACTOME_FBXW7_MUTANTS_AND_NOTCH1_IN_CANCER')
    assert (entity.entity_type, entity.namespace, entity.identifier, entity.label) == (
        'ontology_class', GENE_SET,
        'REACTOME_FBXW7_MUTANTS_AND_NOTCH1_IN_CANCER', 'REACTOME_FBXW7_MUTANTS_AND_NOTCH1_IN_CANCER',
    )
    assert edges(out) == [
        ('REACTOME_FBXW7_MUTANTS_AND_NOTCH1_IN_CANCER', 'subclass_of', COLLECTION, 'C2:CP:REACTOME')
    ]
    assert annotation_values(entity, 'msigdb:url') == [
        'https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/REACTOME_FBXW7_MUTANTS_AND_NOTCH1_IN_CANCER'
    ]


def test_collections_are_nested_as_the_site_nests_them():
    by_code = {r['id']: r for r in terms() if r['kind'] == 'collection'}
    assert by_code['C2:CP:REACTOME']['parent'] == 'C2:CP'
    assert by_code['C2:CP']['parent'] == 'C2'
    assert by_code['C2']['parent'] is None
    assert by_code['C5:GO:CC']['parent'] == 'C5:GO'
    assert by_code['C3:MIR:MIRDB']['parent'] == 'C3:MIR'
    out = extract(by_code['C2:CP:REACTOME'], 'terms')
    assert edges(out) == [('C2:CP:REACTOME', 'subclass_of', COLLECTION, 'C2:CP')]
    (entity,) = (e for e in out.entities.values() if e.label)
    assert entity.label == 'Reactome subset of CP'
    assert extract(by_code['C2'], 'terms').relations == []


def test_the_systematic_name_is_a_second_identifier_and_the_first_is_unchanged():
    # The first identifier is the identity of the term and the object of every
    # membership; the accession of MSigDB is a further identifier of the same set.
    for name, fields in ENTRIES:
        entity, out = gene_set(name)
        assert identifiers(entity)[0] == (GENE_SET, name), name
        assert sorted(identifiers(entity)) == sorted(
            [(GENE_SET, name), ('name', name), ('msigdb', fields['systematicName'])]
        ), name
        assert (entity.namespace, entity.identifier) == (GENE_SET, name)
        assert [i['is_canonical'] for i in entity.identifiers] == [True, False, False]
        assert edges(out) == [(name, 'subclass_of', COLLECTION, fields['collection'])], name


def test_a_go_set_has_its_go_term_as_a_cross_reference_and_no_relation_to_it():
    entity, out = gene_set('GOBP_SINGLE_STRAND_BREAK_REPAIR')
    assert identifiers(entity)[0] == (GENE_SET, entity.identifier)
    assert ('msigdb', 'M22128') in identifiers(entity)
    assert annotation_values(entity, 'xref') == ['GO:0000012']
    assert annotation_values(entity, 'publications') == []
    # Only the collection is related to the set; the GO term is not an object.
    assert {r.predicate for r in out.relations} == {'subclass_of'}
    assert {out.entities[r.object_entity_key].namespace for r in out.relations} == {COLLECTION}
    assert len(out.entities) == 2
    entity, _ = gene_set('GOCC_CUL7_RING_UBIQUITIN_LIGASE_COMPLEX')
    assert annotation_values(entity, 'xref') == ['GO:0031467']


def test_an_hpo_set_has_its_hp_term_as_a_cross_reference():
    entity, _ = gene_set('HP_2_4_FINGER_CUTANEOUS_SYNDACTYLY')
    assert identifiers(entity)[0] == (GENE_SET, entity.identifier)
    assert ('msigdb', 'M46181') in identifiers(entity)
    assert annotation_values(entity, 'xref') == ['HP:0010709']
    assert annotation_values(entity, 'publications') == []


def test_a_reactome_set_has_the_stable_id_as_a_cross_reference():
    entity, _ = gene_set('REACTOME_FBXW7_MUTANTS_AND_NOTCH1_IN_CANCER')
    assert identifiers(entity)[0] == (GENE_SET, entity.identifier)
    assert ('msigdb', 'M27197') in identifiers(entity)
    assert annotation_values(entity, 'xref') == ['R-HSA-2644605']


def test_a_wikipathways_set_has_the_pathway_id_as_a_cross_reference():
    entity, _ = gene_set('WP_TCA_CYCLE_NUTRIENT_USE_AND_INVASIVENESS_OF_OVARIAN_CANCER')
    assert identifiers(entity)[0] == (GENE_SET, entity.identifier)
    assert ('msigdb', 'M39621') in identifiers(entity)
    assert annotation_values(entity, 'xref') == ['WP2868']


def test_a_cgp_set_has_its_publication_and_the_table_of_the_paper_is_no_cross_reference():
    entity, _ = gene_set('SHARMA_ASTROCYTOMA_WITH_NF1_SYNDROM')
    assert identifiers(entity)[0] == (GENE_SET, entity.identifier)
    assert ('msigdb', 'M216') in identifiers(entity)
    assert annotation_values(entity, 'publications') == ['PMID:17283119']
    # "Fig. 2A" is the figure of the publication the genes come from.
    assert term('SHARMA_ASTROCYTOMA_WITH_NF1_SYNDROM')['exact_source'] == 'Fig. 2A'
    assert annotation_values(entity, 'xref') == []


def test_a_hallmark_set_has_its_publication_and_no_exact_source():
    entity, _ = gene_set('HALLMARK_NOTCH_SIGNALING')
    assert identifiers(entity)[0] == (GENE_SET, entity.identifier)
    assert ('msigdb', 'M5903') in identifiers(entity)
    assert annotation_values(entity, 'publications') == ['PMID:26771021']
    assert annotation_values(entity, 'xref') == []
    assert term('HALLMARK_NOTCH_SIGNALING')['exact_source'] is None


def test_an_exact_source_that_is_not_an_accession_of_these_vocabularies_is_left_out():
    # The exact source of a GTRD set is the UniProt accession of the factor, of a
    # KEGG MEDICUS set the KEGG network, of a positional set a description.
    for name, exact in (
        ('SOD1_TARGET_GENES', 'P00441'),
        ('KEGG_MEDICUS_ENV_FACTOR_NNK_TO_DNA_ADDUCTS', 'N01381'),
        ('chr6q11', 'Ensembl 115 human gene annotation on the GRCh38.p14 reference chromosomes only.'),
    ):
        assert term(name)['exact_source'] == exact
        entity, _ = gene_set(name)
        assert annotation_values(entity, 'xref') == [], name
    entity, _ = gene_set('SOD1_TARGET_GENES')
    assert annotation_values(entity, 'publications') == ['PMID:30445619']


def test_a_cross_reference_needs_the_exact_form_of_an_accession():
    assert source_xref('GO:0006915') == 'GO:0006915'
    assert source_xref('HP:0001507') == 'HP:0001507'
    assert source_xref('R-HSA-109581') == 'R-HSA-109581'
    assert source_xref('WP100') == 'WP100'
    for value in (
        None, '', ' GO:0006915', 'GO:0006915 ', 'GO:006915', 'GO:00069150', 'go:0006915', 'GO0006915',
        'HP:001507', 'HP:0001507; GO:0006915', 'HP:00010709', 'R-HSA-', 'R-MMU-109581', 'R-HSA-0',
        'WP', 'WP0', 'WP100_r1', 'wp100', 'N00001', 'P00441', 'Table 1', 'module_1',
    ):
        assert source_xref(value) is None, value


def test_a_publication_needs_a_pubmed_identifier():
    for pmid, expected in (('26771021', ['PMID:26771021']), (None, [])):
        row = dict(term('HALLMARK_NOTCH_SIGNALING'), pmid=pmid)
        out = extract(row, 'terms')
        (entity,) = (e for e in out.entities.values() if e.label)
        assert annotation_values(entity, 'publications') == expected


def test_gene_is_associated_with_every_set_it_is_in():
    members = genes()
    assert members['NOTCH1'] == [
        'HALLMARK_NOTCH_SIGNALING', 'REACTOME_FBXW7_MUTANTS_AND_NOTCH1_IN_CANCER',
    ]
    assert members['RBX1'] == [
        'HALLMARK_NOTCH_SIGNALING',
        'REACTOME_FBXW7_MUTANTS_AND_NOTCH1_IN_CANCER',
        'GOCC_CUL7_RING_UBIQUITIN_LIGASE_COMPLEX',
    ]
    assert members['TWIST1'] == ['HP_2_4_FINGER_CUTANEOUS_SYNDACTYLY']
    row = {'gene': 'RBX1', 'gene_sets': members['RBX1']}
    out = extract(row, 'genes')
    assert edges(out) == [
        ('RBX1', 'associated_with', GENE_SET, 'GOCC_CUL7_RING_UBIQUITIN_LIGASE_COMPLEX'),
        ('RBX1', 'associated_with', GENE_SET, 'HALLMARK_NOTCH_SIGNALING'),
        ('RBX1', 'associated_with', GENE_SET, 'REACTOME_FBXW7_MUTANTS_AND_NOTCH1_IN_CANCER'),
    ]
    (gene,) = (e for e in out.entities.values() if e.entity_type == 'gene')
    assert (gene.namespace, gene.identifier, gene.taxon) == ('genesymbol', 'RBX1', '9606')


def test_membership_statements_equal_the_symbols_of_the_sets():
    total = sum(len(fields['geneSymbols']) for _, fields in ENTRIES)
    assert sum(len(sets) for sets in genes().values()) == total
    # A gene is one record however many sets name it (CUL1 and NOTCH1 are in two
    # of the sets of the sample, RBX1 and SKP1 in three).
    assert len(genes()) == len({g for _, fields in ENTRIES for g in fields['geneSymbols']})
    assert len(genes()) < total


def test_statement_object_is_the_entity_of_the_terms_dataset():
    out = extract({'gene': 'TWIST1', 'gene_sets': genes()['TWIST1']}, 'genes')
    (relation,) = out.relations
    target = out.entities[relation.object_entity_key]
    entity, _ = gene_set('HP_2_4_FINGER_CUTANEOUS_SYNDACTYLY')
    assert (target.entity_type, target.namespace, target.identifier) == (
        entity.entity_type, entity.namespace, entity.identifier,
    )


def test_the_second_identifier_changes_neither_the_identity_nor_the_statement_objects():
    # Before the systematic name was added the term had the identifiers
    # (gene set, name); the statement object of a membership is built from the
    # gene set name alone. Both are the same as with the third identifier.
    row = term('GOBP_SINGLE_STRAND_BREAK_REPAIR')
    record = resource.datasets()['terms'].mapper(row)
    older = record._replace(identifiers=[i for i in record.identifiers if i.type != 'msigdb'])
    before = SilverExtractor('msigdb', 'terms')
    before.process_record(older, row, 'fixture:1', 1)
    after = extract(row, 'terms')
    by_name = lambda out: next(e for e in out.entities.values() if e.identifier == row['id'])
    assert (by_name(before).entity_type, by_name(before).namespace, by_name(before).identifier) == (
        by_name(after).entity_type, by_name(after).namespace, by_name(after).identifier,
    )
    collection_keys = lambda out: sorted(
        out.entities[r.object_entity_key].identifier for r in out.relations
    )
    assert collection_keys(before) == collection_keys(after) == ['C5:GO:BP']
    members = extract({'gene': 'XRCC1', 'gene_sets': ['GOBP_SINGLE_STRAND_BREAK_REPAIR']}, 'genes')
    (relation,) = members.relations
    assert members.entities[relation.object_entity_key].identifier == row['id']


def test_gene_symbols_are_kept_as_written():
    # DHFRP5 is a pseudogene symbol of the positional set chr6q11; the symbol is not
    # translated or checked against a gene list.
    assert 'DHFRP5' in genes()


def test_a_membership_row_has_only_the_gene_and_its_sets():
    assert all(set(row) == {'gene', 'gene_sets'} for row in iter_genes(collections=SAMPLE))


def test_biocarta_and_kegg_legacy_are_left_out_for_their_license():
    assert 'C2:CP:BIOCARTA' not in COLLECTIONS
    assert 'C2:CP:KEGG_LEGACY' not in COLLECTIONS
    assert not any(
        'biocarta' in str(c['download'].url) or 'kegg_legacy' in str(c['download'].url)
        for c in COLLECTIONS.values() if c['download']
    )
    # Their sets are not read, because their files are not: a BioCarta entry stays an
    # entry of a file that no collection names.
    assert BIOCARTA[1]['collection'] == 'C2:CP:BIOCARTA'
    assert BIOCARTA[0] not in {r['id'] for r in terms()}


def test_aggregated_files_are_not_read():
    stems = {c['download'].filename for c in COLLECTIONS.values() if c['download']}
    assert not any(s.startswith(('c2.all', 'c2.cp.v', 'c3.all', 'c3.mir.v', 'c3.tft.v', 'c4.all',
                                 'c5.all', 'c5.go.v', 'c7.all', 'msigdb.')) for s in stems)
    assert len(stems) == 23


def test_the_json_files_of_the_release_are_read_and_not_the_gmt_files():
    for collection in COLLECTIONS.values():
        download = collection['download']
        if download is None:
            continue
        assert download.url.startswith(
            'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/'
        )
        assert download.url.endswith('.v2026.1.Hs.json')
        assert download.filename.endswith('.v2026.1.Hs.json')


def test_every_parent_exists_and_every_leaf_has_a_file():
    parents = {c['parent'] for c in COLLECTIONS.values() if c['parent']}
    assert parents <= set(COLLECTIONS)
    for code, collection in COLLECTIONS.items():
        assert collection['download'] is not None or code in parents, code


def test_the_sets_of_a_file_are_read_in_the_order_of_their_names():
    # A JSON object has no order; the GMT files are in the order of the names.
    reversed_file = Json(H_NOTCH, ('HALLMARK_APOPTOSIS', dict(H_NOTCH[1], systematicName='M5902')))
    sample = collections(H=reversed_file)
    assert [r['id'] for r in iter_terms(collections=sample) if r['kind'] == 'gene_set'] == [
        'HALLMARK_APOPTOSIS', 'HALLMARK_NOTCH_SIGNALING',
    ]
    assert genes(sample)['APH1A'] == ['HALLMARK_APOPTOSIS', 'HALLMARK_NOTCH_SIGNALING']


def test_a_gene_set_in_two_files_raises():
    sample = collections(H=Json(H_NOTCH), C1=Json((H_NOTCH[0], dict(H_NOTCH[1], collection='C1'))))
    with pytest.raises(ValueError, match='is in the collections'):
        terms(sample)


def test_a_systematic_name_of_two_sets_raises():
    sample = collections(H=Json(H_NOTCH), C1=Json((C1[0], dict(C1[1], systematicName='M5903'))))
    with pytest.raises(ValueError, match='M5903 is of the gene sets'):
        terms(sample)


def test_a_set_in_the_file_of_another_collection_raises():
    with pytest.raises(ValueError, match='names the collection'):
        terms(collections(C1=Json(H_NOTCH)))


def test_a_malformed_set_raises():
    for change, message in (
        ({'systematicName': '5903'}, 'malformed systematic name'),
        ({'systematicName': ''}, 'malformed systematic name'),
        ({'pmid': '26771021, 25000000'}, 'malformed PubMed id'),
        ({'pmid': '026771021'}, 'malformed PubMed id'),
        ({'geneSymbols': 'NOTCH1'}, 'no list of gene symbols'),
    ):
        entry = (H_NOTCH[0], dict(H_NOTCH[1], **change))
        with pytest.raises(ValueError, match=message):
            terms(collections(H=Json(entry)))
    for field in ('systematicName', 'pmid', 'exactSource', 'geneSymbols', 'msigdbURL', 'collection'):
        fields = {k: v for k, v in H_NOTCH[1].items() if k != field}
        with pytest.raises(ValueError, match='lacks one of'):
            terms(collections(H=Json((H_NOTCH[0], fields))))


def test_a_gene_listed_twice_in_a_set_raises():
    doubled = (H_NOTCH[0], dict(H_NOTCH[1], geneSymbols=H_NOTCH[1]['geneSymbols'] + ['NOTCH1']))
    with pytest.raises(ValueError, match='lists a gene twice'):
        genes(collections(H=Json(doubled)))


def test_an_empty_file_and_an_html_page_raise():
    with pytest.raises(ValueError, match='no gene sets'):
        terms(collections(H=Json()))
    with pytest.raises(ValueError, match='expected a JSON file'):
        terms(collections(H=Json(text='<!DOCTYPE html>\n<html>')))
    with pytest.raises(ValueError, match='one JSON object'):
        terms(collections(H=Json(text='[]')))


def test_unknown_term_kind_and_collection_raise():
    with pytest.raises(ValueError, match='term kind'):
        msigdb.map_term({'kind': 'tissue', 'id': 'x'})
    with pytest.raises(ValueError, match='Unknown MSigDB collection'):
        msigdb.map_term({'kind': 'gene_set', 'id': 'x', 'url': None, 'collection': 'C99'})


def test_no_annotation_other_than_the_taxon_on_a_membership():
    out = extract({'gene': 'TWIST1', 'gene_sets': ['HP_2_4_FINGER_CUTANEOUS_SYNDACTYLY']}, 'genes')
    assert {a['term'] for r in out.relations for a in r.annotations} <= {'in_taxon'}


def test_resource_declares_its_datasets_and_license():
    assert msigdb.config.id == ResourceCv.MSIGDB
    assert msigdb.config.license == LicenseCV.COMPOSITE
    assert set(resource.datasets()) == {'terms', 'genes'}
