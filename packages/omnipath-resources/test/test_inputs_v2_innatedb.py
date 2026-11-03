"""InnateDB curated interactions: MITAB helpers, mapper and extractor on source lines."""

import io
from types import SimpleNamespace

import pytest
from biolink_model.datamodel.model import slots
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import _mitab, innatedb

HEADER = '\t'.join([
    '#unique_identifier_A', 'unique_identifier_B', 'alt_identifier_A',
    'alt_identifier_B', 'alias_A', 'alias_B', 'interaction_detection_method',
    'author', 'pmid', 'ncbi_taxid_A', 'ncbi_taxid_B', 'interaction_type',
    'source_database', 'idinteraction_in_source_db', 'confidence_score',
    'expansion_method', 'biological_role_A', 'biological_role_B',
    'exp_role_A', 'exp_role_B', 'interactor_type_A', 'interactor_type_B',
    'xrefs_A', 'xrefs_B', 'xrefs_interaction', 'annotations_A',
    'annotations_B', 'annotations_interaction', 'ncbi_taxid_host_organism',
    'parameters_interaction', 'creation_date', 'update_date', 'checksum_A',
    'checksum_B', 'checksum_interaction', 'negative', 'features_A',
    'features_B', 'stoichiometry_A', 'stoichiometry_B',
    'participant_identification_method_A',
    'participant_identification_method_B',
]) + '\n'

UNSPECIFIED = 'psi-mi:"MI:0499"(unspecified role)'
NEUTRAL = 'psi-mi:"MI:0497"(neutral component)'
BAIT = 'psi-mi:"MI:0496"(bait)'
PREY = 'psi-mi:"MI:0498"(prey)'
PROTEIN = 'psi-mi:"MI:0326"(protein)'
INFERRED = 'psi-mi:"MI:0363"(inferred by author)'
COIP = 'psi-mi:"MI:0007"(anti tag coimmunoprecipitation)'
PHYSICAL = 'psi-mi:"MI:0915"(physical association)'


def line(
    gene_a, gene_b, ensembl_a, ensembl_b, alias_a, alias_b, method, author,
    pmid, interaction, host, taxon_a='9606(Human)', taxon_b='9606(Human)',
    interaction_type=PHYSICAL, source='MI:0974(innatedb)', expansion='-',
    biological=(UNSPECIFIED, UNSPECIFIED), experimental=(NEUTRAL, NEUTRAL),
    types=(PROTEIN, PROTEIN), comment='-', negative='false',
    identification=(INFERRED, INFERRED),
):
    """A line of the file; the columns without arguments are always ``-``."""
    return '\t'.join([
        f'innatedb:{gene_a}', f'innatedb:{gene_b}',
        f'ensembl:{ensembl_a}', f'ensembl:{ensembl_b}', alias_a, alias_b,
        method, author, f'pubmed:{pmid}', f'taxid:{taxon_a}',
        f'taxid:{taxon_b}', interaction_type, source, interaction,
        'lpr:3|hpr:3|np:1|', expansion, *biological, *experimental, *types,
        '-', '-', '-', '-', '-', comment, f'taxid:{host}', '-', '2008/03/30',
        '2008/03/30', '-', '-', '-', negative, '-', '-', '-', '-',
        *identification,
    ]) + '\n'


LY96 = (
    'uniprotkb:LY96_HUMAN|refseq:NP_056179|uniprotkb:Q9Y6Y9|'
    'refseq:NP_001182726|hgnc:LY96(display_short)'
)
TLR4 = (
    'refseq:NP_612564|refseq:NP_612567|uniprotkb:O00206|'
    'uniprotkb:TLR4_HUMAN|refseq:NP_003257|hgnc:TLR4(display_short)'
)
IRAK1 = (
    'uniprotkb:IRAK1_HUMAN|refseq:NP_001560|uniprotkb:P51617|'
    'refseq:NP_001020414|refseq:NP_001020413|hgnc:IRAK1(display_short)'
)
# innatedb_all.mitab.gz line 2: human proteins tested in mouse cells.
PHYSICAL_LINE = line(
    'IDBG-25842', 'IDBG-82738', 'ENSG00000154589', 'ENSG00000136869',
    LY96, TLR4, COIP, 'Shimazu et al. (1999)', '10359581',
    'innatedb:IDB-113240', '10090', experimental=(PREY, BAIT),
)
# Line 243: IRAK1 phosphorylates itself; the roles name enzyme and target.
PHOSPHORYLATION_LINE = line(
    'IDBG-90782', 'IDBG-90782', 'ENSG00000184216', 'ENSG00000184216',
    IRAK1, IRAK1, 'psi-mi:"MI:0423"(in-gel kinase assay)',
    'Zhang et al. (2002)', '11751856', 'innatedb:IDB-113479', '9606',
    interaction_type='psi-mi:"MI:0217"(phosphorylation reaction)',
    biological=(
        'psi-mi:"MI:0501"(enzyme)', 'psi-mi:"MI:0502"(enzyme target)',
    ),
)
# Line 980: IKBKB protein on the DNA of the CXCL8 gene; no host organism.
DNA_LINE = line(
    'IDBG-19987', 'IDBG-23954', 'ENSG00000104365', 'ENSG00000169429',
    'uniprotkb:IKKB_HUMAN|refseq:NP_001229707|uniprotkb:O14920|'
    'refseq:NP_001547|refseq:NP_001177649|hgnc:IKBKB(display_short)',
    'uniprotkb:IL8_HUMAN|refseq:NP_000575|uniprotkb:P10145|'
    'hgnc:CXCL8(display_short)',
    COIP, 'Hoberg et al. (2004)', '15494311', 'innatedb:IDB-114221', '0',
    types=(PROTEIN, 'psi-mi:"MI:0319"(deoxyribonucleic acid)'),
)
# Line 8406: a mouse protein and the RNA of a mouse gene.
RNA_LINE = line(
    'IDBG-193424', 'IDBG-205243', 'ENSMUSG00000007670', 'ENSMUSG00000006932',
    'uniprotkb:FUBP2_MOUSE|refseq:NP_034743|uniprotkb:Q3U0V1|'
    'mgi:Khsrp(display_short)',
    'uniprotkb:CTNB1_MOUSE|refseq:NP_031640|uniprotkb:Q02248|'
    'refseq:NP_001159374|mgi:Ctnnb1(display_short)',
    COIP, 'Gherzi et al. (2006)', '17177604', 'innatedb:IDB-155988',
    '10090', taxon_a='10090(Mouse)', taxon_b='10090(Mouse)',
    experimental=(BAIT, PREY),
    types=(PROTEIN, 'psi-mi:"MI:0320"(ribonucleic acid)'),
    identification=(
        'psi-mi:"MI:0421"(identification by antibody)',
        'psi-mi:"MI:0088"(primer specific pcr)',
    ),
)
# Line 19214: a record of the allergy and asthma curation project.
ALLERGY_LINE = line(
    'IDBG-26087', 'IDBG-65194', 'ENSG00000127445', 'ENSG00000154229',
    'refseq:NP_006212|uniprotkb:PIN1_HUMAN|uniprotkb:Q13526|'
    'hgnc:PIN1(display_short)',
    'uniprotkb:KPCA_HUMAN|refseq:NP_002728|uniprotkb:P17252|'
    'hgnc:PRKCA(display_short)',
    COIP, 'Shen et al. (2008)', '18188456', 'InnateDB Allergy:IDB-648953',
    '9606', experimental=(BAIT, PREY),
)
# Line 1349: the interaction type column holds a detection method.
UNTYPED_LINE = line(
    'IDBG-46584', 'IDBG-85142', 'ENSG00000126247', 'ENSG00000101966',
    'uniprotkb:CPNS1_HUMAN|refseq:NP_001003962|refseq:NP_001740|'
    'uniprotkb:P04632|hgnc:CAPNS1(display_short)',
    'uniprotkb:P98170|uniprotkb:XIAP_HUMAN|refseq:NP_001158|'
    'refseq:NP_001191330|hgnc:XIAP(display_short)',
    'psi-mi:"MI:0435"(protease assay)', 'Kobayashi et al. (2002)',
    '12121983', 'innatedb:IDB-114606', '0',
    interaction_type='psi-mi:"MI:0401"(biochemical)',
)
# Line 13721: the second gene has no symbol and no protein.
PLACEHOLDER_LINE = line(
    'IDBG-14629', 'IDBG-303636', 'ENSG00000108953', 'ENSG00000237708',
    'uniprotkb:1433E_HUMAN|refseq:NP_006752|uniprotkb:P62258|'
    'hgnc:YWHAE(display_short)',
    'innatedb:IDBG=303636(display_short)',
    COIP, 'Zuo et al. (2010)', '20462248', 'innatedb:IDB-188540', '9606',
    experimental=(BAIT, PREY),
)
# all.mitab.gz line 5: a BioGRID record in InnateDB's all-interactions file.
IMPORTED_LINE = line(
    'IDBG-96746', 'IDBG-73725', 'ENSG00000084072', 'ENSG00000118058',
    'refseq:NP_982281|uniprotkb:PPIE_HUMAN|uniprotkb:Q9UNP9|'
    'refseq:NP_006103|refseq:NP_001181936|hgnc:PPIE(display_short)',
    'uniprotkb:MLL1_HUMAN|refseq:NP_005924|uniprotkb:Q03164|'
    'refseq:NP_001184033|hgnc:KMT2A(display_short)',
    'psi-mi:"MI:0018"(two hybrid)', 'Fair et al.(2001)', '11313484',
    'BIOGRID:10022', '32644', source='MI:0463(biogrid)',
    experimental=(PREY, BAIT),
    comment='comment:"BioGRID Evidence Code: Two-hybrid; "',
)


def read(*lines):
    opener = SimpleNamespace(result=io.StringIO(HEADER + ''.join(lines)))
    return list(innatedb.resource.datasets()['interactions']._raw_parser(opener))


def extract(rows):
    out = SilverExtractor('innatedb', 'interactions')
    for index, row in enumerate(rows):
        record = innatedb.resource.interactions.mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def identifiers(entity):
    return [(i.type, i.value) for i in entity.identifiers]


def values(relation, term, scope='relation'):
    return [
        a['value'] for a in relation.annotations
        if a['term'] == term and a['scope'] == scope
    ]


# Shared MITAB helpers ----------------------------------------------------


def test_mitab_item_has_database_accession_and_description():
    item = _mitab.parse_mitab_item('psi-mi:"MI:0018"(two hybrid)')
    assert item == ('psi-mi', 'MI:0018', 'two hybrid')
    assert _mitab.parse_mitab_item('taxid:9606(human-293f)') == (
        'taxid', '9606', 'human-293f',
    )
    assert _mitab.parse_mitab_item('uniprotkb:P01730') == (
        'uniprotkb', 'P01730', None,
    )
    # The database name is case-insensitive and may contain spaces.
    assert _mitab.parse_mitab_item('InnateDB Allergy:IDB-648953') == (
        'innatedb allergy', 'IDB-648953', None,
    )
    assert _mitab.parse_mitab_item('-') is None
    assert _mitab.parse_mitab_item('') is None


def test_mitab_item_keeps_parentheses_and_quotes_that_belong_to_values():
    # The description is the parenthesis that ends the item.
    assert _mitab.parse_mitab_item(
        'entrez gene/locuslink:p38(MAPK)(gene name synonym)'
    ) == ('entrez gene/locuslink', 'p38(MAPK)', 'gene name synonym')
    assert _mitab.parse_mitab_item(
        'go:"GO:0035401"("histone kinase activity (H3-Y41 specific)")'
    ) == ('go', 'GO:0035401', 'histone kinase activity (H3-Y41 specific)')
    # DIP writes its own identifiers without a database name.
    assert _mitab.parse_mitab_item('DIP-617N') == (None, 'DIP-617N', None)
    # BioGRID quotes the first author, whose year is not a description.
    assert _mitab.parse_mitab_item('"Marti A (1997)"') == (
        None, 'Marti A (1997)', None,
    )
    assert _mitab.parse_mitab_item('comment:"vector: pDEST | 2u"') == (
        'comment', 'vector: pDEST | 2u', None,
    )


def test_mitab_cell_splits_outside_quotes_only():
    cell = 'comment:"a|b"|pubmed:1|-'
    assert [item.accession for item in _mitab.parse_mitab_cell(cell)] == [
        'a|b', '1',
    ]
    assert _mitab.parse_mitab_cell('-') == []
    assert _mitab.parse_mitab_cell(None) == []


def test_mitab_terms_publications_and_taxa():
    assert _mitab.mi_terms(
        'MI:0114(x-ray crystallography)|psi-mi:MI:0496(bait)|'
        'psi-mi:"MI:0114"(x-ray crystallography)|psi-mi:jak2_mouse(display_long)'
    ) == ['MI:0114', 'MI:0496']
    # DIP repeats its own record number under the pubmed database name.
    assert _mitab.pmid_curies(
        'pubmed:9168119|pubmed:DIP-17838S|imex:IM-22293|pubmed:9168119'
    ) == ['PMID:9168119']
    assert _mitab.pmid_curies('unassigned1304') == []
    assert _mitab.taxon_curies(
        'taxid:9606(human)|taxid:9606(Homo sapiens)'
    ) == ['NCBITaxon:9606']
    # In vitro, unknown and the zero of InnateDB are not organisms.
    assert _mitab.taxon_curies('taxid:-1(in vitro)') == []
    assert _mitab.taxon_curies('taxid:0') == []


def test_mitab_identifiers_follow_the_database_name_and_accession_syntax():
    cell = (
        'biogrid:112315|entrez gene/locuslink:MAP2K4|'
        'entrez gene/locuslink:6416|uniprot/swiss-prot:P45985|'
        'refseq:NP_003001|refseq:NM_003010|ensembl:ENSP00000262445|'
        'ensembl:ZAK|unknown database:X1'
    )
    assert _mitab.identifier_pairs(cell) == [
        (Namespace.BIOGRID, '112315'),
        (Namespace.ENTREZ, '6416'),
        (Namespace.UNIPROT, 'P45985'),
        (Namespace.REFSEQ_PROTEIN, 'NP_003001'),
        (Namespace.REFSEQ, 'NM_003010'),
        (Namespace.ENSP, 'ENSP00000262445'),
    ]


def test_mitab_predicate_is_direct_only_for_direct_interaction_types():
    direct = slots.directly_physically_interacts_with
    assert _mitab.physical_predicate(['MI:0407']) is direct
    # A phosphorylation reaction is an enzymatic reaction, a direct interaction.
    assert _mitab.physical_predicate(['MI:0217']) is direct
    assert _mitab.physical_predicate(['MI:0915']) is slots.interacts_with
    assert _mitab.physical_predicate(['MI:0914']) is slots.interacts_with
    # Colocalization, genetic interaction and unknown terms assert neither.
    for other in ('MI:0403', 'MI:0208', 'MI:0401', 'MI:9999'):
        assert _mitab.physical_predicate([other]) is None
    assert _mitab.physical_predicate([]) is None


# InnateDB ----------------------------------------------------------------


def test_module_import_defines_the_curated_dataset_without_downloading():
    assert list(innatedb.resource.datasets()) == ['interactions']
    assert innatedb.config.id.name == 'INNATEDB'
    url = innatedb.resource.interactions.download.url
    assert url.endswith('/innatedb_all.mitab.gz')


def test_reader_keeps_the_columns_of_the_header_line():
    rows = read(PHYSICAL_LINE, DNA_LINE)
    assert len(rows) == 2
    assert len(rows[0]) == 42
    assert rows[0]['#unique_identifier_A'] == 'innatedb:IDBG-25842'
    assert rows[1]['interactor_type_B'] == (
        'psi-mi:"MI:0319"(deoxyribonucleic acid)'
    )


def test_reader_rejects_a_file_with_another_layout():
    def parse(text):
        parser = innatedb.resource.datasets()['interactions']._raw_parser
        return list(parser(SimpleNamespace(result=io.StringIO(text))))

    with pytest.raises(ValueError, match='received an HTML document'):
        parse('<!DOCTYPE html>\n<html>\n')
    # The IMEx file InnateDB also offers has no header line.
    with pytest.raises(ValueError, match="missing .*'unique_identifier_B'"):
        parse('uniprotkb:P20491\tuniprotkb:Q62120\nuniprotkb:Q62120\tuniprotkb:P20491\n')
    with pytest.raises(ValueError, match='line 2: .*fields beyond the header'):
        parse(HEADER + PHYSICAL_LINE.rstrip('\n') + '\t-\n')
    assert parse('') == []
    assert parse(HEADER) == []


def test_physical_association_keeps_identifiers_method_roles_and_reference():
    (row,) = read(PHYSICAL_LINE)
    relation = innatedb.map_interaction(row)
    assert slot_name(relation.predicate) == 'interacts_with'
    assert identifiers(relation.subject) == [
        ('innatedb_gene', 'IDBG-25842'),
        (Namespace.ENSG, 'ENSG00000154589'),
        (Namespace.UNIPROT_ENTRY, 'LY96_HUMAN'),
        (Namespace.REFSEQ_PROTEIN, 'NP_056179'),
        (Namespace.UNIPROT, 'Q9Y6Y9'),
        (Namespace.REFSEQ_PROTEIN, 'NP_001182726'),
        (Namespace.GENESYMBOL, 'LY96'),
    ]
    assert identifiers(relation.object)[0] == ('innatedb_gene', 'IDBG-82738')
    assert [(i.type, i.value) for i in relation.identifiers] == [
        ('innatedb', 'IDB-113240'),
    ]
    # The participant is the gene record; no protein sequence is asserted.
    assert relation.subject.molecular_form is None

    out = extract([row])
    (observed,) = out.relations
    assert observed.predicate == 'interacts_with'
    assert values(observed, 'has_evidence_of_type') == ['MI:0007']
    assert values(observed, 'original_predicate') == ['MI:0915']
    assert values(observed, 'publications') == ['PMID:10359581']
    assert values(observed, 'in_taxon') == ['NCBITaxon:10090']
    assert values(observed, 'psi_mi:experimental_role', 'subject') == ['MI:0498']
    assert values(observed, 'psi_mi:experimental_role', 'object') == ['MI:0496']
    assert values(
        observed, 'psi_mi:participant_identification_method', 'subject',
    ) == ['MI:0363']
    # An unspecified role is not a role.
    assert not values(observed, 'psi_mi:biological_role', 'subject')
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert {e.taxon for e in out.entities.values()} == {'9606'}


def test_enzymatic_reaction_is_direct_and_gets_no_invented_effect():
    (observed,) = extract(read(PHOSPHORYLATION_LINE)).relations
    assert observed.predicate == 'directly_physically_interacts_with'
    assert values(observed, 'original_predicate') == ['MI:0217']
    assert values(observed, 'psi_mi:biological_role', 'subject') == ['MI:0501']
    assert values(observed, 'psi_mi:biological_role', 'object') == ['MI:0502']
    terms = {a['term'] for a in observed.annotations}
    assert not {term for term in terms if term.endswith('_qualifier')}
    # A self interaction stays one entity on both sides.
    assert observed.subject_entity_key == observed.object_entity_key


def test_dna_and_rna_participants_keep_the_type_the_source_states():
    dna, rna = (innatedb.map_interaction(row) for row in read(DNA_LINE, RNA_LINE))
    assert (dna.subject.type, dna.object.type) == ('protein', 'gene')
    assert (rna.subject.type, rna.object.type) == ('protein', 'rna_product')
    assert identifiers(rna.subject) == [
        ('innatedb_gene', 'IDBG-193424'),
        (Namespace.ENSG, 'ENSMUSG00000007670'),
        (Namespace.UNIPROT_ENTRY, 'FUBP2_MOUSE'),
        (Namespace.REFSEQ_PROTEIN, 'NP_034743'),
        (Namespace.UNIPROT, 'Q3U0V1'),
        (Namespace.GENESYMBOL, 'Khsrp'),
    ]
    out = extract(read(DNA_LINE, RNA_LINE))
    observed_dna, observed_rna = out.relations
    # The file writes taxid:0 where no host organism was recorded.
    assert not values(observed_dna, 'in_taxon')
    assert values(observed_rna, 'in_taxon') == ['NCBITaxon:10090']
    assert values(
        observed_rna, 'psi_mi:participant_identification_method', 'object',
    ) == ['MI:0088']
    assert sorted({e.taxon for e in out.entities.values()}) == ['10090', '9606']


def test_record_numbers_and_placeholder_names():
    allergy, placeholder = (
        innatedb.map_interaction(row)
        for row in read(ALLERGY_LINE, PLACEHOLDER_LINE)
    )
    assert [(i.type, i.value) for i in allergy.identifiers] == [
        ('innatedb', 'IDB-648953'),
    ]
    # A gene without a symbol has its record number as display name.
    assert identifiers(placeholder.object) == [
        ('innatedb_gene', 'IDBG-303636'),
        (Namespace.ENSG, 'ENSG00000237708'),
    ]


def test_symbol_that_comes_from_ensembl_is_a_gene_symbol():
    row = read(PLACEHOLDER_LINE)[0]
    row['alias_B'] = 'uniprotkb:MLTK_HUMAN|uniprotkb:Q9NYL2|ensembl:ZAK(display_short)'
    assert identifiers(innatedb.map_interaction(row).object)[2:] == [
        (Namespace.UNIPROT_ENTRY, 'MLTK_HUMAN'),
        (Namespace.UNIPROT, 'Q9NYL2'),
        (Namespace.GENESYMBOL, 'ZAK'),
    ]


def test_term_that_is_not_an_interaction_type_gives_the_broad_predicate():
    (observed,) = extract(read(UNTYPED_LINE)).relations
    assert observed.predicate == 'interacts_with'
    assert values(observed, 'original_predicate') == ['MI:0401']
    assert values(observed, 'has_evidence_of_type') == ['MI:0435']


def test_no_relation_from_imported_negative_or_unknown_records():
    (imported,) = read(IMPORTED_LINE)
    assert innatedb.exclusion_reason(imported) == 'source database is not InnateDB'
    assert innatedb.map_interaction(imported) is None

    def changed(**fields):
        row = dict(read(PHYSICAL_LINE)[0])
        row.update(fields)
        return row

    cases = {
        'negative result': changed(negative='true'),
        'unknown interaction type': changed(
            interaction_type='psi-mi:"MI:0403"(colocalization)',
        ),
        'unknown interactor type': changed(
            interactor_type_B='psi-mi:"MI:0328"(small molecule)',
        ),
        'participant without an InnateDB gene record': changed(
            unique_identifier_B='-',
        ),
    }
    for reason, row in cases.items():
        assert innatedb.exclusion_reason(row) == reason
        assert innatedb.map_interaction(row) is None
    assert not extract([imported, *cases.values()]).relations
    assert innatedb.exclusion_reason(read(PHYSICAL_LINE)[0]) is None
