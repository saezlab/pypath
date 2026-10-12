"""Source contract of the DOMINO interaction input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.domino import exclusion_reason, modified_residue, resource
from pypath.inputs_v2.parsers.domino import (
    COLUMNS,
    FIELDS,
    iter_domino,
    parse_binding_site,
    parse_ptms,
)


def row(**cells):
    """A row of the file; cells that are not given are ``-``, as in the file."""
    values = dict.fromkeys(FIELDS, '-')
    values.update(
        parameters='',
        source_databases='MI:0471(mint)',
        biological_role_a='MI:0499(unspecified role)',
        biological_role_b='MI:0499(unspecified role)',
        experimental_role_a='MI:0496(bait)',
        experimental_role_b='MI:0498(prey)',
        interactor_type_a='protein',
        interactor_type_b='protein',
    )
    values.update(cells)
    return values


# Rows of 2009-10-22-domino-full-binary.mitab26. Long alias cells are
# shortened and the confidence, dataset and figure legend cells left out.
TWO_RANGES_ROW = row(
    id_a='uniprotkb:O00213',
    id_b='uniprotkb:P00533',
    aliases_a=(
        'uniprotkb:Fe65 protein(gene name synonym)|uniprotkb:RIR(gene name synonym)'
        '|uniprotkb:APBB1(gene name)'
    ),
    aliases_b='uniprotkb:ERBB1(gene name synonym)|uniprotkb:EGFR(gene name)',
    detection_methods='MI:0089(protein array)',
    publications='pubmed:16273093',
    taxid_a='taxid:9606(Homo sapiens)',
    taxid_b='taxid:9606(Homo sapiens)',
    interaction_types='MI:0407(direct interaction)',
    interaction_ids='mint:MINT-74468',
    annotations_interaction='kd:1500 nanomolar|figure-legend:ST4C',
    host_organism='taxid:-1(in vitro)',
    binding_site_a='370-509,543-664;-;interpro:IPR006020(PTB_PID)',
    binding_site_b='989-1005;-;-',
    ptm_b=(
        '998;MI:0170(phosphorylated residue);range MINT-1497958: original '
        'sequence: 998-998 : HLPSPTDSNFYRALMDEEDM;-'
    ),
)
TWO_DOMAINS_ROW = row(
    id_a='uniprotkb:P06239',
    id_b='uniprotkb:P20963',
    aliases_a='uniprotkb:LSK(gene name synonym)|uniprotkb:LCK(gene name)',
    aliases_b='uniprotkb:CD3Z(gene name synonym)|uniprotkb:CD247(gene name)',
    detection_methods='MI:0096(pull down)',
    publications='pubmed:8626561',
    taxid_a='taxid:9606(Homo sapiens)',
    taxid_b='taxid:9606(Homo sapiens)',
    interaction_types='MI:0915(physical association)',
    interaction_ids='mint:MINT-60782',
    host_organism='taxid:-1(in vitro)',
    binding_site_a='63-208;-;interpro:IPR000980(SH2),interpro:IPR001452(SH3)',
)
PREREQUISITE_PTM_ROW = row(
    id_a='uniprotkb:O15117',
    id_b='uniprotkb:Q86WV1',
    aliases_a='uniprotkb:SLAP-130(gene name synonym)|uniprotkb:FYB(gene name)',
    aliases_b='uniprotkb:SKAP55(gene name synonym)|uniprotkb:SKAP1(gene name)',
    detection_methods='MI:0107(surface plasmon resonance)',
    publications='doi:10.1074/jbc.M508774200',
    taxid_a='taxid:9606(Homo sapiens)',
    taxid_b='taxid:9606(Homo sapiens)',
    interaction_types='MI:0407(direct interaction)',
    interaction_ids='mint:MINT-7147926',
    host_organism='taxid:-1(in vitro)',
    caution='curation-level:MIMIx',
    binding_site_a='514-570;-;interpro:IPR001452(SH3)',
    binding_site_b='288-300;-;-',
    ptm_b=(
        '294;MI:0170(phosphorylated residue);prerequisite-ptm;-'
        '|298;MI:0170(phosphorylated residue);prerequisite-ptm;-'
    ),
)
ISOFORM_ROW = row(
    id_a='uniprotkb:O15350-1',
    id_b='uniprotkb:Q96J02',
    aliases_a='uniprotkb:P73(gene name synonym)|uniprotkb:TP73(gene name)',
    aliases_b='uniprotkb:ITCH(gene name)',
    detection_methods='MI:0084(phage display)',
    publications='pubmed:15678106',
    taxid_a='taxid:9606(Homo sapiens)',
    taxid_b='taxid:9606(Homo sapiens)',
    interaction_types='MI:0915(physical association)',
    interaction_ids='mint:MINT-3389044',
    host_organism='taxid:-1(in vitro)',
    binding_site_a='452-489;-;-',
    binding_site_b='250-520;-;interpro:IPR002349(WW)',
)
ENZYME_ROW = row(
    id_a='uniprotkb:O13016',
    id_b='uniprotkb:P46108',
    aliases_a='uniprotkb:PTPN1(gene name)|uniprotkb:PTP1B(gene name synonym)',
    aliases_b='uniprotkb:CRK(gene name)|uniprotkb:p38(gene name synonym)',
    detection_methods='MI:0434(phosphatase assay)',
    publications='pubmed:12799422',
    taxid_a='taxid:9031(Gallus gallus)',
    taxid_b='taxid:9606(Homo sapiens)',
    interaction_types='MI:0203(dephosphorylation reaction)',
    interaction_ids='mint:MINT-6763908',
    biological_role_a='MI:0501(enzyme)',
    biological_role_b='MI:0502(enzyme target)',
    experimental_role_a='MI:0497(neutral component)',
    experimental_role_b='MI:0497(neutral component)',
    xrefs_a='uniprotkb:O13016(O13016)',
    xrefs_b='uniprotkb:P46108(P46108)',
    host_organism='taxid:9606(Homo sapiens)',
    binding_site_a='0;-;interpro:IPR012265(Tyr_Phos_no_rcpt)',
    mutations_a='215;MI:0573(mutation disrupting interaction);MI:0074(mutation analysis);-',
    mutations_b='221;MI:0119(mutation decreasing interaction);MI:0074(mutation analysis);-',
)
YEAST_ORF_ROW = row(
    id_a='uniprotkb:P22696',
    id_b='cygd:YPR197C',
    aliases_a=(
        'uniprotkb:J1452(orf name)|uniprotkb:YJR017C(locus name)'
        '|uniprotkb:PIN1(gene name synonym)|uniprotkb:ESS1(gene name)'
    ),
    detection_methods='MI:0089(protein array)',
    publications='pubmed:16606443',
    taxid_a='taxid:4932(Saccharomyces cerevisiae)',
    taxid_b='taxid:4932(Saccharomyces cerevisiae)',
    interaction_types='MI:0915(physical association)',
    interaction_ids='mint:MINT-2851384',
    xrefs_a='uniprotkb:YJR017C(ESS1)',
    host_organism='taxid:-1(in vitro)',
    binding_site_a='5-45;-;interpro:IPR001202(WW_Rsp5_WWP)',
)
COLOCALIZATION_ROW = row(
    id_a='uniprotkb:O00401',
    id_b='uniprotkb:Q15811',
    aliases_a='uniprotkb:WASL(gene name)',
    aliases_b='uniprotkb:ITSN(gene name synonym)|uniprotkb:ITSN1(gene name)',
    detection_methods='MI:0023(colocalization/visualisation technologies)',
    publications='pubmed:11584276',
    taxid_a='taxid:9606(Homo sapiens)',
    taxid_b='taxid:9606(Homo sapiens)',
    interaction_types='MI:0403(colocalization)',
    interaction_ids='mint:MINT-14607',
    experimental_role_a='MI:0499(unspecified role)',
    experimental_role_b='MI:0499(unspecified role)',
    host_organism='taxid:-4(in vivo)',
    binding_site_a='265-391;MI:0033(deletion analysis);-',
    binding_site_b='740-806;MI:0033(deletion analysis);interpro:IPR001452(SH3)',
)
NEGATIVE_ROW = row(
    id_a='uniprotkb:O54864',
    id_b='uniprotkb:Q61686',
    aliases_a='uniprotkb:Suv39h(gene name synonym)|uniprotkb:Suv39h1(gene name)',
    aliases_b='uniprotkb:Hp1a(gene name synonym)|uniprotkb:Cbx5(gene name)',
    detection_methods='MI:0096(pull down)',
    publications='pubmed:12565857',
    taxid_a='taxid:10090(Mus musculus)',
    taxid_b='taxid:10090(Mus musculus)',
    interaction_types='MI:0407(direct interaction)',
    interaction_ids='mint:MINT-20117',
    host_organism='taxid:-1(in vitro)',
    caution='negative:the interactions are not observed under these experimental conditions',
    binding_site_b='121-179;-;interpro:IPR008251(Chromo_shadow)',
    negative='the interactions are not observed under these experimental conditions',
)
PEPTIDE_ROW = row(
    id_a='',
    id_b='uniprotkb:O00459',
    aliases_b='uniprotkb:PI3-kinase p85 subunit beta(gene name synonym)|uniprotkb:PIK3R2(gene name)',
    detection_methods='MI:0016(circular dichroism)',
    publications='doi:10.1002/psc.1132',
    taxid_a='taxid:-2(chemical synthesis)',
    taxid_b='taxid:9606(Homo sapiens)',
    interaction_types='MI:0407(direct interaction)',
    interaction_ids='mint:MINT-7226016',
    experimental_role_a='MI:0497(neutral component)',
    experimental_role_b='MI:0497(neutral component)',
    interactor_type_a='peptide',
    host_organism='taxid:-1(in vitro)',
    caution='curation-level:MIMIx',
    binding_site_a='0;-;interpro:IPR006568',
    binding_site_b='4-80;-;interpro:IPR001452',
)

PROTEIN = {'identifier': None, 'coordinate_system': 'protein', 'position_base': 1}


def mapper(source_row):
    return resource.datasets()['interactions'].mapper(source_row)


def extract(source_row):
    record = mapper(source_row)
    assert record is not None
    out = SilverExtractor('domino', 'interactions')
    out.process_record(record, source_row, 'fixture:1', 1)
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
        (a['term'], a['value'])
        for a in relation.annotations
        if a['scope'] == scope
    )


def regions(entity):
    return [
        (
            region['type'],
            region['identifier'],
            region['position'],
            region['end_position'],
            region['coordinate_reference'],
        )
        for region in (entity.molecular_form or {}).get('regions') or []
    ]


def test_direct_interaction_with_method_publication_and_record():
    relation, first, second = extract(TWO_RANGES_ROW)
    assert relation.predicate == 'directly_physically_interacts_with'
    assert relation.upstream_id == 'MINT-74468'
    assert annotations(relation) == [
        ('has_evidence_of_type', 'MI:0089'),
        ('original_predicate', 'MI:0407'),
        ('publications', 'PMID:16273093'),
        ('supporting_data_source', 'MI:0471'),
    ]
    assert identifiers(first) == {
        (Namespace.UNIPROT, 'O00213'),
        (Namespace.GENESYMBOL, 'APBB1'),
        (Namespace.SYNONYM, 'Fe65 protein'),
        (Namespace.SYNONYM, 'RIR'),
    }
    assert identifiers(second) == {
        (Namespace.UNIPROT, 'P00533'),
        (Namespace.GENESYMBOL, 'EGFR'),
        (Namespace.SYNONYM, 'ERBB1'),
    }
    assert (first.namespace, first.identifier) == ('uniprot', 'O00213')
    assert {first.entity_type, second.entity_type} == {'protein'}
    assert {first.taxon, second.taxon} == {'9606'}
    assert annotations(relation, 'subject') == [
        ('in_taxon', 'NCBITaxon:9606'),
        ('psi_mi:experimental_role', 'MI:0496'),
    ]
    assert annotations(relation, 'object') == [
        ('in_taxon', 'NCBITaxon:9606'),
        ('psi_mi:experimental_role', 'MI:0498'),
    ]


def test_binding_site_ranges_are_regions_named_by_the_interpro_entry():
    _, first, second = extract(TWO_RANGES_ROW)
    entry = {'ns': 'interpro', 'id': 'IPR006020'}
    assert regions(first) == [
        ('binding site', entry, 370, 509, PROTEIN),
        ('binding site', entry, 543, 664, PROTEIN),
    ]
    # No domain is named for the site of the second protein.
    assert regions(second) == [('binding site', None, 989, 1005, PROTEIN)]


def test_modified_residue_is_a_modification_of_the_participant():
    _, first, second = extract(TWO_RANGES_ROW)
    assert first.molecular_form['modifications'] is None
    (modification,) = second.molecular_form['modifications']
    assert modification == {
        'term': 'phosphorylated residue',
        # The residue is read from the window of the note: 998 is the eleventh
        # residue of HLPSPTDSNFYRALMDEEDM.
        'residue': 'Y',
        'position': 998,
        'end_position': 998,
        'coordinate_reference': PROTEIN,
        'description': None,
    }
    _, _, second = extract(PREREQUISITE_PTM_ROW)
    assert [
        (m['term'], m['position'])
        for m in second.molecular_form['modifications']
    ] == [('phosphorylated residue', 294), ('phosphorylated residue', 298)]


def ptm(cell):
    (entry,) = parse_ptms(cell)
    return entry


def test_note_window_is_read_with_or_without_the_prerequisite_suffix():
    window = ptm(
        '998;MI:0170(phosphorylated residue);range MINT-1497979: original '
        'sequence: 998-998 : HLPSPTDSNFYRALMDEEDM;-'
    )['window']
    assert window == {'start': 998, 'end': 998, 'sequence': 'HLPSPTDSNFYRALMDEEDM'}
    # The source writes the suffix straight after the residues.
    glued = ptm(
        '5;MI:0165(n6-methyl-lysine);range MINT-6542511: original sequence: '
        '5-5 : MARTKQTARKSTGGprerequisite-ptm;-'
    )
    assert glued['window']['sequence'] == 'MARTKQTARKSTGG'
    assert glued['note'].endswith('prerequisite-ptm')
    assert ptm('294;MI:0170(phosphorylated residue);prerequisite-ptm;-')['window'] is None


@pytest.mark.parametrize(
    ('cell', 'residue'),
    [
        # The modified residue is the eleventh of the window of twenty.
        (
            '1222;MI:0170(phosphorylated residue);range MINT-1491479: original '
            'sequence: 1222-1222 : AFSPAFDNLYYWDQDPPERG;-',
            'Y',
        ),
        (
            '367;MI:0170(phosphorylated residue);range MINT-2789388: original '
            'sequence: 367-367 : NDVPDCRRTISAPVVRPKDA;-',
            'S',
        ),
        # Among the first ten residues the window starts at the first one.
        (
            '3;MI:0176(o-phospho-serine);range MINT-721509: original sequence: '
            '3-3 : MASGVAVSDGVI;-',
            'S',
        ),
        (
            '5;MI:0165(n6-methyl-lysine);range MINT-6542511: original sequence: '
            '5-5 : MARTKQTARKSTGGprerequisite-ptm;-',
            'K',
        ),
        # No window, no residue.
        ('294;MI:0170(phosphorylated residue);prerequisite-ptm;-', None),
        # The residue the layout points to is not one the modification can be on.
        (
            '1272;MI:0176(o-phospho-serine);range MINT-1523137: original sequence: '
            '1272-1272 : IISSIKEEMEPGFREVSFYY;-',
            None,
        ),
        (
            '20;MI:0165(n6-methyl-lysine);range MINT-2830336: original sequence: '
            '20-20 : GLGKGGAKRHRKVLRDNIQG;-',
            None,
        ),
        # The note gives a range for a modification that has one position, or the
        # modification has a range.
        (
            '681-685;MI:0170(phosphorylated residue);range MINT-6793956: original '
            'sequence: 681-685 : MSRDIYSTDYYRVGGRTMLPIRWMprerequisite-ptm;-',
            None,
        ),
        # A window of ten residues does not have the layout.
        (
            '198;MI:0177(o-phospho-threonine);range MINT-716016: original sequence: '
            '198-198 : PKKPGLRRRQ;-',
            None,
        ),
        # A modification the crosswalk does not list gets no residue.
        (
            '14;MI:0121(acetylated residue);range MINT-6491909: original sequence: '
            '14-14 : TKQTARKSTGGKAPRKQLATprerequisite-ptm;-',
            None,
        ),
    ],
)
def test_residue_is_read_from_the_window_only_where_everything_agrees(cell, residue):
    assert modified_residue(ptm(cell)) == residue


def test_note_that_names_another_position_than_the_modification_gives_no_residue():
    entry = ptm(
        '998;MI:0170(phosphorylated residue);range MINT-1497979: original '
        'sequence: 998-998 : HLPSPTDSNFYRALMDEEDM;-'
    )
    assert modified_residue(entry) == 'Y'
    assert modified_residue({**entry, 'start': 999, 'end': 999}) is None
    assert modified_residue({**entry, 'start': None, 'end': None}) is None


def test_range_with_two_domains_is_not_assigned_to_either():
    relation, first, _ = extract(TWO_DOMAINS_ROW)
    assert relation.predicate == 'interacts_with'
    assert regions(first) == [
        ('binding site', None, 63, 208, PROTEIN),
        ('binding site', {'ns': 'interpro', 'id': 'IPR000980'}, None, None, None),
        ('binding site', {'ns': 'interpro', 'id': 'IPR001452'}, None, None, None),
    ]


def test_isoform_is_the_reference_of_its_positions():
    _, first, second = extract(ISOFORM_ROW)
    isoform = {'ns': 'uniprot', 'id': 'O15350-1'}
    assert (Namespace.UNIPROT, 'O15350-1') in identifiers(first)
    assert first.molecular_form['isoform_identifier'] == isoform
    assert regions(first) == [
        ('binding site', None, 452, 489, {**PROTEIN, 'identifier': isoform}),
    ]
    assert regions(second) == [
        ('binding site', {'ns': 'interpro', 'id': 'IPR002349'}, 250, 520, PROTEIN),
    ]


def test_enzymatic_reaction_is_a_direct_interaction_with_roles():
    relation, enzyme, target = extract(ENZYME_ROW)
    assert relation.predicate == 'directly_physically_interacts_with'
    assert ('original_predicate', 'MI:0203') in annotations(relation)
    assert ('in_taxon', 'NCBITaxon:9606') in annotations(relation)
    assert annotations(relation, 'subject') == [
        ('in_taxon', 'NCBITaxon:9031'),
        ('psi_mi:biological_role', 'MI:0501'),
        ('psi_mi:experimental_role', 'MI:0497'),
    ]
    assert ('psi_mi:biological_role', 'MI:0502') in annotations(relation, 'object')
    # A site of unknown position (0) keeps its domain and gets no position.
    assert regions(enzyme) == [
        ('binding site', {'ns': 'interpro', 'id': 'IPR012265'}, None, None, None),
    ]
    assert (enzyme.taxon, target.taxon) == ('9031', '9606')


def test_mutations_are_not_features_of_the_interacting_proteins():
    # A mutation that disrupts the interaction is not a property of the
    # proteins that interact.
    relation, enzyme, target = extract(ENZYME_ROW)
    assert enzyme.molecular_form['variants'] is None
    assert target.molecular_form is None
    assert relation.predicate != 'affects'
    assert not any('qualifier' in a['term'] for a in relation.annotations)


def test_doi_is_the_publication_of_an_article_cited_by_doi():
    relation, _, _ = extract(PREREQUISITE_PTM_ROW)
    assert ('publications', 'doi:10.1074/jbc.M508774200') in annotations(relation)
    relation, _, _ = extract(dict(TWO_RANGES_ROW, publications='-'))
    assert all(a['term'] != 'publications' for a in relation.annotations)


def test_yeast_orf_name_identifies_a_participant():
    relation, _, second = extract(YEAST_ORF_ROW)
    assert relation.predicate == 'interacts_with'
    assert identifiers(second) == {('cygd', 'YPR197C')}
    assert second.taxon == '4932'


@pytest.mark.parametrize(
    'source_row,reason',
    [
        (NEGATIVE_ROW, 'negative result'),
        (PEPTIDE_ROW, 'peptide without an identifier'),
        (COLOCALIZATION_ROW, 'colocalization'),
        (
            dict(TWO_RANGES_ROW, interaction_types='MI:0208(genetic interaction)'),
            'unknown interaction type',
        ),
        (dict(TWO_RANGES_ROW, interaction_types='-'), 'unknown interaction type'),
        (
            dict(TWO_RANGES_ROW, id_b='genbank_protein_gi:AAH62987'),
            'participant without a known identifier',
        ),
        (dict(TWO_RANGES_ROW, interactor_type_b='gene'), 'unknown interactor type'),
    ],
)
def test_rows_that_give_no_relation(source_row, reason):
    assert exclusion_reason(source_row) == reason
    assert mapper(source_row) is None


def test_negative_result_is_excluded_whichever_column_says_so():
    assert mapper(dict(NEGATIVE_ROW, caution='-')) is None
    assert mapper(dict(NEGATIVE_ROW, negative='-')) is None


def test_obsolete_physical_interaction_type_is_the_broad_predicate():
    relation, _, _ = extract(
        dict(TWO_RANGES_ROW, interaction_types='MI:0218(physical interaction)')
    )
    assert relation.predicate == 'interacts_with'


def test_domino_record_number_is_not_a_mint_accession():
    relation, _, _ = extract(
        dict(TWO_RANGES_ROW, interaction_ids='mint:DOMINO-1234')
    )
    assert relation.upstream_id == 'fixture:1'
    relation, _, _ = extract(
        dict(TWO_RANGES_ROW, interaction_ids='mint:MINT-74468|imex:IM-11788-3')
    )
    assert relation.upstream_id == 'MINT-74468'


def test_binding_site_cell_grammar():
    assert parse_binding_site('-') is None
    assert parse_binding_site('15-92') is None
    assert parse_binding_site(
        '431-545,763-880;MI:0033(deletion analysis);interpro:IPR001478(PDZ)'
    ) == {
        'ranges': [(431, 545), (763, 880)],
        'unknown_ranges': [],
        'method': 'MI:0033',
        'domains': [('IPR001478', 'PDZ')],
        'unknown_domains': [],
    }
    # 0 stands for an unknown position; an end before the start is not a range.
    assert parse_binding_site('0;-;interpro:IPR001452')['unknown_ranges'] == ['0']
    assert parse_binding_site('351-349;-;-') == {
        'ranges': [],
        'unknown_ranges': ['351-349'],
        'method': None,
        'domains': [],
        'unknown_domains': [],
    }
    site = parse_binding_site('58-128;MI:0033(deletion analysis);interpro:CSD(CSD)')
    assert site['domains'] == [] and site['unknown_domains'] == ['interpro:CSD(CSD)']


def test_ptm_cell_grammar():
    assert parse_ptms('-') == []
    assert parse_ptms(
        '20;MI:0165(n6-methyl-lysine);prerequisite-ptm;MI:0114(x-ray crystallography)'
        '|0;MI:0166("n6,n6-dimethyl-lysine");-;-'
    ) == [
        {
            'start': 20,
            'end': 20,
            'term': 'MI:0165',
            'name': 'n6-methyl-lysine',
            'note': 'prerequisite-ptm',
            'method': 'MI:0114',
            'window': None,
        },
        {
            'start': None,
            'end': None,
            'term': 'MI:0166',
            'name': 'n6,n6-dimethyl-lysine',
            'note': None,
            'method': None,
            'window': None,
        },
    ]
    assert parse_ptms('137;o-phospho-serine;-;-') == []


HEADER = '\t'.join(header for header, _ in COLUMNS) + '\t'


def line(source_row):
    return '\t'.join(source_row[field] for field in FIELDS)


def table(lines):
    return SimpleNamespace(result=io.StringIO('\n'.join(lines) + '\n'))


def test_reader_names_both_interactor_type_columns():
    rows = list(iter_domino(table([HEADER, line(PEPTIDE_ROW), line(ISOFORM_ROW)])))
    assert rows == [PEPTIDE_ROW, ISOFORM_ROW]
    assert (rows[0]['interactor_type_a'], rows[0]['interactor_type_b']) == (
        'peptide',
        'protein',
    )


def test_reader_reads_bytes():
    data = '\n'.join([HEADER, line(ISOFORM_ROW)]).encode('utf-8')
    assert list(iter_domino(SimpleNamespace(result=io.BytesIO(data)))) == [
        ISOFORM_ROW
    ]


def test_reader_rejects_a_changed_header():
    header = HEADER.replace('binding site A', 'Feature(s) interactor A')
    with pytest.raises(ValueError, match='expected header'):
        list(iter_domino(table([header, line(ISOFORM_ROW)])))


def test_reader_rejects_a_short_row_and_an_empty_file():
    with pytest.raises(ValueError, match='row 2 has 39 cells'):
        list(iter_domino(table([HEADER, line(ISOFORM_ROW).rsplit('\t', 1)[0]])))
    with pytest.raises(ValueError, match='no rows'):
        list(iter_domino(table([HEADER])))
