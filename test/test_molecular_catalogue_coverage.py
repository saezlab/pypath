"""Source feature observations stay separate from alternative catalogue states."""

import csv
import io
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from biolink_model.datamodel import model
from omnipath_build.silver import SilverExtractor
from pypath.inputs_v2 import brenda, chembl, mirbase, recon3d, uniprot
from pypath.inputs_v2.parsers import brenda as bp, recon3d as rp
from pypath.internals.tabular_builder import (
    CV,
    EntityBuilder,
    IdentifiersBuilder,
)


def opener(text):
    return SimpleNamespace(result={'data': io.StringIO(text)})


def extract(record):
    result = SilverExtractor('fixture', 'molecular')
    result.process_record(record, raw_payload={}, row_id='1', row_number=1)
    return result


@pytest.mark.parametrize(
    'entity_type,namespace,value,expected',
    [
        (model.Protein, 'uniprot', 'P04637-2', True),
        (model.Protein, 'refseq', 'NP_000537.3', True),
        (model.Protein, 'ensp', 'ENSP00000269305.4', True),
        (model.RNAProduct, 'enst', 'ENST00000269305.9', True),
        (model.Protein, 'uniprot', 'P04637', False),
        (model.Gene, 'entrez', '7157', False),
        (model.ChemicalEntity, 'uniprot', 'P04637-2', False),
    ],
)
def test_generic_primary_product_capture(
    entity_type, namespace, value, expected
):
    record = EntityBuilder(
        entity_type=entity_type,
        identifiers=IdentifiersBuilder(
            CV(term=namespace, value=value),
            CV(term='uniprot', value='P04637-3'),
        ),
    )({})
    assert bool(record.molecular_form) == expected
    if record.molecular_form:
        ids = record.molecular_form['sequence_identifiers'] or []
        assert {'ns': 'uniprot', 'id': 'P04637-3'} not in ids


def tsv(row):
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=list(row), delimiter='\t')
    writer.writeheader()
    writer.writerow(row)
    return stream.getvalue()


def test_uniprot_catalogue_alternatives_are_individual_observations():
    row = {
        'Entry': 'P04637',
        'Sequence': 'MASG',
        'Sequence version': '4',
        'Mutagenesis': 'MUTAGEN 2; /note="A->G,V: two alternatives; distinct mutants"; /evidence="ECO:0000269|PubMed:123"; MUTAGEN <3; /note="Missing: uncertain region"',
        'Modified residue': 'MOD_RES 3; /note="Phosphoserine"',
    }
    features = list(uniprot._catalogue_feature_rows(opener(tsv(row))))
    assert len(features) == 4
    forms = [
        uniprot.catalogue_features_schema(row).molecular_form
        for row in features
    ]
    variants = [form['variants'][0] for form in forms if form['variants']]
    assert [variant['alternate'] for variant in variants] == ['G', 'V', None]
    assert variants[2]['position'] is None
    assert 'distinct mutants' in variants[0]['description']
    assert not forms[0]['modifications']
    assert forms[-1]['modifications'][0]['term'] == 'Phosphoserine'
    assert (
        variants[0]['coordinate_reference']['identifier']['ns']
        == 'protein_sequence_sha256'
    )
    assert features[0]['feature_publications'] == ['123']
    assert uniprot.proteins_schema(row).molecular_form['sequence_identifiers']


def test_uniprot_real_tsv_features_keep_ids_and_positions():
    text = (
        Path(__file__).parent / 'data/molecular/uniprot-P04637.tsv'
    ).read_text()
    rows = list(uniprot._catalogue_feature_rows(opener(text)))
    assert len(rows) > 100
    assert any(
        row['molecular_form'] and row['molecular_form']['variants']
        for row in rows
    )
    assert any(
        row['molecular_form'] and row['molecular_form']['modifications']
        for row in rows
    )
    assert any(
        row['molecular_form']
        and any(
            '-PRO_' in identifier['id']
            for identifier in row['molecular_form']['sequence_identifiers']
            or []
        )
        for row in rows
    )
    assert all(
        extract(uniprot.catalogue_features_schema(row)).entities for row in rows
    )


def test_brenda_protein_records_and_alternatives_do_not_broadcast_state():
    text = """ID\t1.1.1.1
PR\t#1# Homo sapiens {P04637; source: UniProt} <1>
PR\t#2# Homo sapiens {P16455; source: UniProt} <2>
EN\t#1# A2G/S3D (#1# one double mutant <1>) <1>
EN\t#1# A2V (#1# a separate mutant
\twith a continued note <2>) <2>
PM\t#2# no glycoprotein <2>
PM\t#2# glycoprotein <1>
RF\t<1> Citation {Pubmed:123}
RF\t<2> Citation {Pubmed:456}
///
"""
    rows = list(bp.iter_molecular_forms(opener(text)))
    assert [row['UniProt'] for row in rows] == [
        'P04637',
        'P04637',
        'P16455',
        'P16455',
    ]
    assert rows[1]['Refs'] == ['456']
    forms = [brenda.molecular_forms_schema(row).molecular_form for row in rows]
    assert [variant['position'] for variant in forms[0]['variants']] == [2, 3]
    assert len(forms[1]['variants']) == 1
    assert forms[2] is None
    assert forms[3]['modifications'][0]['term'] == 'glycoprotein'
    assert (
        brenda.schema({'UniProt': ['P04637'], 'EC': '1.1.1.1'}).molecular_form
        is None
    )


def test_brenda_ambiguous_accessions_do_not_vote_as_one_protein():
    text = 'ID\t1.1.1.1\nPR\t#1# Human {P04637 AND P16455; source: UniProt} <1>\nEN\t#1# A2G <1>\n///\n'
    (row,) = bp.iter_molecular_forms(opener(text))
    assert row['UniProt'] is None
    record = brenda.molecular_forms_schema(row)
    assert record.identifiers[0].type == 'brenda_protein_record'
    assert row['source_accessions'] == ['P04637', 'P16455']


def test_chembl_assay_variants_remain_per_assay_and_representative():
    base = {
        'target_chembl_id': 'CHEMBLT1',
        'target_type': 'SINGLE PROTEIN',
        'target_component_uniprot_accessions': 'P16455',
        'variant_accession': 'P16455',
        'variant_isoform': 1,
        'variant_version': 2,
    }
    mutant = chembl.target_builder(
        {
            **base,
            'variant_id': 12,
            'variant_mutation': 'A2G/S3D',
            'variant_sequence': 'MGDG',
        }
    )
    other = chembl.target_builder(
        {**base, 'variant_id': 13, 'variant_mutation': 'A2V'}
    )
    assert len(mutant.molecular_form['variants']) == 2
    assert other.molecular_form['variants'][0]['alternate'] == 'V'
    assert (
        mutant.molecular_form['sequence_identifiers'][-1]['ns']
        == 'chembl_representative_sequence_sha256'
    )
    assert (
        mutant.molecular_form['variants'][0]['coordinate_reference'][
            'identifier'
        ]['id']
        == 'P16455-1.2'
    )
    assert chembl.target_builder(base).molecular_form is None
    assert (
        chembl.target_builder(
            {
                **base,
                'target_type': 'PROTEIN COMPLEX',
                'variant_id': 12,
                'variant_mutation': 'A2G',
            }
        ).molecular_form
        is None
    )
    unknown = chembl.target_builder({**base, 'variant_id': -1})
    assert unknown.molecular_form['variants'][0]['position'] is None
    assert (
        chembl.target_builder(
            {**base, 'variant_accession': 'P04637', 'variant_id': 12}
        ).molecular_form['isoform_identifier']
        is None
    )


def test_chembl_component_records_preserve_commas_and_sequence_alignment():
    components = [
        {
            'component_id': 1,
            'component_type': 'PROTEIN',
            'db_source': 'SWISS-PROT',
            'accession': 'P04637-2',
            'description': 'first, with comma',
            'sequence': 'MAAK',
        },
        {
            'component_id': 2,
            'component_type': 'PROTEIN',
            'db_source': 'SWISS-PROT',
            'accession': 'P04637-3',
            'description': 'second',
            'sequence': 'MAAT',
        },
    ]
    record = chembl.targets_schema(
        {
            'target_type': 'PROTEIN COMPLEX',
            'chembl_id': 'CHEMBLT1',
            'component_records': json.dumps(components),
        }
    )
    first, second = [member.member for member in record.membership]
    assert first.annotations[0].value == 'first, with comma'
    assert first.molecular_form['isoform_identifier']['id'] == 'P04637-2'
    assert (
        first.molecular_form['sequence_identifiers']
        != second.molecular_form['sequence_identifiers']
    )


def embl(precursor, location):
    return f"""ID   example standard; RNA;
AC   {precursor};
FT   miRNA           {location}
FT                   /accession="MIMAT1"
FT                   /product="example-5p"
SQ   Sequence 8 BP;
     acguacgu 8
//
"""


def test_mirbase_products_keep_precursor_ranges_and_actual_sequences():
    rows = list(
        mirbase._matures_raw(opener(embl('MI1', '1..4') + embl('MI2', '5..8')))
    )
    assert len(rows) == 1
    assert rows[0]['precursors'] == ['MI1', 'MI2']
    assert rows[0]['sequence'] == 'ACGU'
    assert [region['position'] for region in rows[0]['precursor_regions']] == [
        1,
        5,
    ]
    form = mirbase.matures_schema(rows[0]).molecular_form
    assert form['sequence_identifiers'][0]['ns'] == 'transcript_sequence_sha256'
    assert form['modifications'] is None
    assert {
        relation.predicate
        for relation in extract(mirbase.matures_schema(rows[0])).relations
    } == {'derives_from'}


def test_mirbase_different_sequences_and_fuzzy_positions_are_not_collapsed():
    rows = list(
        mirbase._matures_raw(
            opener(
                embl('MI1', '1..3') + embl('MI2', '2..4') + embl('MI3', '<2..4')
            )
        )
    )
    assert {row['sequence'] for row in rows} == {'ACG', 'CGU', None}
    fuzzy = next(row for row in rows if row['sequence'] is None)
    assert fuzzy['precursor_regions'][0]['location'] == '<2..4'
    assert fuzzy['precursor_regions'][0]['position'] is None
    assert mirbase.matures_schema(fuzzy).molecular_form is None


def test_recon3d_selectors_remain_source_context_without_guessed_forms():
    data = {
        'genes': [{'id': '1_AT1', 'name': 'A'}, {'id': '1_AT2', 'name': 'A'}],
        'metabolites': [],
        'reactions': [
            {
                'id': 'r1',
                'metabolites': {},
                'gene_reaction_rule': '1_AT1 or 1_AT2',
                'lower_bound': 0,
                'upper_bound': 1,
            }
        ],
    }
    (gene,) = rp._parse_genes(data)
    assert gene['source_selectors'] == ['1_AT1', '1_AT2']
    assert recon3d.genes_schema(gene).molecular_form is None
    (row,) = rp._parse_reactions(data)
    assert row['gene_rule_clauses'] == [['1']]
    assert row['source_gene_rule_clauses'] == [['1_AT1'], ['1_AT2']]
    (relation,) = extract(recon3d.reactions_schema(row)).relations
    assert json.loads(
        next(
            annotation['value']
            for annotation in relation.annotations
            if annotation['term'] == 'recon3d:source_product_clauses'
        )
    ) == [['1_AT1'], ['1_AT2']]


@pytest.mark.parametrize('accession', ['P04637-2', 'P04637-PRO_0000185703'])
def test_other_resource_member_and_endpoint_suffixes_survive(accession):
    from pypath.inputs_v2 import cellphonedb, corum, tcdb

    cpdb = cellphonedb.partner_a_builder({'partner_a': accession})
    assert cpdb.type == 'protein'
    assert cpdb.molecular_form['sequence_identifiers'][0]['id'] == accession
    complex_ = corum.complexes_schema(
        {'ComplexID': '1', 'subunits(UniProt IDs)': accession}
    )
    assert (
        complex_.membership[0].member.molecular_form['sequence_identifiers'][0][
            'id'
        ]
        == accession
    )
    transporter = tcdb._transporters_schema({'uniprot': accession})
    assert (
        transporter.molecular_form['sequence_identifiers'][0]['id'] == accession
    )


def test_recon3d_complex_catalogue_keeps_separate_selector_alternatives():
    data = {
        'reactions': [
            {'gene_reaction_rule': '(1_AT1 and 2_AT1) or (1_AT2 and 2_AT1)'}
        ]
    }
    (row,) = rp._parse_enzyme_complexes(data)
    assert len(row['source_gene_rule_clauses']) == 2
    record = recon3d.enzyme_complexes_schema(row)
    assert json.loads(record.annotations[0].value) == [
        ['1_AT1', '2_AT1'],
        ['1_AT2', '2_AT1'],
    ]
    assert all(
        member.member.molecular_form is None for member in record.membership
    )



def test_conflicting_same_site_substitutions_keep_ambiguous_source_description():
    from pypath.inputs_v2._molecular_forms import protein_variants
    variant, = protein_variants('A2G/A2V')
    assert variant['description'] == 'A2G/A2V'
    assert variant['position'] is None
    assert variant['alternate'] is None
