"""Source contract of the LncRNADisease interaction input."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.naming import Namespace

from pypath.inputs_v2.lncdisease import resource
from pypath.inputs_v2.parsers.lncdisease import iter_interactions


def line(*values):
    assert len(values) == 10
    return '\t'.join(values)


# Rows of LncRNADisease_data2.txt (the sentences are shortened).
BINDING = line(
    '5', '7SK', 'HEXIM1', 'RNA-Protein', 'binding',
    'The Hexim1 binding site is located in the 24-87 region of 7SK RNA.',
    'Human', 'RN7SK; 7SK', 'NR_001445.2', '19244621',
)
COEXPRESSION = line(
    '45', 'AF087999', 'FLJ21870', 'RNA-RNA', 'Co-expression',
    'Half-sbsRNA1 siRNA increased the levels of SERPINE1 and FLJ21870 mRNAs.',
    'Human', 'RBM4; LARK; RBM4A', 'AF087999', '21307942',
)
TF_REGULATION = line(
    '24', '7SK', 'OCT1', 'DNA-TF', 'Regulatory',
    'PTF and Oct-1 enhance transcription from the 7SK promoter.',
    'Human', 'RN7SK; 7SK', 'NR_001445.2', '1535687',
)
REGULATION = line(
    '56', 'Airn', 'Igf2r', 'RNA-DNA', 'Regulatory',
    'Controls Igf2r imprinting cluster in mouse', 'Mouse', 'Airn; Air; AI256653',
    'NR_002853.2', '16874305',
)
BINDING_AND_REGULATION = line(
    '40', '7SL', 'MnSOD', 'RNA-RNA', 'binding/Regulatory',
    'The 7SL cytoplasmic non-coding RNA has been shown to bind MnSOD mRNA.',
    'Human', 'RN7SL1; 7SL; 7L1a', 'NR_002715.1', '20951849',
)
GENOMIC_ACCESSION = line(
    '93', 'BCYRN1P2', 'FMR1', 'RNA-Protein', 'binding',
    'The neural phenotype of BC1 (BCYRN1P2) is consistent with its interaction with FMRP.',
    'Human', 'BCYRN1P2; BC200g', 'NG_001085.2?from=101&to=289&report=genbank', '19239885',
)
GI_ACCESSION = line(
    '61', 'AK017619', 'Dab2', 'RNA-RNA', 'Co-expression',
    'The expression profiles of Dab2 and the intronic ncRNA were correlated.',
    'Mouse', 'Dab2; p96; Doc2', '12856954', '18562676',
)
GENE_LEVEL = line(
    '54', 'Airn', 'Igf2r', 'DNA-DNA', 'Regulatory',
    'Deletion of the Airn promoter on the paternal allele results in activation of Igf2r.',
    'Mouse', 'Airn; Air; AI256653', 'NR_002853.2', '19727951',
)
NO_SPECIES = line(
    '355', 'NRSE', 'REST', 'RNA-Protein', 'Regulatory',
    'The ncRNAs NRSE and HSR have been shown to regulate the activity of REST and HSF-1.',
    'N/A', 'N/A', 'N/A', '19191220',
)
NO_LEVEL = line(
    '429', 'SRG1', 'Spt2', 'N/A', 'N/A',
    'We find that Spt2 and Spt6 are required for the repression of SER3 by SRG1 transcription.',
    'Saccharomyces cerevisiae S288c', 'SRG1', 'N/A', '21220514',
)
YEAST = line(
    '166', 'GAL10', 'H3', 'RNA-Protein', 'Regulatory',
    'Cryptic ncRNA transcription can direct H3 K36me3 formation.',
    'Saccharomyces cerevisiae S288c', 'GAL10', 'NM_001178367.1', '19061643',
)


def table(*lines):
    return SimpleNamespace(result=io.StringIO('\r\n'.join(lines) + '\r\n'))


def parse(text):
    (row,) = list(iter_interactions(table(text)))
    return row


def extract(text):
    row = parse(text)
    record = resource.datasets()['interactions'].mapper(row)
    assert record is not None
    out = SilverExtractor('lncdisease', 'interactions')
    out.process_record(record, row, 'fixture:1', 1)
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


def test_binding_is_a_physical_interaction_from_lncrna_to_protein():
    relation, lncrna, protein = extract(BINDING)
    assert relation.predicate == 'interacts_with'
    assert lncrna.entity_type == 'noncoding_rna_product'
    assert (lncrna.namespace, lncrna.identifier) == ('refseq', 'NR_001445.2')
    assert identifiers(lncrna) == {
        (Namespace.REFSEQ, 'NR_001445.2'),
        (Namespace.NAME, '7SK'),
        (Namespace.SYNONYM, 'RN7SK'),
        (Namespace.SYNONYM, '7SK'),
    }
    assert protein.entity_type == 'protein'
    assert identifiers(protein) == {(Namespace.NAME, 'HEXIM1')}
    assert {lncrna.taxon, protein.taxon} == {'9606'}
    assert annotations(relation) == [
        ('description', 'The Hexim1 binding site is located in the 24-87 region of 7SK RNA.'),
        ('publications', 'PMID:19244621'),
    ]


def test_coexpression_is_not_a_physical_interaction():
    relation, _, rna = extract(COEXPRESSION)
    assert relation.predicate == 'coexpressed_with'
    assert rna.entity_type == 'rna_product'
    assert (Namespace.GENBANK, 'AF087999') in identifiers(extract(COEXPRESSION)[1])


def test_regulation_without_a_stated_direction_is_a_plain_association():
    relation, lncrna, gene = extract(REGULATION)
    assert relation.predicate == 'associated_with'
    assert gene.entity_type == 'gene'
    assert lncrna.taxon == gene.taxon == '10090'


def test_transcription_factor_regulates_the_gene_of_the_lncrna():
    relation, tf, gene = extract(TF_REGULATION)
    # The source's level is DNA-TF: the TF acts on the lncRNA gene, so the relation runs from the TF.
    assert relation.predicate == 'affects'
    assert tf.entity_type == 'protein'
    assert identifiers(tf) == {(Namespace.NAME, 'OCT1')}
    assert gene.entity_type == 'gene'
    assert identifiers(gene) == {
        (Namespace.NAME, '7SK'),
        (Namespace.SYNONYM, 'RN7SK'),
        (Namespace.SYNONYM, '7SK'),
    }
    # No effect direction or aspect is stated, so none is invented.
    assert {a['term'] for a in relation.annotations} <= {
        'description', 'publications', 'in_taxon',
    }


def test_gene_level_has_no_transcript_accession():
    _, lncrna_gene, _ = extract(GENE_LEVEL)
    assert lncrna_gene.entity_type == 'gene'
    assert all(i['ns'] != 'refseq' for i in lncrna_gene.identifiers)


def test_binding_and_regulation_is_carried_as_binding_only():
    relation, _, _ = extract(BINDING_AND_REGULATION)
    assert relation.predicate == 'interacts_with'


def test_genomic_region_accession_is_not_an_rna_identifier():
    _, lncrna, _ = extract(GENOMIC_ACCESSION)
    assert identifiers(lncrna) == {
        (Namespace.NAME, 'BCYRN1P2'),
        (Namespace.SYNONYM, 'BCYRN1P2'),
        (Namespace.SYNONYM, 'BC200g'),
    }


def test_gi_number_is_a_genbank_gi_identifier():
    _, lncrna, _ = extract(GI_ACCESSION)
    assert (Namespace.GENBANK_GI, '12856954') in identifiers(lncrna)


def test_missing_species_gives_no_taxon_and_strain_gives_the_strain_taxon():
    _, lncrna, protein = extract(NO_SPECIES)
    assert lncrna.taxon is None and protein.taxon is None
    assert identifiers(lncrna) == {(Namespace.NAME, 'NRSE')}
    _, lncrna, _ = extract(YEAST)
    assert lncrna.taxon == '559292'


def test_row_without_an_interaction_level_gives_no_relation():
    row = parse(NO_LEVEL)
    assert (row['level'], row['first_type'], row['second_type']) == (None, None, None)
    assert resource.datasets()['interactions'].mapper(row) is None


def test_unknown_values_are_errors():
    mapper = resource.datasets()['interactions'].mapper
    with pytest.raises(ValueError, match='Unknown interaction class'):
        mapper(dict(parse(BINDING), interaction_class='Inhibition'))
    with pytest.raises(ValueError, match='Unknown species'):
        mapper(dict(parse(BINDING), species='Zebrafish'))
    with pytest.raises(ValueError, match='Unknown second partner type'):
        mapper(dict(parse(BINDING), level='RNA-Lipid', second_type='Lipid'))
    with pytest.raises(ValueError, match='Unknown sequence accession'):
        mapper(dict(parse(BINDING), accession='XYZ'))


@pytest.mark.parametrize('number', ['0123456', '123456789'])
def test_number_with_a_leading_zero_or_nine_digits_is_not_a_publication(number):
    text = TF_REGULATION.replace('1535687', number)
    assert parse(text)['pmids'] == []
    relation, _, _ = extract(text)
    assert not [1 for a in relation.annotations if a['term'] == 'publications']


def test_reader_splits_the_aliases_and_reads_na_as_missing():
    row = parse(NO_SPECIES)
    assert row['aliases'] == [] and row['species'] is None and row['accession'] is None
    assert parse(BINDING)['aliases'] == ['RN7SK', '7SK']
    assert parse(BINDING)['pmids'] == ['19244621']


def test_reader_rejects_a_changed_layout_and_an_empty_download():
    with pytest.raises(ValueError, match='fields'):
        list(iter_interactions(table('1\t7SK\tABO')))
    with pytest.raises(ValueError, match='no rows'):
        list(iter_interactions(table()))
