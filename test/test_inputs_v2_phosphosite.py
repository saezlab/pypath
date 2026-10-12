"""PhosphoSitePlus kinase-substrate sites: readers, join, mapper and extractor on source rows."""

import io
from types import SimpleNamespace

import pytest
from omnipath_build.silver import SilverExtractor
from omnipath_core.biolink import slot_name
from omnipath_core.naming import Namespace

from pypath.inputs_v2 import phosphosite
from pypath.inputs_v2.parsers import phosphosite as phosphosite_parser

# The lines above the table, as in the source file (license statement shortened
# after its first sentences).
PREAMBLE = (
    '050522\n'
    'Data extracted from PhosphoSitePlus(R), created by Cell Signaling '
    'Technology Inc. PhosphoSitePlus is licensed under a Creative Commons '
    'Attribution-NonCommercial-ShareAlike 3.0 Unported License. Attribution '
    'must be given in written, oral and digital presentations to '
    'PhosphoSitePlus, www.phosphosite.org.\n'
    '\n'
)
HEADER = '\t'.join(phosphosite_parser.COLUMNS) + '\n'
# Rows copied from the source file.
ROWS = (
    # 0: the first row of the file, an in vitro reaction of mouse proteins
    'Dyrk2\tDYRK2\tQ5U4C9\tmouse\tNDEL1\t83431\tQ9ERR1\tNdel1\tmouse\tS336\t1869686801\tLGSsRPSsAPGMLPL\t\t \tX\t\n'
    # 1: seen in vivo and in vitro; a human kinase on a rat substrate
    'PRKCD\tPKCD\tQ05655\thuman\tsyndecan-4\t24771\tP34901\tSdc4\trat\tS183\t449863\tMKKKDEGsYDLGKKP\tSyndecan\tX\tX\t\n'
    # 2: in vivo only; a window that starts before the sequence does
    'PRKCD\tPKCD\tQ05655\thuman\tL-plastin\t3936\tP13796\tLCP1\thuman\tS5\t450852\t___MARGsVsDEEMM\t\tX\t \t\n'
    # 3: a substrate isoform
    'CAMK2A\tCAMK2A\tQ9UQM7\thuman\tTau iso8\t4137\tP10636-8\tMAPT\thuman\tS262\t447930\tNVKsKIGstENLkHQ\tTubulin-binding\t \tX\t\n'
    # 4: a kinase isoform, with antibody catalogue numbers
    'VRK2\tVRK2 iso2\tQ86Y07-2\thuman\tp53\t22059\tP02340\tTp53\tmouse\tS18\t447490\tIsLELPLsQEtFsGL\tP53_TAD\tX\t \t9284; 9286; 12571\n'
    # 5: a kinase with a RefSeq protein accession
    'Rps6ka3\tRSK2\tNP_001178933\trat\tBSN\t29138\tO88778\tBsn\trat\tS2845\t456576\tKTLQRSLsDPKPLSP\t\t \tX\t\n'
    # 6: a substrate with a GenBank protein accession and no gene ID
    'CAMK2A\tCAMK2A\tQ9UQM7\thuman\tCACNA1C\t\tBAA34185.2\tCACNA1C\tguinea pig\tT1603\t7437803\tRTALRIKtEGNLEQA\tCa_chan_IQ\t \tX\t\n'
    # 7: a histidine site
    'NME1\tNME1\tP15531\thuman\tNME1\t4830\tP15531\tNME1\thuman\tH118\t474528\tQVGRNIIhGsDsVEs\tNDK\tX\tX\t\n'
    # 8, 9: organism labels that do not name a species
    'CAMK2A\tCAMK2A\tQ9UQM7\thuman\tPEA-15\t100736554\tQ9Z297\tPea15\thamster\tS116\t450281\tKDIIRQPsEEEIIKL\t\t \tX\t\n'
    'UL97\tUL97\tP16788\therpesvirus\tRb\t5925\tP06400\tRB1\thuman\tS780\t447978\tstrPPtLsPIPHIPr\tRb_C\tX\t \t9307; 3590; 8180\n'
    # 10: a record on a UniProt sequence variant
    'AKT1\tAkt1\tP31749\thuman\tKv11.1 iso5\t3757\tQ12809_VAR_014388\tKCNH2\thuman\tT897\t4839200\tSFRRRtDtDtEQPGE\t\tX\t \t\n'
    # 11: a record on a sequence UniProt has withdrawn
    'CSNK2A1\tCK2A1\tP68400\thuman\tTR-alpha iso5\t21833\tP63058_OBS\tThra\tmouse\tS474\t449309\tASSLSSSssDEDTEV\t\t \tX\t\n'
    # 12: a substrate accession of no recognised kind
    'PRKACA\tPKACA\tP17612\thuman\tH2B\t525512\tHSBO22\tH2BC3\tcow\tS36\t449903\tRKRsRKEsYSVYVYK\tHistone\t \tX\t\n'
    # 13: a kinase with only a name
    '\tCDK5\t\tno information\tHMGB1\t25459\tP63159\tHmgb1\trat\tS181\t20506700\tGVVkAEKsKKKKEEE\t\t \tX\t\n'
)


def parse(text=PREAMBLE + HEADER + ROWS):
    opener = SimpleNamespace(result=io.StringIO(text))
    return list(
        phosphosite_parser.iter_kinase_substrate(
            opener, organisms=phosphosite.ORGANISM_LABELS
        )
    )


def extract(rows):
    out = SilverExtractor('phosphosite', 'enzyme_substrate')
    for index, row in enumerate(rows):
        record = phosphosite.resource.enzyme_substrate.mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def annotations(relation, term, scope='relation'):
    return [
        a['value']
        for a in relation.annotations
        if a['term'] == term and a['scope'] == scope
    ]


def site(out, relation):
    return out.entities[relation.object_entity_key].molecular_form['modifications'][0]


def test_reader_skips_the_preamble_and_reads_site_and_reaction_columns():
    rows = parse()
    assert len(rows) == 14
    first = rows[0]
    assert first['release'] == '050522'
    assert (first['KIN_ACC_ID'], first['SUB_ACC_ID']) == ('Q5U4C9', 'Q9ERR1')
    assert (first['residue'], first['position']) == ('S', 336)
    assert first['SITE_+/-7_AA'] == 'LGSsRPSsAPGMLPL'
    assert [(row['in_vivo'], row['in_vitro']) for row in rows[:3]] == [
        (False, True),
        (True, True),
        (True, False),
    ]
    assert (rows[7]['residue'], rows[7]['position']) == ('H', 118)


def test_reader_rejects_a_changed_file():
    # Another table, an error page, a file without the license statement and
    # an empty table are errors, not empty datasets.
    with pytest.raises(ValueError, match='other columns'):
        parse(PREAMBLE + 'GENE\tKINASE\tKIN_ACC_ID\n')
    with pytest.raises(ValueError, match='column names'):
        parse('<html>\n<head><title>404 Not Found</title></head>\n')
    with pytest.raises(ValueError, match='CC BY-NC-SA'):
        parse('050522\n\n' + HEADER + ROWS)
    with pytest.raises(ValueError, match='no rows'):
        parse(PREAMBLE + HEADER)
    with pytest.raises(ValueError, match='IN_VIVO_RXN'):
        parse(PREAMBLE + HEADER + ROWS.split('\n')[0].replace('\t \tX\t', '\t?\tX\t') + '\n')
    with pytest.raises(ValueError, match='window'):
        parse(PREAMBLE + HEADER + ROWS.split('\n')[0].replace('LGSsRPSsAPGMLPL', 'LGSsRPStAPGMLPL') + '\n')


def test_kinase_increases_substrate_phosphorylation_at_the_site():
    record = phosphosite.resource.enzyme_substrate.mapper(parse()[0])
    assert slot_name(record.predicate) == 'affects'
    assert record.subject.identifiers == [
        (Namespace.UNIPROT, 'Q5U4C9'),
        (Namespace.GENESYMBOL, 'Dyrk2'),
        (Namespace.NAME, 'DYRK2'),
    ]
    assert record.object.identifiers == [
        (Namespace.UNIPROT, 'Q9ERR1'),
        (Namespace.GENESYMBOL, 'Ndel1'),
        (Namespace.ENTREZ, '83431'),
        (Namespace.NAME, 'NDEL1'),
    ]
    assert record.subject.molecular_form is None
    assert record.object.molecular_form['modifications'] == [
        {
            'term': 'phosphorylated residue',
            'residue': 'S',
            'position': 336,
            'end_position': 336,
            'coordinate_reference': {
                'identifier': {'ns': 'uniprot', 'id': 'Q9ERR1'},
                'coordinate_system': 'protein',
                'position_base': 1,
            },
            'description': 'S336 (LGSsRPSsAPGMLPL)',
        }
    ]

    out = extract(parse()[:1])
    (relation,) = out.relations
    assert relation.predicate == 'affects'
    assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']
    assert annotations(relation, 'object_direction_qualifier') == ['increased']
    assert annotations(relation, 'has_evidence_of_type') == ['ECO:0000181']
    assert {e.entity_type for e in out.entities.values()} == {'protein'}
    assert {e.taxon for e in out.entities.values()} == {'10090'}


def test_reaction_columns_become_evidence_types():
    out = extract(parse()[:3])
    assert [annotations(r, 'has_evidence_of_type') for r in out.relations] == [
        ['ECO:0000181'],
        ['ECO:0000178', 'ECO:0000181'],
        ['ECO:0000178'],
    ]


def test_species_are_kept_as_written_without_translation():
    rows = parse()
    # A human kinase and a rat substrate stay a human and a rat protein.
    out = extract([rows[1]])
    (relation,) = out.relations
    assert out.entities[relation.subject_entity_key].taxon == '9606'
    assert out.entities[relation.object_entity_key].taxon == '10116'
    assert out.entities[relation.subject_entity_key].identifier == 'Q05655'
    assert out.entities[relation.object_entity_key].identifier == 'P34901'
    # "hamster" and "herpesvirus" name no single species: no taxon is given.
    hamster = extract([rows[8]])
    assert hamster.entities[hamster.relations[0].object_entity_key].taxon is None
    virus = extract([rows[9]])
    assert virus.entities[virus.relations[0].subject_entity_key].taxon is None
    assert virus.entities[virus.relations[0].object_entity_key].taxon == '9606'
    # A label the crosswalk does not know stops the reader.
    with pytest.raises(ValueError, match='Unknown PhosphoSite organism'):
        parse(PREAMBLE + HEADER + ROWS.split('\n')[0].replace('\tmouse\tS336', '\tzebrafish\tS336') + '\n')


def test_isoforms_and_other_accession_kinds_are_kept():
    rows = parse()
    substrate = phosphosite.resource.enzyme_substrate.mapper(rows[3]).object
    assert substrate.identifiers[0] == (Namespace.UNIPROT, 'P10636-8')
    assert substrate.molecular_form['isoform_identifier'] == {'ns': 'uniprot', 'id': 'P10636-8'}
    # The position is in the isoform's sequence.
    assert substrate.molecular_form['modifications'][0]['coordinate_reference'][
        'identifier'
    ] == {'ns': 'uniprot', 'id': 'P10636-8'}

    kinase = phosphosite.resource.enzyme_substrate.mapper(rows[4]).subject
    assert kinase.identifiers[0] == (Namespace.UNIPROT, 'Q86Y07-2')
    assert kinase.molecular_form['isoform_identifier'] == {'ns': 'uniprot', 'id': 'Q86Y07-2'}
    assert kinase.molecular_form['modifications'] is None

    refseq = phosphosite.resource.enzyme_substrate.mapper(rows[5]).subject
    assert refseq.identifiers[0] == (Namespace.REFSEQ_PROTEIN, 'NP_001178933')

    genbank = phosphosite.resource.enzyme_substrate.mapper(rows[6]).object
    assert genbank.identifiers == [
        (Namespace.GENBANK, 'BAA34185.2'),
        (Namespace.GENESYMBOL, 'CACNA1C'),
        (Namespace.NAME, 'CACNA1C'),
    ]
    assert genbank.molecular_form['modifications'][0]['coordinate_reference'][
        'identifier'
    ] == {'ns': 'genbank', 'id': 'BAA34185.2'}


def test_variant_withdrawn_and_unrecognised_accessions_have_no_reference_sequence():
    rows = parse()
    variant = phosphosite.resource.enzyme_substrate.mapper(rows[10]).object
    assert variant.identifiers[0] == (Namespace.UNIPROT, 'Q12809')
    assert variant.molecular_form['isoform_identifier'] is None
    assert [v['description'] for v in variant.molecular_form['variants']] == ['VAR_014388']
    modification = variant.molecular_form['modifications'][0]
    assert (modification['residue'], modification['position']) == ('T', 897)
    assert modification['coordinate_reference']['identifier'] is None

    withdrawn = phosphosite.resource.enzyme_substrate.mapper(rows[11]).object
    assert withdrawn.identifiers[0] == (Namespace.UNIPROT, 'P63058')
    assert withdrawn.molecular_form['variants'] is None
    assert withdrawn.molecular_form['modifications'][0]['coordinate_reference']['identifier'] is None

    # "HSBO22" is not an accession of a kind the table otherwise uses: the
    # substrate keeps its gene symbol and gene ID, and gets no accession.
    other = phosphosite.resource.enzyme_substrate.mapper(rows[12]).object
    assert other.identifiers == [
        (Namespace.GENESYMBOL, 'H2BC3'),
        (Namespace.ENTREZ, '525512'),
        (Namespace.NAME, 'H2B'),
    ]
    assert other.molecular_form['modifications'][0]['coordinate_reference']['identifier'] is None


def test_sites_of_one_pair_stay_distinct():
    rows = parse()
    second = dict(rows[0], SUB_MOD_RSD='S198', residue='S', position=198)
    out = extract([rows[0], second])
    assert len(out.relations) == 2
    assert out.relations[0].relation_key != out.relations[1].relation_key
    assert {site(out, r)['position'] for r in out.relations} == {336, 198}


def test_no_relation_without_an_identified_kinase_and_no_invented_values():
    rows = parse()
    # The kinase of row 13 has a name but neither accession nor gene symbol.
    assert phosphosite.resource.enzyme_substrate.mapper(rows[13]) is None

    out = extract(rows)
    assert len(out.relations) == 13
    for relation in out.relations:
        # Phosphorylation is not an effect on activity, and rows read from the
        # table alone have no publications.
        assert relation.predicate == 'affects'
        assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']
        assert annotations(relation, 'object_direction_qualifier') == ['increased']
        assert annotations(relation, 'publications') == []
    # The residue letter is the source's, also where it is not S, T or Y.
    assert site(out, out.relations[7])['residue'] == 'H'


def test_source_is_the_rescued_copy_under_its_stated_license():
    download = phosphosite.resource.enzyme_substrate.download
    assert download.url.startswith('https://rescued.omnipathdb.org/phosphosite/')
    assert phosphosite.config.license.name == 'CC_BY_NC_SA_3_0'
    # The BioPAX file of the evidence is the rescued copy of the same release.
    assert phosphosite.evidence_download.url == (
        'https://rescued.omnipathdb.org/phosphosite/Kinase_substrates.owl.gz'
    )


# ---------------------------------------------------------------------------
# BioPAX evidence file. Catalysis elements copied from Kinase_substrates.owl
# (May 2022), trimmed: titles, years, URLs, protein names and the taxonomy
# cross-references are left out, and the objects the file defines once are
# defined at their first use here and referred to by rdf:resource afterwards,
# as in the file. They are:
#   Catalysis_43850_457     PKCZ -> MARK3 T564: four PubMed records and a
#                           PhosphoSitePlus record without PubMed ID; antibody
#                           and mass spectrometry methods
#   Catalysis_6062646_637   PKACA -> rat CACNA1D S121: one PubMed record, no method
#   Catalysis_17538100_832, Catalysis_17538100_595
#                           ERK2 and ERK1 -> Raptor S8: one reaction, so the
#                           same evidence; the second refers to the reaction
#   Catalysis_2566108_7240  PKACA -> MCAM S85 on a RefSeq accession (the table
#                           has the rat protein under a UniProt accession)
#   Catalysis_9730604_1875247300
#                           a kinase without accession (PKCE)
# ---------------------------------------------------------------------------

OWL_HEAD = (
    '<?xml version="1.0"?>\n'
    '<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" '
    'xmlns:bp="http://www.biopax.org/release/biopax-level3.owl#">\n'
)
OWL_CATALYSES = """\
 <bp:Catalysis rdf:ID="Catalysis_43850_457">
  <bp:controller>
   <bp:Protein rdf:ID="Protein_457">
    <bp:entityReference>
     <bp:ProteinReference rdf:ID="po_457">
      <bp:xref>
       <bp:UnificationXref rdf:ID="xref_457">
        <bp:db>UniProtKB</bp:db>
        <bp:id>Q05513</bp:id>
       </bp:UnificationXref>
      </bp:xref>
      <bp:organism>
       <bp:BioSource rdf:ID="organism_human">
        <bp:displayName>human</bp:displayName>
       </bp:BioSource>
      </bp:organism>
     </bp:ProteinReference>
    </bp:entityReference>
   </bp:Protein>
  </bp:controller>
  <bp:controlled>
   <bp:BiochemicalReaction rdf:ID="BiochemicalReaction_43850">
    <bp:evidence>
     <bp:Evidence rdf:ID="evidence_4143_EvidenceCodeVocabulary_0421">
      <bp:xref>
       <bp:PublicationXref rdf:ID="pmid_15084291">
        <bp:id>15084291</bp:id>
        <bp:db>pubmed</bp:db>
       </bp:PublicationXref>
      </bp:xref>
      <bp:evidenceCode>
       <bp:EvidenceCodeVocabulary rdf:ID="EvidenceCodeVocabulary_0421">
        <bp:xref>
         <bp:UnificationXref rdf:ID="method_0421">
          <bp:id>MI:0421</bp:id>
          <bp:db>PSI-MI</bp:db>
         </bp:UnificationXref>
        </bp:xref>
        <bp:term>identification by antibody</bp:term>
       </bp:EvidenceCodeVocabulary>
      </bp:evidenceCode>
     </bp:Evidence>
    </bp:evidence>
    <bp:evidence>
     <bp:Evidence rdf:ID="evidence_21390221_null">
      <bp:xref>
       <bp:PublicationXref rdf:ID="pmid_21983960">
        <bp:id>21983960</bp:id>
        <bp:db>pubmed</bp:db>
       </bp:PublicationXref>
      </bp:xref>
     </bp:Evidence>
    </bp:evidence>
    <bp:evidence>
     <bp:Evidence rdf:ID="evidence_27303810_EvidenceCodeVocabulary_0427">
      <bp:xref>
       <bp:PublicationXref rdf:ID="recid_27303810">
                </bp:PublicationXref>
      </bp:xref>
      <bp:evidenceCode>
       <bp:EvidenceCodeVocabulary rdf:ID="EvidenceCodeVocabulary_0427">
        <bp:xref>
         <bp:UnificationXref rdf:ID="method_0427">
          <bp:id>MI:0427</bp:id>
          <bp:db>PSI-MI</bp:db>
         </bp:UnificationXref>
        </bp:xref>
        <bp:term>Identification by mass spectrometry</bp:term>
       </bp:EvidenceCodeVocabulary>
      </bp:evidenceCode>
     </bp:Evidence>
    </bp:evidence>
    <bp:left>
     <bp:Protein rdf:ID="Protein_721">
      <bp:entityReference>
       <bp:ProteinReference rdf:ID="po_721">
        <bp:xref>
         <bp:UnificationXref rdf:ID="xref_721">
          <bp:db>UniProtKB</bp:db>
          <bp:id>P27448</bp:id>
         </bp:UnificationXref>
        </bp:xref>
        <bp:organism rdf:resource="#organism_human"/>
       </bp:ProteinReference>
      </bp:entityReference>
     </bp:Protein>
    </bp:left>
    <bp:right>
     <bp:Protein rdf:ID="ModifiedProtein_43850">
      <bp:feature>
       <bp:ModificationFeature rdf:ID="ModificationFeature_T564p">
        <bp:featureLocation>
         <bp:SequenceSite rdf:ID="SequenceSite_564">
          <bp:sequencePosition>564</bp:sequencePosition>
          <bp:positionStatus>EQUAL</bp:positionStatus>
         </bp:SequenceSite>
        </bp:featureLocation>
        <bp:modificationType>
         <bp:SequenceModificationVocabulary rdf:ID="SequenceModificationVocabulary_T-p">
          <bp:xref>
           <bp:UnificationXref rdf:ID="UnificationXref_T-p">
            <bp:id>MOD:00047</bp:id>
            <bp:db>PSI-MI</bp:db>
           </bp:UnificationXref>
          </bp:xref>
          <bp:term>O-phospho-L-threonine</bp:term>
         </bp:SequenceModificationVocabulary>
        </bp:modificationType>
       </bp:ModificationFeature>
      </bp:feature>
      <bp:entityReference rdf:resource="#po_721"/>
     </bp:Protein>
    </bp:right>
    <bp:evidence>
     <bp:Evidence rdf:ID="evidence_34666515_EvidenceCodeVocabulary_0427">
      <bp:xref>
       <bp:PublicationXref rdf:ID="pmid_24719451">
        <bp:id>24719451</bp:id>
        <bp:db>pubmed</bp:db>
       </bp:PublicationXref>
      </bp:xref>
      <bp:evidenceCode rdf:resource="#EvidenceCodeVocabulary_0427"/>
     </bp:Evidence>
    </bp:evidence>
    <bp:evidence>
     <bp:Evidence rdf:ID="evidence_27423500_null">
      <bp:xref>
       <bp:PublicationXref rdf:ID="pmid_23312004">
        <bp:id>23312004</bp:id>
        <bp:db>pubmed</bp:db>
       </bp:PublicationXref>
      </bp:xref>
     </bp:Evidence>
    </bp:evidence>
   </bp:BiochemicalReaction>
  </bp:controlled>
 </bp:Catalysis>
 <bp:Catalysis rdf:ID="Catalysis_6062646_637">
  <bp:controller>
   <bp:Protein rdf:ID="Protein_637">
    <bp:entityReference>
     <bp:ProteinReference rdf:ID="po_637">
      <bp:xref>
       <bp:UnificationXref rdf:ID="xref_637">
        <bp:db>UniProtKB</bp:db>
        <bp:id>P17612</bp:id>
       </bp:UnificationXref>
      </bp:xref>
      <bp:organism rdf:resource="#organism_human"/>
     </bp:ProteinReference>
    </bp:entityReference>
   </bp:Protein>
  </bp:controller>
  <bp:controlled>
   <bp:BiochemicalReaction rdf:ID="BiochemicalReaction_6062646">
    <bp:evidence>
     <bp:Evidence rdf:ID="evidence_5930025_null">
      <bp:xref>
       <bp:PublicationXref rdf:ID="pmid_19074150">
        <bp:id>19074150</bp:id>
        <bp:db>pubmed</bp:db>
       </bp:PublicationXref>
      </bp:xref>
     </bp:Evidence>
    </bp:evidence>
    <bp:right>
     <bp:Protein rdf:ID="ModifiedProtein_6062646">
      <bp:feature>
       <bp:ModificationFeature rdf:ID="ModificationFeature_S121p">
        <bp:featureLocation>
         <bp:SequenceSite rdf:ID="SequenceSite_121">
          <bp:sequencePosition>121</bp:sequencePosition>
          <bp:positionStatus>EQUAL</bp:positionStatus>
         </bp:SequenceSite>
        </bp:featureLocation>
        <bp:modificationType>
         <bp:SequenceModificationVocabulary rdf:ID="SequenceModificationVocabulary_S-p">
          <bp:xref>
           <bp:UnificationXref rdf:ID="UnificationXref_S-p">
            <bp:id>MOD:00046</bp:id>
            <bp:db>PSI-MI</bp:db>
           </bp:UnificationXref>
          </bp:xref>
          <bp:term>O-phospho-L-serine</bp:term>
         </bp:SequenceModificationVocabulary>
        </bp:modificationType>
       </bp:ModificationFeature>
      </bp:feature>
      <bp:entityReference>
       <bp:ProteinReference rdf:ID="po_15477">
        <bp:xref>
         <bp:UnificationXref rdf:ID="xref_15477">
          <bp:db>UniProtKB</bp:db>
          <bp:id>P27732</bp:id>
         </bp:UnificationXref>
        </bp:xref>
        <bp:organism>
         <bp:BioSource rdf:ID="organism_rat">
          <bp:displayName>rat</bp:displayName>
         </bp:BioSource>
        </bp:organism>
       </bp:ProteinReference>
      </bp:entityReference>
     </bp:Protein>
    </bp:right>
    <bp:left>
     <bp:Protein rdf:ID="Protein_15477">
      <bp:entityReference rdf:resource="#po_15477"/>
     </bp:Protein>
    </bp:left>
   </bp:BiochemicalReaction>
  </bp:controlled>
 </bp:Catalysis>
 <bp:Catalysis rdf:ID="Catalysis_17538100_832">
  <bp:controller>
   <bp:Protein rdf:ID="Protein_832">
    <bp:entityReference>
     <bp:ProteinReference rdf:ID="po_832">
      <bp:xref>
       <bp:UnificationXref rdf:ID="xref_832">
        <bp:db>UniProtKB</bp:db>
        <bp:id>P28482</bp:id>
       </bp:UnificationXref>
      </bp:xref>
      <bp:organism rdf:resource="#organism_human"/>
     </bp:ProteinReference>
    </bp:entityReference>
   </bp:Protein>
  </bp:controller>
  <bp:controlled>
   <bp:BiochemicalReaction rdf:ID="BiochemicalReaction_17538100">
    <bp:evidence>
     <bp:Evidence rdf:ID="evidence_17431523_EvidenceCodeVocabulary_0074">
      <bp:xref>
       <bp:PublicationXref rdf:ID="pmid_21071439">
        <bp:id>21071439</bp:id>
        <bp:db>pubmed</bp:db>
       </bp:PublicationXref>
      </bp:xref>
      <bp:evidenceCode>
       <bp:EvidenceCodeVocabulary rdf:ID="EvidenceCodeVocabulary_0074">
        <bp:xref>
         <bp:UnificationXref rdf:ID="method_0074">
          <bp:id>MI:0074</bp:id>
          <bp:db>PSI-MI</bp:db>
         </bp:UnificationXref>
        </bp:xref>
        <bp:term>mutation analysis</bp:term>
       </bp:EvidenceCodeVocabulary>
      </bp:evidenceCode>
     </bp:Evidence>
    </bp:evidence>
    <bp:right>
     <bp:Protein rdf:ID="ModifiedProtein_17538100">
      <bp:feature>
       <bp:ModificationFeature rdf:ID="ModificationFeature_S8p">
        <bp:featureLocation>
         <bp:SequenceSite rdf:ID="SequenceSite_8">
          <bp:sequencePosition>8</bp:sequencePosition>
          <bp:positionStatus>EQUAL</bp:positionStatus>
         </bp:SequenceSite>
        </bp:featureLocation>
        <bp:modificationType rdf:resource="#SequenceModificationVocabulary_S-p"/>
       </bp:ModificationFeature>
      </bp:feature>
      <bp:entityReference>
       <bp:ProteinReference rdf:ID="po_6245">
        <bp:xref>
         <bp:UnificationXref rdf:ID="xref_6245">
          <bp:db>UniProtKB</bp:db>
          <bp:id>Q8N122</bp:id>
         </bp:UnificationXref>
        </bp:xref>
        <bp:organism rdf:resource="#organism_human"/>
       </bp:ProteinReference>
      </bp:entityReference>
     </bp:Protein>
    </bp:right>
    <bp:left>
     <bp:Protein rdf:ID="Protein_6245">
      <bp:entityReference rdf:resource="#po_6245"/>
     </bp:Protein>
    </bp:left>
   </bp:BiochemicalReaction>
  </bp:controlled>
 </bp:Catalysis>
 <bp:Catalysis rdf:ID="Catalysis_17538100_595">
  <bp:controller>
   <bp:Protein rdf:ID="Protein_595">
    <bp:entityReference>
     <bp:ProteinReference rdf:ID="po_595">
      <bp:xref>
       <bp:UnificationXref rdf:ID="xref_595">
        <bp:db>UniProtKB</bp:db>
        <bp:id>P27361</bp:id>
       </bp:UnificationXref>
      </bp:xref>
      <bp:organism rdf:resource="#organism_human"/>
     </bp:ProteinReference>
    </bp:entityReference>
   </bp:Protein>
  </bp:controller>
  <bp:controlled rdf:resource="#BiochemicalReaction_17538100"/>
 </bp:Catalysis>
 <bp:Catalysis rdf:ID="Catalysis_2566108_7240">
  <bp:controller>
   <bp:Protein rdf:ID="Protein_7240">
    <bp:entityReference>
     <bp:ProteinReference rdf:ID="po_7240">
      <bp:xref>
       <bp:UnificationXref rdf:ID="xref_7240">
        <bp:db>UniProtKB</bp:db>
        <bp:id>P00517</bp:id>
       </bp:UnificationXref>
      </bp:xref>
      <bp:organism>
       <bp:BioSource rdf:ID="organism_cow">
        <bp:displayName>cow</bp:displayName>
       </bp:BioSource>
      </bp:organism>
     </bp:ProteinReference>
    </bp:entityReference>
   </bp:Protein>
  </bp:controller>
  <bp:controlled>
   <bp:BiochemicalReaction rdf:ID="BiochemicalReaction_2566108">
    <bp:evidence>
     <bp:Evidence rdf:ID="evidence_11906_EvidenceCodeVocabulary_0427">
      <bp:xref>
       <bp:PublicationXref rdf:ID="pmid_17564427">
        <bp:id>17564427</bp:id>
        <bp:db>pubmed</bp:db>
       </bp:PublicationXref>
      </bp:xref>
      <bp:evidenceCode rdf:resource="#EvidenceCodeVocabulary_0427"/>
     </bp:Evidence>
    </bp:evidence>
    <bp:right>
     <bp:Protein rdf:ID="ModifiedProtein_2566108">
      <bp:feature>
       <bp:ModificationFeature rdf:ID="ModificationFeature_S85p">
        <bp:featureLocation>
         <bp:SequenceSite rdf:ID="SequenceSite_85">
          <bp:sequencePosition>85</bp:sequencePosition>
          <bp:positionStatus>EQUAL</bp:positionStatus>
         </bp:SequenceSite>
        </bp:featureLocation>
        <bp:modificationType rdf:resource="#SequenceModificationVocabulary_S-p"/>
       </bp:ModificationFeature>
      </bp:feature>
      <bp:entityReference>
       <bp:ProteinReference rdf:ID="po_23452">
        <bp:xref>
         <bp:UnificationXref rdf:ID="xref_23452">
          <bp:db>RefSeq</bp:db>
          <bp:id>NP_076473</bp:id>
         </bp:UnificationXref>
        </bp:xref>
        <bp:organism rdf:resource="#organism_rat"/>
       </bp:ProteinReference>
      </bp:entityReference>
     </bp:Protein>
    </bp:right>
    <bp:left>
     <bp:Protein rdf:ID="Protein_23452">
      <bp:entityReference rdf:resource="#po_23452"/>
     </bp:Protein>
    </bp:left>
   </bp:BiochemicalReaction>
  </bp:controlled>
 </bp:Catalysis>
 <bp:Catalysis rdf:ID="Catalysis_9730604_1875247300">
  <bp:controller>
   <bp:Protein rdf:ID="Protein_1875247300">
    <bp:entityReference>
     <bp:ProteinReference rdf:ID="po_1875247300">
      <bp:xref>
       <bp:UnificationXref rdf:ID="xref_1875247300">
        <bp:db/>
        <bp:id/>
       </bp:UnificationXref>
      </bp:xref>
      <bp:organism>
       <bp:BioSource rdf:ID="organism_eukaryotes">
        <bp:displayName>eukaryotes</bp:displayName>
       </bp:BioSource>
      </bp:organism>
     </bp:ProteinReference>
    </bp:entityReference>
   </bp:Protein>
  </bp:controller>
  <bp:controlled>
   <bp:BiochemicalReaction rdf:ID="BiochemicalReaction_9730604">
    <bp:evidence>
     <bp:Evidence rdf:ID="evidence_28118862_null">
      <bp:xref>
       <bp:PublicationXref rdf:ID="pmid_23804652">
        <bp:id>23804652</bp:id>
        <bp:db>pubmed</bp:db>
       </bp:PublicationXref>
      </bp:xref>
     </bp:Evidence>
    </bp:evidence>
    <bp:evidence>
     <bp:Evidence rdf:ID="evidence_9689605_EvidenceCodeVocabulary_0074">
      <bp:xref>
       <bp:PublicationXref rdf:ID="pmid_19625451">
        <bp:id>19625451</bp:id>
        <bp:db>pubmed</bp:db>
       </bp:PublicationXref>
      </bp:xref>
      <bp:evidenceCode rdf:resource="#EvidenceCodeVocabulary_0074"/>
     </bp:Evidence>
    </bp:evidence>
    <bp:evidence>
     <bp:Evidence rdf:ID="evidence_9689605_EvidenceCodeVocabulary_0421">
      <bp:xref rdf:resource="#pmid_19625451"/>
      <bp:evidenceCode rdf:resource="#EvidenceCodeVocabulary_0421"/>
     </bp:Evidence>
    </bp:evidence>
    <bp:right>
     <bp:Protein rdf:ID="ModifiedProtein_9730604">
      <bp:feature>
       <bp:ModificationFeature rdf:ID="ModificationFeature_S369p">
        <bp:featureLocation>
         <bp:SequenceSite rdf:ID="SequenceSite_369">
          <bp:sequencePosition>369</bp:sequencePosition>
          <bp:positionStatus>EQUAL</bp:positionStatus>
         </bp:SequenceSite>
        </bp:featureLocation>
        <bp:modificationType rdf:resource="#SequenceModificationVocabulary_S-p"/>
       </bp:ModificationFeature>
      </bp:feature>
      <bp:entityReference>
       <bp:ProteinReference rdf:ID="po_9731400">
        <bp:xref>
         <bp:UnificationXref rdf:ID="xref_9731400">
          <bp:db>UniProtKB</bp:db>
          <bp:id>Q95168</bp:id>
         </bp:UnificationXref>
        </bp:xref>
        <bp:organism>
         <bp:BioSource rdf:ID="organism_dog">
          <bp:displayName>dog</bp:displayName>
         </bp:BioSource>
        </bp:organism>
       </bp:ProteinReference>
      </bp:entityReference>
     </bp:Protein>
    </bp:right>
    <bp:left>
     <bp:Protein rdf:ID="Protein_9731400">
      <bp:entityReference rdf:resource="#po_9731400"/>
     </bp:Protein>
    </bp:left>
   </bp:BiochemicalReaction>
  </bp:controlled>
 </bp:Catalysis>

"""
OWL = OWL_HEAD + OWL_CATALYSES + '</rdf:RDF>\n'

# Rows of the table for the same sites, copied from the source file.
JOIN_ROWS = (
    # 0: in the BioPAX file, four PubMed records, two methods
    'PRKCZ\tPKCZ\tQ05513\thuman\tMARK3\t4140\tP27448\tMARK3\thuman\tT564\t466001\tRGTASRStFHGQPRE\t\tX\tX\t\n'
    # 1: one PubMed record, no method
    'PRKACA\tPKACA\tP17612\thuman\tCACNA1D\t29716\tP27732\tCacna1d\trat\tS121\t6062715\tPIRRACIsIVDWKPF\t\t \tX\t\n'
    # 2, 3: two kinases of one site
    'MAPK3\tERK1\tP27361\thuman\tRaptor\t57521\tQ8N122\tRPTOR\thuman\tS8\t17538200\tMESEMLQsPLLGLGE\t\tX\tX\t\n'
    'MAPK1\tERK2\tP28482\thuman\tRaptor\t57521\tQ8N122\tRPTOR\thuman\tS8\t17538200\tMESEMLQsPLLGLGE\t\tX\tX\t\n'
    # 4: the substrate has another accession than in the BioPAX file
    'PRKACA\tPKACA\tP00517\tcow\tMCAM\t78967\tQ9EPF2\tMcam\trat\tS85\t2566107\tVHQGKGQsEPGEYEH\tV-set\t \tX\t\n'
    # 5: a kinase without accession, as in the BioPAX file
    '\tPKCE\t\tno information\tZO2\t403854\tQ95168\tTJP2\tdog\tS369\t9730701\tQLVVLRDsKQTLINI\t\tX\t \t\n'
    # 6: no catalysis in the BioPAX file
    'PRKCD\tPKCD\tQ05655\thuman\tCDK5\t1020\tQ00535\tCDK5\thuman\tT77\t900975700\tLHsDkKLtLVFEFCD\tPkinase\tX\t \t\n'
)


def opener_of(text):
    return SimpleNamespace(result=io.StringIO(text))


def catalyses(text=OWL):
    return {
        record['catalysis']: record
        for record in phosphosite_parser.iter_catalysis_evidence(
            opener_of(text), methods=phosphosite.METHODS
        )
    }


def joined(min_matched=0.5, owl=OWL):
    evidence = SimpleNamespace(open=lambda force_refresh=False: opener_of(owl))
    return list(
        phosphosite_parser.iter_kinase_substrate_evidence(
            opener_of(PREAMBLE + HEADER + JOIN_ROWS),
            organisms=phosphosite.ORGANISM_LABELS,
            methods=phosphosite.METHODS,
            evidence_download=evidence,
            min_matched=min_matched,
        )
    )


def test_biopax_reader_follows_references_to_publications_and_methods():
    records = catalyses()
    assert len(records) == 6
    several = records['Catalysis_43850_457']
    assert (several['kinase'], several['substrate']) == ('Q05513', 'P27448')
    assert (several['residue'], several['position']) == ('T', 564)
    # Four PubMed records in file order; the PhosphoSitePlus record without a
    # PubMed ID ("recid_27303810") is counted and not carried. The mass
    # spectrometry evidence code of the last two is a reference to a term
    # defined in the first evidence.
    assert several['publications'] == ['15084291', '21983960', '24719451', '23312004']
    assert several['methods'] == ['MI:0421', 'MI:0427']
    assert several['records_without_pmid'] == 1
    assert (several['kinase_organism'], several['substrate_organism']) == ('human', 'human')

    # A PubMed record without an evidence code gives a publication, no method.
    no_method = records['Catalysis_6062646_637']
    assert no_method['publications'] == ['19074150']
    assert no_method['methods'] == []
    assert no_method['substrate_organism'] == 'rat'

    # The accession is as the file writes it, with its database.
    refseq = records['Catalysis_2566108_7240']
    assert (refseq['substrate'], refseq['substrate_database']) == ('NP_076473', 'RefSeq')
    assert refseq['kinase_database'] == 'UniProtKB'

    # A kinase reference without an accession.
    nameless = records['Catalysis_9730604_1875247300']
    assert (nameless['kinase'], nameless['kinase_database']) == (None, None)
    assert nameless['substrate'] == 'Q95168'


def test_catalyses_of_one_site_share_the_evidence_of_its_reaction():
    records = catalyses()
    first = records['Catalysis_17538100_832']
    second = records['Catalysis_17538100_595']
    # The second catalysis refers to the reaction defined in the first.
    assert (first['kinase'], second['kinase']) == ('P28482', 'P27361')
    assert first['publications'] == second['publications'] == ['21071439']
    assert first['methods'] == second['methods'] == ['MI:0074']
    assert (second['substrate'], second['residue'], second['position']) == ('Q8N122', 'S', 8)


def test_biopax_reader_rejects_what_it_does_not_know():
    # A method outside the crosswalk, or with another label.
    with pytest.raises(ValueError, match='method not understood'):
        catalyses(OWL.replace('MI:0421', 'MI:0018'))
    with pytest.raises(ValueError, match='method not understood'):
        catalyses(OWL.replace('identification by antibody', 'two hybrid'))
    # A publication from another database.
    with pytest.raises(ValueError, match='database'):
        catalyses(OWL.replace('<bp:db>pubmed</bp:db>', '<bp:db>doi</bp:db>', 1))
    # A site whose residue, position and modification term disagree.
    with pytest.raises(ValueError, match='disagree'):
        catalyses(OWL.replace('<bp:sequencePosition>564<', '<bp:sequencePosition>565<'))
    with pytest.raises(ValueError, match='disagree'):
        catalyses(OWL.replace('MOD:00047', 'MOD:00046'))
    # A reference to an object the file does not define.
    with pytest.raises(ValueError, match='unresolved'):
        catalyses(OWL.replace('#po_721', '#po_0'))
    # No catalysis at all.
    with pytest.raises(ValueError, match='no catalysis'):
        catalyses(OWL_HEAD + '</rdf:RDF>\n')


def test_only_well_formed_pubmed_ids_are_carried():
    text = OWL.replace('<bp:id>19074150<', '<bp:id>0019074150<')
    assert catalyses(text)['Catalysis_6062646_637']['publications'] == []
    text = OWL.replace('<bp:id>19074150<', '<bp:id>190741500<')
    assert catalyses(text)['Catalysis_6062646_637']['publications'] == []
    # PubMed has IDs of five and six digits, and the file holds some.
    text = OWL.replace('<bp:id>19074150<', '<bp:id>208844<')
    assert catalyses(text)['Catalysis_6062646_637']['publications'] == ['208844']


def test_rows_take_the_publications_and_methods_of_their_catalysis():
    rows = joined()
    assert len(rows) == 7
    by_site = {(r['KIN_ACC_ID'], r['SUB_ACC_ID'], r['position']): r for r in rows}
    several = by_site[('Q05513', 'P27448', 564)]
    assert several['evidence_record'] is True
    assert several['publications'] == ['15084291', '21983960', '24719451', '23312004']
    assert several['evidence_methods'] == ['MI:0421', 'MI:0427']
    # The columns of the table are untouched.
    assert (several['in_vivo'], several['in_vitro']) == (True, True)
    assert several['SITE_+/-7_AA'] == 'RGTASRStFHGQPRE'
    no_method = by_site[('P17612', 'P27732', 121)]
    assert (no_method['publications'], no_method['evidence_methods']) == (['19074150'], [])
    # Two kinases of a site carry the same evidence.
    erk1, erk2 = by_site[('P27361', 'Q8N122', 8)], by_site[('P28482', 'Q8N122', 8)]
    assert erk1['publications'] == erk2['publications'] == ['21071439']


def test_rows_without_a_catalysis_and_catalyses_without_a_row_are_left_out():
    rows = joined()
    by_site = {(r['KIN_ACC_ID'], r['SUB_ACC_ID'], r['position']): r for r in rows}
    # The table has the rat MCAM under a UniProt accession, the BioPAX file
    # under a RefSeq accession: no match is made through the gene name.
    mcam = by_site[('P00517', 'Q9EPF2', 85)]
    assert mcam['evidence_record'] is False
    assert (mcam['publications'], mcam['evidence_methods']) == ([], [])
    # A row and a catalysis that both lack the kinase accession do not match.
    nameless = by_site[('', 'Q95168', 369)]
    assert nameless['evidence_record'] is False
    assert nameless['publications'] == []
    # A site the BioPAX file does not have.
    assert by_site[('Q05655', 'Q00535', 77)]['evidence_record'] is False
    # The catalysis on NP_076473 and the one without kinase accession add no row.
    assert {r['SUB_ACC_ID'] for r in rows} == {
        'P27448', 'P27732', 'Q8N122', 'Q9EPF2', 'Q95168', 'Q00535'
    }
    assert sum(r['evidence_record'] for r in rows) == 4


def test_join_stops_when_the_files_do_not_belong_together():
    # Of the six catalyses of the fixture four match: below the default share.
    with pytest.raises(ValueError, match='different releases'):
        joined(min_matched=0.95)
    # A site that two catalyses of the file claim.
    twice = list(catalyses().values())
    twice.append(dict(twice[0], catalysis='Catalysis_x_x'))
    with pytest.raises(ValueError, match='two catalysis records'):
        list(phosphosite_parser.join_evidence([], twice))


def mapped(rows):
    out = SilverExtractor('phosphosite', 'enzyme_substrate')
    for index, row in enumerate(rows):
        record = phosphosite.resource.enzyme_substrate.mapper(row)
        if record is not None:
            out.process_record(record, row, f'fixture:{index}', index)
    return out


def test_relations_carry_publications_and_methods_next_to_the_table_evidence():
    out = mapped(joined())
    # Row 5 has no kinase identity and gives no relation.
    assert len(out.relations) == 6
    several, no_method, erk1, erk2, mcam, cdk5 = out.relations
    # The ECO terms of the table come first, then the PSI-MI methods.
    assert annotations(several, 'has_evidence_of_type') == [
        'ECO:0000178', 'ECO:0000181', 'MI:0421', 'MI:0427'
    ]
    assert annotations(several, 'publications') == [
        'PMID:15084291', 'PMID:21983960', 'PMID:24719451', 'PMID:23312004'
    ]
    assert annotations(no_method, 'has_evidence_of_type') == ['ECO:0000181']
    assert annotations(no_method, 'publications') == ['PMID:19074150']
    assert annotations(erk1, 'publications') == annotations(erk2, 'publications') == ['PMID:21071439']
    assert annotations(erk1, 'has_evidence_of_type') == ['ECO:0000178', 'ECO:0000181', 'MI:0074']
    for relation in (several, no_method, erk1, erk2):
        assert relation.predicate == 'affects'
        assert annotations(relation, 'object_aspect_qualifier') == ['phosphorylation']
        assert annotations(relation, 'object_direction_qualifier') == ['increased']


def test_relation_without_a_catalysis_keeps_what_the_table_gave_and_gets_no_evidence():
    out = mapped(joined())
    mcam, cdk5 = out.relations[4], out.relations[5]
    # No publication and no PSI-MI term invented for a row without catalysis.
    assert annotations(mcam, 'publications') == []
    assert annotations(mcam, 'has_evidence_of_type') == ['ECO:0000181']
    assert annotations(cdk5, 'publications') == []
    assert annotations(cdk5, 'has_evidence_of_type') == ['ECO:0000178']
    # The site is as before.
    assert site(out, mcam)['position'] == 85


def test_joined_relations_equal_table_only_relations_apart_from_the_evidence():
    table_only = mapped(parse(PREAMBLE + HEADER + JOIN_ROWS))
    with_evidence = mapped(joined())
    assert len(table_only.relations) == len(with_evidence.relations)
    for before, after in zip(table_only.relations, with_evidence.relations):
        assert before.relation_key == after.relation_key
        assert (before.subject_entity_key, before.object_entity_key) == (
            after.subject_entity_key,
            after.object_entity_key,
        )
        assert [a for a in after.annotations if a['term'] not in ('publications', 'has_evidence_of_type')] == [
            a for a in before.annotations if a['term'] != 'has_evidence_of_type'
        ]
        table_terms = annotations(before, 'has_evidence_of_type')
        assert annotations(after, 'has_evidence_of_type')[: len(table_terms)] == table_terms
    assert table_only.entities.keys() == with_evidence.entities.keys()
