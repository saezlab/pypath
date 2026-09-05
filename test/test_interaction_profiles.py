"""Small interaction vocabulary: roles, qualifiers and observation preservation."""
import pytest
from test.test_inputs_v2_chemical_migration import project
from omnipath_core.interaction_profiles import interaction_label
from omnipath_core.biolink import annotation_term, annotation_value, predicate


def test_binding_measurement_retains_comparator_without_activity_threshold():
    out = project('bindingdb', 'interactions_schema', {'BindingDB MonomerID':'1', 'UniProt (SwissProt) Primary ID of Target Chain 1':'P00533', 'Ki (nM)':'>10000'})
    r = out.relations[-1]
    assert r.predicate == 'interacts_with'
    assert interaction_label(r.predicate, r.annotations) == 'Binds'
    q = next(a['quantity'] for a in r.annotations if a.get('quantity'))
    assert q['has_numeric_value'] == 10000
    assert '>' in str(q)


@pytest.mark.parametrize('assay,action,label,predicate', [('B','', 'Binds','interacts_with'), ('F','',None,'interacts_with'), ('F','INHIBITOR',None,'affects')])
def test_chembl_distinguishes_binding_from_functional_effect(assay,action,label,predicate):
    out=project('chembl','activities_schema',{'molecule_chembl_id':'CHEMBL1','target_chembl_id':'CHEMBL2','target_type':'SINGLE PROTEIN','assay_type':assay,'action_type':action})
    r=out.relations[-1]
    assert r.predicate == predicate
    assert interaction_label(r.predicate,r.annotations) == label


@pytest.mark.parametrize('role,expected', [('Transporter','Transports'),('Receptor; Transporter','Transports'),('Receptor',None)])
def test_mebocost_transport_is_protein_to_chemical_and_keeps_roles(role,expected):
    from pypath.inputs_v2.mebocost import get_interactions_schema
    r=get_interactions_schema('9606')({'ID':'1','HMDB_ID':'HMDB0000001','Gene_name':'ABCA1','Annotation':role})
    annotations=[{'term':annotation_term(a.term),'value':annotation_value(a.term,a.value)} for a in r.annotations]
    assert interaction_label(r.predicate,annotations)==expected
    assert {a['value'] for a in annotations if a['term']=='original_predicate'} == set(role.split('; '))
    if expected:
        assert r.subject.type == 'protein'
        assert r.object.type == 'chemical_entity'
    assert not any(annotation_term(a.term).endswith('_qualifier') for e in [r.subject,r.object] for a in e.annotations or [])


@pytest.mark.parametrize('module,mapper,row', [
 ('tcdb','_transport_schema',{'transporter_uniprot':'P00533','substrate_id':'CHEBI:1','tcid':'1.A.1.1.1'}),
 ('mrclinksdb','human_transporters_schema',{'uniprot_id':'P00533','hmdb_id':'HMDB0000001','taxon_id':'9606'}),
])
def test_transport_qualifiers_survive_extraction(module,mapper,row):
    out=project(module,mapper,row);r=out.relations[-1]
    assert interaction_label(r.predicate,r.annotations)=='Transports'
    assert out.entities[r.subject_entity_key].entity_type=='protein'
    assert out.entities[r.object_entity_key].entity_type=='chemical_entity'


def test_cellinker_enzyme_role_is_not_transport():
    from pypath.inputs_v2.cellinker import _make_metabolite_protein_mapper
    row={'hmdb_id':'HMDB0000001','HMDB_ID':'HMDB0000001','HMDB ID':'HMDB0000001','UNIPROT_ID':'P00533','uniprot_id':'P00533','METABOLITE_NAME':'X','ENZYME_NAME':'Y'}
    for role in ['enzyme','transporter']:
        r=_make_metabolite_protein_mapper('human',role)(row)
        assert r is not None
        assert predicate(r.predicate) == ('affects' if role=='transporter' else 'associated_with')
