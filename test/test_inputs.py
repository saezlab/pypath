import pypath.inputs.cancerdrugs_db as cancerdrugs_db
import pypath.inputs.cancerdrugsdb as cancerdrugsdb
import pypath.inputs.biomodels as biomodels
from collections import OrderedDict

def test_cancerdrugs_db():
    # test annotations
    a = cancerdrugs_db.cancerdrugs_db_annotations()
    b = a.get('46220502')
    c = next(b.__iter__()).label

    assert 100 < len(a) < 1000
    assert c == 'Abemaciclib'

    # test interactions
    d = cancerdrugs_db.cancerdrugs_db_interactions()
    e = next(d.__iter__())

    assert 500 < len(d) < 20000
    # how to test for specific interactions when they are not 1 to 1?

def test_cancerdrugsdb():
    # test annotations
    a = cancerdrugsdb.cancerdrugsdb_annotations()
    b = a.get('46220502')
    c = next(b.__iter__()).drug_label

    assert 100 < len(a) < 1000
    assert c == 'Abemaciclib'

    # test interactions
    d = cancerdrugsdb.cancerdrugsdb_interactions()

    assert 500 < len(d) < 20000
    # how to test for specific interactions when they are not 1 to 1?

def test_biomodels_one_model():
    model = biomodels.get_single_model_information("BIOMD0000000299")

    assert model["name"] == 'Leloup1999_CircadianRhythms_Neurospora'

def test_biomodels_all_models():
    models = biomodels.get_all_models()

    assert (
        len(models) > 2000 
        and isinstance(models[0], dict) 
        and 'format' in models[0].keys()
        )

def test_biomodels_download():
    model = biomodels.get_single_model_information("BIOMD0000000299")

    dl = biomodels.get_single_model_main_file(
        model["publicationId"], model["files"]["main"][0]["name"])

    od = dl["sbml"]
    md = od["model"]

    assert (
        isinstance(od, OrderedDict)
        and isinstance(md, OrderedDict)
        and all(key in md.keys() 
                for key in [
                    "listOfCompartments", 
                    "listOfSpecies", 
                    "listOfParameters", 
                    "listOfRules", 
                    "annotation"
                ])
    )

    