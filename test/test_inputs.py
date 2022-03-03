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
    model = biomodels.get_single_model_information(
        "MODEL1508040001", invalidate = False)

    assert model["name"] == 'Tomàs-Gamisans2016 - Genome-Scale Metabolic Model of Pichia pastoris (version 2)'

def test_biomodels_all_models():
    models = biomodels.get_all_models(invalidate = False)

    assert (
        len(models) > 2000 
        and isinstance(models[0], dict) 
        and 'format' in models[0].keys()
        )

def test_biomodels_download():
    """
    MODEL1508040001: is v2 on webpage, download fails because ".2" is missing from URL
    https://www.ebi.ac.uk/biomodels/model/download/MODEL1508040001.2?filename=MODEL1508040001_url.xml
    """
    model = biomodels.get_single_model_information("MODEL1508040001")

    model_id = model.get("publicationId") or model.get("submissionId")
    version = model.get("history").get("revisions")[-1].get("version")

    dl = biomodels.get_single_model_main_file(
            model_id, 
            model["files"]["main"][0]["name"], 
            version = version,
            invalidate = False
        )

    od = dl["sbml"]
    md = od["model"]

    assert (
        isinstance(od, dict)
        and isinstance(md, dict)
        and all(key in md.keys() 
                for key in [
                    "listOfCompartments", 
                    "listOfSpecies", 
                    "annotation"
                ])
    )

def test_biomodels_download2():
    """
    MODEL1509220012: doesn't fail after deleting cache..??
    https://www.ebi.ac.uk/biomodels/model/download/MODEL1508040001.2?filename=MODEL1508040001_url.xml
    """
    model = biomodels.get_single_model_information("MODEL1509220012")

    model_id = model.get("publicationId") or model.get("submissionId")
    version = model.get("history").get("revisions")[-1].get("version")

    dl = biomodels.get_single_model_main_file(
            model_id, 
            model["files"]["main"][0]["name"], 
            version = version,
            invalidate = False
        )

    od = dl["sbml"]
    md = od["model"]

    assert (
        isinstance(od, dict)
        and isinstance(md, dict)
        and all(key in md.keys() 
                for key in [
                    "listOfCompartments", 
                    "listOfSpecies", 
                    "annotation"
                ])
    )

def test_biomodels_parse():
    parse = biomodels.parse_biomodels(invalidate = False)
    