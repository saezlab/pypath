import pypath.inputs.cancerdrugs_db as cancerdrugs_db

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