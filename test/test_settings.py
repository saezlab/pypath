import pypath.share.settings as settings


def test_get():
    res = settings.get("module_name")
    assert res == "pypath"
