import pypath.share.settings as settings


def test_get():
    res = settings.get("module_name")
    assert res == "pypath"

def test_compare_to_old():
    old_settings = settings._defaults_old.keys()
    new_settings = settings._defaults.keys()

    assert all(elem in old_settings for elem in new_settings)
