"""Importable mapper fixture shared by production multiprocessing tests."""

import importlib
import os
import shutil
import sys
import uuid

import pytest

from library_fixture import SHARED_KEY_ENV, shared_root


def pytest_configure(config):
    """Under xdist, let every worker share one reference template per variant.

    Exported before workers are spawned so they inherit it; see library_fixture.
    """
    if (
        not hasattr(config, "workerinput")
        and getattr(config.option, "numprocesses", None)
        and SHARED_KEY_ENV not in os.environ
    ):
        os.environ[SHARED_KEY_ENV] = uuid.uuid4().hex
        config._omnipath_shared_key = True


def pytest_unconfigure(config):
    if getattr(config, "_omnipath_shared_key", False):
        root = shared_root()
        os.environ.pop(SHARED_KEY_ENV, None)
        if root is not None:
            shutil.rmtree(root, ignore_errors=True)


@pytest.fixture
def mapper_module(tmp_path, monkeypatch):
    name = "pipeline_mapper_fixture"
    (tmp_path / f"{name}.py").write_text("""
def mapper(row):
    if 'original' in row:
        row['original'] = 'mutated'
    if 'id' in row:
        return {'type':'named_thing','identifiers':[{'type':'fixture','value':row['id']}]}
    def entity(identifier):
        return {'type':'protein','identifiers':[{'type':'uniprot','value':identifier}]}
    return {'subject':entity('P1'),'predicate':'affects','object':entity('P2')}
class Input:
    mapper = staticmethod(mapper)
items = first = second = repeated = Input()
""")
    monkeypatch.syspath_prepend(str(tmp_path))
    sys.modules.pop(name, None)
    return importlib.import_module(name)
