"""Build dispatch honors caller network configuration without importing a builder."""

import os
import sys
from types import ModuleType

from omnipath_api.jobs import build_ops


def test_worker_dispatch_preserves_proxy_environment(monkeypatch):
    monkeypatch.setenv("__CURSOR_SANDBOX_ENV_RESTORE", "fixture")
    monkeypatch.setenv("HTTPS_PROXY", "http://fixture.invalid:1234")
    package = ModuleType("omnipath_build")
    package.__path__ = []
    pipeline = ModuleType("omnipath_build.pipeline")
    cachedir = ModuleType("omnipath_build.cachedir_compat")
    pipeline.build_all = lambda **kwargs: dict(proxy=os.environ["HTTPS_PROXY"], **kwargs)
    cachedir.patch_cachedir_opener = lambda: None
    monkeypatch.setitem(sys.modules, "omnipath_build", package)
    monkeypatch.setitem(sys.modules, "omnipath_build.pipeline", pipeline)
    monkeypatch.setitem(sys.modules, "omnipath_build.cachedir_compat", cachedir)
    monkeypatch.setattr(build_ops, "_import_build", lambda: package)
    result = build_ops.BuildOps().build_all(max_records=20)
    assert result == {"proxy": "http://fixture.invalid:1234", "max_records": 20}
    assert os.environ["HTTPS_PROXY"] == "http://fixture.invalid:1234"
