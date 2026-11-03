"""Installed-file provenance identifies same-version runtime/vocabulary changes."""

import shutil

import pytest

from omnipath_build import runtime_files
from omnipath_build.provenance import runtime_provenance


@pytest.fixture
def installed_packages(tmp_path, monkeypatch):
    roots = {}
    # A wheel-like site-packages layout with no Git or builder checkout.
    for module in runtime_files._SOFTWARE.values():
        root = tmp_path / "site-packages" / module
        root.mkdir(parents=True)
        (root / "__init__.py").write_text("# installed runtime\n")
        roots[module] = root
    for module, filename in (
        ("omnipath_core", "vocab/organisms.yaml"),
        ("pypath", "internals/cv_terms/identifiers.py"),
        ("biolink_model", "schema/biolink_model.yaml"),
    ):
        path = roots[module] / filename
        path.parent.mkdir(parents=True)
        path.write_text("version: exact-vocabulary\n")
    (roots["omnipath_resolver"] / "_omnipath_resolver.abi3.so").write_bytes(b"compiled-native-v1")
    monkeypatch.setattr(runtime_files, "_package_roots", lambda module: (roots[module],))
    monkeypatch.setattr(runtime_files.importlib.metadata, "version", lambda distribution: "0.1.0")
    return roots


def test_installed_only_provenance_covers_code_native_and_all_vocabularies(installed_packages):
    value = runtime_provenance()
    assert value["code_revision"] is None
    assert value["software"]["format"] == "installed-package-content-v1"
    packages = value["software"]["packages"]
    assert set(packages) == set(runtime_files._SOFTWARE)
    assert (
        packages["omnipath-resolver"]["native_files"]["_omnipath_resolver.abi3.so"]["size_bytes"]
        == 18
    )
    for name in ("omnipath-core", "pypath-omnipath", "biolink-model"):
        assert len(packages[name]["vocabulary"]["content_sha256"]) == 64
    assert all(package["git"] is None for package in packages.values())
    assert value["build_code_sha256"] == packages["omnipath-build"]["content_sha256"]


@pytest.mark.parametrize(
    "module,relative",
    [
        ("omnipath_core", "__init__.py"),
        ("omnipath_core", "vocab/organisms.yaml"),
        ("omnipath_resolver", "_omnipath_resolver.abi3.so"),
        ("pypath", "__init__.py"),
        ("pypath", "internals/cv_terms/identifiers.py"),
        ("biolink_model", "schema/biolink_model.yaml"),
    ],
)
def test_same_version_runtime_edits_change_exact_identity(installed_packages, module, relative):
    first = runtime_files.software_provenance()
    path = installed_packages[module] / relative
    path.write_bytes(path.read_bytes() + b"changed-runtime")
    second = runtime_files.software_provenance()
    assert first["sha256"] != second["sha256"]
    distribution = next(
        name for name, target in runtime_files._SOFTWARE.items() if target == module
    )
    assert (
        first["packages"][distribution]["content_sha256"]
        != second["packages"][distribution]["content_sha256"]
    )
    assert second["packages"][distribution]["version"] == "0.1.0"


def test_installation_path_and_bytecode_do_not_change_runtime_identity(
    installed_packages, tmp_path, monkeypatch
):
    first = runtime_files.software_provenance()
    moved = tmp_path / "other-install"
    roots = {}
    for module, path in installed_packages.items():
        roots[module] = moved / module
        shutil.copytree(path, roots[module])
        cache = roots[module] / "__pycache__"
        cache.mkdir()
        (cache / "__init__.cpython-312.pyc").write_bytes(b"machine-local-bytecode")
    monkeypatch.setattr(runtime_files, "_package_roots", lambda module: (roots[module],))
    assert runtime_files.software_provenance()["sha256"] == first["sha256"]


def test_file_cache_reuses_unchanged_content_hashes(installed_packages):
    runtime_files._content_hash.cache_clear()
    first = runtime_files.software_provenance()
    before = runtime_files._content_hash.cache_info()
    assert runtime_files.software_provenance() == first
    after = runtime_files._content_hash.cache_info()
    assert after.misses == before.misses
    assert after.hits > before.hits


def test_required_native_component_cannot_silently_be_unfingerprinted(installed_packages):
    (installed_packages["omnipath_resolver"] / "_omnipath_resolver.abi3.so").unlink()
    with pytest.raises(RuntimeError, match="native resolver"):
        runtime_files.software_provenance()


def test_ambient_ancestor_git_repo_is_not_claimed_for_installed_package(
    installed_packages, monkeypatch
):
    class Result:
        stdout = "/unrelated/repo\n"

    monkeypatch.setattr(runtime_files.subprocess, "run", lambda *args, **kwargs: Result())
    assert runtime_files._git_provenance(installed_packages["omnipath_core"]) is None
