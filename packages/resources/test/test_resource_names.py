"""Every resource has a slug, a short and a full name that follow the naming rules."""

import importlib
import pkgutil

import pytest

import omnipath_resources.inputs_v2 as inputs
from omnipath_build.discovery import resolve_source_name
from omnipath_resources.inputs_v2.resource_names import resolve_names, validate_resource_name

# Modules whose slug is not their module name; the module keeps its file name.
RENAMED = {
    "rampdb": "ramp",
    "omnipath_ontology": "omnipathontology",
    "phenol_explorer": "phenolexplorer",
    "psi_mi": "psimi",
    "switches_elm": "switcheselm",
}


@pytest.mark.parametrize("module, slug", sorted(RENAMED.items()))
def test_a_module_with_an_underscore_or_an_old_name_has_the_slug_without_it(module, slug):
    assert resolve_source_name(module) == slug


def test_other_modules_keep_their_name_as_the_slug():
    assert resolve_source_name("signor") == "signor"
    assert resolve_source_name("connectomedb") == "connectomedb2025"


def configs():
    for module in pkgutil.iter_modules(inputs.__path__):
        if module.ispkg or module.name.startswith("_"):
            continue
        try:
            loaded = importlib.import_module(f"omnipath_resources.inputs_v2.{module.name}")
        except ImportError:  # an input that needs the resolver's native extension
            continue
        config = getattr(loaded, "config", None)
        if config is not None:
            yield module.name, config


def test_every_resource_follows_the_naming_rules():
    broken = {}
    for module, config in configs():
        slug = resolve_source_name(module)
        names = resolve_names(name=config.name, slug=slug)
        problems = validate_resource_name(names.slug, names.short, names.full)
        if problems:
            broken[module] = problems
    assert not broken, broken


def test_the_short_and_full_names_are_carried_for_every_resource():
    for module, config in configs():
        names = resolve_names(name=config.name, slug=resolve_source_name(module))
        assert names.short and names.full, module
