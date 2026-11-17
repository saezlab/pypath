"""The interaction-class rules name the predicates the build stores: Biolink's."""

from omnipath_core import biolink
from omnipath_postgres.relational.classify import interaction_class


def rules():
    return interaction_class._load_rules()


def test_every_rule_names_a_biolink_predicate():
    # A rule keyed on a name the build never stores matches nothing and hides a class.
    for name in rules()["predicates"]:
        assert biolink.predicate(name) == name


def test_every_rule_names_a_class_of_the_vocabulary():
    known = {entry["name"] for entry in rules()["classes"]}
    assert set(rules()["predicates"].values()) <= known
    assert rules()["default"] in known


def test_the_miRBase_maturation_edge_is_a_maturation():
    # Each mature miRNA derives_from its precursor.
    assert rules()["predicates"]["derives_from"] == "maturation"


def test_causal_regulation_stays_signaling():
    mapped = rules()["predicates"]
    for name in ("affects", "regulates", "enabled_by", "catalyzes"):
        assert mapped[name] == "signaling"
