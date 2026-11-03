"""The complete scientific oracle must work in a checkout without Git history."""

from oracles import (
    CV_ROOT,
    SOURCE_HASHES,
    read_checked,
    read_deployed_cosmos,
    read_frozen,
    read_legacy,
)
import pytest


def test_all_frozen_sources_are_tracked_and_hash_checked():
    for relative in SOURCE_HASHES["legacy"]:
        read_legacy(relative)
    for relative in SOURCE_HASHES["cv"]:
        read_frozen(CV_ROOT / relative)
    assert read_deployed_cosmos()


def test_changed_oracle_fails_instead_of_skipping(tmp_path):
    source = tmp_path / "oracle.sql"
    source.write_text("changed SQL")
    with pytest.raises(AssertionError, match="Immutable oracle changed"):
        read_checked(source, "0" * 64)
