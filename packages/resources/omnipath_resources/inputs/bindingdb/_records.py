from __future__ import annotations

from typing import NamedTuple

__all__ = [
    "BindingdbInteraction",
    "BindingdbLigand",
    "BindingdbTarget",
    "AllostericRegulation",
    "ReactionConstant",
]


class BindingdbLigand(NamedTuple):
    """BindingDB ligand record."""
    name: str | None = None
    smiles: str | None = None
    inchi: str | None = None
    inchi_key: str | None = None
    pubchem: str | None = None


class BindingdbTarget(NamedTuple):
    """BindingDB target record."""
    name: str | None = None
    organism: str | None = None
    ncbi_tax_id: int | None = None
    uniprot: str | None = None
    regions_mutations: list[list[str]] | None = None


class BindingdbInteraction(NamedTuple):
    """BindingDB interaction record."""
    ligand: BindingdbLigand
    target: BindingdbTarget


class ReactionConstant(NamedTuple):
    """Reaction constant for allosteric regulation."""
    type: str | None = None
    value: str | None = None
    conditions: str | None = None
    pubmeds: list[str] | None = None


class AllostericRegulation(NamedTuple):
    """Allosteric regulation record."""
    action: str | None = None
    compound: str | None = None
    organism: str | None = None
    protein: str | None = None
    id_type: str | None = None
    wrong_ec: str | None = None
    pubmeds: list[str] | None = None
    reaction_constants: list[ReactionConstant] | None = None