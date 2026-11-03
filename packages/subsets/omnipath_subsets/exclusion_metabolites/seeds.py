"""Read and validate the seed files of the exclusion list.

Three files in ``seed/`` hold everything a person edits:

* ``metaproviz_seed.yaml``: the MetaProViz lists, vendored with their origin.
* ``project_cofactors.yaml``: the reviewed core of cofactors.
* ``overrides.yaml``: members a curator adds, and members a curator removes.

A malformed file stops the build with the reason. A seed that is merely absent
from the release is not an error here: resolution counts it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

SEED_DIR = Path(__file__).with_name("seed")

METAPROVIZ_FILE = "metaproviz_seed.yaml"
COFACTOR_FILE = "project_cofactors.yaml"
OVERRIDE_FILE = "overrides.yaml"

CLASSES = ("ions", "atoms", "small_molecules", "xenobiotics", "cofactors", "cofactors_extended")

# The chemical namespaces the release carries as identifier types.
NAMESPACES = (
    "chebi",
    "kegg",
    "hmdb",
    "pubchem",
    "lipidmaps",
    "swisslipids",
    "chembl",
    "drugbank",
)


class SeedError(ValueError):
    """A seed file is malformed."""


@dataclass(frozen=True)
class Seed:
    """One identifier a curator listed for a class."""

    exclusion_class: str
    namespace: str
    identifier: str
    assignment_source: str
    source_file: str
    name: str | None = None
    note: str | None = None

    @property
    def line(self) -> str:
        """The seed as a reviewer would read it in the file."""
        return f"{self.source_file}: {self.exclusion_class} {self.namespace} {self.identifier}"


@dataclass(frozen=True)
class SeedSet:
    seeds: list[Seed]
    vetoes: list[Seed]
    paths: tuple[Path, ...]
    source: dict = field(default_factory=dict)


def _load(path: Path) -> dict:
    try:
        document = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as error:
        raise SeedError(f"{path.name}: not valid YAML: {error}") from error
    if document is None:
        return {}
    if not isinstance(document, dict):
        raise SeedError(f"{path.name}: expected a mapping at the top level")
    return document


def _class(name: str, value) -> str:
    if value not in CLASSES:
        raise SeedError(f"{name}: unknown class {value!r}; valid classes are {', '.join(CLASSES)}")
    return value


def _namespace(name: str, value) -> str:
    if value not in NAMESPACES:
        raise SeedError(
            f"{name}: unknown namespace {value!r}; valid namespaces are {', '.join(NAMESPACES)}"
        )
    return value


def _identifier(name: str, value) -> str:
    text = "" if value is None else str(value).strip()
    if not text:
        raise SeedError(f"{name}: empty identifier")
    return text


def _unique(seeds: list[Seed], what: str) -> None:
    seen = set()
    for seed in seeds:
        key = (seed.exclusion_class, seed.namespace, seed.identifier)
        if key in seen:
            raise SeedError(f"{seed.source_file}: duplicate {what} {' '.join(key)}")
        seen.add(key)


def _metaproviz(path: Path) -> tuple[list[Seed], dict]:
    document = _load(path)
    seeds = []
    for group in document.get("seeds") or []:
        exclusion_class = _class(path.name, group.get("class"))
        namespace = _namespace(path.name, group.get("namespace"))
        for identifier in group.get("ids") or []:
            seeds.append(
                Seed(
                    exclusion_class,
                    namespace,
                    _identifier(path.name, identifier),
                    "metaproviz_seed",
                    path.name,
                )
            )
    _unique(seeds, "seed")
    return seeds, dict(document.get("source") or {})


def _cofactors(path: Path) -> list[Seed]:
    document = _load(path)
    members = document.get("members") or []
    if not members:
        return []
    exclusion_class = _class(path.name, document.get("class"))
    namespace = _namespace(path.name, document.get("namespace"))
    seeds = [
        Seed(
            exclusion_class,
            namespace,
            _identifier(path.name, member.get("id")),
            "project_seed",
            path.name,
            name=member.get("name"),
        )
        for member in members
    ]
    _unique(seeds, "seed")
    return seeds


def _overrides(path: Path) -> tuple[list[Seed], list[Seed]]:
    document = _load(path)

    def read(key: str, source: str) -> list[Seed]:
        listed = [
            Seed(
                _class(path.name, entry.get("class")),
                _namespace(path.name, entry.get("namespace")),
                _identifier(path.name, entry.get("id")),
                source,
                path.name,
                note=entry.get("note"),
            )
            for entry in document.get(key) or []
        ]
        _unique(listed, key)
        return listed

    return read("add", "override_add"), read("veto", "override_veto")


def read_seeds(directory: Path | str = SEED_DIR) -> SeedSet:
    """Read the three files of ``directory``, in the order a reviewer reads them."""
    directory = Path(directory)
    paths = tuple(directory / name for name in (METAPROVIZ_FILE, COFACTOR_FILE, OVERRIDE_FILE))
    for path in paths:
        if not path.is_file():
            raise SeedError(f"{path.name}: missing from {directory}")
    curated, source = _metaproviz(paths[0])
    added, vetoes = _overrides(paths[2])
    return SeedSet(
        seeds=curated + _cofactors(paths[1]) + added,
        vetoes=vetoes,
        paths=paths,
        source=source,
    )
