"""
The rule file of the metabolite classes: the classes, and the sources that claim molecules.

The file is read and checked before the build touches the release, so a mistake in it
stops the build with the reason instead of publishing a wrong table.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

RULE_FILE = Path(__file__).with_name('metabolite_class.yaml')

# The keys that name a signal, and the other keys each signal takes.
_SIGNALS = {
    'resource': set(),
    'chebi_role_root': set(),
    'predicates': {'resources'},
}
_COMMON = {'priority', 'class'}


class RuleError(ValueError):
    """
    The rule file cannot be used.
    """


@dataclass(frozen = True)
class MetaboliteClass:
    """
    A class of the rule file.

    Args:
        name: The name a row holds and a gate admits.
        precedence: The rank of the class; the lowest is the best.
    """

    name: str
    precedence: int


@dataclass(frozen = True)
class Source:
    """
    One source of claims.

    Args:
        signal: ``resource``, ``chebi_role`` or ``network``.
        priority: The order the source runs in; the lowest runs first.
        metabolite_class: The class the source claims molecules for.
        resource: For a ``resource`` signal, the release resource.
        chebi_role_root: For a ``chebi_role`` signal, the ChEBI role term.
        predicates: For a ``network`` signal, the predicates of the reactions.
        resources: For a ``network`` signal, the resources of the reactions.
    """

    signal: str
    priority: int
    metabolite_class: str
    resource: str | None = None
    chebi_role_root: str | None = None
    predicates: tuple[str, ...] = ()
    resources: tuple[str, ...] = ()


@dataclass(frozen = True)
class Rules:
    """
    The checked content of a rule file.

    Args:
        classes: The classes, best first.
        default: The class of a molecule that no source claims; it holds no rows.
        sources: The sources, in the order they run.
        path: The file the rules were read from.
    """

    classes: tuple[MetaboliteClass, ...]
    default: str
    sources: tuple[Source, ...]
    path: Path


    def precedence(self, name: str) -> int:
        """
        The precedence of the class ``name``.
        """

        return next(c.precedence for c in self.classes if c.name == name)


    @property
    def claimable(self) -> tuple[str, ...]:
        """
        The classes a row can hold: every class but the default, best first.
        """

        return tuple(c.name for c in self.classes if c.name != self.default)


def _source(entry, classes: dict[str, int], default: str, index: int) -> Source:
    """
    One source of the rule file, checked; raises ``RuleError``.
    """

    where = f'source {index + 1}'

    if not isinstance(entry, dict):

        raise RuleError(f'{where} is not a mapping')

    named = [key for key in _SIGNALS if key in entry]

    if len(named) != 1:

        raise RuleError(f'{where} must name exactly one of {", ".join(_SIGNALS)}')

    signal = named[0]
    unknown = set(entry) - _COMMON - {signal} - _SIGNALS[signal]

    if unknown:

        raise RuleError(f'{where} has unknown keys: {", ".join(sorted(unknown))}')

    missing = _COMMON - set(entry)

    if missing:

        raise RuleError(f'{where} needs {", ".join(sorted(missing))}')

    name = entry['class']

    if name not in classes:

        raise RuleError(f'{where} names the unknown class {name!r}')

    if name == default:

        raise RuleError(f'{where} claims the default class {name!r}, which holds no rows')

    priority = entry['priority']

    if not isinstance(priority, int) or isinstance(priority, bool):

        raise RuleError(f'{where} has a priority that is not an integer')

    if signal == 'resource':

        return Source('resource', priority, name, resource = str(entry['resource']))

    if signal == 'chebi_role_root':

        return Source(
            'chebi_role',
            priority,
            name,
            chebi_role_root = str(entry['chebi_role_root']),
        )

    predicates, resources = entry['predicates'], entry.get('resources')

    if not predicates or not isinstance(predicates, list):

        raise RuleError(f'{where} needs a list of predicates')

    if not resources or not isinstance(resources, list):

        raise RuleError(f'{where} needs a list of resources')

    return Source(
        'network',
        priority,
        name,
        predicates = tuple(sorted(map(str, predicates))),
        resources = tuple(sorted(map(str, resources))),
    )


def read_rules(path: Path | str | None = None) -> Rules:
    """
    The rules of a rule file; raises ``RuleError``.

    Args:
        path: The rule file; the packaged one by default.

    Returns:
        The classes, the default class and the sources in the order they run.
    """

    path = Path(path) if path is not None else RULE_FILE
    data = yaml.safe_load(path.read_text(encoding = 'utf-8')) or {}
    classes = {}

    for entry in data.get('classes') or []:

        name, precedence = entry.get('name'), entry.get('precedence')

        if not name or not isinstance(precedence, int) or name in classes:

            raise RuleError(f'class {entry!r} needs a unique name and an integer precedence')

        classes[name] = precedence

    if len(set(classes.values())) != len(classes):

        raise RuleError('two classes share a precedence')

    default = data.get('default')

    if default not in classes:

        raise RuleError(f'the default class {default!r} is not a class')

    sources = tuple(
        _source(entry, classes, default, index)
        for index, entry in enumerate(data.get('sources') or [])
    )

    if len({s.priority for s in sources}) != len(sources):

        raise RuleError('two sources share a priority')

    # A source that runs first must not give a worse class than one that runs later:
    # the first claim and the group's best class then name the same class.
    ordered = sorted(sources, key = lambda s: s.priority)

    for earlier, later in zip(ordered, ordered[1:]):

        if classes[earlier.metabolite_class] > classes[later.metabolite_class]:

            raise RuleError(
                f'priority {earlier.priority} gives {earlier.metabolite_class!r}, which '
                f'ranks below {later.metabolite_class!r} of priority {later.priority}'
            )

    return Rules(
        classes = tuple(
            MetaboliteClass(name, precedence)
            for name, precedence in sorted(classes.items(), key = lambda item: item[1])
        ),
        default = default,
        sources = ordered,
        path = path,
    )
