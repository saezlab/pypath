"""
Parser of the MIPS mammalian protein-protein interaction file (MPPI).

MPPI is one PSI-MI XML level 1 file. Each ``<interaction>`` is one experiment
of one publication on its participants.
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from pypath.inputs_v2.parsers.base import iter_psimi_xml, read_bytes


# What a row states about a participant; two rows are mirrors only when all of
# it is alike with the participants exchanged.
_PARTICIPANT_FIELDS = ('db', 'id', 'label', 'name', 'taxon', 'confidence')


def iter_interactions(
    opener, **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Yield one flat row per ``<interaction>`` of the MPPI file.

    Fields: ``index`` (position of the interaction in the file), ``pubmed``
    (the ``id`` of the experiment's PubMed reference, as written),
    ``method_label`` and ``method_mi`` (the detection method's label and
    PSI-MI accession), ``participants`` (their number) and, for the first two
    participants ``a`` and ``b``: ``<x>_db`` and ``<x>_id`` (the primary
    reference), ``<x>_label`` (the MIPS protein number), ``<x>_name``,
    ``<x>_taxon`` and ``<x>_confidence`` (the participant's boolean flag), and
    ``mirror``.

    The file lists most interactions of two different proteins twice, once in
    each order, as two experiments with the same publication and method and
    the same two participants with every field alike. ``mirror`` is True on
    the second of such rows (the row whose participants, taken the other way
    round, were read before) and False on every other row, including a row
    that repeats an earlier one in the same order. The experiment carries no
    number or other field that would tell two such rows apart.
    """
    content = read_bytes(opener)
    if not content:
        raise ValueError('The MPPI download is empty')
    count = 0
    read: set[tuple[Any, ...]] = set()
    for interaction in iter_psimi_xml(content, level='1'):
        experiments = interaction['experiments']
        if len(experiments) != 1:
            raise ValueError(
                f'MPPI interaction {interaction["index"]} has '
                f'{len(experiments)} experiments; the file has one for each'
            )
        (experiment,) = experiments
        publication = next(
            (
                ref
                for ref in experiment['bibliography']
                if ref['primary'] and (ref['db'] or '').lower() == 'pubmed'
            ),
            None,
        )
        method = experiment['detection_method'] or {}
        row = {
            'index': interaction['index'],
            'pubmed': publication['id'] if publication else None,
            'method_label': method.get('label'),
            'method_mi': method.get('mi'),
            'participants': len(interaction['participants']),
        }
        for side, participant in zip('ab', interaction['participants']):
            interactor = participant['interactor'] or {}
            primary = next(
                (ref for ref in interactor.get('xrefs', ()) if ref['primary']),
                {},
            )
            row.update(
                {
                    f'{side}_db': primary.get('db'),
                    f'{side}_id': primary.get('id'),
                    f'{side}_label': interactor.get('label'),
                    f'{side}_name': interactor.get('name'),
                    f'{side}_taxon': (interactor.get('organism') or {}).get(
                        'taxon'
                    ),
                    f'{side}_confidence': (
                        participant['confidence'] or {}
                    ).get('value'),
                }
            )
        row['mirror'] = False
        if row['participants'] == 2:
            first, second = (
                tuple(row[f'{side}_{field}'] for field in _PARTICIPANT_FIELDS)
                for side in 'ab'
            )
            experiment_key = (row['pubmed'], row['method_mi'])
            row['mirror'] = first != second and (second, first, *experiment_key) in read
            if not row['mirror']:
                read.add((first, second, *experiment_key))
        count += 1
        yield row
    if not count:
        raise ValueError('The MPPI file has no interactions')
