"""
Parser of the NetPath PSI-MI export and of the NetPath pathway list.

NetPath distributed one PSI-MI XML 2.5 file per signaling pathway, zipped into
one archive. The pathway files are read with the shared PSI-MI XML reader
(:func:`iter_psimi_xml`); this module collects the interaction, its
experiments and its distinct, identified interactors into one row per
interaction and experiment. The names of the pathways are on a copy of the
site's "browse" page, read with the shared link reader.
"""

from __future__ import annotations

from collections.abc import Generator
import re
from typing import Any

from pypath.inputs_v2.parsers.base import (
    content_text,
    iter_html_links,
    iter_psimi_xml,
    read_text,
)

# Database names of interactor cross-references -> parsed row field.
_XREF_FIELDS = {
    'HPRD': 'hprd',
    'Entrez gene': 'entrez',
    'Refseq': 'refseq',
    'GenBank': 'genbank',
    'OMIM': 'omim',
}
# The cross-reference fields of a participant record.
_ACCESSION_FIELDS = (*_XREF_FIELDS.values(), 'uniprot')
# The form of a UniProtKB accession (the pattern of the UniProt manual).
# NetPath files UniProt accessions under the "GenBank" label. The two forms
# cannot be confused: a GenBank protein accession is three letters and five
# to seven digits (``AAC24716``), a UniProtKB accession starts with one
# letter and a digit. A value of the "GenBank" field in the UniProt form is
# read as a UniProt accession (``uniprot``); the other values stay ``genbank``.
UNIPROT_RE = re.compile(
    r'[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}'
)

_PATHWAY_MEMBER_RE = re.compile(r'NetPath_(\d+)_PSIMI\.xml$')
_PATHWAY_LINK_RE = re.compile(r'^pathways\?path_id=(NetPath_\d+)$')


def _archive_members(opener) -> list[tuple[str, bytes]]:
    """Pathway files of the archive as (name, content), in pathway order."""
    result = getattr(opener, 'result', None) if opener else None
    if not isinstance(result, dict):
        raise ValueError('NetPath: expected a zip archive of PSI-MI files')
    members = []
    for name, content in result.items():
        match = _PATHWAY_MEMBER_RE.search(name)
        if match:
            members.append((int(match[1]), name, content))
    if not members:
        raise ValueError('NetPath: the archive has no NetPath_*_PSIMI.xml file')
    return [
        (name, content_text(content).encode('utf-8'))
        for _, name, content in sorted(members)
    ]


def _unique(values: list[str]) -> list[str]:
    return list(dict.fromkeys(value for value in values if value))


def _participant(participant: dict[str, Any]) -> dict[str, Any]:
    """Flat record of one PSI-MI participant and its interactor."""
    interactor = participant['interactor'] or {}
    record: dict[str, Any] = {field: [] for field in _ACCESSION_FIELDS}
    for xref in interactor.get('xrefs') or []:
        if xref['db'] not in _XREF_FIELDS:
            raise ValueError(
                f'NetPath: unknown interactor database {xref["db"]!r}'
            )
        field = _XREF_FIELDS[xref['db']]
        if field == 'genbank' and UNIPROT_RE.fullmatch(xref['id'].strip()):
            field = 'uniprot'
        record[field].append(xref['id'])
    interactor_type = interactor.get('type') or {}
    organism = interactor.get('organism') or {}
    roles = participant['experimental_roles']
    record.update(
        interactor_id=interactor.get('id'),
        symbol=interactor.get('label'),
        name=interactor.get('name'),
        synonyms=_unique(
            [alias['value'] for alias in interactor.get('aliases') or []]
        ),
        type=interactor_type.get('label'),
        type_mi=interactor_type.get('mi'),
        taxon=organism.get('taxon'),
        organism=organism.get('label'),
        experimental_role=[role['label'] for role in roles],
        experimental_role_mi=[role['mi'] for role in roles],
    )
    return record


def _identified(participant: dict[str, Any]) -> bool:
    return any(participant[field] for field in _ACCESSION_FIELDS)


def _distinct_interactors(
    participants: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """One record per interactor, in source order, with all its roles.

    NetPath lists the subunit of a homo-oligomer once per copy.
    """
    distinct: dict[Any, dict[str, Any]] = {}
    for participant in participants:
        record = distinct.setdefault(
            participant['interactor_id'],
            {**participant, 'experimental_role': [], 'experimental_role_mi': []},
        )
        for field in ('experimental_role', 'experimental_role_mi'):
            record[field].extend(
                role for role in participant[field] if role not in record[field]
            )
    return list(distinct.values())


def iter_netpath_observations(
    opener, **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
    """Yield one dictionary per NetPath interaction and experiment.

    Output fields:

    - ``pathway_id``: ``NetPath_<n>``, from the entry's source reference.
    - ``entry_index``: position of the interaction in the pathway file.
    - ``interaction_name``: the interaction's full name.
    - ``interaction_type``, ``interaction_type_mi``: label and PSI-MI term.
    - ``experiment_id``: the experiment's ID within the entry.
    - ``pubmed``: PubMed IDs of the experiment's ``<bibref>``.
    - ``detection_method``, ``detection_method_mi``: label and PSI-MI term.
    - ``host_taxon``, ``host_organism``: NCBI Taxonomy ID and label of the
      experiment's host organism.
    - ``participant_count``: number of participants the source lists.
    - ``participants``: the distinct interactors that have a database
      cross-reference, in source order. Each has ``interactor_id``, ``symbol``
      (short label), ``name``, ``synonyms``, ``type``, ``type_mi``, ``taxon``,
      ``organism``, ``hprd``, ``entrez``, ``refseq``, ``genbank``,
      ``uniprot``, ``omim`` (lists, as written), ``experimental_role`` and
      ``experimental_role_mi`` (lists). ``uniprot`` holds the values of the
      source's "GenBank" field that are in the form of a UniProtKB accession
      (see ``UNIPROT_RE``); ``genbank`` holds the other values.
    - ``unidentified_participants``: interactors without a cross-reference.
      The export has two, both a repeated label typed "small molecule" in a
      ligand-receptor interaction; they are not mapped.
    """
    for name, content in _archive_members(opener):
        for interaction in iter_psimi_xml(content, level='2.5'):
            pathways = [
                xref['id']
                for xref in interaction['source']['xrefs']
                if xref['db'] == 'NetPath'
            ]
            types = interaction['interaction_types']
            if len(pathways) != 1 or len(types) != 1:
                raise ValueError(
                    f'NetPath: {name} entry {interaction["entry_index"]} '
                    'needs one pathway reference and one interaction type'
                )
            participants = _distinct_interactors(
                list(map(_participant, interaction['participants']))
            )
            base = {
                'pathway_id': f'NetPath_{pathways[0]}',
                'entry_index': interaction['entry_index'],
                'interaction_name': interaction['name'],
                'interaction_type': types[0]['label'],
                'interaction_type_mi': types[0]['mi'],
                'participant_count': len(interaction['participants']),
                'participants': [p for p in participants if _identified(p)],
                'unidentified_participants': [
                    p for p in participants if not _identified(p)
                ],
            }
            if not interaction['experiments']:
                raise ValueError(
                    f'NetPath: {name} entry {interaction["entry_index"]} '
                    'has no experiment'
                )
            for experiment in interaction['experiments']:
                method = experiment['detection_method'] or {}
                hosts = experiment['host_organisms']
                if len(hosts) > 1:
                    raise ValueError(
                        f'NetPath: {name} experiment with several hosts'
                    )
                host = hosts[0] if hosts else {}
                yield {
                    **base,
                    'experiment_id': experiment['id'],
                    'pubmed': _unique(
                        [
                            xref['id']
                            for xref in experiment['bibliography']
                            if xref['db'].lower() == 'pubmed'
                        ]
                    ),
                    'detection_method': method.get('label'),
                    'detection_method_mi': method.get('mi'),
                    'host_taxon': host.get('taxon'),
                    'host_organism': host.get('label'),
                }


def iter_netpath_interactions(opener, **kwargs: Any):
    """Observations with one or two distinct identified interactors."""
    for row in iter_netpath_observations(opener, **kwargs):
        if 1 <= len(row['participants']) <= 2:
            yield row


def iter_netpath_complexes(opener, **kwargs: Any):
    """Observations with three or more distinct identified interactors."""
    for row in iter_netpath_observations(opener, **kwargs):
        if len(row['participants']) >= 3:
            yield row


def iter_netpath_pathways(opener, **_kwargs: Any):
    """Yield ``pathway_id`` and ``name`` for each pathway of the browse page."""
    found = False
    for link in iter_html_links(read_text(opener)):
        match = _PATHWAY_LINK_RE.match(link['href'])
        if match and link['text']:
            found = True
            yield {'pathway_id': match[1], 'name': link['text']}
    if not found:
        raise ValueError('NetPath: the browse page lists no pathway')
