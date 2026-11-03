"""Parser parity: the legacy pypath input and its inputs_v2 port, side by side.

For each resource, run the legacy function and the inputs_v2 datasets, reduce
both to comparable keys in the source's own identifiers, and report what only
one side has. Both sides read the same cached file and neither translates
identifiers, so a difference is a difference in parsing (each result lists
the files each side read, with their checksums).

Keys by kind:
- relation: (subject, object) with the sign (+, -, 0 or none) compared on the
  pairs both sides have;
- complex: the set of member identifiers;
- keyed_complex: the complex's own identifier, then its members where they
  are plain identifiers on both sides (legacy expands nested complexes and
  families into every combination of their members, inputs_v2 keeps them as
  references);
- annotation: (entity, value).

Legacy identifier translation (pypath.utils.mapping) is replaced by identity,
so legacy keeps the source's identifiers too.

Nothing is downloaded: both sides read the files in PYPATH_DOWNLOAD_DATADIR
(the build's download directory), and a legacy input that asks for a file the
inputs_v2 module does not have fails with NotCached. Each resource runs in its
own process. The legacy inputs need pycurl, which the project environment does
not install. On nicesrv, from the build checkout:

    PYPATH_DOWNLOAD_DATADIR=/root/data/pypath-cache \
        uv run --frozen --with pycurl python scripts/legacy_parity.py \
        trrust corum --out /root/data/scratch/parity/out
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import traceback
from typing import Any, Literal

SAMPLE = 25

Kind = Literal['relation', 'complex', 'keyed_complex', 'annotation']


# Setup ---------------------------------------------------------------------
# Nothing is downloaded. inputs_v2 opens the files already in the build's
# download directory; a legacy download is served the inputs_v2 file of the
# same request (URL and POST data), and any other legacy download fails.

_files: list[dict[str, Any]] = []
_served: dict[tuple, str] = {}


class NotCached(RuntimeError):
    """A download the cache cannot serve."""


def _sha256(path: str | os.PathLike) -> str | None:
    try:
        digest = hashlib.sha256()
        with open(path, 'rb') as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b''):
                digest.update(chunk)
        return digest.hexdigest()
    except OSError:
        return None


def _record_file(side: str, url: str, path: str | None) -> None:
    entry = {'side': side, 'url': url, 'path': path}
    if path and os.path.isfile(path):
        entry['size'] = os.path.getsize(path)
        entry['sha256'] = _sha256(path)
        entry['mtime'] = time.strftime(
            '%Y-%m-%d', time.localtime(os.path.getmtime(path))
        )
    if entry not in _files:
        _files.append(entry)


def _text(value: Any) -> str:
    return value.decode() if isinstance(value, bytes) else str(value)


def _request(url: str, post: Any = None) -> tuple:
    """A download request as both generations can spell it."""
    url = re.sub(r'^https?://', '', url.strip()).rstrip('/')
    if isinstance(post, dict):
        post = post.items()
    params = frozenset((_text(k), _text(v)) for k, v in post or ())
    return url, params


def serve(module: str, aliases: dict[str, str] | None = None) -> None:
    """Let legacy downloads read the files of an inputs_v2 module."""
    import importlib

    from pypath.share.downloads import get_data_dir

    data_dir = get_data_dir()
    resource = importlib.import_module(f'pypath.inputs_v2.{module}').resource
    for dataset in resource.datasets().values():
        download = dataset.download
        if download is None or callable(download.url) or callable(download.filename):
            continue
        kwargs = download.download_kwargs or {}
        post = kwargs.get('query') if kwargs.get('post') else None
        path = str(data_dir / download.subfolder / download.filename)
        _served[_request(download.url, post)] = path
    for legacy_url, v2_url in (aliases or {}).items():
        for (url, params), path in list(_served.items()):
            if url == _request(v2_url)[0]:
                _served[_request(legacy_url)[0], params] = path


def setup(scratch: Path) -> None:
    """Make both input generations importable, offline, logging the files read."""
    from omnipath_build.cachedir_compat import patch_cachedir_opener

    patch_cachedir_opener()

    # The legacy package writes its own caches; keep them out of the build's
    # download directory.
    import pypath.share.cache as legacy_cache

    legacy_dir = scratch / 'legacy-cache'
    legacy_dir.mkdir(parents=True, exist_ok=True)
    legacy_cache.get_cachedir = lambda *a, **k: legacy_dir

    # Module-level lookups of the legacy package that call services whose
    # formats changed; nothing compared here depends on them.
    import pypath.inputs.unichem as unichem
    import pypath.inputs.mirbase as mirbase

    unichem.unichem_sources = lambda *a, **k: {}
    mirbase.mirbase_organisms = lambda *a, **k: {}

    # Legacy inputs translate identifiers with tables downloaded from UniProt;
    # here every name translates to itself, so that legacy keeps the source's
    # identifiers as inputs_v2 does.
    import pypath.utils.mapping as mapping
    import pypath.inputs.uniprot_db as uniprot_db

    mapping.map_name = lambda name, *a, **k: {name} if name else set()
    mapping.map_name0 = lambda name, *a, **k: name
    mapping.map_names = lambda names, *a, **k: {n for n in names if n}
    mapping.label = lambda name, *a, **k: name
    mapping.id_from_label = lambda label, *a, **k: {label} if label else set()
    mapping.id_from_label0 = lambda label, *a, **k: label
    uniprot_db.is_uniprot = lambda name, *a, **k: bool(_ISOFORM.match(f'{name}-1'))
    uniprot_db.is_swissprot = uniprot_db.is_uniprot

    # Organism names, without the Ensembl and NCBI tables the legacy lookups
    # download.
    import pypath.utils.taxonomy as taxonomy

    organisms = {
        9606: ('human', 'Homo sapiens', 'hs'),
        10090: ('mouse', 'Mus musculus', 'mm'),
        10116: ('rat', 'Rattus norvegicus', 'rn'),
    }
    by_name = {
        name.lower(): taxid for taxid, names in organisms.items() for name in names
    }
    by_name |= {'murine': 10090}

    def ensure_ncbi_tax_id(taxon):
        if taxon is None:
            return None
        if isinstance(taxon, int) or str(taxon).isdigit():
            return int(taxon)
        return by_name.get(str(taxon).strip().lower())

    def ensure_common_name(taxon, lower=False):
        taxid = ensure_ncbi_tax_id(taxon)
        name = organisms[taxid][0] if taxid in organisms else None
        return name if lower or name is None else name.capitalize()

    taxonomy.ensure_ncbi_tax_id = ensure_ncbi_tax_id
    taxonomy.ensure_common_name = ensure_common_name
    taxonomy.ensure_latin_name = lambda t: organisms.get(ensure_ncbi_tax_id(t), (None, None))[1]

    from cachedir._open import Opener
    import pypath.inputs_v2.base as v2_base
    from pypath.share.downloads import get_data_dir
    import pypath.share.curl as curl

    def cached_download_and_open(
        url, filename, subfolder, large=True, encoding='utf-8',
        default_mode='r', ext=None, needed=None, **_download_kwargs,
    ):
        path = get_data_dir() / subfolder / filename
        if not path.is_file() or not path.stat().st_size:
            raise NotCached(f'inputs_v2 file not in the cache: {path} ({url})')
        _record_file('v2', url, str(path))
        return Opener(
            path=str(path), ext=ext, needed=needed, large=large,
            default_mode=default_mode, encoding=encoding,
        )

    v2_base.download_and_open = cached_download_and_open

    # Legacy inputs ported to the download manager call the shared
    # download_and_open: served the inputs_v2 file of the same request.
    import pypath.share.downloads as downloads

    def legacy_download_and_open(url, filename, subfolder, large=True, encoding='utf-8',
                                 default_mode='r', ext=None, needed=None, **kwargs):
        post = kwargs.get('query') if kwargs.get('post') else None
        path = _served.get(_request(url, post))
        if path is None or not os.path.isfile(path):
            raise NotCached(f'legacy download with no inputs_v2 file: {url}')
        _record_file('legacy', url, path)
        return Opener(
            path=path, ext=ext, needed=needed, large=large,
            default_mode=default_mode, encoding=encoding,
        )

    original_download_and_open = downloads.download_and_open
    downloads.download_and_open = legacy_download_and_open
    # Modules that imported the function by name before this point.
    for module in list(sys.modules.values()):
        if (
            module is not None
            and module is not v2_base
            and getattr(module, 'download_and_open', None) is original_download_and_open
        ):
            module.download_and_open = legacy_download_and_open

    # ... and some call the download manager itself for the file's path.
    class CachedDownloads:
        def download(self, url, *args, **kwargs):
            post = kwargs.get('query') if kwargs.get('post') else None
            path = _served.get(_request(url, post))
            if path is None or not os.path.isfile(path):
                raise NotCached(f'legacy download with no inputs_v2 file: {url}')
            _record_file('legacy', url, path)
            return path

        def __getattr__(self, name):
            raise NotCached(f'download manager: {name}')

    downloads.get_download_manager = lambda *a, **k: CachedDownloads()

    curl_init = curl.Curl.__init__

    def cached_curl_init(self, url, *args, **kwargs):
        if kwargs.get('setup', True) is False or os.path.exists(url):
            return curl_init(self, url, *args, **kwargs)
        post = kwargs.get('binary_data') or kwargs.get('post')
        path = _served.get(_request(url, post))
        if path is None or not os.path.isfile(path):
            raise NotCached(
                f'legacy download with no inputs_v2 file: {url} {post or ""}'.strip()
            )
        _record_file('legacy', url, path)
        kwargs.pop('binary_data', None)
        kwargs.pop('post', None)
        curl_init(self, path, *args, **kwargs)

    curl.Curl.__init__ = cached_curl_init

    # Some legacy inputs download with requests or urllib instead of Curl.
    import requests.adapters
    import urllib.request

    def no_download(*args, **_kwargs):
        # HTTPAdapter.send(self, request, ...) or urlopen(url_or_request, ...)
        request = args[1] if len(args) > 1 else args[0] if args else None
        url = getattr(request, 'url', None) or getattr(request, 'full_url', None) or request
        raise NotCached(f'legacy download outside Curl: {url}')

    requests.adapters.HTTPAdapter.send = no_download
    urllib.request.urlopen = no_download


# Keys ----------------------------------------------------------------------

_ISOFORM = re.compile(r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})-\d+$')


def tokens(value: Any) -> frozenset[str]:
    """Spellings of one identifier that count as the same: without a CURIE
    prefix, upper case, and a UniProt isoform also as its canonical entry."""
    if value is None:
        return frozenset()
    value = str(value).strip()
    if not value:
        return frozenset()
    if value.upper().startswith('INCHIKEY='):
        value = value[9:]
    out = {value.upper()}
    # HMDB accessions of five digits are the old spelling of seven.
    if re.fullmatch(r'HMDB\d{5}', value.upper()):
        out.add(f'HMDB00{value[4:]}')
    if ':' in value:
        out.add(value.split(':', 1)[1].strip().upper())
    if re.fullmatch(r'\d+(?:_\d+)+', value):
        # Legacy SPIKE writes a complex as the Entrez IDs of its members.
        value = value.replace('_', ',')
        out.add(value)
    if ',' in value:
        # A list of complex members, in any order.
        out.add(','.join(sorted(v.strip().upper() for v in value.split(',') if v.strip())))
    for v in list(out):
        if match := _ISOFORM.match(v):
            out.add(match.group(1))
    return frozenset(out)


def id_tokens(values: Iterable[Any]) -> frozenset[str]:
    """Spellings of several identifiers; a list among them is its items."""
    out = set()
    for value in values or ():
        if isinstance(value, (list, tuple, set, frozenset)):
            out |= id_tokens(value)
        else:
            out |= tokens(value)
    return frozenset(out)


def sign_of(value: Any) -> str | None:
    """+, -, 0 or None from the many ways inputs say it."""
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return '+' if value > 0 else '-' if value < 0 else '0'
    text = str(value).strip().lower()
    if not text or text in {'nan', 'none', 'unknown', '?'}:
        return None
    if text in {'1', '+', '+1'}:
        return '+'
    if text in {'-1', '-'}:
        return '-'
    if text in {'0'}:
        return '0'
    negative = ('down', 'inhib', 'decreas', 'repress', 'negativ', 'antagon', 'destab')
    positive = ('up', 'activ', 'stimul', 'increas', 'induc', 'positiv', 'agon')
    has_negative = any(word in text for word in negative)
    has_positive = any(word in text for word in positive) and not has_negative
    if has_negative and not any(text.startswith(w) for w in positive):
        return '-'
    if has_positive:
        return '+'
    if 'unknown' in text or 'unspecified' in text or 'neutral' in text:
        return '0'
    return None


@dataclass
class Rel:
    subject: frozenset[str]
    object: frozenset[str]
    sign: str | None = None
    label: str = ''


@dataclass
class Spec:
    name: str
    kind: Kind
    legacy: Callable[[], Iterable[Any]]
    v2: Callable[[], Iterable[Any]]
    # The inputs_v2 module whose files the legacy downloads are served.
    serve: str = ''
    aliases: dict[str, str] = field(default_factory=dict)
    note: str = ''
    # Only the first records of each inputs_v2 dataset (for the largest
    # resources): every v2 relation should be in legacy, but legacy pairs
    # missing in v2 then say nothing.
    v2_subset: int | None = None


# inputs_v2 records -----------------------------------------------------------

def name(term: Any) -> str:
    """A Biolink slot, enum member or CV term as its plain name."""
    for attr in ('name', 'text'):
        value = getattr(term, attr, None)
        if isinstance(value, str):
            return value.replace(' ', '_')
    return str(getattr(term, 'value', term))


def _ref_ids(ref) -> list[str]:
    if hasattr(ref, 'identifiers'):
        return [i.value for i in ref.identifiers or ()]
    return [getattr(ref, 'identifier', None)]


def _member_tokens(ref) -> frozenset[str]:
    """A complex as its members: for each namespace all members have, their
    identifiers of it joined in order."""
    members = [m.member for m in getattr(ref, 'membership', None) or () if not m.is_parent]
    if len(members) < 2:
        return frozenset()
    by_type = [
        {name(i.type): str(i.value).strip().upper() for i in member.identifiers or ()}
        for member in members
    ]
    shared = set.intersection(*(set(ids) for ids in by_type))
    return frozenset(','.join(sorted(ids[t] for ids in by_type)) for t in shared)


def _ref_taxa(ref) -> set[str]:
    return {
        name(a.value).removeprefix('NCBITaxon:')
        for a in getattr(ref, 'annotations', None) or ()
        if name(a.term).endswith('in_taxon')
    }


# Aspects whose direction is a regulatory sign. A decreased phosphorylation
# (a phosphatase) says nothing about activation or inhibition, and the legacy
# inputs give no sign for it.
_SIGNED_ASPECTS = {None, 'activity', 'abundance', 'activity_or_abundance', 'expression', 'stability'}


def _v2_sign(relation) -> str | None:
    annotations = {
        name(a.term): name(a.value) for a in relation.annotations or ()
    }
    if annotations.get('object_aspect_qualifier') not in _SIGNED_ASPECTS:
        return None
    if 'object_direction_qualifier' in annotations:
        return sign_of(annotations['object_direction_qualifier'])
    return sign_of(name(relation.predicate).rsplit(':', 1)[-1].replace('_', ' '))


def v2_relations(records, taxon: str | None = None) -> Iterable[Rel]:
    for relation in records:
        if relation is None:  # a row the mapper drops
            continue
        if taxon and any(
            (taxa := _ref_taxa(ref)) and taxon not in taxa
            for ref in (relation.subject, relation.object)
        ):
            continue
        yield Rel(
            id_tokens(_ref_ids(relation.subject)) | _member_tokens(relation.subject),
            id_tokens(_ref_ids(relation.object)) | _member_tokens(relation.object),
            _v2_sign(relation),
            label=f'{_ref_ids(relation.subject)} -[{name(relation.predicate)}]-> {_ref_ids(relation.object)}',
        )


def _primary_id(entity, prefer: str) -> str | None:
    identifiers = list(getattr(entity, 'identifiers', None) or ())
    for identifier in identifiers:
        if prefer in name(identifier.type).lower():
            return identifier.value
    return identifiers[0].value if identifiers else getattr(entity, 'identifier', None)


def v2_complexes(records, prefer: str = 'uniprot', taxon: str | None = None) -> Iterable[frozenset[str]]:
    for entity in records:
        if entity is None:
            continue
        members = [m.member for m in entity.membership or () if not m.is_parent]
        if taxon and any((t := _ref_taxa(m)) and taxon not in t for m in members):
            continue
        ids = {_primary_id(m, prefer) for m in members}
        ids.discard(None)
        if ids:
            yield frozenset(min(tokens(i), key=len) for i in ids)


_V2_SUBSET: int | None = None


def datasets(module: str, *names: str, **kwargs) -> Callable[[], Iterable[Any]]:
    def run():
        import importlib

        resource = importlib.import_module(f'pypath.inputs_v2.{module}').resource
        if _V2_SUBSET:
            kwargs.setdefault('max_records', _V2_SUBSET)
        for name in names:
            yield from getattr(resource, name)(**kwargs)

    return run


# Comparison ----------------------------------------------------------------

def compare_relations(legacy: Iterable[Rel], v2: Iterable[Rel]) -> dict[str, Any]:
    # Each legacy record with every spelling of its pair (prefix, isoform,
    # abbreviation); records with the same spellings are one pair.
    legacy_records: dict[frozenset, set] = defaultdict(set)
    for rel in legacy:
        pairs = frozenset((a, b) for a in rel.subject for b in rel.object)
        if pairs:
            legacy_records[pairs].add(rel.sign)
    legacy_signs = {pair for pairs in legacy_records for pair in pairs}

    v2_total = 0
    v2_unmatched: list[str] = []
    v2_reversed: list[str] = []
    v2_ids: set[str] = set()
    matched: set[tuple[str, str]] = set()
    sign_pairs: dict[tuple[str, str], set] = defaultdict(set)
    v2_signs = Counter()
    v2_all: set[tuple[str, str]] = set()
    for rel in v2:
        v2_total += 1
        v2_ids |= rel.subject | rel.object
        v2_all.update((a, b) for a in rel.subject for b in rel.object)
        v2_signs[rel.sign] += 1
        hits = {(a, b) for a in rel.subject for b in rel.object if (a, b) in legacy_signs}
        if hits:
            matched |= hits
            for hit in hits:
                sign_pairs[hit].add(rel.sign)
        elif any((b, a) in legacy_signs for a in rel.subject for b in rel.object):
            v2_reversed.append(rel.label)
        else:
            v2_unmatched.append(rel.label)

    # A legacy record is found when any spelling of its pair matched.
    found = [rec for rec in legacy_records if rec & matched]
    unmatched = [rec for rec in legacy_records if not rec & matched]
    # An undirected pair that v2 has the other way round.
    reversed_ = [rec for rec in unmatched if any((b, a) in v2_all for a, b in rec)]
    missing = [rec for rec in unmatched if not any((b, a) in v2_all for a, b in rec)]

    sign_agree = sign_disagree = 0
    disagreements = []
    for rec in found:
        pair = next(iter(rec & matched))
        # A sign of 0 states no direction, as no sign does.
        legacy_sign = legacy_records[rec] - {None, '0'}
        v2_sign = {s for p in rec & matched for s in sign_pairs[p]} - {None, '0'}
        if not legacy_sign or not v2_sign:
            continue
        if legacy_sign == v2_sign:
            sign_agree += 1
        else:
            sign_disagree += 1
            if len(disagreements) < SAMPLE:
                disagreements.append(
                    {'pair': list(pair), 'legacy': sorted(legacy_sign), 'v2': sorted(v2_sign)}
                )

    def known(rec):
        return max((a in v2_ids) + (b in v2_ids) for a, b in rec)

    missing_known = Counter(known(rec) for rec in missing)
    return {
        'legacy_pairs': len(legacy_records),
        'legacy_pairs_found_in_v2': len(found),
        'legacy_pairs_reversed_in_v2': len(reversed_),
        'legacy_pairs_missing_in_v2': len(missing),
        'missing_by_ids_known_to_v2': {
            'both': missing_known[2], 'one': missing_known[1], 'neither': missing_known[0],
        },
        'v2_relations': v2_total,
        'v2_relations_without_legacy_pair': len(v2_unmatched),
        'v2_relations_reversed_in_legacy': len(v2_reversed),
        'v2_signs': {str(k): v for k, v in v2_signs.most_common()},
        'signs_compared': sign_agree + sign_disagree,
        'signs_agree': sign_agree,
        'signs_disagree': sign_disagree,
        'samples': {
            'legacy_missing_in_v2': [list(min(r)) for r in missing[:SAMPLE]],
            'v2_without_legacy_pair': v2_unmatched[:SAMPLE],
            'v2_reversed_in_legacy': v2_reversed[:SAMPLE],
            'sign_disagreements': disagreements,
        },
    }


def compare_sets(legacy: Iterable[Any], v2: Iterable[Any], what: str) -> dict[str, Any]:
    legacy_set = set(legacy)
    v2_set = set(v2)
    both = legacy_set & v2_set
    legacy_only = sorted(legacy_set - v2_set, key=str)
    v2_only = sorted(v2_set - legacy_set, key=str)
    show = sorted if what == 'complex' else list
    result = {
        'legacy_keys': len(legacy_set),
        'v2_keys': len(v2_set),
        'both': len(both),
        'legacy_only': len(legacy_only),
        'v2_only': len(v2_only),
        'samples': {
            'legacy_only': [show(x) for x in legacy_only[:SAMPLE]],
            'v2_only': [show(x) for x in v2_only[:SAMPLE]],
        },
    }
    if what == 'complex':
        # A complex whose member set differs slightly is not a missing complex.
        v2_members = set().union(*v2_set) if v2_set else set()
        result['legacy_only_with_all_members_in_v2'] = sum(
            1 for c in legacy_only if c <= v2_members
        )
    return result


def compare_keyed_complexes(legacy: Iterable[tuple], v2: Iterable[tuple]) -> dict[str, Any]:
    legacy_by_id: dict[str, set] = defaultdict(set)
    for key, members in legacy:
        legacy_by_id[key].add(members)
    v2_by_id = dict(v2)
    shared = sorted(legacy_by_id.keys() & v2_by_id.keys())
    nested = [k for k in shared if any(m.startswith('SIGNOR-') for m in v2_by_id[k])]
    flat = [k for k in shared if k not in set(nested)]
    differ = [k for k in flat if v2_by_id[k] not in legacy_by_id[k]]
    return {
        'legacy_ids': len(legacy_by_id),
        'v2_ids': len(v2_by_id),
        'both': len(shared),
        'legacy_only': len(legacy_by_id.keys() - v2_by_id.keys()),
        'v2_only': len(v2_by_id.keys() - legacy_by_id.keys()),
        'both_with_nested_members': len(nested),
        'both_flat_members_equal': len(flat) - len(differ),
        'both_flat_members_differ': len(differ),
        'samples': {
            'legacy_only': sorted(legacy_by_id.keys() - v2_by_id.keys())[:SAMPLE],
            'v2_only': sorted(v2_by_id.keys() - legacy_by_id.keys())[:SAMPLE],
            'members_differ': [
                {
                    'id': k,
                    'legacy': [sorted(m) for m in legacy_by_id[k]][:3],
                    'v2': sorted(v2_by_id[k]),
                }
                for k in differ[:SAMPLE]
            ],
        },
    }


def v2_keyed_complexes(records, id_namespace: str, prefer: str = 'uniprot') -> Iterable[tuple]:
    for entity in records:
        if entity is None:
            continue
        key = _primary_id(entity, id_namespace)
        members = [m.member for m in entity.membership or () if not m.is_parent]
        ids = {_primary_id(m, prefer) for m in members} - {None}
        yield key, frozenset(min(tokens(i), key=len) for i in ids)


# Resources -----------------------------------------------------------------
# Each legacy adapter reads the legacy function's least processed output, so
# that neither side has translated identifiers.

def _complex_tokens(cplex) -> frozenset[str]:
    return frozenset(min(tokens(c), key=len) for c in cplex.components)


def _signor_interactions():
    from pypath.inputs import signor

    for row in signor.signor_interactions(raw_records=True):
        yield Rel(tokens(row['IDA']), tokens(row['IDB']), sign_of(row['EFFECT']))


def _signor_complexes():
    from pypath.inputs import signor

    for cplex in signor.signor_complexes().values():
        for key in cplex.ids['SIGNOR']:
            yield key, _complex_tokens(cplex)


def _signalink():
    from pypath.inputs import signalink

    for rec in signalink.signalink_interactions(exclude_secondary=False):
        yield Rel(tokens(rec.id_a), tokens(rec.id_b), sign_of(rec.effect))


def _collectri():
    from pypath.inputs import collectri

    for rec in collectri.collectri_raw(mirna=True):
        yield Rel(tokens(rec.tf), tokens(rec.target), sign_of(rec.effect))


def _dorothea():
    from pypath.inputs import dorothea

    for rec in dorothea.dorothea_full_raw().itertuples():
        yield Rel(tokens(rec.tf), tokens(rec.target), sign_of(int(rec.mor)))


def _trrust():
    from pypath.inputs import trrust

    for rec in trrust.trrust_interactions('human'):
        yield Rel(
            tokens(rec.source_genesymbol), tokens(rec.target_genesymbol),
            sign_of(rec.effect),
        )


def _corum():
    from pypath.inputs import corum

    for cplex in corum.corum_complexes(9606).values():
        yield _complex_tokens(cplex)


def _phosphosite():
    from pypath.inputs import phosphosite

    for rec in phosphosite.phosphosite_enzyme_substrate(raw=True, organism='human'):
        yield Rel(tokens(rec['kinase']), tokens(rec['substrate']))


def _complexportal():
    from pypath.inputs import complexportal

    for cplex in complexportal.complexportal_complexes(9606).values():
        yield _complex_tokens(cplex)


def _cellphonedb():
    from pypath.inputs import cellphonedb

    def ids(partner):
        names = (
            [i for ids in partner.ids.values() for i in ids]
            if hasattr(partner, 'ids') else [partner]
        )
        # A metabolite system ("Acetylcholine_byCHAT") is its metabolite in v2.
        return id_tokens(names) | id_tokens(str(n).split('_by')[0] for n in names if '_by' in str(n))

    for rec in cellphonedb.cellphonedb_interactions():
        yield Rel(ids(rec.id_a), ids(rec.id_b))


def _cellphonedb_roles():
    from pypath.inputs import cellphonedb
    from pypath.resources.urls import urls

    for table, key in (('proteins', 'uniprot'), ('complexes', 'complex_name')):
        annotations = cellphonedb._cellphonedb_annotations(
            url=urls['cellphonedb_git'][table],
            # The source's own names: no translation to UniProt.
            name_method=lambda rec, key=key: rec[key],
        )
        for name, annotation in annotations.items():
            for role in ('receptor', 'secreted', 'transmembrane', 'peripheral', 'integrin'):
                if getattr(annotation, role):
                    yield (name, role)


def _v2_cellphonedb_roles():
    from pypath.inputs_v2 import cellphonedb

    for dataset, key in (('protein_roles', 'uniprot'), ('complex_roles', 'complex_name')):
        for row in getattr(cellphonedb.resource, dataset).raw():
            yield (row[key], row['role'])


def _legacy_rel(module: str, function: str, subject, object_, sign=None, **kwargs):
    """A legacy function's records as pairs: subject, object and sign are
    attribute names (or callables of the record)."""

    def get(rec, field):
        if callable(field):
            return field(rec)
        return rec[field] if isinstance(rec, dict) else getattr(rec, field)

    def run():
        import importlib

        records = getattr(importlib.import_module(f'pypath.inputs.{module}'), function)(**kwargs)
        for rec in records:
            yield Rel(
                id_tokens(get(rec, f) for f in _fields(subject)),
                id_tokens(get(rec, f) for f in _fields(object_)),
                sign_of(get(rec, sign)) if sign else None,
            )

    return run


def _fields(spec) -> tuple:
    return spec if isinstance(spec, tuple) else (spec,)


def _legacy_annotations(module: str, function: str, value, key=None, **kwargs):
    """A legacy annotation function (key: set of records) as pairs of the key
    and the value; value (and key) are attribute names or callables of a
    record, and a callable can give several values."""

    def run():
        import importlib

        records = getattr(importlib.import_module(f'pypath.inputs.{module}'), function)(**kwargs)
        for k, annotations in records.items():
            for a in annotations:
                values = value(a) if callable(value) else getattr(a, value)
                if not isinstance(values, (list, tuple, set)):
                    values = [values]
                subject = tokens(k) if key is None else id_tokens(_fields(getattr(a, key)))
                for v in values:
                    if v not in (None, ''):
                        yield Rel(subject, tokens(v))

    return run


def v2_pairs(records, taxon: str | None = None, transitive: bool = False) -> Iterable[Rel]:
    """Relations as pairs, and the associations and members of entities: an
    entity with its associated term, a member with its group, through nested
    groups (with ``transitive``, a member also with the groups above its own)."""

    def members(entity, above):
        ids = id_tokens(_ref_ids(entity))
        groups = [ids, *above] if transitive else [ids]
        for membership in entity.membership or ():
            member = id_tokens(_ref_ids(membership.member))
            for group in groups:
                yield Rel(member, group)
            yield from members(membership.member, groups)

    for record in records:
        if record is None:
            continue
        if hasattr(record, 'subject'):
            yield from v2_relations([record], taxon)
            continue
        ids = id_tokens(_ref_ids(record))
        for association in record.associations or ():
            yield Rel(ids, tokens(association.object.identifier))
        for relation in record.ontology_relations or ():
            yield Rel(ids, tokens(relation.object.identifier))
        yield from members(record, [])


def v2_id_pairs(records, key: str, values: tuple[str, ...]) -> Iterable[Rel]:
    """An entity's identifier of one namespace with each of its identifiers
    of others (``key`` and ``values`` are parts of namespace names)."""
    for record in records:
        if record is None:
            continue
        identifiers = [(name(i.type).lower(), i.value) for i in record.identifiers or ()]
        keys = id_tokens(v for t, v in identifiers if key in t)
        for t, v in identifiers:
            if any(n in t for n in values) and keys:
                yield Rel(keys, tokens(v))


def v2_annotation_pairs(records, term: str) -> Iterable[Rel]:
    """An entity with the values of one of its annotations."""
    for record in records:
        if record is None:
            continue
        ids = id_tokens(_ref_ids(record))
        for annotation in record.annotations or ():
            if term in name(annotation.term):
                yield Rel(ids, tokens(name(annotation.value)))


def _legacy_id_pairs(module: str, function: str, key, values: dict, limit: int | None = None, **kwargs):
    """Legacy metabolite records as pairs of their ID and each cross-reference
    (``values``: field to a prefix to strip, or None)."""

    def get(rec, field):
        return rec.get(field) if isinstance(rec, dict) else getattr(rec, field, None)

    def run():
        import importlib

        records = getattr(importlib.import_module(f'pypath.inputs.{module}'), function)(**kwargs)
        for rec in itertools.islice(records, limit):
            subject = tokens(key(rec) if callable(key) else get(rec, key))
            for field in values:
                found = get(rec, field)
                for v in found if isinstance(found, (list, tuple, set)) else [found]:
                    if v not in (None, '', 'CHEBI:'):
                        yield Rel(subject, tokens(v))

    return run


def _complex_or_id(value) -> str:
    """A legacy partner: its identifier, or a complex as its members."""
    if hasattr(value, 'components'):
        return ','.join(sorted(str(c) for c in value.components))
    return value


def datasets_starting(module: str, *prefixes: str) -> Callable[[], Iterable[Any]]:
    """The datasets of a module whose names start with one of the prefixes."""

    def run():
        import importlib

        resource = importlib.import_module(f'pypath.inputs_v2.{module}').resource
        names = [n for n in resource.datasets() if n.startswith(prefixes)]
        yield from datasets(module, *names)()

    return run


def _cellchat():
    from pypath.inputs import cellchatdb
    from pypath.inputs_v2.cellchat import metabolite

    frame = cellchatdb.cellchatdb_download(9606)['interaction']
    for rec in frame.itertuples():
        # A metabolite ligand ("PGE2-PTGES") is the metabolite in inputs_v2,
        # which keeps CellChat's abbreviation as a synonym.
        ligand = tokens(rec.ligand)
        if rec.annotation == 'Non-protein Signaling' and (abbreviation := metabolite(rec.ligand)):
            ligand |= tokens(abbreviation)
        yield Rel(ligand, tokens(rec.receptor))


def _matrisome():
    from pypath.inputs import matrisome

    for uniprot, annotations in matrisome.matrisome_annotations(9606).items():
        for annotation in annotations:
            yield (uniprot, annotation.mainclass)


def _v2_matrisome():
    for relation in datasets('matrisome', 'human')():
        if relation is None:
            continue
        for gene in _ref_ids(relation.subject)[:2]:
            yield (gene, relation.object.identifier)


SPECS: dict[str, Spec] = {}


def spec(*args, **kwargs):
    s = Spec(*args, **kwargs)
    SPECS[s.name] = s
    return s


spec(
    'signor', 'relation', _signor_interactions,
    lambda: v2_relations(datasets('signor', 'interactions')(), taxon='9606'),
    serve='signor',
    note='Legacy reads the CSV export of human interactions; inputs_v2 reads the causalTab export.',
)
spec(
    'signor_complexes', 'keyed_complex', _signor_complexes,
    lambda: v2_keyed_complexes(datasets('signor', 'complexes')(), 'signor'),
    serve='signor',
)
spec(
    'signalink', 'relation', _signalink,
    lambda: v2_relations(datasets('signalink', 'interactions')(), taxon='9606'),
    serve='signalink',
    note='Legacy with secondary sources included, as inputs_v2 has them; human pairs only.',
)
spec(
    'collectri', 'relation', _collectri,
    lambda: v2_relations(datasets('collectri', 'interactions')()),
    serve='collectri',
    note='Legacy collectri_raw with miRNA targets, before translation to UniProt.',
)
spec(
    'dorothea', 'relation', _dorothea,
    lambda: v2_relations(datasets('dorothea', 'level_a', 'level_b', 'level_c', 'level_d', 'level_e')()),
    serve='dorothea',
    note='All confidence levels on both sides.',
)
spec(
    'trrust', 'relation', _trrust,
    lambda: v2_relations(datasets('trrust', 'human_interactions')()),
    serve='trrust',
)
spec(
    'corum', 'complex', _corum,
    lambda: v2_complexes(datasets('corum', 'complexes')(), taxon='9606'),
    serve='corum',
)
spec(
    'phosphosite', 'relation', _phosphosite,
    lambda: v2_relations(datasets('phosphosite', 'enzyme_substrate')(), taxon='9606'),
    serve='phosphosite',
    note='Kinase-substrate pairs with both proteins human; legacy states no sign.',
)
spec(
    'complexportal', 'complex', _complexportal,
    lambda: v2_complexes(datasets('complexportal', 'human_complexes')()),
    serve='complexportal',
    note='Legacy reads the PSI-MI XML archive; inputs_v2 reads ComplexTAB.',
)
spec(
    'cellphonedb', 'relation', _cellphonedb,
    lambda: v2_relations(datasets('cellphonedb', 'interactions')()),
    serve='cellphonedb',
)
spec(
    'cellphonedb_roles', 'annotation', _cellphonedb_roles, _v2_cellphonedb_roles,
    serve='cellphonedb',
    note=(
        'Role flags of the proteins and complexes (a class counts as its role), by '
        "the source's names; the interactions need legacy's UniProt translation."
    ),
)
spec(
    'cellchat', 'relation', _cellchat,
    lambda: v2_relations(datasets('cellchat', 'human_interactions')()),
    serve='cellchat',
    aliases={
        'https://raw.githubusercontent.com/jinworks/CellChat/master/data/CellChatDB.human.rda':
        'https://raw.githubusercontent.com/jinworks/CellChat/main/data/CellChatDB.human.rda',
    },
    note=(
        'Legacy side is the interaction table of the R data file (its processing needs '
        'UniProt mapping); a metabolite ligand label also matches its abbreviation.'
    ),
)
spec(
    'matrisome', 'annotation', _matrisome, _v2_matrisome,
    serve='matrisome',
    note='Legacy reads the rescued XLS master list; inputs_v2 reads the current XLSX.',
)


spec(
    'ca1', 'relation',
    _legacy_rel('ca1', 'ca1_interactions', ('genesymbol_source', 'uniprot_source'),
                ('genesymbol_target', 'uniprot_target'), 'effect'),
    lambda: v2_relations(datasets('ca1', 'interactions')()),
    serve='ca1',
)
spec(
    'arn', 'relation',
    _legacy_rel('netbiol', 'arn_interactions', 'source_uniprot', 'target_uniprot', 'effect'),
    lambda: v2_relations(datasets('arn', 'interactions')()),
    serve='arn',
)
spec(
    'nrf2ome', 'relation',
    _legacy_rel('netbiol', 'nrf2ome_interactions', 'source_uniprot', 'target_uniprot', 'effect'),
    lambda: v2_relations(datasets('nrf2ome', 'interactions')()),
    serve='nrf2ome',
)
spec(
    'macrophage', 'relation',
    _legacy_rel('macrophage', 'macrophage_interactions', 'source_genesymbol', 'target_genesymbol'),
    lambda: v2_relations(datasets('macrophage', 'interactions')()),
    serve='macrophage',
)
spec(
    'hprd', 'relation',
    _legacy_rel('hprd', 'hprd_enzyme_substrate', 'kinase', ('substrate', 'substrate_refseqp'),
                in_vivo=False),
    lambda: v2_relations(datasets('hprd', 'enzyme_substrate')()),
    serve='hprd',
    note='Enzyme-substrate pairs, legacy with the in vitro ones as inputs_v2; legacy gives no sign.',
)
spec(
    'li2012', 'relation',
    _legacy_rel('li2012', 'li2012_enzyme_substrate', 'kinase', 'substrate'),
    lambda: v2_relations(datasets('li2012', 'kinase_substrate')()),
    serve='li2012',
)
spec(
    'phosphonetworks', 'relation',
    _legacy_rel('phosphonetworks', 'phosphonetworks_enzyme_substrate', 'kinase', 'substrate'),
    lambda: v2_relations(datasets('phosphonetworks', 'enzyme_substrate')()),
    serve='phosphonetworks',
)
spec(
    'phosphopoint', 'relation',
    _legacy_rel('phosphopoint', 'phosphopoint_interactions',
                ('source_genesymbol', 'source_entrez'), ('target_genesymbol', 'target_entrez')),
    lambda: v2_relations(datasets('phosphopoint', 'interactions', 'enzyme_substrate')()),
    serve='phosphopoint',
    note='inputs_v2 splits the pairs into interactions and enzyme-substrate by category.',
)
spec(
    'dbptm', 'relation',
    lambda: (
        r for r in _legacy_rel('dbptm', 'dbptm_enzyme_substrate', 'kinase', 'substrate')()
        if r.subject
    ),
    lambda: v2_relations(datasets('dbptm', 'enzyme_substrate')()),
    serve='dbptm',
    note='Legacy records with an enzyme; legacy also lists the PTM sites without one.',
)
spec(
    'guidetopharma', 'relation',
    _legacy_rel('guide2pharma', 'guide2pharma_interactions', 'ligand', 'target', 'effect',
                endogenous=False),
    lambda: v2_relations(datasets('guidetopharma', 'interactions')()),
    serve='guidetopharma',
    note='Legacy with all ligands, not only endogenous ones; human.',
)
spec(
    'ramilowski2015', 'relation',
    _legacy_rel('ramilowski2015', 'ramilowski_interactions', 'ligand', 'receptor'),
    lambda: v2_relations(datasets('ramilowski2015', 'interactions')()),
    serve='ramilowski2015',
)
spec(
    'lrdb', 'relation',
    _legacy_rel('lrdb', 'lrdb_interactions', 'ligand_genesymbol', 'receptor_genesymbol'),
    lambda: v2_relations(datasets('lrdb', 'interactions')()),
    serve='lrdb',
)
spec(
    'celltalkdb', 'relation',
    _legacy_rel('celltalkdb', 'celltalkdb_interactions', 'ligand_genesymbol', 'receptor_genesymbol'),
    lambda: v2_relations(datasets('celltalkdb', 'human_interactions')()),
    serve='celltalkdb',
)

# Third batch.
spec('adhesome', 'relation', _legacy_rel('adhesome', 'adhesome_interactions', 'source', 'target', 'effect'),
     lambda: v2_relations(datasets('adhesome', 'interactions', 'website_interactions')()), serve='adhesome')
spec('cancerdrugsdb', 'relation',
     _legacy_rel('cancerdrugsdb', 'cancerdrugsdb_interactions',
                 ('drug_pubchem', 'drug_chembl', 'drug_drugbank', 'drug_label'), 'target_uniprot'),
     lambda: v2_relations(datasets('cancerdrugsdb', 'drug_targets')()), serve='cancerdrugsdb')
spec('cellcall', 'relation', _legacy_rel('cellcall', 'cellcall_interactions', 'ligand_uniprot', 'receptor_uniprot'),
     lambda: v2_relations(datasets('cellcall', 'human_core_interactions', 'human_extended_interactions')()),
     serve='cellcall')
spec('cellcellinteractions', 'relation',
     _legacy_annotations('cellcellinteractions', 'cellcellinteractions_annotations', 'mainclass'),
     lambda: v2_pairs(datasets('cellcellinteractions', 'genes')()), serve='cellcellinteractions')
spec('celltypist', 'relation', _legacy_annotations('celltypist', 'celltypist_annotations', 'cell_ontology'),
     lambda: v2_pairs(datasets('celltypist', 'markers')()), serve='celltypist')
spec('compleat', 'complex',
     lambda: (_complex_tokens(c) for c in __import__('pypath.inputs.compleat', fromlist=['x']).compleat_complexes().values()),
     lambda: v2_complexes(datasets('compleat', 'literature_complexes', 'predicted_complexes')()), serve='compleat')
spec('cpad', 'relation', _legacy_annotations('cpad', 'cpad_annotations', 'cancer', key='pathway'),
     lambda: v2_pairs(datasets('cpad', 'entries')()), serve='cpad', note='Pathway and cancer.')
spec('cpdb', 'relation', _legacy_rel('cpdb', 'cpdb_interactions', lambda r: r[0], lambda r: r[1]),
     lambda: v2_relations(datasets('cpdb', 'interactions')()), serve='cpdb')
spec('cytosig', 'relation', _legacy_annotations('cytosig', 'cytosig_annotations', 'cytokine'),
     lambda: v2_pairs(datasets('cytosig', 'signature')()), serve='cytosig')
spec('domino', 'relation', _legacy_rel('domino', 'domino_interactions', 'uniprot_A', 'uniprot_B'),
     lambda: v2_relations(datasets('domino', 'interactions')()), serve='domino')
spec('elm', 'relation', _legacy_rel('elm', 'elm_interactions', 'uniprot_motif', 'uniprot_domain'),
     lambda: v2_relations(datasets('elm', 'interactions')()), serve='elm')
spec('encode', 'relation', _legacy_rel('encode', 'encode_tf_mirna_interactions', 'tf_genesymbol', 'mirna'),
     lambda: v2_relations(datasets('encode', 'tf_mirna')()), serve='encode',
     note='The TF-miRNA table; legacy has no function for the TF-gene tables.')
spec('gpcrdb', 'relation', _legacy_annotations('gpcrdb', 'gpcrdb_annotations', 'subfamily'),
     lambda: v2_pairs(datasets('gpcrdb', 'families')()), serve='gpcrdb', note='Receptor and its subfamily.')
spec('havugimana', 'complex',
     lambda: (_complex_tokens(c) for c in __import__('pypath.inputs.havugimana', fromlist=['x']).havugimana_complexes().values()),
     lambda: v2_complexes(datasets('havugimana', 'complexes')()), serve='havugimana')
def _score(relation) -> float:
    for annotation in relation.annotations or ():
        if 'confidence_score' in name(annotation.term):
            value = annotation.value
            value = getattr(value, 'has_numeric_value', None) or (
                json.loads(value).get('has_numeric_value') if isinstance(value, str) and value.startswith('{') else value
            )
            try:
                return float(value)
            except (TypeError, ValueError):
                return 0.0
    return 0.0


spec('hippie', 'relation', _legacy_rel('hippie', 'hippie_interactions', 'id_a', 'id_b'),
     lambda: v2_relations(r for r in datasets('hippie', 'interactions')() if r and _score(r) >= .75),
     serve='hippie',
     note=('Legacy keeps the interactions of score 0.75 or more (it fails without the threshold); '
           'v2 keeps all with their scores and is compared on the same threshold.'))
spec('hpo', 'relation', _legacy_annotations('hpo', 'hpo_annotations', 'hpo_id'),
     lambda: v2_pairs(datasets('hpo', 'annotations')()), serve='hpo', note='Gene and phenotype term.')
spec('htri', 'relation',
     _legacy_rel('htri', 'htri_interactions', ('entrez_tf', 'genesymbol_tf'), ('entrez_target', 'genesymbol_target')),
     lambda: v2_relations(datasets('htri', 'interactions')()), serve='htri')
spec('humap', 'complex',
     lambda: (_complex_tokens(c) for c in __import__('pypath.inputs.humap', fromlist=['x']).humap2and3_complexes().values()),
     lambda: v2_complexes(datasets('humap', 'humap2_complexes', 'humap3_complexes')()), serve='humap')
spec('lambert2018', 'relation',
     _legacy_annotations('lambert2018', 'lambert2018_annotations',
                         lambda a: (['human_transcription_factors'] if a.is_tf else []) + [a.binding_domain]),
     lambda: v2_pairs(datasets('lambert2018', 'transcription_factors', 'families')()), serve='lambert2018',
     note='Gene with the TF list and its DNA-binding domain family.')
spec('lncdisease', 'relation', _legacy_rel('lncdisease', 'lncdisease_interactions', 'source', 'target'),
     lambda: v2_relations(datasets('lncdisease', 'interactions')()), serve='lncdisease')
spec('lncrnadb', 'relation', _legacy_rel('lncrnadb', 'lncrnadb_interactions', 'lncrna', 'partner'),
     lambda: v2_relations(datasets('lncrnadb', 'associations')()), serve='lncrnadb')
spec('mcam', 'relation',
     lambda: (Rel(tokens(u), tokens('human_cell_adhesion_molecule'))
              for u in __import__('pypath.inputs.mcam', fromlist=['x']).mcam_cell_adhesion_molecules()),
     lambda: v2_pairs(datasets('mcam', 'proteins')()), serve='mcam')
spec('mirdeathdb', 'relation', _legacy_rel('mirdeathdb', 'mirdeathdb_interactions', 'mirna', 'target_entrez'),
     lambda: v2_relations(datasets('mirdeathdb', 'interactions')()), serve='mirdeathdb')
spec('mirecords', 'relation',
     _legacy_rel('mirecords', 'mirecords_interactions', 'mirna_name', ('target_refseq', 'target_genesymbol')),
     lambda: v2_relations(datasets('mirecords', 'interactions')()), serve='mirecords')
spec('mppi', 'relation', _legacy_rel('mppi', 'mppi_interactions', lambda r: r[2], lambda r: r[6]),
     lambda: v2_relations(datasets('mppi', 'interactions')()), serve='mppi')
spec('ncrdeathdb', 'relation', _legacy_rel('ncrdeathdb', 'ncrdeathdb_interactions', 'ncrna', 'target_gene'),
     lambda: v2_relations(datasets('ncrdeathdb', 'interactions', 'homology_predictions')()), serve='ncrdeathdb')
spec('netpath', 'relation',
     _legacy_rel('netpath', 'netpath_interactions', lambda r: (r[0], r[1]), lambda r: (r[2], r[3])),
     lambda: v2_relations(datasets('netpath', 'interactions')()), serve='netpath')
spec('oreganno', 'relation', _legacy_rel('oreganno', 'oreganno_interactions', 'tf', 'target'),
     lambda: v2_relations(datasets('oreganno', 'interactions')()), serve='oreganno')
spec('panglaodb', 'relation', _legacy_annotations('panglaodb', 'panglaodb_annotations', 'cell_type'),
     lambda: v2_pairs(datasets('panglaodb', 'markers')()), serve='panglaodb', note='Gene and cell type.')
spec('pazar', 'relation', _legacy_rel('pazar', 'pazar_interactions', 'tf', 'target'),
     lambda: v2_relations(datasets('pazar', 'interactions')()), serve='pazar')
spec('pdzbase', 'relation', _legacy_rel('pdzbase', 'pdzbase_interactions', 'uniprot_pdz', 'uniprot_ligand'),
     lambda: v2_relations(datasets('pdzbase', 'interactions')()), serve='pdzbase')
spec('phobius', 'relation',
     _legacy_annotations('phobius', 'phobius_annotations',
                         lambda a: (['transmembrane'] if a.tm_helices else [])
                         + (['signal_peptide'] if a.signal_peptide else [])),
     lambda: v2_pairs(datasets('phobius', 'predictions')()), serve='phobius')
spec('phosphatome', 'relation',
     _legacy_annotations('phosphatome', 'phosphatome_annotations',
                         lambda a: [a.fold, a.family, a.subfamily]
                         + (['Protein substrates'] if a.has_protein_substrates else [])
                         + (['Non-protein substrates'] if a.has_non_protein_substrates else [])),
     lambda: v2_pairs(datasets('phosphatome', 'class_annotations', 'families')(), transitive=True),
     serve='phosphatome')
spec('tfcensus', 'relation',
     _legacy_annotations('tfcensus', 'tfcensus_annotations', lambda a: ['human_transcription_factors']),
     lambda: v2_pairs(datasets('tfcensus', 'transcription_factors')()), serve='tfcensus')
spec('topdb', 'relation',
     _legacy_annotations('topdb', 'topdb_annotations', lambda a: ['protein']),
     lambda: (Rel(r.subject, tokens('protein')) for r in v2_pairs(datasets('topdb', 'proteins')())),
     serve='topdb', note='The proteins only: legacy reads the regions, v2 the protein classes.')
spec('wang', 'relation',
     _legacy_rel('wang', 'wang_interactions', ('entrez_source', 'genesymbol_source'),
                 ('entrez_target', 'genesymbol_target'), 'effect'),
     lambda: v2_relations(datasets('wang', 'interactions')()), serve='wang')
spec('zhong2015', 'relation',
     _legacy_annotations('zhong2015', 'zhong2015_annotations',
                         lambda a: {v: k for k, v in __import__('pypath.inputs_v2.zhong2015', fromlist=['x']).CAM_TYPES.items()}.get(a.type, a.type)),
     lambda: v2_pairs(datasets('zhong2015', 'genes')()), serve='zhong2015')
spec('laudanna_directions', 'relation',
     _legacy_rel('laudanna', 'laudanna_directions', 'source_genesymbol', 'target_genesymbol'),
     lambda: v2_relations(datasets('laudanna', 'signaling_flow')()), serve='laudanna')
spec('laudanna_effects', 'relation',
     _legacy_rel('laudanna', 'laudanna_effects', 'source_genesymbol', 'target_genesymbol', 'effect'),
     lambda: v2_relations(datasets('laudanna', 'signaling_direction')()), serve='laudanna')
spec('instruct', 'relation',
     _legacy_rel('instruct', 'get_instruct', lambda r: r['uniprots'][0], lambda r: r['uniprots'][1]),
     lambda: v2_relations(datasets('instruct', 'interactions')()), serve='instruct')
spec('tcdb', 'relation',
     lambda: (Rel(tokens(u), tokens(c[0])) for u, c in __import__('pypath.inputs.tcdb', fromlist=['x']).tcdb_classes().items()),
     lambda: v2_annotation_pairs(datasets('tcdb', 'transporters')(), 'topic'), serve='tcdb',
     note='Protein and TC number.')
spec('recon3d_metabolites', 'relation',
     _legacy_id_pairs('recon3d', 'recon3d_metabolites', 'base_id', ('chebi', 'hmdb', 'kegg', 'metanetx')),
     lambda: v2_id_pairs(
         (m.member for r in datasets('recon3d', 'reactions')() if r for m in r.membership or ()),
         'bigg', ('chebi', 'hmdb', 'kegg', 'metanetx'),
     ),
     serve='recon3d',
     note='BiGG metabolite and its cross-references; v2 has the metabolites of reactions as their members.')
spec('recon3d_reactions', 'relation',
     lambda: (Rel(tokens(r['id']), tokens(m.rsplit('_', 1)[0]))
              for r in __import__('pypath.inputs.recon3d', fromlist=['x']).recon3d_reactions() for m in r['metabolites']),
     lambda: (Rel(r.object, r.subject) for r in v2_pairs(datasets('recon3d', 'reactions')())),
     serve='recon3d', note='Reaction and its metabolites (BiGG).')
spec('imm1415_metabolites', 'relation',
     _legacy_id_pairs('imm1415', 'imm1415_metabolites', 'base_id', ('chebi', 'hmdb', 'kegg', 'metanetx')),
     lambda: v2_id_pairs(datasets('imm1415', 'metabolites')(), 'bigg', ('chebi', 'hmdb', 'kegg', 'metanetx')),
     serve='imm1415', note='BiGG metabolite and its cross-references.')
spec('lipidmaps', 'relation',
     _legacy_id_pairs('lipidmaps', 'lipidmaps_lipids', 'id', ('chebi', 'inchikey', 'pubchem', 'hmdb', 'swisslipids')),
     lambda: v2_id_pairs(datasets('lipidmaps', 'lipids')(), 'lipidmaps', ('chebi', 'inchikey', 'pubchem', 'hmdb', 'swisslipids')),
     serve='lipidmaps', note='LIPID MAPS ID and its cross-references.')
spec('swisslipids', 'relation',
     _legacy_id_pairs('swisslipids', 'swisslipids_lipids', 'id', ('chebi', 'inchikey', 'lipidmaps', 'hmdb', 'metanetx')),
     lambda: v2_id_pairs(datasets('swisslipids', 'lipids')(), 'swisslipids', ('chebi', 'inchikey', 'lipidmaps', 'hmdb', 'metanetx')),
     serve='swisslipids', note='SwissLipids ID and its cross-references.')
spec('kea', 'relation', _legacy_rel('kea', 'kea_interactions', 'enzyme', 'substrate'),
     lambda: v2_relations(datasets('kea', 'enzyme_substrate')()), serve='kea')
spec('icellnet', 'relation',
     _legacy_rel('icellnet', 'icellnet_interactions',
                 lambda r: _complex_or_id(r.ligand), lambda r: _complex_or_id(r.receptor)),
     lambda: v2_relations(datasets('icellnet', 'interactions')()), serve='icellnet')
spec('spike', 'relation',
     _legacy_rel('spike', 'spike_interactions', ('entrez_a', 'genesymbol_a'), ('entrez_b', 'genesymbol_b')),
     lambda: v2_relations(datasets_starting('spike', 'regulations', 'interactions')()), serve='spike',
     note='Regulations and interactions of all quality levels.')
spec('spike_complexes', 'complex',
     lambda: (_complex_tokens(c) for c in __import__('pypath.inputs.spike', fromlist=['x']).spike_complexes().values()),
     lambda: v2_complexes(datasets('spike', 'complexes')()), serve='spike')
spec('psi_mi', 'relation',
     lambda: (Rel(tokens(r.id), tokens(parent))
              for r in __import__('pypath.inputs.psimi', fromlist=['x']).psimi_ontology()
              for parent in r.parent_ids or ()),
     lambda: v2_pairs(datasets('psi_mi', 'terms')()), serve='psi_mi',
     aliases={
         'https://raw.githubusercontent.com/HUPO-PSI/psi-mi-CV/refs/heads/master/psi-mi.obo':
         'https://raw.githubusercontent.com/HUPO-PSI/psi-mi-CV/master/psi-mi.obo',
     },
     note='Terms and their is_a parents.')
spec('italk', 'relation', _legacy_rel('italk', 'italk_interactions', 'ligand', 'receptor'),
     lambda: v2_relations(datasets('italk', 'interactions')()), serve='italk',
     aliases={'https://github.com/Coolgenome/iTALK/blob/master/data/LR_database.rda?raw=true':
              'https://raw.githubusercontent.com/Coolgenome/iTALK/master/data/LR_database.rda'})
spec('talklr', 'relation', _legacy_rel('talklr', 'talklr_interactions', 'ligand', 'receptor'),
     lambda: v2_relations(datasets('talklr', 'interactions', 'putative_interactions')()), serve='talklr',
     aliases={'https://github.com/yuliangwang/talklr/raw/master/data/receptor_ligand.rda':
              'https://raw.githubusercontent.com/yuliangwang/talklr/master/data/receptor_ligand.rda'})
spec('progeny', 'relation', _legacy_annotations('progeny', 'progeny_annotations', 'pathway'),
     lambda: v2_pairs(datasets('progeny', 'responses')()), serve='progeny',
     aliases={'https://github.com/saezlab/progeny/raw/master/data/model_human_full.rda':
              'https://raw.githubusercontent.com/saezlab/progeny/master/data/model_human_full.rda'},
     note='Gene and pathway.')
spec('kinasedotcom', 'relation',
     _legacy_annotations('kinasedotcom', 'kinasedotcom_annotations', lambda a: [a.group, a.family, a.subfamily]),
     lambda: v2_pairs(datasets('kinasedotcom', 'families')(), transitive=True), serve='kinasedotcom',
     aliases={'https://rescued.omnipathdb.org/Kincat_Hsap.08.02.xls':
              'http://kinase.com/static/colt/data/human/kinome/tables/Kincat_Hsap.08.02.xls'},
     note="Kinase and its group, family and subfamily; legacy reads OmniPath's copy of the file.")
spec('mir2disease', 'relation', _legacy_rel('mir2disease', 'mir2disease_interactions', 'mirna', 'target_genesymbol'),
     lambda: v2_relations(datasets('mir2disease', 'mirna_targets')()), serve='mir2disease',
     aliases={'http://rescued.omnipathdb.org/miR2Disease_miRtar.txt':
              'http://watson.compbio.iupui.edu:8080/miR2Disease/download/miRtar.txt'},
     note="Legacy reads OmniPath's copy of the file.")
spec('cancercellmap', 'relation', _legacy_rel('cancercellmap', 'ccmap_interactions', 'source_uniprot', 'target_uniprot'),
     lambda: v2_relations(datasets('cancercellmap', 'interactions')()), serve='cancercellmap')
spec('dip', 'relation', _legacy_rel('dip', 'dip_interactions', lambda r: r[0], lambda r: r[1]),
     lambda: v2_relations(datasets('dip', 'interactions')()), serve='dip',
     aliases={'https://rescued.omnipathdb.org/dip/Hsapi20170205CR.txt.gz':
              'https://rescued.omnipathdb.org/dip/Hsapi20170205.txt.gz'},
     note='Legacy reads the core subset of the release by default; here it reads the full file as v2.')
spec('hmdb', 'relation',
     _legacy_id_pairs('hmdb.metabolites', 'raw', 'accession', ('chebi_id', 'inchikey', 'kegg_id', 'pubchem_compound_id'),
                      limit=20000),
     lambda: v2_id_pairs(datasets('hmdb', 'metabolites')(), 'hmdb', ('chebi', 'inchikey', 'kegg', 'pubchem')),
     serve='hmdb', v2_subset=20000,
     note='HMDB accession and its cross-references, of the first 20,000 metabolites of the file on both sides.')


# Driver --------------------------------------------------------------------

def run(spec: Spec) -> dict[str, Any]:
    """Compare one resource. inputs_v2 runs first: some legacy modules patch
    shared packages (pypath.inputs.rdata changes the rdata converter)."""
    _files.clear()
    _served.clear()
    serve(spec.serve, spec.aliases)
    global _V2_SUBSET
    _V2_SUBSET = spec.v2_subset
    result: dict[str, Any] = {'resource': spec.name, 'kind': spec.kind, 'note': spec.note}
    if spec.v2_subset:
        result['v2_subset'] = spec.v2_subset
    outputs = {}
    for side, produce in (('v2', spec.v2), ('legacy', spec.legacy)):
        started = time.time()
        try:
            outputs[side] = list(produce())
        except Exception as error:
            result['error'] = f'{side}: {error!r}'
            result['traceback'] = traceback.format_exc()
            result['files'] = list(_files)
            return result
        result[f'{side}_seconds'] = round(time.time() - started, 1)
    legacy, v2 = outputs['legacy'], outputs['v2']
    if spec.kind == 'relation':
        result |= compare_relations(legacy, v2)
    elif spec.kind == 'keyed_complex':
        result |= compare_keyed_complexes(legacy, v2)
    else:
        result |= compare_sets(legacy, v2, spec.kind)
    result['files'] = list(_files)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('resources', nargs='*', help=f'default: all of {", ".join(SPECS)}')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--in-process', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    names = args.resources or list(SPECS)
    if not args.in_process:
        # One process per resource, so that no legacy module's patches reach
        # another resource.
        for name in names:
            subprocess.run(
                [sys.executable, __file__, name, '--out', str(args.out), '--in-process'],
                check=False,
            )
        return

    setup(args.out)
    for name in names:
        print(f'== {name}', file=sys.stderr, flush=True)
        result = run(SPECS[name])
        (args.out / f'{name}.json').write_text(json.dumps(result, indent=2, default=str) + '\n')
        summary = {k: v for k, v in result.items() if k not in {'samples', 'files', 'traceback'}}
        print(json.dumps(summary, indent=2, default=str), file=sys.stderr, flush=True)


if __name__ == '__main__':
    main()
