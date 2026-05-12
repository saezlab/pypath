"""
Shared building blocks for inputs_v2 datasets.
"""

from __future__ import annotations

from collections.abc import Callable, Generator
import csv
from dataclasses import dataclass
import functools
import hashlib
import itertools
import json
import os
from pathlib import Path
from typing import Any, Literal, Protocol

from pypath.inputs_v2.raw_records import (
    ProvenancedRecord,
    RawRecordProvenance,
    RawSnapshot,
    accept_raw_snapshot,
    default_raw_records_root,
    iter_changed_raw_record_dicts,
    iter_raw_record_dicts,
    materialize_raw_records,
    reuse_raw_snapshot_if_unchanged,
)
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    OntologyAnnotationCv,
    OntologyCv,
    ResourceAnnotationCv,
    ResourceCv,
    UpdateCategoryCV,
)
from pypath.internals.ontology_schema import (
    OntologyDocument,
    OntologyTerm,
    OntologyTypedef,
)
from pypath.internals.silver_schema import Annotation, Entity, Identifier
from pypath.share.downloads import download_and_open


class _Resolver(Protocol):
    def __call__(self, **kwargs: Any) -> str: ...


def _resolve(value: str | _Resolver, **kwargs: Any) -> str:
    if callable(value):
        return value(**kwargs)
    return value


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in sorted(value.items())}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(item) for item in value]
    try:
        json.dumps(value, sort_keys=True)
    except TypeError:
        return repr(value)
    return value


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _download_fingerprint(
    download: Download | None,
    opener: Any,
    kwargs: dict[str, Any],
) -> dict[str, Any] | None:
    if download is None or opener is None:
        return None
    path_value = getattr(opener, 'path', None)
    if not path_value:
        return None
    path = Path(path_value)
    if not path.exists():
        return None
    stat = path.stat()
    download_kwargs = dict(download.download_kwargs or {})
    return {
        'path': str(path),
        'sha256': _file_sha256(path),
        'size': stat.st_size,
        'mtime_ns': stat.st_mtime_ns,
        'url': _resolve(download.url, **kwargs),
        'filename': _resolve(download.filename, **kwargs),
        'subfolder': download.subfolder,
        'needed': list(download.needed or []),
        'download_kwargs': _jsonable(download_kwargs),
    }


def _download_fingerprint_from_cache(
    download: Download | None,
    kwargs: dict[str, Any],
) -> dict[str, Any] | None:
    if download is None:
        return None
    filename = _resolve(download.filename, **kwargs)
    cache_root = Path(os.environ.get('PYPATH_DOWNLOAD_DATADIR', 'pypath-data'))
    path = cache_root / download.subfolder / filename
    if not path.exists():
        return None
    stat = path.stat()
    download_kwargs = dict(download.download_kwargs or {})
    return {
        'path': str(path),
        'sha256': _file_sha256(path),
        'size': stat.st_size,
        'mtime_ns': stat.st_mtime_ns,
        'url': _resolve(download.url, **kwargs),
        'filename': filename,
        'subfolder': download.subfolder,
        'needed': list(download.needed or []),
        'download_kwargs': _jsonable(download_kwargs),
    }


def _callable_name(parser: Callable[..., Any]) -> Any:
    if isinstance(parser, functools.partial):
        return {
            'callable': _callable_name(parser.func),
            'args': _jsonable(parser.args),
            'keywords': _jsonable(parser.keywords or {}),
        }
    return f'{getattr(parser, "__module__", "")}.{getattr(parser, "__qualname__", repr(parser))}'


def _parser_contract(parser: Callable[..., Any], kwargs: dict[str, Any]) -> dict[str, Any]:
    return {
        'raw_parser': _callable_name(parser),
        'kwargs': _jsonable(kwargs),
    }


@dataclass(frozen=True)
class ResourceConfig:
    id: ResourceCv
    name: str
    url: str
    license: LicenseCV
    update_category: UpdateCategoryCV
    description: str
    pubmed: str | None = None
    primary_category: str | None = None
    annotation_ontologies: tuple[OntologyCv, ...] = ()
    resource_kind: str = 'data_resource'

    def metadata(self) -> Entity:
        annotations = [
            Annotation(term=ResourceAnnotationCv.LICENSE, value=str(self.license)),
            Annotation(term=ResourceAnnotationCv.UPDATE_CATEGORY, value=str(self.update_category)),
        ]
        if self.pubmed:
            annotations.append(Annotation(term=IdentifierNamespaceCv.PUBMED, value=self.pubmed))
        annotations.extend([
            Annotation(term=ResourceAnnotationCv.URL, value=self.url),
            Annotation(term=ResourceAnnotationCv.DESCRIPTION, value=self.description),
        ])

        return Entity(
            type=EntityTypeCv.CV_TERM,
            identifiers=[
                Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=self.id),
                Identifier(type=IdentifierNamespaceCv.NAME, value=self.name),
            ],
            annotations=annotations,
        )


@dataclass(frozen=True)
class Download:
    url: str | _Resolver
    filename: str | _Resolver
    subfolder: str
    large: bool = True
    encoding: str | None = 'utf-8'
    default_mode: str = 'r'
    ext: str | None = None
    needed: list[str] | None = None
    download_kwargs: dict[str, Any] | None = None

    def open(self, *, force_refresh: bool = False, **kwargs: Any):
        download_kwargs = dict(self.download_kwargs or {})
        url = _resolve(self.url, **kwargs)
        filename = _resolve(self.filename, **kwargs)
        return download_and_open(
            url=url,
            filename=filename,
            subfolder=self.subfolder,
            large=self.large,
            encoding=self.encoding,
            default_mode=self.default_mode,
            ext=self.ext,
            needed=self.needed,
            force_download=force_refresh,
            **download_kwargs,
        )


DatasetKind = Literal['ontology', 'id_translation']


class Dataset:
    """A Lego brick: download + raw parsing + mapping to Entities."""

    def __init__(
        self,
        download: Download | None,
        mapper: Callable[[dict[str, Any]], Entity],
        raw_parser: Callable[..., Generator[dict[str, Any], None, None]],
        *,
        kind: DatasetKind | None = None,
    ) -> None:
        self.download = download
        self.mapper = mapper
        self._raw_parser = raw_parser
        self.kind = kind
        self._last_raw_snapshot: RawSnapshot | None = None

    def preparse(
        self,
        *,
        source: str,
        dataset: str,
        raw_records_root: str | Path | None = None,
        force_refresh: bool = False,
        **kwargs: Any,
    ) -> RawSnapshot:
        max_lines = kwargs.pop('max_lines', None)
        contract_kwargs = dict(kwargs)
        if max_lines is not None:
            contract_kwargs['max_lines'] = int(max_lines)
        parser_contract = _parser_contract(self._raw_parser, contract_kwargs)
        output_root = Path(raw_records_root) if raw_records_root else default_raw_records_root()
        if not force_refresh:
            cached_fingerprint = _download_fingerprint_from_cache(self.download, kwargs)
            reusable = reuse_raw_snapshot_if_unchanged(
                source=source,
                dataset=dataset,
                output_root=output_root,
                download_fingerprint=cached_fingerprint,
                parser_contract=parser_contract,
            )
            if reusable is not None:
                self._last_raw_snapshot = reusable
                return reusable
        opener = self.download.open(force_refresh=force_refresh, **kwargs) if self.download else None
        download_fingerprint = _download_fingerprint(self.download, opener, kwargs)
        if not force_refresh:
            reusable = reuse_raw_snapshot_if_unchanged(
                source=source,
                dataset=dataset,
                output_root=output_root,
                download_fingerprint=download_fingerprint,
                parser_contract=parser_contract,
            )
            if reusable is not None:
                self._last_raw_snapshot = reusable
                return reusable
        records = self._raw_parser(opener, force_refresh=force_refresh, **kwargs)
        if max_lines is not None:
            records = itertools.islice(records, int(max_lines))
        self._last_raw_snapshot = materialize_raw_records(
            records=records,
            source=source,
            dataset=dataset,
            output_root=output_root,
            download_fingerprint=download_fingerprint,
            parser_contract=parser_contract,
        )
        return self._last_raw_snapshot

    def raw(
        self,
        force_refresh: bool = False,
        *,
        source: str | None = None,
        dataset: str | None = None,
        raw_records_root: str | Path | None = None,
        use_preparse: bool = False,
        changed_only: bool = False,
        raw_snapshot: RawSnapshot | None = None,
        **kwargs: Any,
    ) -> Generator[dict[str, Any], None, None]:
        if use_preparse:
            if not source or not dataset:
                raise ValueError('source and dataset are required when use_preparse=True')
            snapshot = raw_snapshot or self.preparse(
                source=source,
                dataset=dataset,
                raw_records_root=raw_records_root,
                force_refresh=force_refresh,
                **kwargs,
            )
            if raw_snapshot is not None:
                self._last_raw_snapshot = raw_snapshot
            if changed_only:
                yield from iter_changed_raw_record_dicts(snapshot.records_path, snapshot.delta_path)
            else:
                yield from iter_raw_record_dicts(snapshot.records_path)
            return

        opener = self.download.open(force_refresh=force_refresh, **kwargs) if self.download else None
        yield from self._raw_parser(opener, force_refresh=force_refresh, **kwargs)

    def accept_last_preparse(self) -> None:
        if self._last_raw_snapshot is not None:
            accept_raw_snapshot(self._last_raw_snapshot)
            self._last_raw_snapshot = None

    def __call__(self, force_refresh: bool = False, **kwargs: Any) -> Generator[Entity | ProvenancedRecord, None, None]:
        if kwargs.get('use_preparse'):
            source = kwargs.get('source')
            dataset = kwargs.get('dataset')
            if not source or not dataset:
                raise ValueError('source and dataset are required when use_preparse=True')
            snapshot = kwargs.get('raw_snapshot') or self.preparse(
                source=source,
                dataset=dataset,
                raw_records_root=kwargs.get('raw_records_root'),
                force_refresh=force_refresh,
                **{
                    key: value
                    for key, value in kwargs.items()
                    if key not in {'source', 'dataset', 'raw_records_root', 'use_preparse', 'changed_only', 'raw_snapshot'}
                },
            )
            if kwargs.get('raw_snapshot') is not None:
                self._last_raw_snapshot = snapshot
            raw_rows = (
                iter_changed_raw_record_dicts(
                    snapshot.records_path,
                    snapshot.delta_path,
                    include_metadata=True,
                )
                if kwargs.get('changed_only')
                else iter_raw_record_dicts(snapshot.records_path, include_metadata=True)
            )
            for raw_row in raw_rows:
                record = {
                    k: v
                    for k, v in raw_row.items()
                    if k not in {'_source', '_dataset', '_raw_record_key', '_raw_record_id'}
                }
                yield ProvenancedRecord(
                    record=self.mapper(record),
                    provenance=RawRecordProvenance(
                        source=str(raw_row.get('_source') or source),
                        dataset=str(raw_row.get('_dataset') or dataset),
                        snapshot_id=snapshot.snapshot_id,
                        raw_record_key=str(raw_row['_raw_record_key']),
                        raw_record_id=int(raw_row['_raw_record_id']),
                        raw_record_bucket=int(raw_row['raw_record_bucket']),
                        raw_record_part=int(raw_row['raw_record_part']),
                    ),
                )
            return

        for record in self.raw(force_refresh=force_refresh, **kwargs):
            yield self.mapper(record)


class OntologyDataset:
    """Structured ontology dataset rendered as an ontology artifact such as OBO."""

    kind: DatasetKind = 'ontology'

    def __init__(
        self,
        *,
        download: Download | None,
        mapper: Callable[[dict[str, Any]], OntologyTerm | None],
        raw_parser: Callable[..., Generator[dict[str, Any], None, None]],
        ontology_id: str,
        document: OntologyDocument | None = None,
        remark: str | None = None,
        typedefs: list[OntologyTypedef] | None = None,
        extension: str = 'obo',
        file_stem: str | None = None,
    ) -> None:
        self.download = download
        self.mapper = mapper
        self._raw_parser = raw_parser
        self._last_raw_snapshot: RawSnapshot | None = None
        if not ontology_id or not ontology_id.strip():
            raise ValueError('OntologyDataset requires a non-empty ontology_id')
        self.ontology_id = ontology_id.strip()
        self.document = document or OntologyDocument(
            ontology=self.ontology_id,
            remark=remark,
            typedefs=typedefs,
        )
        self.extension = extension
        self.file_stem = file_stem

    def preparse(
        self,
        *,
        source: str,
        dataset: str,
        raw_records_root: str | Path | None = None,
        force_refresh: bool = False,
        **kwargs: Any,
    ) -> RawSnapshot:
        parser_contract = _parser_contract(self._raw_parser, kwargs)
        output_root = Path(raw_records_root) if raw_records_root else default_raw_records_root()
        if not force_refresh:
            cached_fingerprint = _download_fingerprint_from_cache(self.download, kwargs)
            reusable = reuse_raw_snapshot_if_unchanged(
                source=source,
                dataset=dataset,
                output_root=output_root,
                download_fingerprint=cached_fingerprint,
                parser_contract=parser_contract,
            )
            if reusable is not None:
                self._last_raw_snapshot = reusable
                return reusable
        opener = self.download.open(force_refresh=force_refresh, **kwargs) if self.download else None
        download_fingerprint = _download_fingerprint(self.download, opener, kwargs)
        if not force_refresh:
            reusable = reuse_raw_snapshot_if_unchanged(
                source=source,
                dataset=dataset,
                output_root=output_root,
                download_fingerprint=download_fingerprint,
                parser_contract=parser_contract,
            )
            if reusable is not None:
                self._last_raw_snapshot = reusable
                return reusable
        records = self._raw_parser(opener, force_refresh=force_refresh, **kwargs)
        self._last_raw_snapshot = materialize_raw_records(
            records=records,
            source=source,
            dataset=dataset,
            output_root=output_root,
            download_fingerprint=download_fingerprint,
            parser_contract=parser_contract,
        )
        return self._last_raw_snapshot

    def raw(
        self,
        force_refresh: bool = False,
        *,
        source: str | None = None,
        dataset: str | None = None,
        raw_records_root: str | Path | None = None,
        use_preparse: bool = False,
        changed_only: bool = False,
        raw_snapshot: RawSnapshot | None = None,
        **kwargs: Any,
    ) -> Generator[dict[str, Any], None, None]:
        if use_preparse:
            if not source or not dataset:
                raise ValueError('source and dataset are required when use_preparse=True')
            snapshot = raw_snapshot or self.preparse(
                source=source,
                dataset=dataset,
                raw_records_root=raw_records_root,
                force_refresh=force_refresh,
                **kwargs,
            )
            if raw_snapshot is not None:
                self._last_raw_snapshot = raw_snapshot
            if changed_only:
                yield from iter_changed_raw_record_dicts(snapshot.records_path, snapshot.delta_path)
            else:
                yield from iter_raw_record_dicts(snapshot.records_path)
            return

        opener = self.download.open(force_refresh=force_refresh, **kwargs) if self.download else None
        yield from self._raw_parser(opener, force_refresh=force_refresh, **kwargs)

    def accept_last_preparse(self) -> None:
        if self._last_raw_snapshot is not None:
            accept_raw_snapshot(self._last_raw_snapshot)

    def __call__(self, force_refresh: bool = False, **kwargs: Any) -> Generator[OntologyTerm, None, None]:
        for record in self.raw(force_refresh=force_refresh, **kwargs):
            term = self.mapper(record)
            if term is not None:
                yield term


def ontology_term_to_entity(
    term: OntologyTerm,
    *,
    ontology_id: str,
) -> Entity:
    """Convert a structured ontology term into a first-class CV-term entity."""
    identifiers = [Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=term.id)]

    for alt_id in term.alt_ids or []:
        if alt_id and alt_id != term.id:
            identifiers.append(Identifier(type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=alt_id))

    if term.name:
        identifiers.append(Identifier(type=IdentifierNamespaceCv.NAME, value=term.name))

    for synonym in term.synonyms or []:
        if synonym and synonym != term.name:
            identifiers.append(Identifier(type=IdentifierNamespaceCv.SYNONYM, value=synonym))

    annotations: list[Annotation] = [
        Annotation(term=OntologyAnnotationCv.ONTOLOGY_ID, value=ontology_id),
    ]
    if term.definition:
        annotations.append(Annotation(term=OntologyAnnotationCv.DEFINITION, value=term.definition))
    for comment in term.comments or []:
        if comment:
            annotations.append(Annotation(term=OntologyAnnotationCv.COMMENT, value=comment))
    if term.is_obsolete is not None:
        annotations.append(Annotation(term=OntologyAnnotationCv.IS_OBSOLETE, value=str(bool(term.is_obsolete)).lower()))

    return Entity(
        type=EntityTypeCv.CV_TERM,
        identifiers=identifiers,
        annotations=annotations or None,
    )


class ArtifactDataset:
    """Generic non-parquet artifact dataset."""

    def __init__(
        self,
        *,
        renderer: Callable[..., str],
        download: Download | None = None,
        extension: str,
        file_stem: str | None = None,
        kind: DatasetKind | None = None,
    ) -> None:
        self.renderer = renderer
        self.download = download
        self.extension = extension
        self.file_stem = file_stem
        self.kind = kind

    def render(self, force_refresh: bool = False, **kwargs: Any) -> str:
        opener = self.download.open(force_refresh=force_refresh, **kwargs) if self.download else None
        return self.renderer(opener, force_refresh=force_refresh, **kwargs)


class Resource:
    """Container for datasets and metadata."""

    def __init__(self, config: ResourceConfig, **datasets: Dataset) -> None:
        self.config = config
        for name, dataset in datasets.items():
            setattr(self, name, dataset)
        self._datasets = datasets

    def metadata(self) -> Generator[Entity, None, None]:
        yield self.config.metadata()

    def datasets(self) -> dict[str, Dataset]:
        return self._datasets

    def __call__(self) -> Generator[Entity, None, None]:
        return self.metadata()


def _first_handle(opener) -> Any | None:
    if not opener or not opener.result:
        return None
    if isinstance(opener.result, dict):
        return next(iter(opener.result.values()), None)
    return opener.result


def read_opener_text(opener, **_kwargs: Any) -> str:
    """Read text content from a download opener.

    ``Opener.result`` may be a plain file handle, an archive mapping, or a
    line iterator depending on the downloaded file and cachedir settings.
    """
    handle = _first_handle(opener)
    if handle is None:
        return ''
    if hasattr(handle, 'read'):
        content = handle.read()
        return content.decode('utf-8') if isinstance(content, bytes) else str(content)
    return ''.join(
        chunk.decode('utf-8') if isinstance(chunk, bytes) else str(chunk)
        for chunk in handle
    )


def iter_csv(opener, delimiter: str = ',', **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    handle = _first_handle(opener)
    if not handle:
        return
    yield from csv.DictReader(handle, delimiter=delimiter)


def iter_tsv(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    yield from iter_csv(opener, delimiter='\t')


def iter_json(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    handle = _first_handle(opener)
    if not handle:
        return
    data = json.load(handle)
    if isinstance(data, list):
        yield from data
    else:
        yield data


def iter_jsonl(opener, **_kwargs: Any) -> Generator[dict[str, Any], None, None]:
    handle = _first_handle(opener)
    if not handle:
        return
    for line in handle:
        if line.strip():
            yield json.loads(line)
