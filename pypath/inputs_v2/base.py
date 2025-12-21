"""
Shared building blocks for inputs_v2 datasets.
"""

from __future__ import annotations

from collections.abc import Callable, Generator
from dataclasses import dataclass
import csv
import json
from typing import Any, Protocol

from pypath.internals.silver_schema import Entity, Identifier, Annotation
from pypath.internals.cv_terms import (
    EntityTypeCv,
    IdentifierNamespaceCv,
    LicenseCV,
    UpdateCategoryCV,
    ResourceAnnotationCv,
    ResourceCv,
)
from pypath.share.downloads import download_and_open


class _Resolver(Protocol):
    def __call__(self, **kwargs: Any) -> str: ...


def _resolve(value: str | _Resolver, **kwargs: Any) -> str:
    if callable(value):
        return value(**kwargs)
    return value


@dataclass(frozen=True)
class ResourceConfig:
    id: ResourceCv
    name: str
    url: str
    license: LicenseCV
    update_category: UpdateCategoryCV
    description: str
    pubmed: str | None = None

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
    encoding: str = 'utf-8'
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
            force=force_refresh,
            **download_kwargs,
        )


class Dataset:
    """A Lego brick: download + raw parsing + mapping to Entities."""

    def __init__(
        self,
        download: Download | None,
        mapper: Callable[[dict[str, Any]], Entity],
        raw_parser: Callable[..., Generator[dict[str, Any], None, None]],
    ) -> None:
        self.download = download
        self.mapper = mapper
        self._raw_parser = raw_parser

    def raw(self, force_refresh: bool = False, **kwargs: Any) -> Generator[dict[str, Any], None, None]:
        opener = self.download.open(force_refresh=force_refresh, **kwargs) if self.download else None
        yield from self._raw_parser(opener, force_refresh=force_refresh, **kwargs)

    def __call__(self, force_refresh: bool = False, **kwargs: Any) -> Generator[Entity, None, None]:
        for record in self.raw(force_refresh=force_refresh, **kwargs):
            yield self.mapper(record)


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
