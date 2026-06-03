"""
Shared building blocks for inputs_v2 datasets.
"""

from __future__ import annotations

from collections.abc import Callable, Generator
import csv
from dataclasses import dataclass
import functools
import json
from typing import Any, Literal, Protocol

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
from pypath.internals.ontology_schema import OntologyTerm
from pypath.internals.silver_schema import (
    Annotation,
    Entity,
    EntityRef,
    Identifier,
    OntologyRelation,
)
from pypath.share.downloads import download_and_open


class _Resolver(Protocol):
    def __call__(self, **kwargs: Any) -> str: ...


def _resolve(value: str | _Resolver, **kwargs: Any) -> str:
    if callable(value):
        return value(**kwargs)
    return value


def _prepared_cache_available(
    raw_parser: Callable[..., Generator[dict[str, Any], None, None]],
    *,
    force_refresh: bool,
    kwargs: dict[str, Any],
) -> bool:
    if force_refresh:
        return False

    cache_available = getattr(raw_parser, 'prepared_cache_available', None)
    parser_kwargs = dict(kwargs)
    raw_parser_func = getattr(raw_parser, 'func', None)
    partial_keywords = getattr(raw_parser, 'keywords', None)
    if raw_parser_func is not None:
        cache_available = cache_available or getattr(
            raw_parser_func,
            'prepared_cache_available',
            None,
        )
        if partial_keywords:
            parser_kwargs = {**partial_keywords, **parser_kwargs}

    if cache_available is None:
        return False

    return bool(
        cache_available(
            force_refresh=force_refresh,
            **parser_kwargs,
        )
    )


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
    # 3-name model (Milestone M; see pypath.inputs_v2.resource_names):
    #   slug  — all-lowercase canonical id/filter key, no `_`, no spaces
    #   short — the resource's own spelling (display name), no `_`, no spaces
    #   full  — the long name, spaces allowed, no `_`
    # Left None here means "derive": slug from `name`, short = `name`, full from
    # the curated audit registry. `_` is reserved for primary_secondary labels.
    slug: str | None = None
    short: str | None = None
    full: str | None = None
    synonyms: tuple[str, ...] = ()

    def names(self) -> 'ResourceNames':
        """Resolve this resource's (slug, short, full, synonyms)."""
        from pypath.inputs_v2.resource_names import resolve_names

        return resolve_names(
            name=self.name,
            slug=self.slug,
            short=self.short,
            full=self.full,
            synonyms=self.synonyms,
        )

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
        url = _resolve(self.url, force_refresh=force_refresh, **kwargs)
        filename = _resolve(self.filename, force_refresh=force_refresh, **kwargs)
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


DatasetKind = Literal['id_translation']


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

    def raw(
        self,
        force_refresh: bool = False,
        **kwargs: Any,
    ) -> Generator[dict[str, Any], None, None]:
        skip_download_open = bool(kwargs.pop('skip_download_open', False))
        skip_download_open = skip_download_open or _prepared_cache_available(
            self._raw_parser,
            force_refresh=force_refresh,
            kwargs=kwargs,
        )
        opener = (
            None
            if skip_download_open
            else self.download.open(force_refresh=force_refresh, **kwargs)
            if self.download
            else None
        )
        yield from self._raw_parser(opener, force_refresh=force_refresh, **kwargs)

    def __call__(self, force_refresh: bool = False, **kwargs: Any) -> Generator[Entity, None, None]:
        for record in self.raw(force_refresh=force_refresh, **kwargs):
            yield self.mapper(record)


def ontology_term_to_entity(
    term: OntologyTerm,
    *,
    ontology_id: str,
    entity_type: EntityTypeCv | str = EntityTypeCv.CV_TERM,
    identifier_type: IdentifierNamespaceCv | str = (
        IdentifierNamespaceCv.CV_TERM_ACCESSION
    ),
) -> Entity:
    """Convert a structured ontology term into a first-class CV-term entity."""
    identifiers = [Identifier(type=identifier_type, value=term.id)]

    for alt_id in term.alt_ids or []:
        if alt_id and alt_id != term.id:
            identifiers.append(Identifier(type=identifier_type, value=alt_id))

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

    ontology_relations: list[OntologyRelation] = []
    seen_relations: set[tuple[str, str]] = set()

    for parent in term.is_a or []:
        if parent and ('is_a', parent) not in seen_relations:
            seen_relations.add(('is_a', parent))
            ontology_relations.append(
                OntologyRelation(
                    predicate='is_a',
                    object=EntityRef(
                        type=entity_type,
                        identifier_type=identifier_type,
                        identifier=parent,
                    ),
                    ontology_id=ontology_id,
                )
            )

    for relationship in term.relationships or []:
        if not relationship.type or not relationship.target:
            continue
        key = (relationship.type, relationship.target)
        if key in seen_relations:
            continue
        seen_relations.add(key)
        ontology_relations.append(
            OntologyRelation(
                predicate=relationship.type,
                object=EntityRef(
                    type=entity_type,
                    identifier_type=identifier_type,
                    identifier=relationship.target,
                ),
                ontology_id=ontology_id,
            )
        )

    return Entity(
        type=entity_type,
        identifiers=identifiers,
        annotations=annotations or None,
        ontology_relations=ontology_relations or None,
    )


def ontology_entity_mapper(
    term_mapper: Callable[[dict[str, Any]], OntologyTerm | None],
    *,
    ontology_id: str,
    entity_type: EntityTypeCv | str = EntityTypeCv.CV_TERM,
    identifier_type: IdentifierNamespaceCv | str = (
        IdentifierNamespaceCv.CV_TERM_ACCESSION
    ),
) -> Callable[[dict[str, Any]], Entity | None]:
    """Wrap an ontology-term mapper so a regular Dataset emits term entities."""

    def mapper(row: dict[str, Any]) -> Entity | None:
        term = term_mapper(row)
        if term is None:
            return None
        return ontology_term_to_entity(
            term,
            ontology_id=ontology_id,
            entity_type=entity_type,
            identifier_type=identifier_type,
        )

    return mapper


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
