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
    OntologyCv,
    ResourceAnnotationCv,
    ResourceCv,
    UpdateCategoryCV,
)
from biolink_model.datamodel.model import OntologyClass, slots
from omnipath_core.naming import Namespace
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
        """Resolve this resource's (slug, short, full, synonyms).

        The lookup slug is the explicit ``slug`` or the canonical ``ResourceCv``
        member name (``self.id``) — not the inconsistent ``name`` — so the
        authoritative ``resources.json`` entry is found reliably.
        """
        from pypath.inputs_v2.resource_names import resolve_names, slugify

        cv_slug = None
        try:
            cv_slug = slugify(self.id.name)
        except Exception:
            cv_slug = None

        return resolve_names(
            name=self.name,
            slug=self.slug or cv_slug,
            short=self.short,
            full=self.full,
            synonyms=self.synonyms,
        )

    def metadata(self) -> Entity:
        annotations = [
            Annotation(
                term=ResourceAnnotationCv.LICENSE, value=str(self.license)
            ),
            Annotation(
                term=ResourceAnnotationCv.UPDATE_CATEGORY,
                value=str(self.update_category),
            ),
        ]
        if self.pubmed:
            annotations.append(
                Annotation(term=IdentifierNamespaceCv.PUBMED, value=self.pubmed)
            )
        annotations.extend(
            [
                Annotation(term=ResourceAnnotationCv.URL, value=self.url),
                Annotation(
                    term=ResourceAnnotationCv.DESCRIPTION,
                    value=self.description,
                ),
            ]
        )

        return Entity(
            type=EntityTypeCv.CV_TERM,
            identifiers=[
                Identifier(
                    type=IdentifierNamespaceCv.CV_TERM_ACCESSION, value=self.id
                ),
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
        filename = _resolve(
            self.filename, force_refresh=force_refresh, **kwargs
        )
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
        yield from self._raw_parser(
            opener, force_refresh=force_refresh, **kwargs
        )

    def __call__(
        self, force_refresh: bool = False, **kwargs: Any
    ) -> Generator[Entity, None, None]:
        for record in self.raw(force_refresh=force_refresh, **kwargs):
            yield self.mapper(record)


def _ontology_identifier_namespace(identifier, default):
    """Keep imported CURIE namespaces; a document can contain foreign terms."""
    prefix, separator, _ = str(identifier).partition(':')
    if not separator:
        return default
    return {
        'GO': Namespace.GO, 'HP': Namespace.HPO, 'MONDO': Namespace.MONDO,
        'CHEMONTID': Namespace.CHEMONT, 'MI': Namespace.MI,
        'EC': Namespace.EC, 'OM': Namespace.OM,
        'UniProtKB-KW': Namespace.UNIPROT_KEYWORD,
    }.get(prefix, prefix.lower())


def ontology_term_to_entity(
    term: OntologyTerm,
    *,
    ontology_id: str,
    identifier_type: Namespace,
    entity_type: Any = OntologyClass,
    relationship_predicates: dict[str, Any] | None = None,
) -> Entity | None:
    """Serialize ontology structure using the resource's explicit semantic choices.

    Only OBO is_a and source-approved relationship mappings become edges.
    Unmapped published CURIE predicates remain attributes; other fields stay raw.
    """
    if not term.id or term.is_obsolete:
        return None
    identifiers = [
        Identifier(type=_ontology_identifier_namespace(value, identifier_type), value=value)
        for value in dict.fromkeys([term.id, *(term.alt_ids or [])])
        if value
    ]
    identifiers.extend(
        Identifier(type=Namespace.NAME, value=value)
        for value in [term.name]
        if value
    )
    identifiers.extend(
        Identifier(type=Namespace.SYNONYM, value=value)
        for value in term.synonyms or []
        if value
    )
    annotations = []
    for slot, values in (
        (slots.description, [term.definition, *(term.comments or [])]),
        (slots.xref, [value.split()[0] for value in term.xrefs or [] if value]),
    ):
        annotations.extend(
            Annotation(term=slot, value=value) for value in values if value
        )
    relations = []
    seen = set()
    for predicate, target in [
        *((slots.subclass_of, parent) for parent in term.is_a or []),
        *(
            (relationship_predicates.get(rel.type), rel.target)
            for rel in term.relationships or []
            if relationship_predicates and rel.type in relationship_predicates
        ),
    ]:
        key = (str(predicate), target)
        if not target or key in seen:
            continue
        seen.add(key)
        relations.append(
            OntologyRelation(
                predicate=predicate,
                object=EntityRef(entity_type, _ontology_identifier_namespace(target, identifier_type), target),
                ontology_id=ontology_id,
            )
        )
    for rel in term.relationships or []:
        if (
            not relationship_predicates
            or rel.type not in relationship_predicates
        ) and ':' in rel.type:
            annotations.append(Annotation(term=rel.type, value=rel.target))
    return Entity(
        type=entity_type,
        identifiers=identifiers,
        annotations=annotations or None,
        ontology_relations=relations or None,
    )


def ontology_entity_mapper(
    term_mapper: Callable[[dict[str, Any]], OntologyTerm | None],
    *,
    ontology_id: str,
    identifier_type: Namespace,
    entity_type: Any = OntologyClass,
    relationship_predicates: dict[str, Any] | None = None,
) -> Callable[[dict[str, Any]], Entity | None]:
    """Bind source-owned ontology modeling choices to a parsed term mapper."""

    def mapper(row: dict[str, Any]) -> Entity | None:
        term = term_mapper(row)
        return (
            None
            if term is None
            else ontology_term_to_entity(
                term,
                ontology_id=ontology_id,
                identifier_type=identifier_type,
                entity_type=entity_type,
                relationship_predicates=relationship_predicates,
            )
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
        opener = (
            self.download.open(force_refresh=force_refresh, **kwargs)
            if self.download
            else None
        )
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
        return (
            content.decode('utf-8')
            if isinstance(content, bytes)
            else str(content)
        )
    return ''.join(
        chunk.decode('utf-8') if isinstance(chunk, bytes) else str(chunk)
        for chunk in handle
    )


def iter_csv(
    opener, delimiter: str = ',', **_kwargs: Any
) -> Generator[dict[str, Any], None, None]:
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
