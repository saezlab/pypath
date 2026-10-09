"""
Shared building blocks for inputs_v2 datasets.
"""

from __future__ import annotations

from collections.abc import Callable, Generator
import csv
from dataclasses import dataclass
import functools
import json
import re
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
        raw_table: Callable[..., int] | None = None,
    ) -> None:
        self.download = download
        self.mapper = mapper
        self._raw_parser = raw_parser
        self._raw_table = raw_table
        self.kind = kind

    @property
    def has_table(self) -> bool:
        return self._raw_table is not None

    def table(
        self, db, name: str, force_refresh: bool = False, max_records: int | None = None, **kwargs: Any
    ) -> int:
        """Parse into DuckDB table ``name``: ``rid`` (the row's position in
        :meth:`raw`'s output) and one VARCHAR column per field, NULL where a row
        lacks the field. Returns the row count.

        With ``max_records``, the table holds the first rows only, and only a prefix
        of the input is parsed: it grows until enough rows pass the parser's filter.
        """
        kwargs.pop('max_lines', None)
        opener = self.download.open(force_refresh=force_refresh, **kwargs) if self.download else None
        if max_records is None:
            return self._raw_table(db, name, opener, **kwargs)
        # Estimate the filter's pass rate on a probe, then read a quarter more than needed
        # (early lines may pass more often): one parse of the prefix, rarely two.
        lines = 20_000
        count = self._raw_table(db, name, opener, max_lines=lines, **kwargs)
        while count < max_records:
            previous = count
            lines = int(max_records * lines / max(count, 1) * 1.25) + 1000
            count = self._raw_table(db, name, opener, max_lines=lines, **kwargs)
            if count == previous:  # the input has no more lines
                break
        if count > max_records:
            db.execute(f'DELETE FROM {name} WHERE rid >= {int(max_records)}')
            count = int(max_records)
        return count

    def raw(
        self,
        force_refresh: bool = False,
        **kwargs: Any,
    ) -> Generator[dict[str, Any], None, None]:
        """Parsed rows. A dataset with a table parser parses in SQL and streams the
        table's rows; ``row_parser=True`` uses the row parser instead."""
        if self._raw_table is not None and not kwargs.pop('row_parser', False):
            yield from self._table_rows(force_refresh=force_refresh, **kwargs)
            return
        kwargs.pop('row_parser', None)
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

    def _table_rows(
        self, force_refresh: bool = False, batch_rows: int = 8192, **kwargs: Any
    ) -> Generator[dict[str, Any], None, None]:
        """The table parser's rows in ``rid`` order, as the row parser yields them:
        a field the row lacks (NULL) is not a key.

        With ``max_records``, only a prefix of the input is parsed; it grows until
        it holds enough rows (a parser's filter may drop lines).
        """
        import os
        import tempfile

        import duckdb

        max_records = kwargs.pop('max_records', None)
        # Large inputs spill: work on disk next to the downloads, not in /tmp (often RAM).
        parent = (
            os.environ.get('PYPATH_TABLE_TMPDIR')
            or os.environ.get('PYPATH_DOWNLOAD_DATADIR')
            or None
        )
        with tempfile.TemporaryDirectory(prefix='pypath-table-', dir=parent) as directory:
            db = duckdb.connect(os.path.join(directory, 'table.duckdb'))
            try:
                db.execute(f"SET threads={int(os.environ.get('PYPATH_TABLE_THREADS', 4))}")
                db.execute(f"SET memory_limit='{os.environ.get('PYPATH_TABLE_MEMORY', '2GB')}'")
                db.execute(f"SET temp_directory='{directory}/spill'")
                self.table(db, 'raw', force_refresh=force_refresh, max_records=max_records, **kwargs)
                cursor = db.execute('SELECT * EXCLUDE (rid) FROM raw ORDER BY rid')
                for batch in cursor.fetch_record_batch(batch_rows):
                    for row in batch.to_pylist():
                        yield {key: value for key, value in row.items() if value is not None}
            finally:
                db.close()

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


# Relationship types Biolink does not map although a predicate carries their
# meaning: (predicate, exact). Tautomers and enantiomers are symmetric relations.
_RELATIONSHIP_PREDICATES = {
    'RO:0000087': (slots.has_chemical_role, True),  # has role
    'RO:0002211': (slots.regulates, False),  # regulates
    'RO:0002212': (slots.regulates, False),  # negatively regulates
    'RO:0002213': (slots.regulates, False),  # positively regulates
    'RO:0002203': (slots.develops_into, True),  # develops into
    'BFO:0000067': (slots.contains_process, True),  # contains process
    'CL:4030045': (slots.lacks_part, True),  # lacks_part
    'CL:4030046': (slots.lacks_part, False),  # lacks_plasma_membrane_part
    'RO:0018036': (slots.chemically_similar_to, False),  # is tautomer of
    'RO:0018039': (slots.chemically_similar_to, False),  # is enantiomer of
}


@functools.cache
def relationship_predicate(relationship_type: str) -> tuple[Any, bool] | None:
    """Biolink predicate of an ontology relationship CURIE, and whether it is exact.

    Uses Biolink's own exact mappings, else a broader predicate that Biolink
    lists the relationship under (narrow mappings). A directional relationship
    under a symmetric predicate (e.g. ``related_to``) would lose which term is
    the subject, so it has no predicate.
    """
    from omnipath_core.biolink import is_symmetric, predicate, schema

    if relationship_type in _RELATIONSHIP_PREDICATES:
        return _RELATIONSHIP_PREDICATES[relationship_type]
    view = schema()
    for kind in ('exact_mappings', 'narrow_mappings'):
        candidates = []
        for name, slot in view.all_slots().items():
            if relationship_type in (getattr(slot, kind) or ()):
                try:
                    predicate(name)
                except ValueError:  # not a predicate, or deprecated
                    continue
                candidates.append(name)
        if candidates:
            # The most specific of several mapped predicates.
            name = max(sorted(candidates), key=lambda c: len(view.slot_ancestors(c)))
            exact = kind == 'exact_mappings'
            return (predicate(name), exact) if exact or not is_symmetric(name) else None
    return None


_OBO_ESCAPE_RE = re.compile(r'\\(.)')
_IDENTIFIERS_ORG_RE = re.compile(r'https?://identifiers\.org/([^/]+)/(.+)')
# Publication xrefs (ChemOnt cites books by ISBN) and their CURIE prefixes.
_PUBLICATION_PREFIXES = {
    'PMID': 'PMID',
    'PMCID': 'PMC',
    'DOI': 'DOI',
    'ISBN': 'ISBN',
    'ISBN-10': 'ISBN',
    'ISBN-13': 'ISBN',
}


def _obo_unescape(value: str) -> str:
    """OBO escapes (``\\:``, ``\\"``) of an unquoted value."""
    return _OBO_ESCAPE_RE.sub(r'\1', value)


def _obo_xref(value: str) -> tuple[Any, str] | None:
    """Annotation of an OBO ``xref``: a cross-reference, publication or URL.

    Property-style xrefs (``search-url: "…"``, ``id-validation-regexp: "…"``)
    describe how to use a database, not a cross-reference, and are dropped.
    """
    text = _obo_unescape(value).strip()
    prefix, _, local = text.partition(':')
    if prefix.lower() in {'url', 'http', 'https'}:
        url = (text if prefix.lower() != 'url' else local).strip().strip('"')
        return (slots.url, url.split()[0]) if url else None
    if not local or local[0].isspace():
        return None
    # OBO 1.2 writes database and ID separately: ``GO:GO\:0044093``.
    local = local.split()[0].removeprefix(f'{prefix}:')
    if prefix.upper() in _PUBLICATION_PREFIXES:
        return slots.publications, f'{_PUBLICATION_PREFIXES[prefix.upper()]}:{local}'
    return slots.xref, f'{prefix}:{local}'


def _curie(value: str) -> str | None:
    """CURIE of a relationship target: a CURIE, or an identifiers.org URI."""
    match = _IDENTIFIERS_ORG_RE.fullmatch(value)
    if match:
        return f'{match[1].upper()}:{match[2]}'
    return value if re.fullmatch(r'[A-Za-z][\w.-]*:(?!//)\S+', value) else None


def ontology_term_to_entity(
    term: OntologyTerm,
    *,
    ontology_id: str,
    identifier_type: Namespace,
    entity_type: Any = OntologyClass,
    relationship_predicates: dict[str, Any] | None = None,
) -> Entity | None:
    """Serialize ontology structure using the resource's explicit semantic choices.

    OBO is_a, source-approved relationship mappings and CURIE relationships
    with a Biolink predicate (``relationship_predicate``) between two terms
    become edges; an edge under a broader predicate keeps the ontology's own
    as ``original_predicate``. Other CURIE relationships remain attributes.
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
    annotations = [
        Annotation(term=slot, value=value)
        for slot, value in (
            (slots.description, term.definition),
            *(('rdfs:comment', value) for value in term.comments or []),
            *filter(None, map(_obo_xref, filter(None, term.xrefs or []))),
        )
        if value
    ]
    edges = [(slots.subclass_of, parent, None) for parent in term.is_a or []]
    for rel in term.relationships or []:
        if relationship_predicates and rel.type in relationship_predicates:
            edges.append((relationship_predicates[rel.type], rel.target, None))
        elif ':' in rel.type:
            mapped = relationship_predicate(rel.type)
            target = _curie(rel.target)
            if mapped and target == rel.target:
                predicate, exact = mapped
                edges.append((predicate, target, None if exact else rel.type))
            else:
                annotations.append(Annotation(term=rel.type, value=target or rel.target))
    relations = []
    seen = set()
    for predicate, target, original in edges:
        key = (str(predicate), target, original)
        if not target or key in seen:
            continue
        seen.add(key)
        relations.append(
            OntologyRelation(
                predicate=predicate,
                object=EntityRef(entity_type, _ontology_identifier_namespace(target, identifier_type), target),
                ontology_id=ontology_id,
                annotations=[Annotation(term=slots.original_predicate, value=original)]
                if original
                else None,
            )
        )
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
