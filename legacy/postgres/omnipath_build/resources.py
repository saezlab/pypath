"""Discover runnable pypath ``inputs_v2`` resources for the build pipeline.

The build pipeline does not keep a hand-written list of datasets.
Instead, it imports the configured inputs package, walks resource modules, and
collects ``Resource``, ``Dataset``, and ``ArtifactDataset`` objects exposed by
pypath. Only entity datasets with raw dataset access are selected for evidence
ingest; id translation and artifact-only datasets are skipped.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass
import importlib
import inspect
import os
from pathlib import Path
import pkgutil
import time


@dataclass(slots=True)
class ResourceFunction:
    """Discovered pypath inputs_v2 resource or dataset callable."""

    source: str
    function_name: str
    qualified_module: str
    call: Callable[[], Iterable] | Callable[[], object]
    resource_id: str
    output_kind: str = 'entity'


class DiscoveryError(RuntimeError):
    """Raised when inputs_v2 resource discovery fails."""


# A source is named after its ``inputs_v2`` module, which need not name the
# release the module actually loads. Where a module carries a resource
# *family* name while loading one specific *release*, the build renames the
# source here, so the slug written to ``data_source`` names the release. This is
# a rename, not an alias: the left-hand name must not survive anywhere.
SOURCE_NAME_OVERRIDES: dict[str, str] = {
    # pypath.inputs_v2.connectomedb declares name='ConnectomeDB2025' and
    # downloads connectomedb.org Current-Release. The bare slug would collide
    # with a later ConnectomeDB2020 onboarding.
    'connectomedb': 'connectomedb2025',
}


def resolve_source_name(module_name: str) -> str:
    """The build's source slug for an ``inputs_v2`` module name."""

    return SOURCE_NAME_OVERRIDES.get(module_name, module_name)


def configure_pypath_download_dir() -> Path:
    """Ensure pypath downloads use a project-local cache directory."""

    configured = os.environ.get('PYPATH_DOWNLOAD_DATADIR')
    if configured:
        data_dir = Path(configured)
    else:
        project_root = Path(__file__).resolve().parents[1]
        data_dir = project_root / 'pypath-data'
        os.environ['PYPATH_DOWNLOAD_DATADIR'] = str(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    if os.environ.get('OMNIPATH_BUILD_PYPATH_PROGRESS', '').lower() not in {
        '1',
        'true',
        'yes',
    }:
        try:
            from pypath.share import settings  # noqa: PLC0415

            settings.setup(progressbars=False)
        except Exception:  # noqa: BLE001
            pass
    return data_dir


def discover_resources(
    database_name: str,
    inputs_package: str = 'pypath.inputs_v2',
    progress: bool = False,
) -> tuple[dict[str, list[ResourceFunction]], dict[str, 'ResourceConfig']]:
    """Return input dataset callables grouped by source module name.

    The second value maps each source slug to its module's ``ResourceConfig``,
    when the module declares a module-level ``config``. It carries the
    ``mints`` declaration onward to ``populate_identifier_authority`` without
    a second package walk.
    """

    del database_name
    started = time.perf_counter()
    scanned_modules = 0
    configure_pypath_download_dir()
    if progress:
        print(
            f'[discover_resources] importing package={inputs_package}',
            flush=True,
        )

    from pypath.inputs_v2.base import (  # noqa: PLC0415
        ArtifactDataset,
        Dataset,
        Resource,
        ResourceConfig,
    )

    try:
        root_module = importlib.import_module(inputs_package)
    except ImportError as exc:
        raise DiscoveryError(
            f'Unable to import inputs package "{inputs_package}": {exc}'
        ) from exc

    package_paths = getattr(root_module, '__path__', None)
    if package_paths is None:
        raise DiscoveryError(
            f'Inputs package "{inputs_package}" is not a namespace package'
        )

    prefix = f'{inputs_package}.'
    discovered: dict[str, list[ResourceFunction]] = {}
    configs: dict[str, ResourceConfig] = {}
    dataset_types = (Dataset, ArtifactDataset)

    for module_info in pkgutil.walk_packages(package_paths, prefix):
        scanned_modules += 1
        module_name = module_info.name
        relative_name = module_name[len(prefix) :]
        if not relative_name or relative_name.split('.')[-1].startswith('_'):
            continue
        source_slug = resolve_source_name(relative_name)

        try:
            module = importlib.import_module(module_name)
        except Exception as exc:  # noqa: BLE001
            if progress:
                print(
                    f'[discover_resources] skipping {module_name}: '
                    f'{exc.__class__.__name__}: {exc}',
                    flush=True,
                )
            continue

        resource_members = [
            (name, obj)
            for name, obj in inspect.getmembers(module)
            if isinstance(obj, Resource)
        ]
        dataset_members = [
            (name, obj)
            for name, obj in inspect.getmembers(module)
            if isinstance(obj, dataset_types)
            and getattr(obj, 'kind', None) != 'id_translation'
        ]

        seen_dataset_names = {name for name, _ in dataset_members}
        for _, resource_obj in resource_members:
            for dataset_name, dataset_obj in resource_obj.datasets().items():
                if (
                    getattr(dataset_obj, 'kind', None) != 'id_translation'
                    and dataset_name not in seen_dataset_names
                ):
                    dataset_members.append((dataset_name, dataset_obj))
                    seen_dataset_names.add(dataset_name)

        if not resource_members and not dataset_members:
            continue

        module_config = getattr(module, 'config', None)
        if isinstance(module_config, ResourceConfig) and source_slug not in configs:
            configs[source_slug] = module_config

        module_functions: list[ResourceFunction] = []
        for _, resource_obj in resource_members:
            module_functions.append(
                ResourceFunction(
                    source=source_slug,
                    function_name='resource',
                    qualified_module=module_name,
                    call=resource_obj,
                    resource_id=source_slug,
                ),
            )
            break

        for dataset_name, dataset_obj in dataset_members:
            output_kind = 'entity'
            if isinstance(dataset_obj, ArtifactDataset):
                output_kind = 'artifact'

            if output_kind == 'entity':

                def dataset_call(
                    dataset_obj=dataset_obj,
                    source_name=source_slug,
                    dataset_name=dataset_name,
                ):
                    return dataset_obj(
                        source=source_name,
                        dataset=dataset_name,
                    )

                dataset_call._raw_dataset = dataset_obj
            else:
                dataset_call = dataset_obj

            module_functions.append(
                ResourceFunction(
                    source=source_slug,
                    function_name=dataset_name,
                    qualified_module=module_name,
                    call=dataset_call,
                    resource_id=source_slug,
                    output_kind=output_kind,
                ),
            )

        if module_functions:
            discovered[source_slug] = module_functions

    if not discovered:
        raise DiscoveryError(
            f'No resource functions found under package "{inputs_package}"'
        )
    if progress:
        function_count = sum(len(functions) for functions in discovered.values())
        print(
            '[discover_resources] done '
            f'modules={scanned_modules} sources={len(discovered)} '
            f'functions={function_count} '
            f'elapsed={time.perf_counter() - started:.1f}s',
            flush=True,
        )
    return discovered, configs


def populate_identifier_authority(
    cur,
    schema: str,
    configs: dict[str, 'ResourceConfig'],
    *,
    progress: bool = False,
) -> int:
    """Write the ``mints`` declarations from ``discover_resources`` into
    ``identifier_authority``: one row per namespace a resource is the
    authority for.

    It looks up against this database's live ``vocab_identifier_type`` and
    ``data_source`` tables, not the resolver's static namespace list. Each
    build assigns its own identifier type IDs, and that narrower,
    resolver-only list does not cover every namespace a resource mints.

    The walk skips a namespace with no matching ``vocab_identifier_type``
    row, or a source with no ``data_source`` row, rather than raising an
    error. The resource walk runs before ingest populates those tables.

    Returns the number of rows written.
    """
    from psycopg2 import sql

    from pypath.internals.cv_terms import cv_term_label_accession

    # spec 011 T058: a mints declaration for one of the resolver's own
    # structure-bearing chemical namespaces (the same set resolver_chemical
    # supplies a lookup for) makes that resource the *structure* authority
    # for it, not just the identifier authority -- what
    # resolver_candidate_authority in duckdb_load.py checks to arbitrate a
    # genuine cross-skeleton disagreement.
    from omnipath_build.duckdb_load import (  # noqa: PLC0415
        RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE,
    )

    structure_bearing_type_names = frozenset(
        RESOLVER_CHEMICAL_SLUG_TO_IDENTIFIER_TYPE.values()
    )

    schema_id = sql.Identifier(schema)
    written = 0

    for source_slug, config in configs.items():
        if not config.mints:
            continue

        cur.execute(
            sql.SQL(
                'SELECT source_id FROM {}.data_source WHERE name = %s'
            ).format(schema_id),
            [source_slug],
        )
        source_row = cur.fetchone()
        if source_row is None:
            if progress:
                print(
                    f'[identifier_authority] no data_source row for '
                    f'{source_slug!r}, skipping',
                    flush=True,
                )
            continue
        source_id = source_row[0]

        for namespace in config.mints:
            type_name = cv_term_label_accession(namespace)
            cur.execute(
                sql.SQL(
                    'SELECT identifier_type_id FROM {}.vocab_identifier_type'
                    ' WHERE name = %s'
                ).format(schema_id),
                [type_name],
            )
            type_row = cur.fetchone()
            if type_row is None:
                if progress:
                    print(
                        f'[identifier_authority] no vocab_identifier_type '
                        f'row for {type_name!r} ({source_slug}), skipping',
                        flush=True,
                    )
                continue

            is_structure_authority = type_name in structure_bearing_type_names
            cur.execute(
                sql.SQL(
                    """
                    INSERT INTO {}.identifier_authority
                        (identifier_type_id, source_id, is_structure_authority)
                    VALUES (%s, %s, %s)
                    ON CONFLICT (identifier_type_id) DO UPDATE
                    SET source_id = EXCLUDED.source_id,
                        is_structure_authority = EXCLUDED.is_structure_authority
                    """
                ).format(schema_id),
                [type_row[0], source_id, is_structure_authority],
            )
            written += 1

    if progress:
        print(
            f'[identifier_authority] wrote {written} row(s) from '
            f'{len(configs)} resource config(s)',
            flush=True,
        )
    return written
