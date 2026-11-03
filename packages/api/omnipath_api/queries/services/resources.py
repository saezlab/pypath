"""Resources query service, including domain SQL and result shaping."""

from __future__ import annotations

import logging
import json
from pathlib import Path
from typing import Any
from urllib.parse import quote


import pyarrow as pa
import pyarrow.parquet as pq


logger = logging.getLogger(__name__)


class ResourcesQueries:
    """Resources queries over ParquetQueries' storage and shaping contract."""

    def list_resources(self) -> list[dict[str, Any]]:
        """Return list of indexed resources with table statistics."""
        self._refresh_inventory_if_stale()
        out = []
        for key, info in sorted(self._visible_resources().items()):
            tables = info["tables"]
            try:
                sizes = {name: Path(path).stat().st_size for name, path in tables.items()}
                out.append(
                    {
                        "key": key,
                        "resource": info["resource"],
                        "queryable": info["resource"] not in self._disabled_query_resources(),
                        "version": info["version"],
                        "entities_count": pq.read_metadata(tables["entity"]).num_rows,
                        "relations_count": pq.read_metadata(tables["relation"]).num_rows,
                        "entities_size_kb": round(
                            sum(
                                v
                                for k, v in sizes.items()
                                if k.startswith("entity_") or k == "entity"
                            )
                            / 1024,
                            1,
                        ),
                        "relations_size_kb": round(
                            sum(v for k, v in sizes.items() if k.startswith("relation")) / 1024, 1
                        ),
                    }
                )
            except (OSError, pa.ArrowException):
                logger.exception("Cannot inspect indexed resource %s", key)
                continue
        return out

    def resolve_resource_info(self, resource_id: str) -> dict[str, Any] | None:
        """Resolve a resource slug or `resource/version` key to its latest inventory row."""
        token = str(resource_id or "").strip().strip("/")
        if not token:
            return None
        self._refresh_inventory_if_stale()
        visible = self._visible_resources()
        if token in visible:
            return visible[token]
        latest = self._visible_latest().get(token)
        if latest:
            return latest
        for info in self._visible_latest().values():
            if info["resource"] == token:
                return info
        return None

    def _parquet_file_record(
        self, name: str, path: Path | None, *, include_schema: bool = False
    ) -> dict[str, Any] | None:
        if path is None:
            return None
        file_path = Path(path)
        if not file_path.is_file():
            return None
        try:
            meta = pq.read_metadata(file_path)
            rows = meta.num_rows
        except (OSError, pa.ArrowException):
            logger.exception("Cannot read artifact metadata %s", file_path)
            raise
        record: dict[str, Any] = {
            "name": name,
            "size_bytes": file_path.stat().st_size,
            "rows": rows,
        }
        if self.public_data_url:
            relative_path = file_path.relative_to(self.data_root).as_posix()
            record["url"] = f"{self.public_data_url}/{quote(relative_path, safe='/')}"
        if include_schema:
            try:
                schema = pq.read_schema(file_path)
                record["columns"] = [
                    {"name": field.name, "type": str(field.type)} for field in schema
                ]
            except (OSError, pa.ArrowException):
                logger.exception("Cannot read artifact schema %s", file_path)
                raise
        return record

    def _resource_file_map(self, info: dict[str, Any]) -> dict[str, Path | None]:
        return {f"{name}.parquet": path for name, path in info["tables"].items()}

    def _download_file_map(
        self, info: dict[str, Any], *, include_evidence: bool = False
    ) -> dict[str, Path | None]:
        files = self._resource_file_map(info)
        names = self.RESOURCE_PARQUET_FILES if include_evidence else self.PUBLIC_DOWNLOAD_FILES
        return {name: files[name] for name in names if name in files}

    def list_resource_files(
        self, resource_id: str, *, include_evidence: bool = False
    ) -> dict[str, Any] | None:
        info = self.resolve_resource_info(resource_id)
        if not info:
            return None
        files = [
            record
            for name, path in self._download_file_map(
                info, include_evidence=include_evidence
            ).items()
            if (record := self._parquet_file_record(name, path, include_schema=True))
        ]
        return {
            "resource_id": info["resource"],
            "version": info["version"],
            "files": files,
        }

    def get_resource_file_path(self, resource_id: str, filename: str) -> Path | None:
        if filename not in self.RESOURCE_PARQUET_FILES:
            return None
        info = self.resolve_resource_info(resource_id)
        if not info:
            return None
        path = self._resource_file_map(info).get(filename)
        if path and Path(path).is_file():
            return Path(path)
        return None

    def iter_resource_file_paths(
        self, resource_id: str, *, include_evidence: bool = False
    ) -> list[tuple[str, Path]]:
        info = self.resolve_resource_info(resource_id)
        if not info:
            return []
        out: list[tuple[str, Path]] = []
        for name, path in self._download_file_map(info, include_evidence=include_evidence).items():
            if path and Path(path).is_file():
                out.append((name, Path(path)))
        return out

    def list_resource_catalog(self) -> list[dict[str, Any]]:
        """The explorer catalog of the release in scope.

        It reads every resource's files, build manifest and statistics (about 1 s for 171
        resources), so it is kept until the data, the release, the query policy or the public
        address of the files changes.
        """
        self._refresh_inventory_if_stale()
        key = (
            self._inventory_fingerprint,
            json.dumps(self._release_scope.get(), sort_keys=True, default=str),
            tuple(sorted(self._disabled_query_resources())),
            getattr(self, "public_data_url", None),
        )
        return self._catalog_cache.get(key, self._build_resource_catalog)

    def list_resource_topics(self) -> list[dict[str, Any]]:
        """The topic tags of the catalog, each with the resources that have it.

        What the topic blocks of the explorer need, without the rest of the catalog: the
        topics with most resources first, their resources in name order.
        """
        topics: dict[str, dict[str, Any]] = {}
        for resource in self.list_resource_catalog():
            for tag in resource.get("tags") or []:
                if tag.get("dimension") != "topic":
                    continue
                topic = topics.setdefault(
                    tag["id"],
                    {
                        "id": tag["id"],
                        "label": tag["label"],
                        "description": tag.get("description", ""),
                        "resources": [],
                    },
                )
                topic["resources"].append(resource["resource_id"])
        for topic in topics.values():
            topic["resources"].sort()
        return sorted(topics.values(), key=lambda t: (-len(t["resources"]), t["label"]))

    def _build_resource_catalog(self) -> list[dict[str, Any]]:
        catalog: list[dict[str, Any]] = []
        disabled = self._disabled_query_resources()
        for info in self._selected_resource_infos(include_disabled=True):
            files = [
                record
                for name, path in self._download_file_map(info).items()
                if (record := self._parquet_file_record(name, path))
            ]
            if not files:
                continue
            entity_rows = next(
                (item["rows"] for item in files if item["name"] == "entity.parquet"), 0
            )
            relation_rows = next(
                (item["rows"] for item in files if item["name"] == "relation.parquet"), 0
            )
            manifest_path = Path(info["tables"]["entity"]).parent / "build_manifest.json"
            try:
                manifest = json.loads(manifest_path.read_text())
            except FileNotFoundError:
                manifest = {}
            except (OSError, ValueError):
                logger.exception("Cannot read build manifest %s", manifest_path)
                raise
            from omnipath_core.resource_metadata import resource_metadata

            try:
                resolution = json.loads(
                    (manifest_path.parent / "resolution_stats.json").read_text()
                )
                if not isinstance(resolution, dict) or "input_entities" not in resolution:
                    resolution = None
            except FileNotFoundError:
                resolution = None
            except (OSError, ValueError):
                logger.exception("Cannot read resolution statistics for %s", manifest_path.parent)
                raise
            max_records = manifest.get("max_records")
            sample_build = max_records is not None if "max_records" in manifest else None
            catalog.append(
                {
                    "max_records": max_records,
                    "sample_build": sample_build,
                    "dataset_subset": manifest.get("datasets"),
                    "resolution_stats": resolution,
                    **resource_metadata(info["resource"], manifest),
                    "resource_id": info["resource"],
                    "queryable": info["resource"] not in disabled,
                    "resource_name": info["resource"],
                    "version": info["version"],
                    "entity_count": entity_rows,
                    "interaction_count": relation_rows,
                    "total_size_bytes": sum(item["size_bytes"] for item in files),
                    "files": files,
                }
            )
        return catalog
