"""Export identifier hub Parquet files by streaming each source onto its own ids."""

from __future__ import annotations

import concurrent.futures
import logging
from collections.abc import Sequence
from pathlib import Path

from omnipath_build.progress import CancelCheck, ProgressCallback, check_cancel, emit

from .dictionaries import write_backend, write_id_type, write_manifest, write_organism
from .sources import HUB_EMITTERS
from .writer import HubParquetWriter

logger = logging.getLogger(__name__)

HUB_NAMES: tuple[str, ...] = tuple(HUB_EMITTERS)
DICTIONARY_NAMES: tuple[str, ...] = ("id_type", "backend", "organism")


def export_hubs(
    output_dir: str | Path,
    *,
    hubs: Sequence[str] | None = None,
    max_records: int | None = 100,
    include_dictionaries: bool = True,
    parallel: int = 1,
    on_progress: ProgressCallback | None = None,
    should_cancel: CancelCheck | None = None,
) -> dict[str, int]:
    """Stream selected hubs into ``output_dir``. Default cap is 100 output rows."""
    target = Path(output_dir)
    target.mkdir(parents=True, exist_ok=True)
    selected = list(HUB_NAMES) if hubs is None else list(hubs)
    unknown = sorted(set(selected) - set(HUB_EMITTERS))
    if unknown:
        raise ValueError(f"Unknown hub(s): {unknown}. Valid: {sorted(HUB_EMITTERS)}")

    counts: dict[str, int] = {}

    def _export_hub(name: str, index: int | None = None) -> tuple[str, int]:
        check_cancel(should_cancel)
        emit(
            on_progress,
            pipeline="resolver",
            stage=f"emit:{name}",
            status="running",
            current=index,
            total=len(selected),
            message=f"Exporting {name}",
        )
        path = target / f"{name}.parquet"
        writer = HubParquetWriter(
            path,
            max_records=max_records,
            on_progress=on_progress,
            progress_stage=f"emit:{name}",
        )
        try:
            HUB_EMITTERS[name](writer)
            check_cancel(should_cancel)
        except BaseException:
            writer.abort()
            emit(
                on_progress,
                pipeline="resolver",
                stage=f"emit:{name}",
                status="failed",
                current=writer.count,
                message=f"Hub {name} failed",
            )
            logger.exception("Hub %s failed", name)
            raise
        cnt = writer.close()
        emit(
            on_progress,
            pipeline="resolver",
            stage=f"emit:{name}",
            status="done",
            current=cnt,
            total=len(selected),
            message=f"Wrote {name} ({cnt:,} rows)",
        )
        logger.info("Wrote %s (%s rows)", path, cnt)
        return name, cnt

    if parallel <= 1 or len(selected) <= 1:
        for index, name in enumerate(selected):
            n, c = _export_hub(name, index)
            counts[n] = c
    else:
        workers = min(int(parallel), len(selected))
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            future_to_name = {executor.submit(_export_hub, name): name for name in selected}
            for future in concurrent.futures.as_completed(future_to_name):
                check_cancel(should_cancel)
                n, c = future.result()
                counts[n] = c

    if include_dictionaries:
        check_cancel(should_cancel)
        emit(
            on_progress,
            pipeline="resolver",
            stage="dictionaries",
            status="running",
            message="Writing id_type / backend / organism dictionaries",
        )
        counts["id_type"] = write_id_type(target / "id_type.parquet")
        counts["backend"] = write_backend(
            target / "backend.parquet",
            [name for name in HUB_NAMES if name != "taxon_species"],
        )
        counts["organism"] = write_organism(target / "organism.parquet")
        emit(
            on_progress,
            pipeline="resolver",
            stage="dictionaries",
            status="done",
            current=counts["id_type"] + counts["backend"] + counts["organism"],
            message="Wrote dictionaries",
        )

    check_cancel(should_cancel)
    emit(
        on_progress,
        pipeline="resolver",
        stage="manifest",
        status="running",
        message="Writing hub manifest",
    )
    write_manifest(target / "manifest.json", counts)
    emit(
        on_progress,
        pipeline="resolver",
        stage="manifest",
        status="done",
        message="Wrote manifest.json",
    )
    return counts
