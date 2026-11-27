"""Write the Parquet subset products of a published release, in dependency order."""

from __future__ import annotations

from pathlib import Path
from typing import Callable

from omnipath_core import subset_layout as layout


def _writers() -> dict[str, Callable]:
    # Imported on use: the writers read DuckDB and Parquet, which a PostgreSQL-only run never needs.
    from omnipath_subsets.exclusion_metabolites import write_exclusion_metabolites

    from omnipath_subsets.metabolite_class import write_metabolite_class
    from omnipath_subsets.metsigdb.parquet_build import write_metsigdb

    return {
        'exclusion_metabolites': write_exclusion_metabolites,
        'metsigdb': write_metsigdb,
        'metabolite_class': write_metabolite_class,
    }


def build_release_subsets(
    data_root: str | Path,
    release: str,
    products: tuple[str, ...] | list[str] | None = None,
    *,
    seed_directory: Path | str | None = None,
) -> dict[str, dict]:
    """Build ``products`` (all by default) of ``release``; returns each product's manifest.

    A product that already exists for the release stops the run before anything
    is written for it; products built earlier in the same run stay.
    """
    available = _writers()
    wanted = list(products) if products else list(layout.PRODUCTS)
    unknown = [name for name in wanted if name not in layout.PRODUCTS]
    if unknown:
        raise layout.SubsetError(f"unknown product {', '.join(unknown)}")
    manifests = {}
    for product in layout.PRODUCTS:
        if product not in wanted:
            continue
        if product not in available:
            raise layout.SubsetError(f"{product} cannot be built as Parquet yet")
        write = available[product]
        manifests[product] = layout.build_product(
            data_root,
            release,
            product,
            lambda directory, context, write=write: write(
                directory, context, seed_directory=seed_directory
            ),
        )
    return manifests
