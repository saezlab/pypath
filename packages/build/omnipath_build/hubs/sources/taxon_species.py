"""Taxon species hub: identity rows from organisms.yaml."""

from __future__ import annotations

from ..dictionaries import load_organisms
from ..writer import HubParquetWriter


def emit(writer: HubParquetWriter) -> None:
    for tax_id in load_organisms():
        if not writer.add("tax_id", str(tax_id), str(tax_id), "0", "organisms"):
            return
