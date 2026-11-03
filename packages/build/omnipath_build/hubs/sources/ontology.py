"""Ontology hub: the names of ontology terms, so a term a resource cites only by its id is named.

One record per term (its CURIE), with the term's ``name`` and its ``alt_id``s as further claims.
Terms are not matched against this hub: ontology terms keep their accession as their identity.
The identity build turns it into ``term_labels.parquet``, which the resolver reads to name terms
whose observation carries no name (UniProt cites GO terms, HMDB ChemOnt classes by id alone).
"""

from __future__ import annotations

from collections.abc import Iterable

from ..schema import CHEMICAL_TAXON
from ..stream import iter_url_lines
from ..writer import HubParquetWriter

HUB = "ontology"
# Ontology -> OBO file (the files the pypath inputs of these ontologies read).
SOURCES = {
    "go": "https://purl.obolibrary.org/obo/go.obo",
    "chemont": "http://classyfire.wishartlab.com/system/downloads/1_0/chemont/ChemOnt_2_1.obo.zip",
    "hpo": "https://purl.obolibrary.org/obo/hp.obo",
    "mondo": "https://purl.obolibrary.org/obo/mondo.obo",
    "psi_mi": "https://raw.githubusercontent.com/HUPO-PSI/psi-mi-CV/master/psi-mi.obo",
    "bto": "https://purl.obolibrary.org/obo/bto.obo",
    "apo": "https://purl.obolibrary.org/obo/apo.obo",
    "doid": "https://purl.obolibrary.org/obo/doid.obo",
    "pato": "https://purl.obolibrary.org/obo/pato.obo",
    "eco": "https://purl.obolibrary.org/obo/eco.obo",
    "mp": "https://purl.obolibrary.org/obo/mp.obo",
    "wbphenotype": "https://purl.obolibrary.org/obo/wbphenotype.obo",
    # The basic editions hold only the ontology's own terms, none merged in from its imports.
    "cl": "https://purl.obolibrary.org/obo/cl/cl-basic.obo",
    "uberon": "https://purl.obolibrary.org/obo/uberon/uberon-basic.obo",
}
# Id prefix of the ontologies whose files may carry the terms of others (a merged edition).
OWN_PREFIX = {
    "bto": "BTO:",
    "apo": "APO:",
    "cl": "CL:",
    "uberon": "UBERON:",
    "doid": "DOID:",
    "pato": "PATO:",
    "eco": "ECO:",
    "mp": "MP:",
    "wbphenotype": "WBPhenotype:",
}


def _terms(lines: Iterable[str], prefix: str = ""):
    """``(id, name, alt_ids)`` of every non-obsolete ``[Term]`` stanza whose id starts with ``prefix``.

    An ``alt_id`` that is not an accession is left out: APO writes a curator's
    name there for two terms, and one alias would then name both.
    """
    from pypath.inputs_v2.parsers.obo import parse_obo_text

    for term in parse_obo_text("\n".join(line.rstrip("\r\n") for line in lines)):
        if term.get("id") and term["id"].startswith(prefix) and term.get("name") and not term.get("is_obsolete"):
            yield term["id"], term["name"], [alt for alt in term.get("alt_ids") or [] if ":" in alt]


def emit(writer: HubParquetWriter, sources: dict[str, Iterable[str]] | None = None) -> None:
    for backend, url in SOURCES.items():
        lines = sources[backend] if sources is not None else iter_url_lines(url)
        for term, name, alt_ids in _terms(lines, OWN_PREFIX.get(backend, "")):
            if not writer.add_identity(HUB, term, CHEMICAL_TAXON, backend):
                return
            if not writer.add("name", name, term, CHEMICAL_TAXON, backend):
                return
            for alt_id in alt_ids:
                if not writer.add(HUB, alt_id, term, CHEMICAL_TAXON, backend):
                    return
