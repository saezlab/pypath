"""PubChem CID hub from CID-InChI-Key.gz then CID-SMILES.gz."""

from __future__ import annotations

from ..schema import CHEMICAL_TAXON
from ..stream import iter_url_lines
from ..writer import HubParquetWriter

INCHIKEY_URL = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-InChI-Key.gz"
SMILES_URL = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz"


def emit(writer: HubParquetWriter) -> None:
    for line in iter_url_lines(INCHIKEY_URL):
        if writer.full:
            return
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            continue
        cid, inchi, inchikey = fields[0].strip(), fields[1].strip(), fields[2].strip()
        if not cid:
            continue
        if not writer.add_identity("pubchem", cid, CHEMICAL_TAXON, "pubchem"):
            return
        if inchi and not writer.add("inchi", inchi, cid, CHEMICAL_TAXON, "pubchem"):
            return
        if inchikey and not writer.add("inchikey", inchikey, cid, CHEMICAL_TAXON, "pubchem"):
            return
    for line in iter_url_lines(SMILES_URL):
        if writer.full:
            return
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 2:
            continue
        cid, smiles = fields[0].strip(), fields[1].strip()
        if not cid or not smiles:
            continue
        if not writer.add("smiles", smiles, cid, CHEMICAL_TAXON, "pubchem"):
            return
