"""Synthetic identifier hubs + identity library shared by the canonicalization tests.

The hubs mirror the shape of the real ``data/reference/hubs/*.parquet`` files
(``source_type, source_id, hub_id, taxonomy_id, backend``) and cover the
cases the matching rules must get right:

* TP53: reviewed accession, TrEMBL fragment, secondary accession, ENSG with a
  version suffix, bare HGNC number, synonym ``P53``; the same symbol ``TP53``
  is a synonym of mouse Trp53.
* calmodulin ``P0DP23``: one protein linked to three Entrez genes.
* gene 999 with two reviewed accessions (canonical falls back to Entrez).
* a yeast accession outside the taxon scope.
* water (ChEBI + HMDB agreeing on one InChIKey), aspirin (ChEBI + ChEMBL),
  a structure-less ChEBI class, a ChEBI class folded into an HMDB compound
  through HMDB's cross-reference, and a ChEBI mixture with two InChIKeys.
"""

from __future__ import annotations

import fcntl
import os
from pathlib import Path
import tempfile
from functools import lru_cache
import shutil

import pyarrow as pa
import pyarrow.parquet as pq

from omnipath_build.hubs.schema import HUB_SCHEMA
from omnipath_build.identity import (
    build_hub_index,
    build_hub_kv,
    build_identity,
    build_identity_kv,
)
from omnipath_build.identity.common import CHEMICAL, HUBS, PROTEIN

WATER = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
ASPIRIN = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
LINKED = "KKKKKKKKKKKKKK-LLLLLLLLLL-N"
MIX_A = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
MIX_B = "DDDDDDDDDDDDDD-BBBBBBBBBB-C"
# WATER's connectivity with another stereo/isotope layer: a different molecule, not a
# protonation state of WATER (those differ only in the last character).
VARIANT = WATER[:15] + "ZZZZZZZZZZ-N"


def _hub(hubs: Path, name: str, rows, backend: str | None = None) -> None:
    table = pa.Table.from_pylist(
        [
            {
                "source_type": a,
                "source_id": b,
                "hub_id": c,
                "taxonomy_id": d,
                "backend": backend or name,
            }
            for a, b, c, d in rows
        ],
        schema=HUB_SCHEMA,
    )
    pq.write_table(table, hubs / f"{name}.parquet")


def write_hubs(hubs: Path) -> None:
    hubs.mkdir(parents=True, exist_ok=True)
    _hub(
        hubs,
        "uniprot",
        [
            ("uniprot", "P04637", "P04637", "9606"),
            ("uniprot_entry", "P53_HUMAN", "P04637", "0"),
            ("genesymbol", "TP53", "P04637", "0"),
            ("genesymbol-syn", "P53", "P04637", "0"),
            ("entrez", "7157", "P04637", "0"),
            ("ensg", "ENSG00000141510.15", "P04637", "0"),
            ("hgnc", "11998", "P04637", "0"),
            ("uniprot-sec", "Q15086", "P04637", "0"),
            ("uniprot", "A0A0U1RQF1", "A0A0U1RQF1", "9606"),
            ("uniprot_entry", "A0A0U1RQF1_HUMAN", "A0A0U1RQF1", "0"),
            ("entrez", "7157", "A0A0U1RQF1", "0"),
            ("genesymbol", "TP53", "A0A0U1RQF1", "0"),
            ("uniprot", "P0DP23", "P0DP23", "9606"),
            ("uniprot_entry", "CALM1_HUMAN", "P0DP23", "0"),
            ("genesymbol", "CALM1", "P0DP23", "0"),
            ("entrez", "801", "P0DP23", "0"),
            ("entrez", "805", "P0DP23", "0"),
            ("entrez", "808", "P0DP23", "0"),
            ("uniprot", "P02340", "P02340", "10090"),
            ("uniprot_entry", "P53_MOUSE", "P02340", "0"),
            ("genesymbol", "Trp53", "P02340", "0"),
            ("genesymbol-syn", "TP53", "P02340", "0"),
            ("entrez", "22059", "P02340", "0"),
            ("uniprot", "P11111", "P11111", "9606"),
            ("uniprot_entry", "AAA_HUMAN", "P11111", "0"),
            ("entrez", "999", "P11111", "0"),
            ("genesymbol", "AAA", "P11111", "0"),
            ("uniprot", "P22222", "P22222", "9606"),
            ("uniprot_entry", "AAB_HUMAN", "P22222", "0"),
            ("entrez", "999", "P22222", "0"),
            ("genesymbol", "AAA", "P22222", "0"),
            ("uniprot", "Q99999", "Q99999", "4932"),
            ("uniprot_entry", "Q99999_YEAST", "Q99999", "0"),
            ("uniprot", "P12345", "P12345", "11676"),
            ("uniprot_entry", "POL_HV1H2", "P12345", "0"),
            ("genesymbol", "pol", "P12345", "0"),
            # Ensembl/RefSeq protein and transcript aliases of TP53's reviewed accession.
            ("ensp", "ENSP00000269305.4", "P04637", "0"),
            ("enst", "ENST00000269305.8", "P04637", "0"),
            ("refseq_protein", "NP_000537.3", "P04637", "0"),
            # A gene whose only catalogued products are unreviewed (entry-name
            # prefix equals the accession): review status must not expand it.
            ("uniprot", "P33331", "P33331", "9606"),
            ("uniprot_entry", "P33331_HUMAN", "P33331", "0"),
            ("entrez", "4242", "P33331", "0"),
            ("genesymbol", "UNREVGENE", "P33331", "0"),
            ("uniprot", "P33332", "P33332", "9606"),
            ("uniprot_entry", "P33332_HUMAN", "P33332", "0"),
            ("entrez", "4242", "P33332", "0"),
            ("genesymbol", "UNREVGENE", "P33332", "0"),
        ],
    )
    # A gene with many explicitly catalogued products still has one gene identity.
    table = pq.read_table(hubs / "uniprot.parquet")
    extra = [
        {
            "source_type": ns,
            "source_id": value,
            "hub_id": f"P{80000 + i}",
            "taxonomy_id": "9606",
            "backend": "uniprot",
        }
        for i in range(25)
        for ns, value in (
            ("uniprot", f"P{80000 + i}"),
            ("entrez", "777"),
            ("genesymbol", "MANYPRODUCTS"),
            ("ensg", "ENSG00000000777"),
            ("hgnc", "777"),
        )
    ]
    pq.write_table(
        pa.concat_tables([table, pa.Table.from_pylist(extra, schema=HUB_SCHEMA)]),
        hubs / "uniprot.parquet",
    )
    _hub(
        hubs,
        "entrez",
        [
            ("entrez", "7157", "7157", "9606"),
            ("ensg", "ENSG00000141510", "7157", "9606"),
            ("entrez", "777", "777", "9606"),
            ("entrez", "55", "55", "9606"),
            ("ensg", "ENSG00000000055", "55", "9606"),
            ("refseq", "NR_000055.1", "55", "9606"),
            ("entrez", "4242", "4242", "9606"),
            ("entrez", "22059", "22059", "10090"),
            ("genesymbol", "Trp53", "22059", "10090"),
        ],
        backend="gene2ensembl",
    )
    pq.write_table(
        pa.table(
            {
                "ncbi_tax_id": ["9606", "10090"],
                "common_name": ["human", "mouse"],
                "latin_name": ["H", "M"],
            }
        ),
        hubs / "organism.parquet",
    )
    _hub(
        hubs,
        "chebi",
        [
            ("chebi", "CHEBI:15377", "CHEBI:15377", "0"),
            ("inchikey", WATER, "CHEBI:15377", "0"),
            ("name", "water", "CHEBI:15377", "0"),
            ("synonym", "H2O", "CHEBI:15377", "0"),
            ("kegg", "C00001", "CHEBI:15377", "0"),
            ("cas", "7732-18-5", "CHEBI:15377", "0"),
            ("chebi", "CHEBI:15365", "CHEBI:15365", "0"),
            ("inchikey", "InChIKey=" + ASPIRIN, "CHEBI:15365", "0"),
            ("name", "acetylsalicylic acid", "CHEBI:15365", "0"),
            ("chebi", "CHEBI:24431", "CHEBI:24431", "0"),
            ("name", "chemical entity", "CHEBI:24431", "0"),
            ("chebi", "CHEBI:16974", "CHEBI:16974", "0"),
            ("name", "some class without structure", "CHEBI:16974", "0"),
            ("chebi", "CHEBI:1", "CHEBI:1", "0"),
            ("inchikey", MIX_A, "CHEBI:1", "0"),
            ("inchikey", MIX_B, "CHEBI:1", "0"),
        ],
    )
    _hub(
        hubs,
        "hmdb",
        [
            ("hmdb", "HMDB0002111", "HMDB0002111", "0"),
            ("inchikey", WATER, "HMDB0002111", "0"),
            ("name", "Water", "HMDB0002111", "0"),
            ("chebi", "15377", "HMDB0002111", "0"),
            ("pubchem", "962", "HMDB0002111", "0"),
            ("hmdb", "HMDB0002024", "HMDB0002024", "0"),
            ("inchikey", LINKED, "HMDB0002024", "0"),
            ("chebi", "16974", "HMDB0002024", "0"),
            ("name", "Linked compound", "HMDB0002024", "0"),
            ("drugbank", "DB09145", "HMDB0002024", "0"),
        ],
    )
    _hub(
        hubs,
        "chembl",
        [
            ("chembl", "CHEMBL25", "CHEMBL25", "0"),
            ("inchikey", ASPIRIN, "CHEMBL25", "0"),
            ("name", "ASPIRIN", "CHEMBL25", "0"),
        ],
    )

    _hub(
        hubs,
        "bigg",
        [
            ("bigg", "h2o", "h2o", "0"),
            ("bigg", "h2o_c", "h2o", "0"),
            ("chebi", "CHEBI:15377", "h2o", "0"),
        ],
    )
    _hub(
        hubs,
        "metanetx",
        [
            ("metanetx", "MNXM2", "MNXM2", "0"),
            ("chebi", "CHEBI:15377", "MNXM2", "0"),
        ],
    )

    for name in CHEMICAL + PROTEIN:
        if not (hubs / f"{name}.parquet").exists():
            _hub(hubs, name, [])


def _write_external_cid_claim_hubs(hubs: Path) -> None:
    """An external resource claims CID 962 for ASPIRIN while HMDB assigns it to WATER.

    Deliberately conflicting PubChem evidence verifies ownership, not real-world
    chemistry. The pubchem hub is empty in the shared hubs, so this lives in its
    own cached template instead of the base one.
    """
    _hub(hubs, "pubchem", [("pubchem", "962", "962", "0"), ("inchikey", ASPIRIN, "962", "0")])


def _write_exact_structure_hubs(hubs: Path) -> None:
    """Authoritative exact-structure identities for the chemical target tests.

    PubChem CID 1 is WATER, ChEBI:1 a VARIANT of WATER's connectivity and
    CHEMBL1 ASPIRIN, so exact structures stay distinct and conflicting
    structures never project. These replace the shared chebi/chembl/pubchem
    hubs (water and the ChEBI mixture disappear), so they need their own template.
    """
    _hub(hubs, "pubchem", [("pubchem", "1", "1", "0"), ("inchikey", WATER, "1", "0")])
    _hub(
        hubs,
        "chebi",
        [("chebi", "CHEBI:1", "CHEBI:1", "0"), ("inchikey", VARIANT, "CHEBI:1", "0")],
    )
    _hub(
        hubs,
        "chembl",
        [("chembl", "CHEMBL1", "CHEMBL1", "0"), ("inchikey", ASPIRIN, "CHEMBL1", "0")],
    )


# Variant name -> hub mutation applied on top of the shared hubs. Only fixtures
# that contradict the shared hubs need a variant; everything additive belongs in
# write_hubs() so it shares the single base build.
VARIANTS = {
    "base": None,
    "external-cid-claim": _write_external_cid_claim_hubs,
    "exact-structures": _write_exact_structure_hubs,
}


SHARED_KEY_ENV = "OMNIPATH_TEST_SHARED_KEY"


def shared_root() -> Path | None:
    """Directory shared by every pytest-xdist worker of one run, else None.

    The controller (see conftest.py) exports ``OMNIPATH_TEST_SHARED_KEY``;
    workers started without it still carry ``PYTEST_XDIST_TESTRUNUID``. Outside
    xdist there is neither, and templates stay per-process.
    """
    key = os.environ.get(SHARED_KEY_ENV) or os.environ.get("PYTEST_XDIST_TESTRUNUID")
    if not key:
        return None
    return Path(tempfile.gettempdir()) / f"omnipath-test-reference-{key}"


def _build_template(root: Path, variant: str) -> Path:
    write_hubs(root / "hubs")
    if VARIANTS[variant]:
        VARIANTS[variant](root / "hubs")
    hubs, index = root / "hubs", root / "hub-index"
    options = dict(memory="512MB", threads=2, min_free_gib=0)
    for hub in HUBS:
        # The hub index builder needs rows; the empty placeholder hubs are skipped.
        if (hubs / f"{hub}.parquet").exists() and pq.read_metadata(
            hubs / f"{hub}.parquet"
        ).num_rows:
            build_hub_index(hub, hubs, index, goslin_cache=root / "goslin", **options)
            build_hub_kv(hub, index, memory="512MB", threads=2, workers=1, min_free_gib=0)
    snapshot = build_identity(index, root / "identity", **options)
    reference = root / "identity" / snapshot["fingerprint"]
    build_identity_kv(reference, min_free_gib=0)
    # Copies below hard-link these files, and every copy reads the template's hub
    # indexes, so an in-place write by a test would silently corrupt the shared
    # template: make that fail loudly instead.
    for directory in (reference, index):
        for path in directory.rglob("*"):
            if path.is_file() and not path.is_symlink():
                path.chmod(0o444)
    return reference


@lru_cache(maxsize=None)
def _template(variant: str = "base"):
    """Build each variant once per process, or once per xdist run across workers.

    Under xdist the first worker to need a variant builds it into the shared
    directory while holding an exclusive ``flock``; the others block on the lock
    and then reuse the finished build (marked by ``reference.txt``, written last,
    so a build that died half way is discarded and redone).
    """
    shared = shared_root()
    if shared is None:
        directory = tempfile.TemporaryDirectory(prefix=f"parquet-test-reference-{variant}-")
        return directory, _build_template(Path(directory.name), variant)
    shared.mkdir(parents=True, exist_ok=True)
    with open(shared / f"{variant}.lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        root, marker = shared / variant, shared / variant / "reference.txt"
        if not marker.exists():
            shutil.rmtree(root, ignore_errors=True)
            root.mkdir()
            reference = _build_template(root, variant)
            marker.write_text(str(reference))
        return None, Path(marker.read_text())


def _link_or_copy(source, destination, *, follow_symlinks=True):
    try:
        os.link(source, destination, follow_symlinks=follow_symlinks)
    except OSError:  # e.g. a different filesystem
        shutil.copy2(source, destination, follow_symlinks=follow_symlinks)


def build_fixture_library(root: Path, variant: str = "base") -> Path:
    """Independent copies of one production-built identity library per test process.

    Each variant is built once per process (lazily, on first use). Callers get
    their own directory tree whose (read-only) files are hard links to the
    template's, since copying ~4k files / 45 MB per test costs seconds. Tests may
    add, unlink or rename files freely; rewriting a file in place raises
    PermissionError instead of touching the shared template.
    """
    _, reference = _template(variant)
    destination = root / "library"
    shutil.copytree(reference, destination, copy_function=_link_or_copy)
    return destination
