"""On-demand resolution runtime over an ``omnipath-identity-v2`` directory.

``resolve(queries, votes)`` answers a batch for ``LibraryMatcher``. Nothing is
materialized for the whole universe: postings are derived per batch from the per-hub ``id``
store (the ``by_id`` rows), the small identity decisions (exceptions, extra entities, gene
products) and the hub ``rec`` store (``records`` plus ``by_record`` rows); entity records are
assembled on first use and cached per identity fingerprint in SQLite. Every point lookup reads
LMDB only (``<hub index>/kv`` and ``<identity dir>/kv``, see ``identity_kv``); the Parquet files
are build artifacts and are not opened. See ``docs/identity-layer-spec.md`` sections 3 and 4.
"""

from __future__ import annotations

from collections import defaultdict, namedtuple
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import threading
import time

from .identity_kv import SYMBOL_NAMESPACES, DecisionsKv, HubKv
from .observations import CODES

FORMAT = "omnipath-identity-v2"
# Bump when record assembly or label policy changes: cached records are discarded.
RECORD_POLICY = "identity-record-v5"
CACHE_ENV = "OMNIPATH_IDENTITY_CACHE"

KINDS = {"chemical": 1, "protein": 2, "gene": 3, "reaction": 4}
# Reactions (proposal "Reactions"): Rhea master reactions anchor; Rhea's own cross-references
# beat MetaNetX's for the same id; model-reaction ids only exist in MetaNetX.
REACTION_HUBS = ("rhea", "metanetx_reaction")
RHEA_XREF_NAMESPACES = frozenset({"kegg_reaction", "reactome", "metacyc_reaction", "ecocyc_reaction", "macie"})
MNX_XREF_NAMESPACES = frozenset({"bigg_reaction", "vmh_reaction", "seed_reaction", "sabiork_reaction"})
# Coarse chemical cross-references: generic stereo or protonation states, model-specific ids.
# A source's secondary one must not veto its specific identifiers (see ``_fallbacks``).
COARSE_CHEMICAL_NAMESPACES = frozenset({"kegg", "bigg", "metanetx"})
GENE_NAMESPACES = frozenset({"hgnc", "ensg", "enst", "refseq", "genesymbol", "genesymbol-syn"})
# Tags carried by access rows (spec section 3).
NATIVE, FALLBACK, SECONDARY, SYMBOL_SYNONYM = "native", "fallback", "secondary", "symbol_synonym"

_INCHIKEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_INCHIKEY_ANY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z0-9]$")
_ACCESSION = re.compile(r"^[A-Z][A-Z0-9]{5,9}(?:-[0-9]+)?$")
MAX_LABEL = 80
# Goslin levels, least to most specific (pygoslin ``LipidLevel`` order).
GOSLIN_LEVELS = (
    "category",
    "class",
    "species",
    "molecular_species",
    "sn_position",
    "structure_defined",
    "full_structure",
    "complete_structure",
)
# Hub rows kept out of ``identifiers``: attributes, and names that only feed labels.
NOT_IDENTIFIERS = frozenset(
    {"smiles", "inchi", "formula", "synonym", "lipid_shorthand", "systematic_name", "iupac_name"}
)
CHEMICAL_NAME_HUBS = ("chebi", "hmdb", "chembl", "pubchem")
# The lipid databases. Goslin reads any systematically named carboxylic acid as a
# fatty acyl, so in a general compound database an amino acid (glutamine:
# "FA 5:0;2NH2,5NH2;5oxo") or a synonym such as water's "HOH" yields one too: a
# fatty-acyl shorthand labels a compound only when it is a lipid database's.
LIPID_HUBS = frozenset({"swisslipids", "lipidmaps"})


def _label_shorthand(hub: str, shorthand: str) -> bool:
    """Whether a Goslin shorthand of a member record may label the entity.

    A wax ester without its alcohol ("WE 6:0/0:0") is how Goslin misreads an
    anion named "-oate", in every hub: it is never a label.
    """
    if shorthand.startswith("WE ") and shorthand.endswith("/0:0"):
        return False
    return hub in LIPID_HUBS or not shorthand.startswith("FA ")

CHEMICAL_HUBS = (
    "chebi", "pubchem", "chembl", "hmdb", "lipidmaps", "swisslipids",
    "bigg", "metanetx", "refmet", "ramp", "kegg",
)  # fmt: skip
# Hub ids that are native identity keys under the hub's own namespace (brief 2.3).
CHEMICAL_OWN = frozenset(
    {
        "chebi",
        "pubchem",
        "chembl",
        "hmdb",
        "lipidmaps",
        "swisslipids",
        "bigg",
        "metanetx",
        "refmet",
        "ramp",
    }
)
PRODUCT_NAMESPACES = frozenset(
    {"uniprot", "uniprot-sec", "uniprot_entry", "ensp", "refseq_protein", "genbank"}
)
# Levels whose Goslin name is a structure: it keys an InChIKey only when unambiguous (rule 5).
STRUCTURE_LEVELS = frozenset({"full_structure", "complete_structure"})
GENE_ALIAS_TYPES = frozenset({"genesymbol", "genesymbol-syn", "hgnc", "ensg"})
_ISOFORM = re.compile(r"^uniprot:([A-Z0-9]+)-[0-9]+$")
_RNA_ACCESSION = re.compile(r"^(?:NM|NR|XM|XR)_[0-9]+$")
CLAIM_TAGS = frozenset({"claim", "version_stripped"})

_Meta = namedtuple("_Meta", "kind anchor taxon quarantined reviewed gene_ids exists")


def _prefixed_anchor(hub: str, anchor: str) -> str:
    if ":" in anchor:
        return anchor
    return ("uniprot:" if hub == "uniprot" else "inchikey:") + anchor


def _symbol_tag(key_ns, source_ns, tag):
    """Brief 2.3b symbol expansions: exact symbols answer both namespaces, synonyms are weaker."""
    if source_ns == "genesymbol" and tag == "claim":
        return "regular"
    if tag in ("claim", SYMBOL_SYNONYM) and source_ns in SYMBOL_NAMESPACES:
        return SYMBOL_SYNONYM
    return None


def _admit(cls, ns, hub, source_ns, tag, identifier, anchor, anchor_count):
    """Tag a by_id row takes as a lookup row for key namespace ``ns``, or None (brief 2.2-2.3)."""
    if cls == "chem":
        if ns == "inchikey":  # only an anchored structure (anchor_count = 1) is an InChIKey key
            ok = anchor_count == 1 and anchor and anchor.removeprefix("inchikey:") == identifier
            return NATIVE if ok else None
        if ns == "goslin":
            return NATIVE
        if ns in CHEMICAL_OWN:
            if tag == "native":
                return NATIVE if hub == ns else None
            if tag == "claim" and hub == ns and ns in ("chebi", "bigg"):
                return FALLBACK if ns == "bigg" else "regular"
            return None
        if ns == "kegg":
            return FALLBACK
        return "regular" if tag in ("claim", "native") else None  # cas, drugbank
    if cls == "ramp_gene":
        return "regular" if tag == "native" else None
    if cls == "reaction":
        if ns in REACTION_HUBS:  # a hub's own ids (Rhea: master and directional ids)
            return NATIVE if hub == ns and source_ns == ns and tag in ("native", "claim") else None
        if tag in CLAIM_TAGS:
            return NATIVE if hub == "rhea" else FALLBACK
        return None
    if cls == "gene":
        if ns in SYMBOL_NAMESPACES:
            return _symbol_tag(ns, source_ns, tag)
        if ns == "refseq":
            ok = source_ns == "refseq" and _RNA_ACCESSION.fullmatch(identifier)
            return "regular" if ok and tag in CLAIM_TAGS else None
        return "regular" if source_ns == ns and tag in CLAIM_TAGS else None
    if cls == "product":
        if ns in ("uniprot", "uniprot-sec"):
            if source_ns == "uniprot" and tag == "native":
                return NATIVE if ns == "uniprot" else None
            if (source_ns == "uniprot-sec" and tag in ("claim", SECONDARY)) or (
                source_ns == "uniprot" and tag == SECONDARY
            ):
                return SECONDARY
            return None
        return "regular" if source_ns == ns and tag in CLAIM_TAGS else None
    return None


def default_cache_dir(path) -> Path:
    return Path(os.environ.get(CACHE_ENV) or Path(path).resolve().parent / "cache")


def open_runtime(path, *, memory_limit=None, cache_dir=None):
    """The runtime for a pinned identity directory; any other directory is rejected."""
    return IdentityRuntime(path, cache_dir=cache_dir, memory_limit=memory_limit)


@lru_cache(maxsize=1 << 20)
def entity_num(entity_id: str) -> int:
    """First 8 bytes (big-endian) of sha256(entity_id) as an unsigned int, shifted right by one."""
    return int.from_bytes(hashlib.sha256(entity_id.encode()).digest()[:8], "big") >> 1


def decode_key(key: bytes, codes: dict[int, str]):
    """Inverse of ``observations.key``: ``(target, route, ns, scope|None, identifier)``."""
    if (
        len(key) < 7
        or key[0] != 1
        or key[1] not in (1, 2, 3)
        or key[2] not in (1, 2)
        or key[5] not in (0, 1)
    ):
        raise ValueError("Unsupported normalized lookup key")
    ns = codes.get(int.from_bytes(key[3:5], "big"))
    if ns is None:
        raise ValueError("Unsupported identifier namespace")
    offset = 10 if key[5] else 6
    identifier = key[offset:].decode()
    if not identifier:
        raise ValueError("Empty identifier")
    scope = str(int.from_bytes(key[6:10], "big")) if key[5] else None
    return key[1], key[2], ns, scope, identifier


def clean_anchor(kind: int, raw) -> str | None:
    """Anchor in the form the decision kernel accepts, else None.

    Chemical anchors are ``inchikey:<key>``, protein anchors ``uniprot:<acc>``; genes and
    lipid-name entities (``goslin:...``) are unanchored for the kernel. Bare keys and
    accessions are accepted and prefixed.
    """
    if not raw:
        return None
    raw = str(raw)
    if kind == 1:
        bare = raw.removeprefix("inchikey:")
        return "inchikey:" + bare if _INCHIKEY.fullmatch(bare) else None
    if kind == 2:
        bare = raw.removeprefix("uniprot:")
        return "uniprot:" + bare if _ACCESSION.fullmatch(bare) else None
    return None


def _taxon(value) -> str | None:
    text = str(value).strip() if value is not None else ""
    return None if text in ("", "0", "None") else text


def _short_first(values):
    """Deterministic choice: shortest, then alphabetical."""
    return min(values, key=lambda v: (len(v), v), default=None)


def _usable(values, limit):
    out = {str(v).strip() for v in values if v is not None}
    return {v for v in out if v and len(v) <= limit}


def choose_label(kind: int, entity_id: str, rows: list[tuple[str, str, str]]) -> str:
    """Label policy (spec rule 13).

    ``rows`` are ``(hub, source_type, value)`` for every member record. First available:
    gene: NCBI Gene symbol or gene_info name; protein: primary gene name > entry name > accession;
    chemical: lipid Goslin shorthand (a fatty-acyl one only from a lipid database's record,
    a wax ester without its alcohol never), else ChEBI > HMDB > ChEMBL > PubChem > other hub name > systematic name
    (InChIKey-shaped and over-long values skipped); else the id local part.
    """
    local = entity_id.split(":", 1)[1] if ":" in entity_id else entity_id
    by_type = defaultdict(list)
    for hub, source_type, value in rows:
        by_type[source_type].append((hub, value))

    def values(source_type, hubs=None):
        return [v for h, v in by_type.get(source_type, ()) if hubs is None or h in hubs]

    if kind == 3:
        # Genes without a UniProt product (tRNAs, many non-coding genes) carry the NCBI
        # symbol only as a gene_info name row.
        symbol = _short_first(_usable(values("genesymbol", {"entrez"}), 120))
        symbol = symbol or _short_first(_usable(values("name", {"entrez"}), 120))
        return symbol or local
    if kind == 2:
        symbol = _short_first(_usable(values("genesymbol", {"uniprot"}), 120)) or _short_first(
            _usable(values("genesymbol"), 120)
        )
        if symbol:
            return symbol
        entry = _short_first(_usable(values("uniprot_entry", {"uniprot"}), 120))
        return entry or local
    if kind == 4:
        names = by_type.get("name", ())
        pick = _short_first(_usable([v for h, v in names if h == "rhea"], 300)) or _short_first(
            _usable([v for _, v in names], 300)
        )
        return pick or local
    # chemical
    if entity_id.startswith("goslin:"):
        # goslin:<level>:<name> - the name is the most specific level it was anchored at
        parts = entity_id.split(":", 2)
        if len(parts) == 3 and parts[2]:
            return parts[2]
    shorthand = _goslin_shorthand(
        v for h, v in by_type.get("goslin", ()) if _label_shorthand(h, str(v).partition(":")[2])
    )
    if shorthand:
        return shorthand
    names = [(h, v) for h, v in by_type.get("name", ())]
    for hub in CHEMICAL_NAME_HUBS:
        pick = _short_first(_chemical_names(v for h, v in names if h == hub))
        if pick:
            return pick
    pick = _short_first(_chemical_names(v for h, v in names if h not in CHEMICAL_NAME_HUBS))
    if pick:
        return pick
    pick = _short_first(_chemical_names(values("systematic_name") + values("iupac_name")))
    return pick or local


def _chemical_names(values):
    return {v for v in _usable(values, MAX_LABEL) if not _INCHIKEY_ANY.fullmatch(v)}


def _goslin_shorthand(values):
    """Most specific ``<level>:<name>`` value, stripped of its level."""
    best = {}
    for value in values:
        level, _, name = str(value).partition(":")
        if level in GOSLIN_LEVELS and name.strip():
            best.setdefault(GOSLIN_LEVELS.index(level), set()).add(name.strip())
    if not best:
        return None
    return _short_first(best[max(best)])


class IdentityRuntime:
    """Resolution runtime over an identity directory and its hub indexes."""

    # Reference libraries this runtime serves.
    libraries = ("gene_protein", "chemical", "reaction")

    def __init__(self, path, *, cache_dir=None, memory_limit=None, threads=None):
        self.path = Path(path)
        manifest = json.loads((self.path / "manifest.json").read_text())
        if manifest.get("format") != FORMAT:
            raise ValueError("An omnipath-identity-v2 identity directory is required")
        if not manifest.get("fingerprint"):
            raise ValueError("Identity manifest has no fingerprint")
        if not manifest.get("hub_indexes"):
            raise ValueError("Identity manifest names no hub indexes")
        if "namespace_codes" in manifest and manifest["namespace_codes"] != CODES:
            raise ValueError("Namespace encoding differs from the identity manifest")
        self.manifest = manifest
        self.fingerprint = manifest["fingerprint"]
        self.hub_dirs = {}
        for hub, directory in manifest["hub_indexes"].items():
            directory = Path(directory)
            self.hub_dirs[hub] = (
                directory if directory.is_absolute() else self.path / directory
            ).resolve()
        self.chemical_hubs = [h for h in CHEMICAL_HUBS if h in self.hub_dirs]
        self.codes = {code: ns for ns, code in CODES.items()}
        # ``memory_limit`` and ``threads`` are accepted for compatibility: LMDB needs neither.
        self._lipids = None
        self._term_labels = None
        # (hub, (source ns, identifier), taxon) -> the by_id rows of that taxon, for hub kv
        # stores not keyed by taxon (``_rows_filtered``).
        self._scoped_rows = {}
        self.decisions = DecisionsKv(self.path, self.fingerprint)
        try:
            self.hub_kv = {hub: HubKv(hub, directory) for hub, directory in self.hub_dirs.items()}
        except BaseException:
            self.decisions.close()
            raise
        self.cache_dir = Path(cache_dir) if cache_dir is not None else default_cache_dir(self.path)
        self._cache = self._open_cache()

    # ------------------------------------------------------------------ cache

    def _open_cache(self):
        directory = self.cache_dir / str(self.fingerprint)
        try:
            directory.mkdir(parents=True, exist_ok=True)
            target = str(directory / "entities.sqlite")
        except OSError:
            target = ":memory:"  # read-only location: still correct, just not shared
        con = sqlite3.connect(target, timeout=120, check_same_thread=False, isolation_level=None)
        if target != ":memory:":
            # Switching a new cache to WAL needs an exclusive lock, and SQLite answers
            # "locked" at once, without waiting, when several workers try it together.
            deadline = time.monotonic() + 120
            while True:
                try:
                    con.execute("PRAGMA journal_mode=WAL")
                    break
                except sqlite3.OperationalError as exc:
                    if "locked" not in str(exc) or time.monotonic() > deadline:
                        raise
                    time.sleep(0.05)
        con.execute("PRAGMA synchronous=NORMAL")
        con.execute("PRAGMA busy_timeout=120000")
        con.execute("BEGIN IMMEDIATE")
        con.execute("CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY, value TEXT)")
        con.execute(
            "CREATE TABLE IF NOT EXISTS records(entity_id TEXT PRIMARY KEY, record TEXT NOT NULL)"
            " WITHOUT ROWID"
        )
        row = con.execute("SELECT value FROM meta WHERE key='record_policy'").fetchone()
        if row is not None and row[0] != RECORD_POLICY:
            con.execute("DELETE FROM records")
        con.execute("INSERT OR REPLACE INTO meta VALUES('record_policy', ?)", (RECORD_POLICY,))
        con.execute("COMMIT")
        return con

    def _cached(self, entity_ids):
        found = {}
        ids = list(entity_ids)
        for offset in range(0, len(ids), 900):
            chunk = ids[offset : offset + 900]
            marks = ",".join("?" * len(chunk))
            for eid, record in self._cache.execute(
                f"SELECT entity_id, record FROM records WHERE entity_id IN ({marks})", chunk
            ):
                found[eid] = json.loads(record)
        return found

    def _store(self, records):
        if not records:
            return
        self._cache.execute("BEGIN IMMEDIATE")
        try:
            self._cache.executemany(
                "INSERT OR REPLACE INTO records VALUES(?, ?)",
                [(eid, json.dumps(rec, separators=(",", ":"))) for eid, rec in records.items()],
            )
        except BaseException:
            self._cache.execute("ROLLBACK")
            raise
        self._cache.execute("COMMIT")

    # ------------------------------------------------- small indexed lookups

    def _exceptions(self, record_ids):
        """record_id -> (entity_id, decision, quarantined) for records with an exception."""
        ids = set(record_ids)
        if not ids:
            return {}
        return {r: (e, d, bool(q)) for r, (e, d, q) in self.decisions.get("exc", ids).items()}

    def _record_candidates(self, record_ids):
        """record_id -> [entity_id]: the anchors a quarantined or ambiguous record points to."""
        ids = set(record_ids)
        if not ids:
            return {}
        return {r: list(e) for r, e in self.decisions.get("cand", ids).items()}

    def _exception_members(self, entity_ids):
        ids = set(entity_ids)
        found = defaultdict(list)
        if ids:
            for entity_id, record_ids in self.decisions.get("exc_members", ids).items():
                found[entity_id].extend(record_ids)
        return found

    def _extra(self, entity_ids):
        """entity_id -> (kind, taxon, quarantined) for entities not derivable from a record."""
        ids = set(entity_ids)
        if not ids:
            return {}
        return {
            e: (k, _taxon(t), bool(q))
            for e, (k, t, q, _preferred) in self.decisions.get("extra", ids).items()
        }

    def _records(self, hub, local_ids):
        """local_id -> (taxon, anchor, anchor_count, reviewed) from a hub's records."""
        ids = set(local_ids)
        if not ids or hub not in self.hub_kv:
            return {}
        return {
            lid: (_taxon(t), a, int(n or 0), bool(rv))
            for lid, (t, a, n, rv) in self.hub_kv[hub].heads(ids).items()
        }

    def _gene_products(self, name, ids):
        ids = set(ids)
        if not ids:
            return {}
        return {
            str(key): [(str(other), _taxon(taxon)) for other, taxon in links]
            for key, links in self.decisions.get(name, ids).items()
        }

    def _by_protein(self, protein_ids):
        """protein entity id -> [(entrez_id, taxon)]"""
        return self._gene_products("gp_protein", protein_ids)

    def _by_gene(self, entrez_ids):
        """entrez id -> [(protein entity id, taxon)]"""
        return self._gene_products("gp_gene", entrez_ids)

    def _lipid_structures(self):
        """Goslin full-structure name -> inchikey entity id (rule 5; unique names only)."""
        if self._lipids is None:
            self._lipids = {}
            for name, key in self.decisions.items("lipid"):
                if name and key:
                    key = str(key)
                    self._lipids[name] = key if key.startswith("inchikey:") else "inchikey:" + key
        return self._lipids

    # ---------------------------------------------------------------- lookup

    def _plan(self, target, route, ns):
        """(class, source namespaces, hubs) a key reads, or None when it can never match."""
        if target == 3:
            if route != 1:
                return None
            if ns in REACTION_HUBS:
                return "reaction", [ns], [ns]
            if ns in RHEA_XREF_NAMESPACES:
                return "reaction", [ns], [h for h in REACTION_HUBS if h in self.hub_dirs]
            if ns in MNX_XREF_NAMESPACES:
                return "reaction", [ns], ["metanetx_reaction"]
            return None
        if target == 1:
            if route != 1:
                return None
            if ns in CHEMICAL_OWN:
                return "chem", [ns], [ns]
            if ns in ("cas", "drugbank", "kegg", "inchikey", "goslin"):
                return "chem", [ns], self.chemical_hubs
            return None
        if route == 2:
            if ns == "entrez":
                return "entrez", [], []
            if ns == "ramp_gene":
                return "ramp_gene", ["ramp_gene"], ["ramp_gene"]
            return None
        if ns in GENE_NAMESPACES:
            sources = ["genesymbol", "genesymbol-syn"] if ns in SYMBOL_NAMESPACES else [ns]
            return "gene", sources, ["uniprot", "entrez"]
        if ns in PRODUCT_NAMESPACES:
            sources = ["uniprot", "uniprot-sec"] if ns in ("uniprot", "uniprot-sec") else [ns]
            return "product", sources, ["uniprot"]
        return None

    # Pairs read per round when rows are filtered by the taxon of their record (format v2).
    _SCOPED_PAIRS = 256
    # Entries of the memo of taxon-filtered rows (format v2); it is emptied when it is full. A
    # symbol has four entries (exact and synonym, in the gene and in the protein hub), most
    # of them empty: 85,000 keys of real builds made 157,000 entries.
    _SCOPED_MEMO = 500_000

    def _by_id(self, hubs, pairs, scopes=None):
        """by_id rows for ``(source ns, identifier)`` pairs: (pair index, hub, row fields...).

        ``scopes``, where given, has for each pair the taxon its rows are wanted for, or
        None. A row whose record has another taxon is left out here; a record without a
        taxon is kept. A gene symbol is the same word in hundreds of species (``TP53``:
        about 350 gene and 400 protein records, of which 115 are human), and the scope
        filter of ``lookup_many`` would drop those rows only after their entities had been
        looked up: 5,000 symbols of one batch were millions of records and several gigabytes.
        """
        index = defaultdict(list)
        for i, pair in enumerate(pairs):
            index[pair].append(i)
        out = []
        for hub in hubs:
            kv = self.hub_kv.get(hub)
            if kv is None:
                continue
            plain, scoped = [], []
            for pair, owners in index.items():
                if scopes and None not in {scopes[i] for i in owners}:
                    scoped.append(pair)
                else:
                    plain.append(pair)
            if plain:
                for pair, rows in self._rows(hub, kv, plain).items():
                    for local_id, tag, anchor, count, _ in rows:
                        out.extend((i, hub, local_id, tag, anchor, count) for i in index[pair])
            if scoped and kv.by_taxon:
                out.extend(self._rows_of_taxon(hub, kv, scoped, index, scopes))
            elif scoped:
                out.extend(self._rows_filtered(hub, kv, scoped, index, scopes))
        return out

    @staticmethod
    def _rows_of_taxon(hub, kv, scoped, index, scopes):
        """Rows of symbol pairs from a kv keyed by taxon: only their taxon's records are read."""
        found = kv.taxon_ids((*pair, scopes[i]) for pair in scoped for i in index[pair])
        heads = kv.heads({local_id for rows in found.values() for local_id, _ in rows})
        out = []
        for pair in scoped:
            for i in index[pair]:
                for local_id, tag in found.get((*pair, scopes[i]), ()):
                    head = heads.get(local_id)
                    if head is None:
                        raise ValueError(f"Hub {hub} by_id row without a record: {local_id}")
                    out.append((i, hub, local_id, tag, head[1], int(head[2] or 0)))
        return out

    def _rows_filtered(self, hub, kv, scoped, index, scopes):
        """Rows of symbol pairs from a kv not keyed by taxon (format v2).

        The rows of the symbol hold every species, and the taxon of a row is in its record,
        so every record of the symbol is read once to learn it. The rows that remain for a
        taxon are kept (``_scoped_rows``): a resource names the same genes in batch after
        batch, and reading a thousand records of other species from disk for each of them
        again was most of the time of a build. The pairs are read a few at a time so that
        the rows of all of them are never held together.
        """
        memo, out, unread = self._scoped_rows, [], []
        for pair in scoped:
            wanted = {scopes[i] for i in index[pair]}
            if any((hub, pair, scope) not in memo for scope in wanted):
                unread.append(pair)
                continue
            for i in index[pair]:
                out.extend(
                    (i, hub, local_id, tag, anchor, count)
                    for local_id, tag, anchor, count in memo[hub, pair, scopes[i]]
                )
        for start in range(0, len(unread), self._SCOPED_PAIRS):
            chunk = unread[start : start + self._SCOPED_PAIRS]
            found = self._rows(hub, kv, chunk)
            if len(memo) + len(chunk) > self._SCOPED_MEMO:
                memo.clear()
            for pair in chunk:
                rows = found.get(pair, ())
                for scope in {scopes[i] for i in index[pair]}:
                    memo[hub, pair, scope] = tuple(
                        row[:4] for row in rows if row[4] in (None, scope)
                    )
                for i in index[pair]:
                    out.extend(
                        (i, hub, local_id, tag, anchor, count)
                        for local_id, tag, anchor, count, taxon in rows
                        if taxon in (None, scopes[i])
                    )
        return out

    @staticmethod
    def _rows(hub, kv, pairs):
        """pair -> [(local_id, tag, anchor, anchor count, taxon)] of the pairs a hub has."""
        found = kv.ids(pairs)
        if not found:
            return {}
        heads = kv.heads({local_id for rows in found.values() for local_id, _ in rows})
        out = {}
        for pair, rows in found.items():
            out[pair] = full = []
            for local_id, tag in rows:
                head = heads.get(local_id)
                if head is None:
                    raise ValueError(f"Hub {hub} by_id row without a record: {local_id}")
                full.append((local_id, tag, head[1], int(head[2] or 0), _taxon(head[0])))
        return out

    def _entity(self, hub, local_id, anchor, anchor_count, exceptions):
        """Entity of a record: its exception, else its anchor, else entrez id, else record id."""
        record_id = f"{hub}:{local_id}"
        if record_id in exceptions:
            entity_id, _, quarantined = exceptions[record_id]
            return entity_id, quarantined
        if anchor_count == 1 and anchor:
            return _prefixed_anchor(hub, anchor), False
        if hub == "entrez":
            return "entrez:" + local_id, False
        return record_id, False

    def lookup_many(self, keys):
        decoded = {key: decode_key(key, self.codes) for key in dict.fromkeys(keys)}
        plans = {}
        for key, (target, route, ns, scope, identifier) in decoded.items():
            # Symbols exist only taxon-scoped: an unscoped symbol key always misses.
            if scope is None and ns in SYMBOL_NAMESPACES:
                continue
            plan = self._plan(target, route, ns)
            if plan:
                plans[key] = plan

        # 1. by_id rows per key, admitted according to brief 2.2-2.3 (hits: key -> rows)
        hits = defaultdict(list)
        groups = defaultdict(list)
        for key, (cls, sources, hubs) in plans.items():
            if sources:
                groups[tuple(hubs)].append(key)
        for hubs, group in groups.items():
            pairs = [(src, decoded[key][4]) for key in group for src in plans[key][1]]
            owners = [key for key in group for _ in plans[key][1]]
            # A symbol key is taxon-scoped (an unscoped one was left out above): only the
            # records of its taxon can pass the scope filter of step 6.
            scopes = [
                decoded[key][3] if decoded[key][2] in SYMBOL_NAMESPACES else None
                for key in owners
            ]
            if not any(scopes):
                scopes = None
            for i, hub, local_id, tag, anchor, count in self._by_id(hubs, pairs, scopes):
                key, (src, identifier) = owners[i], pairs[i]
                cls, _, _ = plans[key]
                out = _admit(cls, decoded[key][2], hub, src, tag, identifier, anchor, count)
                if out:
                    hits[key].append((hub, local_id, out, anchor, count, identifier))

        # 2. records -> entities (exceptions, anchors, entrez ids, record ids)
        record_ids = {f"{hub}:{local_id}" for rows in hits.values() for hub, local_id, *_ in rows}
        exceptions = self._exceptions(record_ids)
        pointing = self._record_candidates(record_ids)
        candidates = defaultdict(list)  # key -> [(entity_id, tag, source hub)]
        quarantined = set()
        for key, rows in hits.items():
            for hub, local_id, tag, anchor, count, identifier in rows:
                record_id = f"{hub}:{local_id}"
                if record_id in pointing:
                    # A quarantined or ambiguous record votes for each anchor it points to,
                    # instead of for its own exception entity.
                    entities = [(entity_id, False) for entity_id in pointing[record_id]]
                else:
                    entities = [self._entity(hub, local_id, anchor, count, exceptions)]
                for entity_id, flagged in entities:
                    if flagged:
                        quarantined.add(entity_id)
                    if decoded[key][2] == "inchikey" and entity_id != "inchikey:" + identifier:
                        continue  # an overridden record does not create its structure's entity
                    candidates[key].append((entity_id, tag, hub))
        for key, (cls, _, _) in plans.items():
            if cls == "entrez":
                candidates[key].append(("entrez:" + decoded[key][4], NATIVE, "entrez"))

        # 3. rule 5: a structure is a goslin key only if the name maps to exactly that InChIKey
        for key, (cls, _, _) in plans.items():
            if cls == "chem" and decoded[key][2] == "goslin":
                identifier = decoded[key][4]
                if identifier.partition(":")[0] in STRUCTURE_LEVELS:
                    unique = self._lipid_structures().get(identifier)
                    candidates[key] = [
                        c
                        for c in candidates[key]
                        if not c[0].startswith("inchikey:") or c[0] == unique
                    ]

        # 4. isoform -> parent projection (products and gene-level postings)
        projected = {}
        isoforms = {
            c[0]
            for key, cs in candidates.items()
            if plans[key][0] in ("product", "gene")
            for c in cs
            if _ISOFORM.fullmatch(c[0])
        }
        if isoforms:
            parents = self._records("uniprot", {_ISOFORM.fullmatch(e).group(1) for e in isoforms})
            for entity_id in isoforms:
                parent = _ISOFORM.fullmatch(entity_id).group(1)
                if parents.get(parent, (None, None, 0, False))[2] == 1:
                    projected[entity_id] = "uniprot:" + parent

        # 5. entity metadata, gene links
        wanted = {c[0] for cs in candidates.values() for c in cs}
        wanted |= set(projected.values())
        meta = self._meta(wanted - set(projected)) if wanted else {}
        links = {}
        proteins = {
            projected.get(c[0], c[0])
            for key, cs in candidates.items()
            if plans[key][0] == "gene"
            for c in cs
            if c[2] == "uniprot"
        }
        if proteins:
            links = self._by_protein(proteins)
            genes = {"entrez:" + g for p in proteins for g, _ in links.get(p, ())}
            meta.update(self._meta((genes | proteins) - set(meta)))

        # 6. posting per key: scope filter, then precedence
        out, nums = {}, {}
        for key, (target, route, ns, scope, identifier) in decoded.items():
            cls = plans[key][0] if key in plans else None
            rows = []
            for entity_id, tag, hub in candidates.get(key, ()):
                entity_id = projected.get(entity_id, entity_id)
                if cls == "gene":
                    rows.extend(self._gene_rows(entity_id, tag, hub, meta, links, quarantined))
                    continue
                m = meta.get(entity_id)
                if m is None or (cls == "entrez" and not m.exists):
                    continue
                if target == 1 and m.kind != 1:
                    continue
                if target == 3 and m.kind != 4:
                    continue
                rows.append(
                    (
                        entity_id,
                        m.kind,
                        m.anchor,
                        m.taxon,
                        m.quarantined or entity_id in quarantined,
                        m.reviewed,
                        [entity_id[7:]] if entity_id.startswith("entrez:") else m.gene_ids,
                        tag,
                    )
                )
            if scope is not None:
                rows = [r for r in rows if r[3] == scope]
            facts = self._posting(ns, rows)
            for fact in facts:
                if nums.setdefault(fact[0], fact[1]) != fact[1]:
                    raise ValueError(
                        f"Entity number collision between {nums[fact[0]]} and {fact[1]}"
                    )
            out[key] = dict(
                candidates=facts,
                gene=target == 2 and (route == 2 or ns in GENE_NAMESPACES),
                products=False,
            )
        return out

    @staticmethod
    def _gene_rows(entity_id, tag, hub, meta, links, quarantined):
        """Gene-level posting rows (rule 10): genes, or proteins linked to exactly one gene."""
        if hub == "entrez":
            genes = [(entity_id, None, None)] if entity_id.startswith("entrez:") else []
        else:
            link = links.get(entity_id, ())
            protein = meta.get(entity_id)
            ids = {g for g, _ in link}
            link_taxa = {t for _, t in link if t}
            if len(ids) != 1 or len(link_taxa) > 1 or protein is None or protein.quarantined:
                return []
            genes = [("entrez:" + next(iter(ids)), next(iter(link_taxa), None), protein.taxon)]
        rows = []
        for gene, link_taxon, protein_taxon in genes:
            m = meta.get(gene)
            if m is None or not m.exists or m.quarantined or gene in quarantined:
                continue
            taxa = {t for t in (m.taxon, link_taxon, protein_taxon) if t}
            if len(taxa) > 1:
                continue
            taxon = m.taxon or link_taxon or protein_taxon
            rows.append((gene, 3, None, taxon, False, False, [gene[7:]], tag))
        return rows

    def _meta(self, entity_ids):
        """Kind, anchor, taxon, quarantine, review state and gene links of entities."""
        ids = sorted(set(entity_ids))
        extras = self._extra(ids)
        proteins = {e: e[8:] for e in ids if e.startswith("uniprot:")}
        genes = {e: e[7:] for e in ids if e.startswith("entrez:")}
        others = {
            e: tuple(e.split(":", 1))
            for e in ids
            if e not in proteins and e not in genes and ":" in e
        }
        uniprot = self._records("uniprot", proteins.values())
        entrez = self._records("entrez", genes.values())
        by_gene = self._by_gene(set(genes.values()) - set(entrez))
        by_protein = self._by_protein(proteins)
        hub_records = defaultdict(set)
        for hub, local in others.values():
            if hub in self.hub_dirs and hub not in ("uniprot", "entrez"):
                hub_records[hub].add(local)
        records = {hub: self._records(hub, locals_) for hub, locals_ in hub_records.items()}
        result = {}
        for entity_id in ids:
            extra = extras.get(entity_id)
            kind, anchor, taxon, reviewed, gene_ids, exists = 1, None, None, False, [], False
            if entity_id in proteins:
                kind = 2
                record = uniprot.get(proteins[entity_id])
                exists = record is not None
                if record:
                    taxon, raw_anchor, count, reviewed = record
                    anchor = entity_id if count == 1 and raw_anchor else None
                gene_ids = sorted({g for g, _ in by_protein.get(entity_id, ())})
            elif entity_id in genes:
                kind = 3
                record = entrez.get(genes[entity_id])
                link = by_gene.get(genes[entity_id], ())
                exists = record is not None or bool(link)
                taxon = record[0] if record else min((t for _, t in link if t), default=None)
            elif entity_id in others:
                hub, local = others[entity_id]
                kind = (
                    4
                    if hub in REACTION_HUBS
                    else 1
                    if hub in CHEMICAL_HUBS or hub in ("inchikey", "goslin")
                    else 3
                )
                exists = False
                record = records.get(hub, {}).get(local)
                if record:
                    exists, taxon = True, record[0]
                if hub == "inchikey":
                    anchor = entity_id
            quarantined = False
            if extra:
                exists = True
                kind = KINDS.get(extra[0], kind) if isinstance(extra[0], str) else int(extra[0])
                taxon = extra[1] or taxon
                quarantined = extra[2]
            result[entity_id] = _Meta(
                kind, clean_anchor(kind, anchor), taxon, quarantined, reviewed, gene_ids, exists
            )
        return result

    @staticmethod
    def _posting(ns, rows):
        """Apply precedence (rule 8) to one posting and build its candidate facts."""
        if not rows:
            return []
        tags = {r[7] for r in rows}
        dropped = set()
        if NATIVE in tags:
            dropped.add(FALLBACK)  # native identity rows beat fallback rows
        if ns == "uniprot" and tags - {SECONDARY}:
            dropped.add(SECONDARY)  # a primary accession beats secondary-accession claims
        exact_taxa = {r[3] for r in rows if r[7] != SYMBOL_SYNONYM}
        facts = {}
        for entity_id, kind, anchor, taxon, quarantined, reviewed, gene_ids, tag in sorted(
            rows, key=lambda r: (r[0], r[7])
        ):
            if tag in dropped or (tag == SYMBOL_SYNONYM and taxon in exact_taxa):
                continue  # exact symbols beat synonyms within one taxon
            if entity_id in facts:
                continue
            facts[entity_id] = [
                entity_num(entity_id),
                entity_id,
                kind,
                anchor,
                bool(quarantined),
                bool(reviewed),
                [str(g) for g in (gene_ids or ())],
            ]
        return sorted(facts.values(), key=lambda f: f[0])

    def lookup(self, key):
        return self.lookup_many([key])[key]

    # --------------------------------------------------------------- records

    def record_many(self, entity_ids):
        wanted = list(dict.fromkeys(entity_ids))
        found = self._cached(wanted)
        missing = [eid for eid in wanted if eid not in found]
        if missing:
            built = self._build_records(missing)
            absent = [eid for eid in missing if eid not in built]
            if absent:
                raise ValueError("Candidate references absent entity: " + absent[0])
            self._store(built)
            found.update(built)
        return {eid: found[eid] for eid in wanted}

    def record(self, entity_id):
        return self.record_many([entity_id])[entity_id]

    def _members(self, entity_ids, extras):
        """entity_id -> member record ids (by_id anchors, own record, exception_members)."""
        members = defaultdict(list)
        for entity_id, record_ids in self._exception_members(entity_ids).items():
            members[entity_id].extend(record_ids)
        # anchored entities: the records whose by_id anchor row names the entity
        lookups = []  # (entity id, hubs, source ns, identifier)
        for entity_id in entity_ids:
            if entity_id.startswith("inchikey:"):
                lookups.append((entity_id, self.chemical_hubs, "inchikey", entity_id[9:]))
            elif entity_id.startswith("uniprot:"):
                lookups.append((entity_id, ["uniprot"], "uniprot", entity_id[8:]))
        by_hubs = defaultdict(list)
        for item in lookups:
            by_hubs[tuple(item[1])].append(item)
        for hubs, items in by_hubs.items():
            pairs = [(ns, identifier) for _, _, ns, identifier in items]
            rows = self._by_id(hubs, pairs)
            exceptions = self._exceptions(f"{hub}:{local}" for _, hub, local, *_ in rows)
            for i, hub, local_id, tag, anchor, count in rows:
                entity_id = items[i][0]
                if hub == "uniprot" and tag != "native":
                    continue
                if self._entity(hub, local_id, anchor, count, exceptions)[0] == entity_id:
                    members[entity_id].append(f"{hub}:{local_id}")
        # record-id entities (own record) and entrez genes (the entrez record)
        own = {}
        for entity_id in entity_ids:
            hub = entity_id.split(":", 1)[0]
            if entity_id.startswith(("inchikey:", "goslin:")) or hub not in self.hub_dirs:
                continue
            if hub == "uniprot" and entity_id not in members:
                own[entity_id] = entity_id
            elif hub != "uniprot":
                own[entity_id] = entity_id
        if own:
            exceptions = self._exceptions(own.values())
            for entity_id, record_id in own.items():
                if exceptions.get(record_id, (entity_id,))[0] == entity_id:
                    members[entity_id].append(record_id)
        return {e: sorted(set(r)) for e, r in members.items()}

    def _by_record(self, record_ids):
        """record id -> [(source_type, value)] from the hubs' by_record rows."""
        by_hub = defaultdict(set)
        for record_id in set(record_ids):
            hub, _, local = record_id.partition(":")
            if hub in self.hub_kv:
                by_hub[hub].add(local)
        rows = defaultdict(list)
        for hub, locals_ in by_hub.items():
            for local_id, pairs in self.hub_kv[hub].rows(locals_).items():
                if pairs:
                    rows[f"{hub}:{local_id}"].extend(pairs)
        return rows

    def _build_records(self, entity_ids):
        entity_ids = list(entity_ids)
        extras = self._extra(entity_ids)
        members = self._members(entity_ids, extras)
        rows = self._by_record({r for ids in members.values() for r in ids})
        meta = self._meta(entity_ids)

        # gene aliases: genesymbol/-syn, hgnc, ensg of proteins linked to exactly that one gene
        genes = [e for e in entity_ids if e.startswith("entrez:")]
        linked = defaultdict(list)
        if genes:
            by_gene = self._by_gene(e[7:] for e in genes)
            proteins = {p for links in by_gene.values() for p, _ in links}
            by_protein = self._by_protein(proteins)
            alias_rows = self._by_record({p for p in proteins if p.startswith("uniprot:")})
            for gene in genes:
                gene_taxon = meta[gene].taxon
                for protein, link_taxon in by_gene.get(gene[7:], ()):
                    ids = {g for g, _ in by_protein.get(protein, ())}
                    taxa = {t for t in (gene_taxon, link_taxon) if t}
                    if ids == {gene[7:]} and len(taxa) <= 1:
                        linked[gene].extend(
                            ("uniprot", st, v)
                            for st, v in alias_rows.get(protein, ())
                            if st in GENE_ALIAS_TYPES
                        )

        built = {}
        for entity_id in entity_ids:
            m = meta[entity_id]
            member_ids = members.get(entity_id, [])
            hub_rows = []
            for record_id in member_ids:
                hub = record_id.split(":", 1)[0]
                hub_rows.extend((hub, st, v) for st, v in rows.get(record_id, ()))
            if not (hub_rows or m.exists):
                continue
            hub_rows_all = hub_rows + linked.get(entity_id, [])
            identifiers = {
                (st, str(v))
                for _, st, v in hub_rows_all
                if st not in NOT_IDENTIFIERS and v is not None
            }
            identifiers.update(tuple(r.split(":", 1)) for r in member_ids)  # each record's own id
            if ":" in entity_id:
                identifiers.add(tuple(entity_id.split(":", 1)))
            built[entity_id] = dict(
                entity_id=entity_id,
                kind=m.kind,
                anchor=m.anchor,
                taxon=m.taxon,
                label=choose_label(m.kind, entity_id, hub_rows),
                identifiers=[list(p) for p in sorted(identifiers)],
                gene_ids=m.gene_ids,
            )
        return built

    # ------------------------------------------------------------- resolve

    # --------------------------------------------------------------- resolve

    def term_label(self, term: str) -> str | None:
        """The ontology's name of a term id (``term_labels.parquet``; absent in older snapshots)."""
        if self._term_labels is None:
            path = self.path / "term_labels.parquet"
            self._term_labels = {}
            if path.is_file():
                import pyarrow.parquet as pq

                table = pq.read_table(path, columns=["term", "label"]).to_pydict()
                self._term_labels = dict(zip(table["term"], table["label"]))
        return self._term_labels.get(term) or self._term_labels.get(term.upper())

    def close(self):
        if self._cache is not None:
            self._cache.close()
            self._cache = None
        for kv in getattr(self, "hub_kv", {}).values():
            kv.close()
        if getattr(self, "decisions", None) is not None:
            self.decisions.close()
            self.decisions = None

    def resolve(self, queries, votes, *, decision_batch_size=4096):
        start = time.perf_counter()
        keys = sorted({v["lookup_key"] for v in votes})
        selected = self.lookup_many(keys)
        postings, metadata = [], {}
        for key in keys:
            value = selected[key]
            postings.append((key, [c[0] for c in value["candidates"]]))
            metadata.update({c[0]: tuple(c) for c in value["candidates"]})
        grouped = defaultdict(list)
        primary_symbol = {
            v["input_id"]
            for v in votes
            if v["ns"] == "genesymbol" and selected[v["lookup_key"]]["candidates"]
        }
        for v in votes:
            if v["ns"] == "genesymbol-syn" and v["input_id"] in primary_symbol:
                continue
            value = selected[v["lookup_key"]]
            grouped[v["input_id"]].append(
                (
                    v["lookup_key"],
                    v["anchor"],
                    v["route"],
                    v["ordinal"],
                    value["gene"] or v.get("gene_only", False),
                    value["products"],
                )
            )
        batches = [(q["input_id"], q["target"], grouped[q["input_id"]]) for q in queries]
        fetched = time.perf_counter()
        result, molecular = self._decide(batches, postings, metadata, decision_batch_size)
        result = self._fallbacks(result, batches, votes, postings, metadata, decision_batch_size)
        decided = time.perf_counter()
        accepted = {eid for _, _, ids, _ in result for eid in ids}
        accepted.update(
            v["protein_entity_id"] for v in molecular.values() if v["protein_entity_id"]
        )
        records = self.record_many(sorted(accepted))
        finished = time.perf_counter()
        return dict(
            results=[
                dict(
                    input_id=i,
                    outcome=o,
                    entities=sorted(ids),
                    candidate_count=n,
                    **molecular.get(i, {}),
                )
                for i, o, ids, n in result
            ],
            records=records,
        ), dict(
            lookup_seconds=fetched - start,
            decision_seconds=decided - fetched,
            entity_fetch_seconds=finished - decided,
            total_seconds=finished - start,
            unique_keys=len(keys),
            candidate_records=len(metadata),
            entity_records=len(records),
        )

    def _fallbacks(self, result, batches, votes, postings, metadata, decision_batch_size):
        """Decide unaccepted chemicals again on narrower evidence, in this order:

        1. without secondary coarse cross-references (KEGG, BiGG, MetaNetX compounds), which
           must not veto the specific identifiers;
        2. on the source's own structure (a stated InChIKey or one derived from its SMILES),
           which decides when the source's identifiers contradict each other;
        3. when every candidate of the specific votes is one molecule in different protonation
           states (InChIKeys equal but for the last character): on the votes for its neutral
           form (``-N``) if the source names it, else on the primary identifier.

        Accepted results never change, and the narrower evidence must still be unique.
        """
        by_vote = {(v["input_id"], v["ordinal"]): v for v in votes if v["target"] == 1}
        groups = {i: (t, group) for i, t, group in batches if t == 1}
        entities = {key: [metadata[n][1] for n in nums] for key, nums in postings}

        def specific(i, vote):
            v = by_vote[(i, vote[3])]
            return v.get("primary") or v["ns"] not in COARSE_CHEMICAL_NAMESPACES

        def structure(i, vote):
            return by_vote[(i, vote[3])]["ns"] == "inchikey"

        def candidates(i):
            """Candidates of the specific votes: coarse ones do not veto here either."""
            return {
                e for vote in groups[i][1] if specific(i, vote) for e in entities.get(vote[0], ())
            }

        def protonation(i, vote):
            if not specific(i, vote):
                return False
            neutral = [e for e in candidates(i) if e.endswith("-N")]
            if neutral:
                return neutral[0] in entities.get(vote[0], ())
            return bool(by_vote[(i, vote[3])].get("primary"))

        def protonation_states(i):
            ids = candidates(i)
            return (
                len(ids) > 1
                and all(e.startswith("inchikey:") for e in ids)
                and len({e[:-1] for e in ids}) == 1
            )

        for keep, eligible in (
            (specific, None),
            (structure, None),
            (protonation, protonation_states),
        ):
            retry = []
            for i, _, ids, _ in result:
                if ids or i not in groups or (eligible and not eligible(i)):
                    continue
                t, group = groups[i]
                narrowed = [vote for vote in group if keep(i, vote)]
                if narrowed and len(narrowed) < len(group):
                    retry.append((i, t, narrowed))
            if retry:
                second, _ = self._decide(retry, postings, metadata, decision_batch_size)
                better = {row[0]: row for row in second if row[2]}
                result = [better.get(row[0], row) for row in result]
        return result

    @staticmethod
    def _decide(batches, postings, metadata, decision_batch_size):
        from omnipath_resolver._omnipath_resolver import (
            resolve_precomputed_batch,
            resolve_molecular_batch,
        )

        result = []
        molecular = {}
        for offset in range(0, len(batches), decision_batch_size):
            chunk = batches[offset : offset + decision_batch_size]
            # Reactions take the chemical decision path: same intersection, quarantine and
            # unique/ambiguous outcome, with no structure anchors. The kernel only knows
            # chemical candidates there, so reaction queries and entities are passed as such.
            chemical = [(q[0], 1, q[2]) for q in chunk if q[1] in (1, 3)]
            biological = [q for q in chunk if q[1] == 2]
            if chemical:
                result.extend(
                    resolve_precomputed_batch(
                        chemical,
                        postings,
                        [(m[0], m[1], 1 if m[2] == 4 else m[2], *m[3:6]) for m in metadata.values()],
                    )
                )
            if biological:
                facts = [
                    tuple(m[:6]) + (list(m[6]) if len(m) > 6 else [],) for m in metadata.values()
                ]
                for row in resolve_molecular_batch(biological, postings, facts):
                    result.append(row[:4])
                    molecular[row[0]] = dict(
                        protein_entity_id=row[4], gene_candidates=row[5], gene_mapping_status=row[6]
                    )
        return result, molecular
