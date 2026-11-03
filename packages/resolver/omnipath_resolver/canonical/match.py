"""Normalize observation evidence and resolve it through the identity library."""

from __future__ import annotations

import copy
import json
import re
from collections import OrderedDict, defaultdict
from dataclasses import dataclass, field
from typing import Any
from omnipath_core.biolink import entity_type as normalize_entity_type
from omnipath_core.naming import normalize_namespace
from .identifiers import normalize_identifier, normalize_ns
from .label import assign_label
from .library import pin_library
from .policy import EntityPolicy, CHEMICAL, PROTEIN_ENTITY_TYPES, get_policy

# Record kind -> reference library (1 chemical, 2 protein, 3 gene, 4 reaction).
REFERENCE_LIBRARY = {1: "chemical", 4: "reaction"}

RESOLVED_BY_MATCHED = ("parquet",)


@dataclass(frozen=True)
class Vote:
    ns: str
    id: str
    taxon: str | None
    kind: str  # "id" | "symbol"
    primary: bool
    # Computed here (SMILES -> InChIKey), not stated by the source: never a primary anchor.
    derived: bool = False

    @property
    def key(self) -> tuple[str, str, str | None]:
        return (self.ns, self.id, self.taxon if self.kind == "symbol" else None)


@dataclass(frozen=True)
class Node:
    node_id: str
    taxon: str | None
    canonical_ns: str
    canonical_id: str
    label: str | None
    aliases: tuple[tuple[str, str], ...]
    uniprot: str | None = None
    entrez: str | None = None


@dataclass
class Match:
    canonical_namespace: str
    canonical_identifier: str
    label: str
    resolved_by: str
    taxon: str | None = None
    node_id: str | None = None
    aliases: dict[str, list[str]] = field(default_factory=dict)

    entity_type: str | None = None
    reference_library: str | None = None
    protein_namespace: str | None = None
    protein_identifier: str | None = None
    protein_node_id: str | None = None
    protein_label: str | None = None
    protein_aliases: dict[str, list[str]] = field(default_factory=dict)
    protein_taxon: str | None = None
    protein_gene_candidates: tuple[str, ...] = ()
    gene_mapping_status: str | None = None
    gene_candidates: tuple[str, ...] = ()
    transcript_namespace: str | None = None
    transcript_identifier: str | None = None

    @property
    def matched(self) -> bool:
        return self.node_id is not None


def _observation_signature(obs: Any) -> str:
    """Everything resolution reads from an observation; not its key, annotations or evidence."""
    return json.dumps(
        [
            getattr(obs, name, None)
            for name in (
                "entity_type",
                "namespace",
                "identifier",
                "taxon",
                "label",
                "identity_scope",
                "identifiers",
                "molecular_form",
            )
        ],
        sort_keys=True,
        default=str,
    )


def votes_for(obs: Any, policy: EntityPolicy) -> tuple[list[Vote], dict[str, list[str]]]:
    """Votes plus every observed identifier grouped by normalized namespace."""
    taxon = str(getattr(obs, "taxon", None) or "").strip() or None
    observed: dict[str, list[str]] = defaultdict(list)
    votes: list[Vote] = []
    seen: set[tuple[str, str]] = set()
    if hasattr(obs, "structure_derivations"):
        obs.structure_derivations = []

    def add(ns_raw: Any, id_raw: Any, primary: bool) -> None:
        if str(ns_raw) == "ramp" and policy.library == "gene_protein":
            ns_raw = "ramp_gene"
        pairs = normalize_identifier(ns_raw, id_raw)
        if not pairs:
            slug = normalize_ns(ns_raw)
            text = str(id_raw or "").strip()
            if slug and text and text not in observed[slug]:
                observed[slug].append(text)
            return
        for ns, ident in pairs:
            # Generic RefSeq namespaces still carry a typed accession. Use its
            # explicit RNA/protein prefix for lookup, retaining versions in forms.
            if ns in {"refseq", "refseq_protein"}:
                if re.match(r"^(AP|NP|XP|YP|WP|ZP)_", ident):
                    ns = "refseq_protein"
                    ident = re.sub(r"\.[0-9]+$", "", ident)
                elif re.match(r"^(NM|NR|XM|XR)_", ident):
                    ns = "refseq"
                    ident = re.sub(r"\.[0-9]+$", "", ident)
            if ident not in observed[ns]:
                observed[ns].append(ident)
            if (ns, ident) in seen:
                continue
            seen.add((ns, ident))
            if ns in policy.match_namespaces:
                votes.append(Vote(ns, ident, taxon, "id", primary))
            elif ns in policy.symbol_namespaces and taxon:
                votes.append(Vote(ns, ident, taxon, "symbol", primary))

    add(getattr(obs, "namespace", ""), getattr(obs, "identifier", ""), True)
    for item in getattr(obs, "identifiers", None) or []:
        if isinstance(item, dict):
            add(
                item.get("ns") or item.get("type") or "",
                item.get("id") or item.get("identifier") or "",
                False,
            )
    form = getattr(obs, "molecular_form", None) or {}
    if policy.library == "gene_protein" and isinstance(form, dict):
        specific = list(form.get("sequence_identifiers") or [])
        if form.get("isoform_identifier"):
            specific.append(form["isoform_identifier"])
        for item in specific:
            if isinstance(item, dict):
                add(item.get("ns"), item.get("id"), False)
    if policy.library == CHEMICAL and not observed.get("inchikey") and observed.get("smiles"):
        from .structures import cached_derivation

        derivations = [cached_derivation(s) for s in sorted(set(observed["smiles"]))]
        if hasattr(obs, "structure_derivations"):
            obs.structure_derivations = derivations
        for result in derivations:
            if result["status"] == "derived" and ("inchikey", result["inchikey"]) not in seen:
                seen.add(("inchikey", result["inchikey"]))
                votes.append(Vote("inchikey", result["inchikey"], None, "id", False, derived=True))
    return votes, observed


class LibraryMatcher:
    """Batch resolution against one pinned immutable identity library."""

    def __init__(self, library_dir, *, defer_aliases=False, memory_limit=None):
        self.library_dir = pin_library(library_dir)
        self.defer_aliases = defer_aliases
        self.memory_limit = memory_limit
        # Resolution against a pinned reference depends only on the observation, so a
        # worker reuses the result when the same observation recurs in later batches
        # (FooDB repeats ~70k compounds across ~5M content rows).
        self._memo: OrderedDict[str, tuple[list, list]] = OrderedDict()
        self.metrics = {
            "batches": 0,
            "observations": 0,
            "memo_hits": 0,
            "lookup_seconds": 0.0,
            "kernel_seconds": 0.0,
            "enrichment_seconds": 0.0,
        }
        from ..identity_runtime import open_runtime

        self.runtime = (
            open_runtime(self.library_dir, memory_limit=memory_limit)
            if self.library_dir is not None and self.library_dir.exists()
            else None
        )
        self.libraries = list(self.runtime.libraries) if self.runtime else []

    def available(self, library):
        return library in self.libraries

    def targets(self, entities):
        from ..observations import (
            observation_bundle,
            flush_lipid_cache,
        )

        prepared = {}
        queries, votes = [], []
        cached, signatures = {}, {}
        for key, obs in entities.items():
            signature = _observation_signature(obs)
            hit = self._memo.get(signature)
            if hit is not None:
                self._memo.move_to_end(signature)
                matches, derivations = hit
                obs.structure_derivations = copy.deepcopy(derivations)
                cached[key] = copy.deepcopy(matches)
                continue
            signatures[key] = signature
            policy = get_policy(obs.entity_type)
            normalized, observed = votes_for(obs, policy)
            prepared[key] = (obs, policy, observed)
            if not self.available(policy.library):
                continue
            query, rows = observation_bundle(key, obs, normalized, observed, policy.library)
            queries.append(query)
            votes.extend(rows)
        results = {
            key: [self._build(obs, policy, observed, None, "unmatched")]
            for key, (obs, policy, observed) in prepared.items()
        }
        for key, matches in results.items():
            for match in matches:
                match.entity_type = prepared[key][0].entity_type
                if prepared[key][1].library == "gene_protein":
                    match.gene_mapping_status = "missing"
                    _retain_asserted_protein(match, prepared[key][0])
                    transcript = _asserted_transcript(prepared[key][0])
                    if transcript:
                        match.transcript_namespace, match.transcript_identifier = transcript
        if not queries:
            return self._remember(results, cached, signatures, entities)
        resolved, metrics = self.runtime.resolve(queries, votes)
        accepted = {row["input_id"]: row for row in resolved["results"]}
        for name, source in (
            ("lookup_seconds", "lookup_seconds"),
            ("kernel_seconds", "decision_seconds"),
            ("enrichment_seconds", "entity_fetch_seconds"),
        ):
            self.metrics[name] += metrics[source]
        for key, resolution in accepted.items():
            ids = resolution["entities"]
            if ids:
                obs, policy, observed = prepared[key]
                matches = []
                for eid in sorted(set(ids)):
                    row = resolved["records"][eid]
                    ns, identifier = eid.split(":", 1)
                    aliases = [tuple(value) for value in row["identifiers"]]
                    reference_aliases = defaultdict(list)
                    for alias_ns, alias_id in aliases:
                        reference_aliases[alias_ns].append(alias_id)
                    node = Node(
                        eid,
                        row["taxon"] or None,
                        ns,
                        identifier,
                        row["label"] if row["label"] != identifier else None,
                        tuple(aliases),
                    )
                    match = self._build(obs, policy, observed, node, "parquet")
                    match.entity_type = obs.entity_type
                    match.reference_library = REFERENCE_LIBRARY.get(row["kind"], "gene_protein")
                    matches.append(match)
                results[key] = matches
            obs, policy, observed = prepared[key]
            for match in results[key]:
                match.entity_type = obs.entity_type
                if policy.library != "gene_protein":
                    continue
                match.gene_mapping_status = resolution.get("gene_mapping_status", "missing")
                match.gene_candidates = tuple(resolution.get("gene_candidates", []))
                transcript = _asserted_transcript(obs)
                if transcript:
                    match.transcript_namespace, match.transcript_identifier = transcript
                product = resolution.get("protein_entity_id")
                if product:
                    row = resolved["records"][product]
                    match.protein_namespace, match.protein_identifier = product.split(":", 1)
                    match.protein_node_id = product
                    match.protein_label = row["label"]
                    match.protein_taxon = row["taxon"] or None
                    match.protein_gene_candidates = tuple(
                        "entrez:" + g for g in row.get("gene_ids", [])
                    )
                    aliases = defaultdict(list)
                    for ns, value in row["identifiers"]:
                        aliases[ns].append(value)
                    match.protein_aliases = dict(aliases)
                else:
                    _retain_asserted_protein(match, obs)
        flush_lipid_cache()
        self.metrics["batches"] += 1
        self.metrics["observations"] += len(queries)
        return self._remember(results, cached, signatures, entities)

    def _remember(self, results, cached, signatures, entities, limit=200_000):
        for key, signature in signatures.items():
            self._memo[signature] = (
                copy.deepcopy(results[key]),
                copy.deepcopy(entities[key].structure_derivations),
            )
        while len(self._memo) > limit:
            self._memo.popitem(last=False)
        self.metrics["memo_hits"] += len(cached)
        results.update(cached)
        return results

    def match(self, entities):
        """Scalar convenience API; ambiguous fan-out retains the input identity."""
        result = {}
        for key, matches in self.targets(entities).items():
            if len(matches) == 1:
                result[key] = matches[0]
            else:
                obs = entities[key]
                policy = get_policy(obs.entity_type)
                _, observed = votes_for(obs, policy)
                result[key] = self._build(obs, policy, observed, None, "multiple_products")
        return result

    def chemical_targets(self, entities):
        return self.targets(entities)

    def _build(
        self,
        obs: Any,
        policy: EntityPolicy,
        observed: dict[str, list[str]],
        node: Node | None,
        rule: str,
    ) -> Match:
        aliases: dict[str, list[str]] = {ns: list(vals) for ns, vals in observed.items()}
        if node is not None:
            for ns, ident in node.aliases:
                if ns in policy.alias_namespaces:
                    bucket = aliases.setdefault(ns, [])
                    if ident not in bucket:
                        bucket.append(ident)
            canonical_ns, canonical_id = (node.canonical_ns, node.canonical_id)
            for ns, ident in ((node.canonical_ns, node.canonical_id), (canonical_ns, canonical_id)):
                if ident not in aliases.setdefault(ns, []):
                    aliases[ns].append(ident)
            label = node.label or assign_label(canonical_ns, canonical_id, aliases, policy)
            return Match(
                canonical_namespace=canonical_ns,
                canonical_identifier=canonical_id,
                label=label,
                resolved_by=rule,
                taxon=node.taxon or (str(getattr(obs, "taxon", None) or "").strip() or None),
                node_id=node.node_id,
                aliases=observed if self.defer_aliases else aliases,
            )

        # A source-supplied full structure remains its identity even when the
        # catalogue has not seen it yet. Conflicting structures stay unresolved.
        from .policy import INCHIKEY_RE

        structures = {v for v in observed.get("inchikey", []) if INCHIKEY_RE.fullmatch(v)}
        if policy.library == CHEMICAL and len(structures) == 1:
            identifier = next(iter(structures))
            return Match(
                "inchikey",
                identifier,
                assign_label("inchikey", identifier, aliases, policy),
                rule,
                taxon=str(getattr(obs, "taxon", None) or "").strip() or None,
                aliases=aliases,
            )

        raw_ns = str(getattr(obs, "namespace", "") or "")
        raw_id = str(getattr(obs, "identifier", "") or "").strip()
        ns = normalize_namespace(raw_ns) or normalize_ns(raw_ns) or "unknown"
        ident = raw_id
        norm = normalize_identifier(raw_ns, raw_id)
        if norm:
            ns, ident = norm[0]
        for pattern, target_ns in policy.rewrites:
            if re.match(pattern, ident, re.IGNORECASE):
                ns = target_ns
                break
        if ident not in aliases.setdefault(ns, []):
            aliases[ns].append(ident)
        label = assign_label(ns, ident, aliases, policy)
        if policy.entity_class == "cv_term" and label == ident:
            # A term cited by its id alone (UniProt's GO terms) is named by its ontology.
            label = (self.runtime.term_label(ident) if self.runtime else None) or label
        return Match(
            canonical_namespace=ns,
            canonical_identifier=ident,
            label=label,
            resolved_by=rule,
            taxon=str(getattr(obs, "taxon", None) or "").strip() or None,
            node_id=None,
            aliases=aliases,
        )

    def close(self):
        from ..observations import flush_lipid_cache

        flush_lipid_cache()
        if self.runtime is not None:
            self.runtime.close()
            self.runtime = None


def _asserted_transcript(obs):
    """Retain participant transcript assertions on matched and unmatched inputs."""
    candidates = set()
    raw_id = str(getattr(obs, "identifier", "") or "").strip()
    pairs = normalize_identifier(getattr(obs, "namespace", ""), raw_id)
    if pairs:
        ns = pairs[0][0]
        if ns == "enst" or (ns == "refseq" and re.match(r"^(NM|NR|XM|XR)_", raw_id)):
            candidates.add((ns, raw_id))
    form = getattr(obs, "molecular_form", None) or {}
    sequences = (form.get("sequence_identifiers") or []) if isinstance(form, dict) else []
    for item in sequences:
        if not isinstance(item, dict):
            continue
        identifier = str(item.get("id", ""))
        pairs = normalize_identifier(item.get("ns", ""), identifier)
        if not pairs:
            continue
        ns = pairs[0][0]
        if ns == "enst" or (ns == "refseq" and re.match(r"^(NM|NR|XM|XR)_", identifier)):
            candidates.add((ns, identifier))
    return next(iter(candidates)) if len(candidates) == 1 else None


def _asserted_protein(obs):
    """Retain one explicitly reported product, without asserting catalogue status.

    Only a protein participant's principal identifier and explicit form IDs
    establish an occurrence. Other identifier aliases are not participant claims.
    Exact versions, isoform suffixes and chain suffixes remain source identities;
    a parent accession must be supplied or resolved independently.
    Without a catalogue-selected product, an explicit entry plus an isoform
    remains conservatively ambiguous, even when their accession bases agree.
    Both assertions remain in the observation; no primary product is chosen.
    """
    if normalize_entity_type(obs.entity_type) not in PROTEIN_ENTITY_TYPES:
        return None
    candidates = set()

    def add(namespace, identifier):
        raw = str(identifier or "").strip()
        pairs = normalize_identifier(namespace, raw)
        if not pairs:
            return
        ns = pairs[0][0]
        if ns == "uniprot":
            pattern = (
                r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
                r"(?:-\d+)?(?:-PRO_\d+)?"
            )
        elif ns == "ensp":
            pattern = r"ENS[A-Z]*P\d+(?:\.\d+)?"
        elif ns in {"refseq", "refseq_protein"}:
            pattern = r"(?:AP|NP|XP|YP|WP|ZP)_[0-9]+(?:\.\d+)?"
        else:
            return
        if re.fullmatch(pattern, raw):
            candidates.add((ns, raw))

    add(getattr(obs, "namespace", ""), getattr(obs, "identifier", ""))
    form = getattr(obs, "molecular_form", None) or {}
    if isinstance(form, dict):
        identifiers = list(form.get("sequence_identifiers") or [])
        if form.get("isoform_identifier"):
            identifiers.append(form["isoform_identifier"])
        for item in identifiers:
            if isinstance(item, dict):
                add(item.get("ns"), item.get("id"))
    # A chain form also carries its enclosing isoform as specificity metadata.
    # Retain the explicitly asserted chain rather than inventing its parent.
    chain_parents = {
        (ns, identifier.rsplit("-PRO_", 1)[0])
        for ns, identifier in candidates
        if ns == "uniprot" and "-PRO_" in identifier
    }
    candidates -= chain_parents
    return next(iter(candidates)) if len(candidates) == 1 else None


def _retain_asserted_protein(match, obs):
    product = _asserted_protein(obs)
    if product is None:
        return
    match.protein_namespace, match.protein_identifier = product
    match.protein_label = getattr(obs, "label", None) or product[1]
    match.protein_taxon = str(getattr(obs, "taxon", None) or "").strip() or None
    match.protein_aliases = {product[0]: [product[1]]}
