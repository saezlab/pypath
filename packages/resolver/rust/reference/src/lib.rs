//! Batch resolution over an already constructed, immutable identity library.
//! This module implements the serving policy.
//! Namespace-specific normalization and relation routing happen before this API.
use std::collections::{BTreeMap, BTreeSet};
use std::io;

pub type EntityId = u64;
pub type Error = Box<dyn std::error::Error + Send + Sync>;
pub type Result<T> = std::result::Result<T, Error>;
pub(crate) fn invalid(message: &str) -> Error {
    io::Error::new(io::ErrorKind::InvalidData, message.to_owned()).into()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
pub enum Kind {
    Chemical = 1,
    Protein = 2,
    Gene = 3,
    Transcript = 4,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum Anchor {
    InchiKey([u8; 27]),
    UniProt(String),
}
impl Anchor {
    /// Syntax only. Semantic validity (including empty-structure sentinels) must
    /// already have been checked by the normalizer/reference builder.
    pub fn inchikey(value: &str) -> Result<Self> {
        let b = value.as_bytes();
        if b.len() != 27
            || b[14] != b'-'
            || b[25] != b'-'
            || b.iter()
                .enumerate()
                .any(|(i, c)| i != 14 && i != 25 && !c.is_ascii_uppercase())
        {
            return Err(invalid("invalid InChIKey syntax"));
        }
        let mut key = [0u8; 27];
        key.copy_from_slice(b);
        Ok(Self::InchiKey(key))
    }
    pub fn kind(&self) -> Kind {
        match self {
            Self::InchiKey(_) => Kind::Chemical,
            Self::UniProt(_) => Kind::Protein,
        }
    }
    /// The policy means same first hash block, not independently verified graph identity.
    pub fn first_block(&self) -> Option<[u8; 14]> {
        match self {
            Self::InchiKey(key) => {
                let mut result = [0u8; 14];
                result.copy_from_slice(&key[..14]);
                Some(result)
            }
            _ => None,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Role {
    /// Only use when the caller genuinely asserts a full primary anchor.
    /// A UniProt-looking accession is not necessarily a PRIMARY accession.
    PrimaryAnchor(Anchor),
    Identity,
    /// A gene identifier whose agreed protein candidates prefer reviewed entries.
    GeneIdentity {
        multiple_products: bool,
    },
    GeneToProtein,
    Connectivity,
}
impl Role {
    fn route(&self) -> u8 {
        match self {
            Self::PrimaryAnchor(_) | Self::Identity | Self::GeneIdentity { .. } => 1,
            Self::GeneToProtein => 2,
            Self::Connectivity => 3,
        }
    }
}

/// Collision-free key encoding of a normalized identifier and lookup route.
/// Namespace IDs and normalization rules must be pinned in the snapshot manifest.
/// `scope` can encode taxon for identifiers requiring it; None is NOT a wildcard.
pub fn encode_key(
    target: Kind,
    role: &Role,
    namespace: u16,
    scope: Option<u32>,
    local_id: &[u8],
) -> Vec<u8> {
    let mut key = Vec::with_capacity(10 + local_id.len());
    key.extend_from_slice(&[1, target as u8, role.route()]);
    key.extend_from_slice(&namespace.to_be_bytes());
    match scope {
        Some(taxon) => {
            key.push(1);
            key.extend_from_slice(&taxon.to_be_bytes());
        }
        None => key.push(0),
    }
    key.extend_from_slice(local_id);
    key
}

#[derive(Debug, Clone)]
pub struct Evidence {
    /// Index of the original identifier in the caller's input; raw input stays with caller.
    pub input_ordinal: usize,
    pub key: Vec<u8>,
    pub role: Role,
}
#[derive(Debug, Clone)]
pub struct Query {
    pub target: Kind,
    pub evidence: Vec<Evidence>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Lookup {
    Missing,
    /// Strictly increasing, non-empty, complete within the pinned snapshot and route.
    Complete(Vec<EntityId>),
    /// No success may be inferred from an incomplete or deliberately capped lookup.
    Incomplete {
        total: Option<u64>,
        reason: String,
    },
}

pub trait PostingIndex {
    /// Input is sorted and deduplicated. Results MUST align one-for-one with keys.
    /// Unsupported namespaces, unreadable shards and failed I/O are errors, not Missing.
    fn lookup_many(
        &self,
        keys: &[Vec<u8>],
        max_postings: usize,
        max_total_postings: usize,
    ) -> Result<Vec<Lookup>>;
}
#[derive(Debug, Clone)]
pub struct EntityMeta {
    pub id: EntityId,
    pub kind: Kind,
    pub anchor: Option<Anchor>,
    pub quarantined: bool,
    pub reviewed: bool,
}
pub trait MetadataIndex {
    /// Results align with sorted, deduplicated IDs. Missing metadata is index corruption.
    fn metadata_many(&self, ids: &[EntityId]) -> Result<Vec<Option<EntityMeta>>>;
}

#[derive(Debug, Clone)]
pub struct Policy {
    pub allow_same_connectivity_multiple: bool,
    pub require_primary_anchor: bool,
    pub allow_entrez_only_fallback: bool,
    pub allow_gene_to_protein: bool,
    pub allow_gene_to_protein_multiple: bool,
    pub allow_connectivity_queries: bool,
    pub max_postings_per_identifier: usize,
    pub max_final_candidates: usize,
    pub max_total_postings_per_batch: usize,
    pub max_total_candidates_per_batch: usize,
}
impl Default for Policy {
    fn default() -> Self {
        Self {
            allow_same_connectivity_multiple: false,
            require_primary_anchor: false,
            allow_entrez_only_fallback: false,
            allow_gene_to_protein: false,
            allow_gene_to_protein_multiple: false,
            allow_connectivity_queries: false,
            max_postings_per_identifier: 100_000,
            max_final_candidates: 10_000,
            max_total_postings_per_batch: 5_000_000,
            max_total_candidates_per_batch: 1_000_000,
        }
    }
}
impl Policy {
    /// The source-structure-first policy used by OmniPath serving adapters.
    pub fn omnipath() -> Self {
        Self {
            allow_same_connectivity_multiple: false,
            allow_gene_to_protein: false,
            allow_entrez_only_fallback: false,
            allow_gene_to_protein_multiple: false,
            ..Self::default()
        }
    }
    fn enabled(&self, role: &Role) -> bool {
        match role {
            Role::GeneToProtein | Role::GeneIdentity { .. } => self.allow_gene_to_protein,
            Role::Connectivity => self.allow_connectivity_queries,
            _ => true,
        }
    }
}
fn participates(query: &Query, evidence: &Evidence, policy: &Policy) -> bool {
    policy.enabled(&evidence.role)
        && !(query.target == Kind::Chemical
            && query
                .evidence
                .iter()
                .any(|e| matches!(e.role, Role::PrimaryAnchor(_)))
            && !matches!(evidence.role, Role::PrimaryAnchor(_)))
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Outcome {
    Unique,
    MultipleSameConnectivity,
    MultipleGeneProducts,
    EntrezOnlyFallback,
    Ambiguous,
    ConflictingEvidence,
    NotFound,
    UnresolvedPrimaryAnchor,
    Incomplete,
    QuarantinedCandidate,
    UnanchoredCandidate,
    InvalidInput,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EvidenceState {
    Found,
    Unmapped,
    IgnoredByPolicy,
    Incomplete,
}
#[derive(Debug, Clone)]
pub struct Assessment {
    pub input_ordinal: usize,
    pub state: EvidenceState,
    pub count: Option<u64>,
}
#[derive(Debug, Clone)]
pub struct Resolution {
    pub outcome: Outcome,
    /// Surviving candidates, NOT necessarily approved matches. Incomplete/invalid
    /// results intentionally contain no candidate prefix that could look complete.
    pub candidates: Vec<EntityId>,
    pub evidence: Vec<Assessment>,
}
impl Resolution {
    pub fn accepted_ids(&self) -> &[EntityId] {
        match self.outcome {
            Outcome::Unique
            | Outcome::MultipleSameConnectivity
            | Outcome::MultipleGeneProducts
            | Outcome::EntrezOnlyFallback => &self.candidates,
            _ => &[],
        }
    }
}

fn intersect(left: &[EntityId], right: &[EntityId]) -> Vec<EntityId> {
    let (mut a, mut b) = (0, 0);
    let mut out = Vec::with_capacity(left.len().min(right.len()));
    while a < left.len() && b < right.len() {
        match left[a].cmp(&right[b]) {
            std::cmp::Ordering::Less => a += 1,
            std::cmp::Ordering::Greater => b += 1,
            std::cmp::Ordering::Equal => {
                out.push(left[a]);
                a += 1;
                b += 1;
            }
        }
    }
    out
}

fn plan(query: &Query, lookups: &BTreeMap<Vec<u8>, Lookup>, policy: &Policy) -> Result<Resolution> {
    let mut assessments = Vec::with_capacity(query.evidence.len());
    let mut anchors = BTreeSet::new();
    let mut input_bad = false;
    let mut incomplete = false;
    let mut missing_anchor = false;
    let mut sets: Vec<&[EntityId]> = Vec::new();
    for evidence in &query.evidence {
        if evidence.key.is_empty() {
            input_bad = true;
        }
        match &evidence.role {
            Role::PrimaryAnchor(a) => {
                anchors.insert(a);
                if a.kind() != query.target {
                    input_bad = true;
                }
            }
            Role::GeneToProtein | Role::GeneIdentity { .. } if query.target != Kind::Protein => {
                input_bad = true
            }
            Role::Connectivity if query.target != Kind::Chemical => input_bad = true,
            _ => (),
        }
        if !participates(query, evidence, policy) {
            assessments.push(Assessment {
                input_ordinal: evidence.input_ordinal,
                state: EvidenceState::IgnoredByPolicy,
                count: None,
            });
            continue;
        }
        let lookup = lookups
            .get(&evidence.key)
            .ok_or_else(|| invalid("missing batch lookup result"))?;
        let (state, count) = match lookup {
            Lookup::Missing => {
                if matches!(&evidence.role, Role::PrimaryAnchor(_)) {
                    missing_anchor = true;
                }
                (EvidenceState::Unmapped, Some(0))
            }
            Lookup::Incomplete { total, .. } => {
                incomplete = true;
                (EvidenceState::Incomplete, *total)
            }
            Lookup::Complete(ids) => {
                if ids.is_empty() || ids.windows(2).any(|x| x[0] >= x[1]) {
                    return Err(invalid(
                        "posting list must be non-empty, sorted and deduplicated",
                    ));
                }
                if ids.len() > policy.max_postings_per_identifier {
                    incomplete = true;
                    (EvidenceState::Incomplete, Some(ids.len() as u64))
                } else {
                    if matches!(&evidence.role, Role::PrimaryAnchor(_)) && ids.len() != 1 {
                        return Err(invalid(
                            "a primary anchor lookup must identify exactly one entity",
                        ));
                    }
                    sets.push(ids);
                    (EvidenceState::Found, Some(ids.len() as u64))
                }
            }
        };
        assessments.push(Assessment {
            input_ordinal: evidence.input_ordinal,
            state,
            count,
        });
    }
    let compatible_claims = policy.allow_same_connectivity_multiple
        && query.target == Kind::Chemical
        && anchors
            .iter()
            .next()
            .and_then(|a| a.first_block())
            .is_some()
        && anchors
            .iter()
            .all(|a| a.first_block() == anchors.iter().next().unwrap().first_block());
    let stop = if input_bad {
        Some(Outcome::InvalidInput)
    } else if anchors.len() > 1 && !compatible_claims {
        Some(Outcome::ConflictingEvidence)
    } else if incomplete {
        Some(Outcome::Incomplete)
    } else if missing_anchor {
        Some(Outcome::UnresolvedPrimaryAnchor)
    } else if sets.is_empty() {
        Some(Outcome::NotFound)
    } else {
        None
    };
    if let Some(outcome) = stop {
        return Ok(Resolution {
            outcome,
            candidates: Vec::new(),
            evidence: assessments,
        });
    }
    // Smaller lists first. All lookups have already been checked: ordering cannot
    // hide a missing primary anchor, disabled route, or incomplete posting list.
    sets.sort_by_key(|ids| ids.len());
    let mut candidates = sets[0].to_vec();
    if query.target == Kind::Chemical && policy.allow_same_connectivity_multiple {
        // Metadata will verify the whole union has one connectivity block;
        // otherwise restore ordinary intersection below.
        candidates = sets
            .iter()
            .flat_map(|s| s.iter().copied())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
    } else {
        for ids in &sets[1..] {
            candidates = intersect(&candidates, ids);
            if candidates.is_empty() {
                break;
            }
        }
    }
    let outcome = if candidates.is_empty() {
        Outcome::ConflictingEvidence
    } else if candidates.len() > policy.max_final_candidates {
        candidates.clear();
        Outcome::Incomplete
    } else {
        Outcome::Ambiguous
    }; // finalized after one batched metadata read
    Ok(Resolution {
        outcome,
        candidates,
        evidence: assessments,
    })
}

/// Batch-deduplicate identifier lookups, intersect all informative evidence, then
/// batch-deduplicate metadata reads. Reference identifier enrichment is a THIRD
/// batched read performed by the caller, only for accepted_ids().
pub fn resolve_batch<P: PostingIndex, M: MetadataIndex>(
    postings: &P,
    metadata: &M,
    queries: &[Query],
    policy: &Policy,
) -> Result<Vec<Resolution>> {
    if policy.max_postings_per_identifier == 0
        || policy.max_final_candidates == 0
        || policy.max_total_postings_per_batch == 0
        || policy.max_total_candidates_per_batch == 0
    {
        return Err(invalid("candidate limits must be positive"));
    }
    let keys: Vec<Vec<u8>> = queries
        .iter()
        .flat_map(|q| q.evidence.iter().filter(|e| participates(q, e, policy)))
        .map(|e| e.key.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let values = postings.lookup_many(
        &keys,
        policy.max_postings_per_identifier,
        policy.max_total_postings_per_batch,
    )?;
    if values.len() != keys.len() {
        return Err(invalid("posting result count mismatch"));
    }
    let lookup_map: BTreeMap<_, _> = keys.into_iter().zip(values).collect();
    let mut plans = Vec::with_capacity(queries.len());
    let mut expanded = 0usize;
    for q in queries {
        let next = plan(q, &lookup_map, policy)?;
        expanded = expanded
            .checked_add(next.candidates.len())
            .ok_or_else(|| invalid("batch candidate count overflow"))?;
        if expanded > policy.max_total_candidates_per_batch {
            return Err(invalid(
                "batch candidate expansion limit exceeded; retry with smaller query batches",
            ));
        }
        plans.push(next);
    }
    let ids: Vec<_> = plans
        .iter()
        .flat_map(|p| p.candidates.iter().copied())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let metas = metadata.metadata_many(&ids)?;
    if ids.len() != metas.len() {
        return Err(invalid("metadata result count mismatch"));
    }
    let mut by_id = BTreeMap::new();
    for (expected, meta) in ids.into_iter().zip(metas) {
        let meta = meta.ok_or_else(|| invalid("missing entity metadata"))?;
        if meta.id != expected {
            return Err(invalid("metadata entity ID mismatch"));
        }
        if let Some(anchor) = &meta.anchor {
            if anchor.kind() != meta.kind {
                return Err(invalid("entity anchor kind mismatch"));
            }
        }
        by_id.insert(expected, meta);
    }
    for (query, resolution) in queries.iter().zip(plans.iter_mut()) {
        if resolution.candidates.is_empty() {
            continue;
        }
        if query.target == Kind::Chemical
            && policy.allow_same_connectivity_multiple
            && !query
                .evidence
                .iter()
                .any(|e| matches!(e.role, Role::PrimaryAnchor(_)))
        {
            let first = by_id[&resolution.candidates[0]]
                .anchor
                .as_ref()
                .and_then(Anchor::first_block);
            let one_block = first.is_some()
                && resolution
                    .candidates
                    .iter()
                    .all(|id| by_id[id].anchor.as_ref().and_then(Anchor::first_block) == first);
            if !one_block {
                let mut sets: Vec<&[EntityId]> = query
                    .evidence
                    .iter()
                    .filter(|e| participates(query, e, policy))
                    .filter_map(|e| match &lookup_map[&e.key] {
                        Lookup::Complete(ids) => Some(ids.as_slice()),
                        _ => None,
                    })
                    .collect();
                sets.sort_by_key(|ids| ids.len());
                let mut common = sets[0].to_vec();
                for ids in &sets[1..] {
                    common = intersect(&common, ids);
                }
                resolution.candidates = common;
                if resolution.candidates.is_empty() {
                    resolution.outcome = Outcome::ConflictingEvidence;
                    continue;
                }
            }
        }
        let mut candidates: Vec<_> = resolution.candidates.iter().map(|id| &by_id[id]).collect();
        // Only the explicit gene-product route may return a retained gene.
        // The builder guarantees this posting exists only without protein products.
        let gene_fallback = policy.allow_entrez_only_fallback
            && policy.allow_gene_to_protein
            && query.target == Kind::Protein
            && candidates.len() == 1
            && candidates[0].kind == Kind::Gene
            && query.evidence.iter().any(|e| {
                matches!(e.role, Role::GeneToProtein)
                    && matches!(lookup_map.get(&e.key), Some(Lookup::Complete(_)))
            })
            && query.evidence.iter().all(|e| {
                matches!(e.role, Role::GeneToProtein)
                    || !matches!(lookup_map.get(&e.key), Some(Lookup::Complete(_)))
            });
        if !gene_fallback && candidates.iter().any(|m| m.kind != query.target) {
            return Err(invalid("posting points to a different entity kind"));
        }
        if candidates.iter().any(|m| m.quarantined) {
            resolution.outcome = Outcome::QuarantinedCandidate;
            continue;
        }
        if policy.require_primary_anchor && candidates.iter().any(|m| m.anchor.is_none()) {
            resolution.outcome = Outcome::UnanchoredCandidate;
            continue;
        }
        let claimed: Vec<_> = query
            .evidence
            .iter()
            .filter_map(|e| match &e.role {
                Role::PrimaryAnchor(a) => Some(a),
                _ => None,
            })
            .collect();
        if !claimed.is_empty()
            && candidates
                .iter()
                .any(|m| !m.anchor.as_ref().is_some_and(|a| claimed.contains(&a)))
        {
            return Err(invalid("primary anchor index and entity metadata disagree"));
        }
        if gene_fallback {
            resolution.outcome = Outcome::EntrezOnlyFallback;
            continue;
        }
        // Soft preference AFTER intersection: per-identifier filtering can
        // manufacture conflicts when only an unreviewed product is shared.
        if query.evidence.iter().any(|e| {
            matches!(e.role, Role::GeneIdentity { .. } | Role::GeneToProtein)
                && matches!(lookup_map.get(&e.key), Some(Lookup::Complete(_)))
        }) && candidates.iter().any(|m| m.reviewed)
        {
            candidates.retain(|m| m.reviewed);
            resolution.candidates.retain(|id| by_id[id].reviewed);
        }
        if candidates.len() == 1 {
            resolution.outcome = Outcome::Unique;
            continue;
        }
        let block = candidates[0].anchor.as_ref().and_then(Anchor::first_block);
        let same_connectivity = query.target == Kind::Chemical
            && block.is_some()
            && candidates
                .iter()
                .all(|m| m.anchor.as_ref().and_then(Anchor::first_block) == block);
        if policy.allow_same_connectivity_multiple && same_connectivity {
            resolution.outcome = Outcome::MultipleSameConnectivity;
            continue;
        }
        let has_gene_evidence = query.evidence.iter().any(|e| {
            if !matches!(
                &e.role,
                Role::GeneToProtein
                    | Role::GeneIdentity {
                        multiple_products: true
                    }
            ) || !policy.allow_gene_to_protein
            {
                return false;
            }
            matches!(lookup_map.get(&e.key), Some(Lookup::Complete(_)))
        });
        let proteins = query.target == Kind::Protein
            && candidates
                .iter()
                .all(|m| matches!(&m.anchor, Some(Anchor::UniProt(_))));
        if policy.allow_gene_to_protein_multiple && has_gene_evidence && proteins {
            resolution.outcome = Outcome::MultipleGeneProducts;
        }
    }
    Ok(plans)
}

/// Test/demo backends. Do NOT load a 270-million-record reference into these maps.
#[derive(Default)]
pub struct MemoryPostings(pub BTreeMap<Vec<u8>, Lookup>);
impl PostingIndex for MemoryPostings {
    fn lookup_many(&self, keys: &[Vec<u8>], cap: usize, total_cap: usize) -> Result<Vec<Lookup>> {
        let mut remaining = total_cap;
        Ok(keys
            .iter()
            .map(|key| match self.0.get(key) {
                Some(Lookup::Complete(ids)) if ids.len() > cap.min(remaining) => {
                    Lookup::Incomplete {
                        total: Some(ids.len() as u64),
                        reason: "posting or batch decode limit".into(),
                    }
                }
                Some(Lookup::Complete(ids)) => {
                    remaining -= ids.len();
                    Lookup::Complete(ids.clone())
                }
                Some(value) => value.clone(),
                None => Lookup::Missing,
            })
            .collect())
    }
}
#[derive(Default)]
pub struct MemoryMetadata(pub BTreeMap<EntityId, EntityMeta>);
impl MetadataIndex for MemoryMetadata {
    fn metadata_many(&self, ids: &[EntityId]) -> Result<Vec<Option<EntityMeta>>> {
        Ok(ids.iter().map(|id| self.0.get(id).cloned()).collect())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RawIdentifier {
    pub namespace: String,
    pub value: String,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReferenceIdentifierRole {
    Anchor,
    IdentityAlias,
    RelatedGene,
    OtherRelationship,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceIdentifier {
    pub namespace: String,
    pub value: String,
    pub role: ReferenceIdentifierRole,
    pub provenance_ids: Vec<u64>,
}
pub trait ForwardIndex {
    /// One result per ID. Enforce total_identifier_limit BEFORE allocating an
    /// unbounded response, or return an error. Missing payloads are corruption.
    fn identifiers_many(
        &self,
        ids: &[EntityId],
        total_identifier_limit: usize,
    ) -> Result<Vec<Option<Vec<ReferenceIdentifier>>>>;
}
#[derive(Debug, Clone)]
pub struct EnrichedMatch {
    pub input_index: usize,
    pub reference_entity_id: EntityId,
    pub input_identifiers: std::sync::Arc<[RawIdentifier]>,
    pub reference_identifiers: std::sync::Arc<[ReferenceIdentifier]>,
    pub snapshot_id: std::sync::Arc<str>,
    pub policy_id: std::sync::Arc<str>,
}

/// Retain both origins without converting input identifiers into reference facts.
/// Caller MUST also retain/return all original `resolutions`, including failures.
/// Multi-matches get separate output rows, never a union of their reference aliases.
pub fn enrich_accepted<F: ForwardIndex>(
    forward: &F,
    raw_inputs: &[Vec<RawIdentifier>],
    resolutions: &[Resolution],
    snapshot_id: &str,
    policy_id: &str,
    max_identifier_occurrences: usize,
) -> Result<Vec<EnrichedMatch>> {
    use std::sync::Arc;
    if raw_inputs.len() != resolutions.len() {
        return Err(invalid("input/resolution count mismatch"));
    }
    if resolutions
        .iter()
        .zip(raw_inputs)
        .any(|(r, raw)| r.evidence.iter().any(|e| e.input_ordinal >= raw.len()))
    {
        return Err(invalid("evidence ordinal outside original input"));
    }
    let ids: Vec<_> = resolutions
        .iter()
        .flat_map(|r| r.accepted_ids().iter().copied())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let payloads = forward.identifiers_many(&ids, max_identifier_occurrences)?;
    if payloads.len() != ids.len() {
        return Err(invalid("forward result count mismatch"));
    }
    let mut by_id: BTreeMap<EntityId, Arc<[ReferenceIdentifier]>> = BTreeMap::new();
    for (id, payload) in ids.into_iter().zip(payloads) {
        by_id.insert(
            id,
            Arc::from(payload.ok_or_else(|| invalid("missing entity identifier payload"))?),
        );
    }
    let snapshot_id: Arc<str> = Arc::from(snapshot_id);
    let policy_id: Arc<str> = Arc::from(policy_id);
    let mut occurrences = 0usize;
    let mut out = Vec::new();
    for (input_index, (raw, resolution)) in raw_inputs.iter().zip(resolutions).enumerate() {
        if resolution.accepted_ids().is_empty() {
            continue;
        }
        let input_identifiers: Arc<[RawIdentifier]> = Arc::from(raw.clone());
        for id in resolution.accepted_ids() {
            let reference_identifiers = by_id[id].clone();
            let row_count = input_identifiers
                .len()
                .checked_add(reference_identifiers.len())
                .ok_or_else(|| invalid("enrichment size overflow"))?;
            occurrences = occurrences
                .checked_add(row_count)
                .ok_or_else(|| invalid("enrichment size overflow"))?;
            if occurrences > max_identifier_occurrences {
                return Err(invalid(
                    "enrichment expansion limit exceeded; stream or reduce query batch size",
                ));
            }
            out.push(EnrichedMatch {
                input_index,
                reference_entity_id: *id,
                input_identifiers: input_identifiers.clone(),
                reference_identifiers,
                snapshot_id: snapshot_id.clone(),
                policy_id: policy_id.clone(),
            });
        }
    }
    Ok(out)
}
