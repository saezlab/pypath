//! Adapter for candidate values the runtime fetched for a batch.
//! Uses the shared biological decision kernel and does not read Parquet.
use anchor_resolution_serving::*;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::BTreeMap;

// Python supplies only the values retrieved for this bounded batch.
type VoteRow = (Vec<u8>, String, u8, usize, bool, bool);
type QueryRow = (String, u8, Vec<VoteRow>);
type MetaRow = (u64, String, u8, Option<String>, bool, bool);
type OutputRow = (String, String, Vec<String>, usize);

fn kind(value: u8) -> Result<Kind> {
    match value {
        1 => Ok(Kind::Chemical),
        2 => Ok(Kind::Protein),
        3 => Ok(Kind::Gene),
        _ => Err("Invalid entity kind".into()),
    }
}

#[pyfunction]
pub fn resolve_precomputed_batch(
    py: Python<'_>,
    queries: Vec<QueryRow>,
    postings: Vec<(Vec<u8>, Vec<u64>)>,
    metadata: Vec<MetaRow>,
) -> PyResult<Vec<OutputRow>> {
    py.detach(move || -> Result<Vec<OutputRow>> {
        let mut index = MemoryPostings::default();
        for (key, ids) in postings {
            if ids.windows(2).any(|w| w[0] >= w[1]) {
                return Err("Posting IDs must be strictly increasing".into());
            }
            index.0.insert(
                key,
                if ids.is_empty() {
                    Lookup::Missing
                } else {
                    Lookup::Complete(ids)
                },
            );
        }
        let mut meta = MemoryMetadata::default();
        let mut names = BTreeMap::new();
        for (id, name, k, anchor, quarantined, reviewed) in metadata {
            let k = kind(k)?;
            let anchor = match anchor {
                None => None,
                Some(a) if k == Kind::Chemical => Some(Anchor::inchikey(
                    a.strip_prefix("inchikey:").ok_or("Bad chemical anchor")?,
                )?),
                Some(a) if k == Kind::Protein => Some(Anchor::UniProt(
                    a.strip_prefix("uniprot:")
                        .ok_or("Bad protein anchor")?
                        .to_owned(),
                )),
                Some(_) => return Err("Unsupported anchored kind".into()),
            };
            names.insert(id, name);
            meta.0.insert(
                id,
                EntityMeta {
                    id,
                    kind: k,
                    anchor,
                    quarantined,
                    reviewed,
                },
            );
        }
        let mut batches = Vec::new();
        for (_, target, votes) in &queries {
            let mut evidence = Vec::new();
            for (key, anchor, route, ordinal, gene, products) in votes {
                let role = if !anchor.is_empty() {
                    Role::PrimaryAnchor(Anchor::inchikey(anchor)?)
                } else if *route == 2 {
                    Role::GeneToProtein
                } else if *gene {
                    Role::GeneIdentity {
                        multiple_products: *products,
                    }
                } else {
                    Role::Identity
                };
                evidence.push(Evidence {
                    key: key.clone(),
                    role,
                    input_ordinal: *ordinal,
                });
            }
            batches.push(Query {
                target: kind(*target)?,
                evidence,
            });
        }
        let results = resolve_batch(&index, &meta, &batches, &Policy::omnipath())?;
        queries
            .iter()
            .zip(results)
            .map(|((id, _, _), r)| {
                let accepted = r
                    .accepted_ids()
                    .iter()
                    .map(|n| {
                        names
                            .get(n)
                            .cloned()
                            .ok_or_else(|| "Missing entity metadata".into())
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok((
                    id.clone(),
                    format!("{:?}", r.outcome),
                    accepted,
                    r.candidates.len(),
                ))
            })
            .collect()
    })
    .map_err(|e| PyValueError::new_err(e.to_string()))
}

// Explicit gene links are distinct from the identifier/alias collection.
type MolecularMetaRow = (u64, String, u8, Option<String>, bool, bool, Vec<String>);
type MolecularOutputRow = (
    String,
    String,
    Vec<String>,
    usize,
    Option<String>,
    Vec<String>,
    String,
);

/// Resolve asserted products and genes independently, then check their links.
/// Gene identifiers never select a protein, even when only one product exists.
#[pyfunction]
pub fn resolve_molecular_batch(
    py: Python<'_>,
    queries: Vec<QueryRow>,
    postings: Vec<(Vec<u8>, Vec<u64>)>,
    metadata: Vec<MolecularMetaRow>,
) -> PyResult<Vec<MolecularOutputRow>> {
    py.detach(move || -> Result<Vec<MolecularOutputRow>> {
        use std::collections::BTreeSet;
        let policy = Policy::default();
        let total_postings = postings
            .iter()
            .try_fold(0usize, |n, (_, ids)| n.checked_add(ids.len()))
            .ok_or("Molecular posting count overflow")?;
        if total_postings > policy.max_total_postings_per_batch {
            return Err("Molecular batch posting limit exceeded; use smaller batches".into());
        }
        let posting_map: BTreeMap<_, _> = postings.into_iter().collect();
        let meta: BTreeMap<_, _> = metadata.into_iter().map(|m| (m.0, m)).collect();
        let mut out = Vec::new();
        let mut total_candidates = 0usize;
        for (input, target, votes) in queries {
            if target != 2 {
                return Err("Molecular queries require the gene/protein domain".into());
            }
            let mut gene_sets: Vec<BTreeSet<String>> = Vec::new();
            let mut protein_sets: Vec<BTreeSet<u64>> = Vec::new();
            let mut count = BTreeSet::new();
            let mut quarantined = false;
            let mut asserted_genes = BTreeSet::new();
            for (key, asserted_gene, route, _, gene, _) in votes {
                if asserted_gene.starts_with("entrez:") {
                    asserted_genes.insert(asserted_gene);
                }
                let ids = posting_map.get(&key).ok_or("Missing molecular posting")?;
                if ids.windows(2).any(|w| w[0] >= w[1]) {
                    return Err("Posting IDs must be strictly increasing".into());
                }
                if ids.len() > Policy::default().max_postings_per_identifier {
                    return Err("Molecular posting exceeds policy limit".into());
                }
                if ids.is_empty() {
                    continue;
                }
                let is_gene = gene || route == 2;
                let mut genes = BTreeSet::new();
                let mut proteins = BTreeSet::new();
                let mut all_linked = true;
                for id in ids {
                    let m = meta.get(id).ok_or("Missing molecular metadata")?;
                    count.insert(*id);
                    quarantined |= m.4;
                    if m.2 == 3 && m.1.starts_with("entrez:") {
                        genes.insert(m.1.clone());
                    } else {
                        all_linked &= !m.6.is_empty();
                        genes.extend(m.6.iter().map(|g| format!("entrez:{g}")));
                    }
                    if !is_gene
                        && m.2 == 2
                        && m.3
                            .as_ref()
                            .is_some_and(|a| a.starts_with("uniprot:") && !a.contains('-'))
                    {
                        proteins.insert(*id);
                    }
                }
                if !is_gene && !proteins.is_empty() {
                    protein_sets.push(proteins);
                }
                // Missing links cannot establish the gene of ambiguous products.
                if all_linked && !genes.is_empty() {
                    gene_sets.push(genes);
                }
            }
            total_candidates = total_candidates
                .checked_add(count.len())
                .ok_or("Molecular candidate count overflow")?;
            if total_candidates > policy.max_total_candidates_per_batch {
                return Err("Molecular batch candidate limit exceeded; use smaller batches".into());
            }
            let combine_genes = |sets: &[BTreeSet<String>]| -> BTreeSet<String> {
                sets.iter()
                    .skip(1)
                    .fold(sets.first().cloned().unwrap_or_default(), |a, b| {
                        a.intersection(b).cloned().collect()
                    })
            };
            // Source-supplied GeneIDs also constrain a known product mapping,
            // including stale/missing catalogue genes. They cannot create a
            // matched reference identity by themselves when the catalogue misses.
            if !asserted_genes.is_empty() && !gene_sets.is_empty() {
                for gene in asserted_genes {
                    gene_sets.push(BTreeSet::from([gene]));
                }
            }
            let genes = combine_genes(&gene_sets);
            let proteins = protein_sets
                .iter()
                .skip(1)
                .fold(protein_sets.first().cloned().unwrap_or_default(), |a, b| {
                    a.intersection(b).copied().collect()
                });
            let product_conflict = !protein_sets.is_empty() && proteins.is_empty();
            let gene_conflict = !gene_sets.is_empty() && genes.is_empty();
            let protein = if proteins.len() == 1 && !product_conflict && !quarantined {
                Some(meta[proteins.first().unwrap()].1.clone())
            } else {
                None
            };
            let status = if gene_conflict || product_conflict {
                "conflict"
            } else if genes.len() == 1 {
                "resolved"
            } else if genes.len() > 1 {
                "ambiguous"
            } else {
                "missing"
            };
            let mut candidates: BTreeSet<String> =
                gene_sets.iter().flat_map(|s| s.iter().cloned()).collect();
            if !gene_conflict {
                candidates = genes.clone();
            }
            let accepted = if !quarantined && !product_conflict && status == "resolved" {
                genes.into_iter().collect()
            } else if let Some(ref p) = protein {
                vec![p.clone()]
            } else {
                Vec::new()
            };
            let outcome = if quarantined {
                "QuarantinedCandidate"
            } else if gene_conflict || product_conflict {
                "ConflictingEvidence"
            } else if accepted.is_empty() {
                "Ambiguous"
            } else {
                "Unique"
            };
            out.push((
                input,
                outcome.into(),
                accepted,
                count.len(),
                protein,
                candidates.into_iter().collect(),
                status.into(),
            ));
        }
        Ok(out)
    })
    .map_err(|e| PyValueError::new_err(e.to_string()))
}
