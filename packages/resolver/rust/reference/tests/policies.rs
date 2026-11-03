use anchor_resolution_serving::*;

fn chem(block: char, suffix: char) -> Anchor {
    Anchor::inchikey(&format!(
        "{}-{}SA-N",
        block.to_string().repeat(14),
        suffix.to_string().repeat(8)
    ))
    .unwrap()
}
fn fixture() -> (MemoryPostings, MemoryMetadata) {
    let mut p = MemoryPostings::default();
    for (key, ids) in [
        ("a", vec![1, 2]),
        ("b", vec![2, 3]),
        ("one", vec![1]),
        ("two", vec![2]),
        ("three", vec![3]),
        ("gene", vec![4, 5]),
        ("protein", vec![5]),
        ("protein_alias", vec![4, 5]),
    ] {
        p.0.insert(key.as_bytes().to_vec(), Lookup::Complete(ids));
    }
    let mut m = MemoryMetadata::default();
    for (id, anchor) in [
        (1, chem('A', 'B')),
        (2, chem('A', 'C')),
        (3, chem('D', 'C')),
        (4, Anchor::UniProt("P00001".into())),
        (5, Anchor::UniProt("P00002".into())),
    ] {
        m.0.insert(
            id,
            EntityMeta {
                id,
                kind: anchor.kind(),
                anchor: Some(anchor),
                quarantined: false,
                reviewed: false,
            },
        );
    }
    (p, m)
}
fn query(target: Kind, values: Vec<(&str, Role)>) -> Query {
    Query {
        target,
        evidence: values
            .into_iter()
            .enumerate()
            .map(|(i, (s, role))| Evidence {
                input_ordinal: i,
                key: s.as_bytes().to_vec(),
                role,
            })
            .collect(),
    }
}
fn run(target: Kind, values: Vec<(&str, Role)>, policy: Policy) -> Resolution {
    let (p, m) = fixture();
    resolve_batch(&p, &m, &[query(target, values)], &policy)
        .unwrap()
        .remove(0)
}
#[test]
fn unique_identity() {
    let r = run(
        Kind::Chemical,
        vec![("one", Role::Identity)],
        Policy::default(),
    );
    assert_eq!(r.outcome, Outcome::Unique);
    assert_eq!(r.accepted_ids(), &[1]);
}
#[test]
fn intersect_not_union() {
    let r = run(
        Kind::Chemical,
        vec![("a", Role::Identity), ("b", Role::Identity)],
        Policy::default(),
    );
    assert_eq!(r.accepted_ids(), &[2]);
}
#[test]
fn order_independent() {
    let a = run(
        Kind::Chemical,
        vec![("a", Role::Identity), ("b", Role::Identity)],
        Policy::default(),
    );
    let b = run(
        Kind::Chemical,
        vec![("b", Role::Identity), ("a", Role::Identity)],
        Policy::default(),
    );
    assert_eq!(a.outcome, b.outcome);
    assert_eq!(a.candidates, b.candidates);
}
#[test]
fn disjoint_same_connectivity_sets_project_when_enabled() {
    let policy = Policy {
        allow_same_connectivity_multiple: true,
        ..Policy::default()
    };
    let r = run(
        Kind::Chemical,
        vec![("one", Role::Identity), ("two", Role::Identity)],
        policy,
    );
    assert_eq!(r.outcome, Outcome::MultipleSameConnectivity);
    assert_eq!(r.accepted_ids(), &[1, 2]);
}
#[test]
fn unknown_alias_preserved_as_unmapped() {
    let r = run(
        Kind::Chemical,
        vec![("one", Role::Identity), ("missing", Role::Identity)],
        Policy::default(),
    );
    assert_eq!(r.accepted_ids(), &[1]);
    assert_eq!(r.evidence[1].state, EvidenceState::Unmapped);
}
#[test]
fn unknown_primary_anchor_blocks_match() {
    let r = run(
        Kind::Chemical,
        vec![
            ("one", Role::Identity),
            ("missing", Role::PrimaryAnchor(chem('A', 'B'))),
        ],
        Policy::default(),
    );
    assert_eq!(r.outcome, Outcome::UnresolvedPrimaryAnchor);
    assert!(r.accepted_ids().is_empty());
}
#[test]
fn not_found_is_distinct_from_conflict() {
    let r = run(
        Kind::Chemical,
        vec![("missing", Role::Identity)],
        Policy::default(),
    );
    assert_eq!(r.outcome, Outcome::NotFound);
}
#[test]
fn same_connectivity_multi_is_explicit() {
    let strict = run(
        Kind::Chemical,
        vec![("a", Role::Identity)],
        Policy::default(),
    );
    assert_eq!(strict.outcome, Outcome::Ambiguous);
    let r = run(
        Kind::Chemical,
        vec![("a", Role::Identity)],
        Policy {
            allow_same_connectivity_multiple: true,
            ..Policy::default()
        },
    );
    assert_eq!(r.outcome, Outcome::MultipleSameConnectivity);
    assert_eq!(r.accepted_ids(), &[1, 2]);
}
#[test]
fn different_connectivity_not_approved() {
    let r = run(
        Kind::Chemical,
        vec![("b", Role::Identity)],
        Policy {
            allow_same_connectivity_multiple: true,
            ..Policy::default()
        },
    );
    assert_eq!(r.outcome, Outcome::Ambiguous);
}
#[test]
fn explicit_anchor_does_not_expand_to_siblings() {
    let r = run(
        Kind::Chemical,
        vec![
            ("a", Role::Identity),
            ("one", Role::PrimaryAnchor(chem('A', 'B'))),
        ],
        Policy {
            allow_same_connectivity_multiple: true,
            ..Policy::default()
        },
    );
    assert_eq!(r.accepted_ids(), &[1]);
}
#[test]
fn compatible_full_anchor_claims_project_when_enabled() {
    let r = run(
        Kind::Chemical,
        vec![
            ("one", Role::PrimaryAnchor(chem('A', 'B'))),
            ("two", Role::PrimaryAnchor(chem('A', 'C'))),
        ],
        Policy {
            allow_same_connectivity_multiple: true,
            ..Policy::default()
        },
    );
    assert_eq!(r.outcome, Outcome::MultipleSameConnectivity);
    assert_eq!(r.accepted_ids(), &[1, 2]);
}
#[test]
fn gene_projection_requires_policy() {
    let r = run(
        Kind::Protein,
        vec![("gene", Role::GeneToProtein)],
        Policy::default(),
    );
    assert_eq!(r.outcome, Outcome::NotFound);
    assert_eq!(r.evidence[0].state, EvidenceState::IgnoredByPolicy);
}
#[test]
fn gene_products_multi_allowed() {
    let r = run(
        Kind::Protein,
        vec![("gene", Role::GeneToProtein)],
        Policy {
            allow_gene_to_protein: true,
            allow_gene_to_protein_multiple: true,
            ..Policy::default()
        },
    );
    assert_eq!(r.outcome, Outcome::MultipleGeneProducts);
    assert_eq!(r.accepted_ids(), &[4, 5]);
}
#[test]
fn gene_plus_protein_can_be_unique() {
    let r = run(
        Kind::Protein,
        vec![("gene", Role::GeneToProtein), ("protein", Role::Identity)],
        Policy {
            allow_gene_to_protein: true,
            allow_gene_to_protein_multiple: true,
            ..Policy::default()
        },
    );
    assert_eq!(r.accepted_ids(), &[5]);
}
#[test]
fn protein_alias_ambiguity_is_not_gene_projection() {
    let r = run(
        Kind::Protein,
        vec![("protein_alias", Role::Identity)],
        Policy {
            allow_gene_to_protein: true,
            allow_gene_to_protein_multiple: true,
            ..Policy::default()
        },
    );
    assert_eq!(r.outcome, Outcome::Ambiguous);
}
#[test]
fn capped_posting_never_becomes_unique() {
    let r = run(
        Kind::Chemical,
        vec![("a", Role::Identity), ("one", Role::Identity)],
        Policy {
            max_postings_per_identifier: 1,
            ..Policy::default()
        },
    );
    assert_eq!(r.outcome, Outcome::Incomplete);
    assert!(r.accepted_ids().is_empty());
}
#[test]
fn final_candidate_cap_is_after_intersection() {
    let r = run(
        Kind::Chemical,
        vec![("a", Role::Identity), ("b", Role::Identity)],
        Policy {
            max_final_candidates: 1,
            ..Policy::default()
        },
    );
    assert_eq!(r.accepted_ids(), &[2]);
}
#[test]
fn quarantined_candidate_is_not_success() {
    let (p, mut m) = fixture();
    m.0.get_mut(&1).unwrap().quarantined = true;
    let r = resolve_batch(
        &p,
        &m,
        &[query(Kind::Chemical, vec![("one", Role::Identity)])],
        &Policy::default(),
    )
    .unwrap();
    assert_eq!(r[0].outcome, Outcome::QuarantinedCandidate);
}
#[test]
fn metadata_absence_is_an_error() {
    let (p, mut m) = fixture();
    m.0.remove(&1);
    assert!(resolve_batch(
        &p,
        &m,
        &[query(Kind::Chemical, vec![("one", Role::Identity)])],
        &Policy::default()
    )
    .is_err());
}
#[test]
fn routes_and_scope_do_not_collide() {
    assert_ne!(
        encode_key(Kind::Protein, &Role::Identity, 1, None, b"123"),
        encode_key(Kind::Protein, &Role::GeneToProtein, 1, None, b"123")
    );
    assert_ne!(
        encode_key(Kind::Protein, &Role::Identity, 1, None, b"123"),
        encode_key(Kind::Protein, &Role::Identity, 1, Some(9606), b"123")
    );
}
#[test]
fn enrichment_keeps_input_and_reference_origins_separate() {
    struct Forward;
    impl ForwardIndex for Forward {
        fn identifiers_many(
            &self,
            ids: &[EntityId],
            _: usize,
        ) -> Result<Vec<Option<Vec<ReferenceIdentifier>>>> {
            Ok(ids
                .iter()
                .map(|id| {
                    Some(vec![ReferenceIdentifier {
                        namespace: "reference_db".into(),
                        value: id.to_string(),
                        role: ReferenceIdentifierRole::IdentityAlias,
                        provenance_ids: vec![42],
                    }])
                })
                .collect())
        }
    }
    let r = run(
        Kind::Chemical,
        vec![("a", Role::Identity)],
        Policy {
            allow_same_connectivity_multiple: true,
            ..Policy::default()
        },
    );
    let raw = vec![vec![RawIdentifier {
        namespace: "input_db".into(),
        value: "a".into(),
    }]];
    let out = enrich_accepted(&Forward, &raw, &[r], "snapshot-1", "policy-1", 10).unwrap();
    assert_eq!(out.len(), 2);
    assert_eq!(out[0].input_identifiers, out[1].input_identifiers);
    assert_ne!(out[0].reference_identifiers, out[1].reference_identifiers);
    assert_eq!(out[0].reference_identifiers[0].provenance_ids, vec![42]);
}

#[test]
fn supplied_structure_overrides_conflicting_cross_reference() {
    let r = run(
        Kind::Chemical,
        vec![
            ("one", Role::PrimaryAnchor(chem('A', 'B'))),
            ("three", Role::Identity),
        ],
        legacy_projection_policy(),
    );
    assert_eq!(r.outcome, Outcome::Unique);
    assert_eq!(r.accepted_ids(), &[1]);
    assert_eq!(r.evidence[1].state, EvidenceState::IgnoredByPolicy);
}
#[test]
fn different_supplied_connectivity_remains_conflicting() {
    let r = run(
        Kind::Chemical,
        vec![
            ("one", Role::PrimaryAnchor(chem('A', 'B'))),
            ("three", Role::PrimaryAnchor(chem('D', 'C'))),
        ],
        legacy_projection_policy(),
    );
    assert_eq!(r.outcome, Outcome::ConflictingEvidence);
}

#[test]
fn entrez_only_fallback_is_explicit_and_route_restricted() {
    let (mut p, mut m) = fixture();
    p.0.insert(b"orphan".to_vec(), Lookup::Complete(vec![6]));
    m.0.insert(
        6,
        EntityMeta {
            id: 6,
            kind: Kind::Gene,
            anchor: None,
            quarantined: false,
            reviewed: false,
        },
    );
    let q = query(Kind::Protein, vec![("orphan", Role::GeneToProtein)]);
    let r = resolve_batch(
        &p,
        &m,
        std::slice::from_ref(&q),
        &legacy_projection_policy(),
    )
    .unwrap()
    .remove(0);
    assert_eq!(r.outcome, Outcome::EntrezOnlyFallback);
    assert_eq!(r.accepted_ids(), &[6]);
    let strict = Policy {
        allow_entrez_only_fallback: false,
        ..legacy_projection_policy()
    };
    assert!(resolve_batch(&p, &m, std::slice::from_ref(&q), &strict).is_err());
    let wrong_route = query(Kind::Protein, vec![("orphan", Role::Identity)]);
    assert!(resolve_batch(&p, &m, &[wrong_route], &legacy_projection_policy()).is_err());
    let conflicting = query(
        Kind::Protein,
        vec![("orphan", Role::GeneToProtein), ("protein", Role::Identity)],
    );
    let r = resolve_batch(&p, &m, &[conflicting], &legacy_projection_policy())
        .unwrap()
        .remove(0);
    assert_eq!(r.outcome, Outcome::ConflictingEvidence);
    m.0.get_mut(&6).unwrap().quarantined = true;
    let r = resolve_batch(&p, &m, &[q], &legacy_projection_policy())
        .unwrap()
        .remove(0);
    assert_eq!(r.outcome, Outcome::QuarantinedCandidate);
    assert!(r.accepted_ids().is_empty());
}

#[test]
fn reviewed_gene_preference_is_after_intersection() {
    let (mut p, mut m) = fixture();
    m.0.get_mut(&4).unwrap().reviewed = true;
    p.0.insert(b"only_unreviewed".to_vec(), Lookup::Complete(vec![5]));
    let qs = vec![
        query(
            Kind::Protein,
            vec![(
                "gene",
                Role::GeneIdentity {
                    multiple_products: false,
                },
            )],
        ),
        query(
            Kind::Protein,
            vec![
                (
                    "gene",
                    Role::GeneIdentity {
                        multiple_products: false,
                    },
                ),
                (
                    "only_unreviewed",
                    Role::GeneIdentity {
                        multiple_products: false,
                    },
                ),
            ],
        ),
        query(Kind::Protein, vec![("gene", Role::Identity)]),
        query(
            Kind::Protein,
            vec![
                ("gene", Role::Identity),
                (
                    "missing",
                    Role::GeneIdentity {
                        multiple_products: false,
                    },
                ),
            ],
        ),
    ];
    let r = resolve_batch(&p, &m, &qs, &legacy_projection_policy()).unwrap();
    assert_eq!(r[0].accepted_ids(), &[4]);
    assert_eq!(r[1].accepted_ids(), &[5]); // agreement must not be destroyed by preference
    assert_eq!(r[2].outcome, Outcome::Ambiguous); // explicit protein aliases unchanged
    assert_eq!(r[3].outcome, Outcome::Ambiguous); // missing gene evidence cannot rank proteins
    m.0.get_mut(&5).unwrap().reviewed = true;
    let r = resolve_batch(&p, &m, &qs[..1], &legacy_projection_policy()).unwrap();
    assert_eq!(r[0].outcome, Outcome::Ambiguous); // two reviewed proteins are not arbitrarily collapsed
}

#[test]
fn stable_gene_id_can_project_multiple_reviewed_products() {
    let (p, mut m) = fixture();
    m.0.get_mut(&4).unwrap().reviewed = true;
    m.0.get_mut(&5).unwrap().reviewed = true;
    let q = query(
        Kind::Protein,
        vec![(
            "gene",
            Role::GeneIdentity {
                multiple_products: true,
            },
        )],
    );
    let r = resolve_batch(&p, &m, &[q], &legacy_projection_policy()).unwrap();
    assert_eq!(r[0].outcome, Outcome::MultipleGeneProducts);
    assert_eq!(r[0].accepted_ids(), &[4, 5]);
}

#[test]
fn omnipath_requires_unique_chemical_even_when_connectivity_matches() {
    let r = run(
        Kind::Chemical,
        vec![("a", Role::Identity)],
        legacy_projection_policy(),
    );
    assert_eq!(r.outcome, Outcome::Ambiguous);
}

#[test]
fn proven_symbol_gene_preserves_explicit_protein_intersection() {
    let (mut p, mut m) = fixture();
    m.0.get_mut(&4).unwrap().reviewed = true;
    p.0.insert(b"explicit_unreviewed".to_vec(), Lookup::Complete(vec![5]));
    let q = query(
        Kind::Protein,
        vec![
            (
                "gene",
                Role::GeneIdentity {
                    multiple_products: true,
                },
            ),
            ("explicit_unreviewed", Role::Identity),
        ],
    );
    let r = resolve_batch(&p, &m, &[q], &legacy_projection_policy()).unwrap();
    assert_eq!(r[0].outcome, Outcome::Unique);
    assert_eq!(r[0].accepted_ids(), &[5]);
}

fn legacy_projection_policy() -> Policy {
    Policy {
        allow_gene_to_protein: true,
        allow_entrez_only_fallback: true,
        allow_gene_to_protein_multiple: true,
        ..Policy::omnipath()
    }
}

#[test]
fn omnipath_never_asserts_products_from_gene_roles() {
    let (p, m) = fixture();
    let qs = [
        query(Kind::Protein, vec![("gene", Role::GeneToProtein)]),
        query(
            Kind::Protein,
            vec![(
                "gene",
                Role::GeneIdentity {
                    multiple_products: true,
                },
            )],
        ),
    ];
    for r in resolve_batch(&p, &m, &qs, &Policy::omnipath()).unwrap() {
        assert_eq!(r.outcome, Outcome::NotFound);
        assert!(r.accepted_ids().is_empty());
    }
}
