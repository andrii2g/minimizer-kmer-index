//! Substitution-only candidate voting and independent exhaustive verification.
use crate::{dna::Kmer, index::Index};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Candidate {
    pub start: u32,
    pub votes: usize,
    pub mismatches: usize,
}

/// Return all in-bounds seed-supported starts, ordered by votes then position.
/// Repeated input seeds cannot inflate the number of supporting query positions.
pub fn candidates(index: &Index, seeds: &[Kmer], reference: &[u8], query: &[u8]) -> Vec<Candidate> {
    if query.is_empty() || query.len() > reference.len() {
        return Vec::new();
    }
    let mut votes: BTreeMap<u32, BTreeSet<u32>> = BTreeMap::new();
    for seed in seeds {
        for &position in index.lookup(seed.key) {
            if let Some(start) = position.checked_sub(seed.position)
                && start as usize <= reference.len() - query.len()
            {
                votes.entry(start).or_default().insert(seed.position);
            }
        }
    }
    let mut result: Vec<_> = votes
        .into_iter()
        .map(|(start, positions)| Candidate {
            start,
            votes: positions.len(),
            mismatches: hamming(
                &reference[start as usize..start as usize + query.len()],
                query,
            )
            .unwrap(),
        })
        .collect();
    result.sort_unstable_by_key(|c| (std::cmp::Reverse(c.votes), c.start));
    result
}

pub fn hamming(a: &[u8], b: &[u8]) -> Option<usize> {
    (a.len() == b.len()).then(|| a.iter().zip(b).filter(|(x, y)| x != y).count())
}

/// Ground-truth oracle, independent of fingerprints. Empty queries have no matches.
pub fn exhaustive_matches(reference: &[u8], query: &[u8], max_mismatches: usize) -> Vec<u32> {
    if query.is_empty() {
        return Vec::new();
    }
    reference
        .windows(query.len())
        .enumerate()
        .filter(|(_, window)| hamming(window, query).unwrap() <= max_mismatches)
        .map(|(i, _)| i as u32)
        .collect()
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Recall {
    pub shared_selected: usize,
    pub surviving: usize,
}

impl Recall {
    pub fn fraction(self) -> Option<f64> {
        (self.surviving > 0).then(|| self.shared_selected as f64 / self.surviving as f64)
    }
}

/// Positional recall at the known query origin, not accidental repeated-key hits.
pub fn seed_recall(
    reference: &[Kmer],
    reference_selected: &[Kmer],
    query: &[Kmer],
    query_selected: &[Kmer],
    origin: u32,
) -> Recall {
    let all: BTreeMap<_, _> = reference.iter().map(|x| (x.position, x.key)).collect();
    let selected: BTreeSet<_> = reference_selected
        .iter()
        .map(|x| (x.position, x.key))
        .collect();
    let query_positions: BTreeSet<_> = query_selected.iter().map(|x| x.position).collect();
    let mut result = Recall::default();
    for item in query {
        if let Some(position) = origin.checked_add(item.position)
            && all.get(&position) == Some(&item.key)
        {
            result.surviving += 1;
            if query_positions.contains(&item.position) && selected.contains(&(position, item.key))
            {
                result.shared_selected += 1;
            }
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        dna::{Generator, kmers},
        minimizer::deque,
    };

    #[test]
    fn candidates_and_recall_agree_with_independent_oracles() {
        let mut rng = Generator::new(199);
        for _ in 0..40 {
            let reference = rng.dna(80);
            let all = kmers(&reference, 4).unwrap();
            let selected = deque(&all, 5).unwrap().selected;
            let all_index = Index::build(&all);
            let small_index = Index::build(&selected);
            for rate in [0.0, 0.1, 0.5, 1.0] {
                let query = rng.mutate(&reference[17..47], rate).unwrap();
                let q = kmers(&query, 4).unwrap();
                let qs = deque(&q, 5).unwrap().selected;
                let full = candidates(&all_index, &q, &reference, &query);
                let small = candidates(&small_index, &qs, &reference, &query);
                for c in &small {
                    assert!(full.iter().any(|x| x.start == c.start));
                    assert!(
                        exhaustive_matches(&reference, &query, c.mismatches).contains(&c.start)
                    );
                    let expected_votes = qs
                        .iter()
                        .filter(|x| {
                            selected
                                .iter()
                                .any(|r| r.position == c.start + x.position && r.key == x.key)
                        })
                        .count();
                    assert_eq!(c.votes, expected_votes);
                }
                // At this threshold, every valid match must have an exact k-mer.
                let verified: BTreeSet<_> = full
                    .iter()
                    .filter(|c| c.mismatches <= 2)
                    .map(|c| c.start)
                    .collect();
                assert_eq!(
                    verified,
                    exhaustive_matches(&reference, &query, 2)
                        .into_iter()
                        .collect()
                );
                let recall = seed_recall(&all, &selected, &q, &qs, 17);
                let surviving = (0..=26)
                    .filter(|&i| reference[17 + i..17 + i + 4] == query[i..i + 4])
                    .count();
                assert_eq!(recall.surviving, surviving);
                assert!(recall.shared_selected <= surviving);
                if rate == 0.0 {
                    assert!(small.iter().any(|x| x.start == 17));
                }
            }
        }
    }

    #[test]
    fn repeated_seeds_boundaries_and_no_survivors() {
        let reference = b"AAAAAAAAAA";
        let all = kmers(reference, 2).unwrap();
        let index = Index::build(&all);
        let mut q = kmers(b"AAAA", 2).unwrap();
        q.extend(q.clone());
        let result = candidates(&index, &q, reference, b"AAAA");
        assert_eq!(result.len(), 7);
        assert!(result.iter().all(|c| c.votes == 3 && c.mismatches == 0));
        assert!(candidates(&index, &q, reference, b"").is_empty());
        assert!(exhaustive_matches(reference, b"AAAAAAAAAAA", 2).is_empty());
        assert_eq!(hamming(b"A", b""), None);
        assert_eq!(Recall::default().fraction(), None);
    }
}
