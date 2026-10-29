//! Windows contain w k-mers. Equal ranks retain the leftmost occurrence.
use crate::dna::Kmer;
use std::collections::VecDeque;

#[derive(Clone, Copy, Debug, Default)]
pub struct Work {
    pub comparisons: u64,
    pub pushes: u64,
    pub pops: u64,
}

#[derive(Debug)]
pub struct Selection {
    pub selected: Vec<Kmer>,
    pub work: Work,
}

fn rank(item: Kmer) -> (u64, u64) {
    (item.hash, item.key)
}

fn emit(selected: &mut Vec<Kmer>, item: Kmer) {
    if selected
        .last()
        .is_none_or(|last| last.position != item.position)
    {
        selected.push(item);
    }
}

/// Fingerprints must be ordered by increasing, distinct sequence position.
pub fn naive(items: &[Kmer], w: usize) -> Result<Selection, String> {
    if w == 0 {
        return Err("w must be positive".into());
    }
    let mut result = Selection {
        selected: Vec::new(),
        work: Work::default(),
    };
    for window in items.windows(w) {
        let mut best = window[0];
        for &item in &window[1..] {
            result.work.comparisons += 1;
            if rank(item) < rank(best) {
                best = item;
            }
        }
        emit(&mut result.selected, best);
    }
    Ok(result)
}

/// Each entry is pushed and removed at most once; auxiliary space is O(w).
pub fn deque(items: &[Kmer], w: usize) -> Result<Selection, String> {
    if w == 0 {
        return Err("w must be positive".into());
    }
    let mut result = Selection {
        selected: Vec::new(),
        work: Work::default(),
    };
    let mut queue: VecDeque<usize> = VecDeque::new();
    for (i, &item) in items.iter().enumerate() {
        while queue.front().is_some_and(|&j| i >= w && j <= i - w) {
            queue.pop_front();
            result.work.pops += 1;
        }
        while let Some(&j) = queue.back() {
            result.work.comparisons += 1;
            if rank(items[j]) <= rank(item) {
                break;
            }
            queue.pop_back();
            result.work.pops += 1;
        }
        queue.push_back(i);
        result.work.pushes += 1;
        if i + 1 >= w {
            emit(&mut result.selected, items[*queue.front().unwrap()]);
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::kmers;

    #[test]
    fn exhaustive_small_dna_matches_window_oracle() {
        for length in 0..=7 {
            for code in 0..(1_usize << (2 * length)) {
                let dna: Vec<_> = (0..length)
                    .map(|i| b"ACGT"[(code >> (2 * i)) & 3])
                    .collect();
                for k in 1..=3 {
                    let items = kmers(&dna, k).unwrap();
                    for w in 1..=5 {
                        let a = naive(&items, w).unwrap();
                        let b = deque(&items, w).unwrap();
                        assert_eq!(a.selected, b.selected);
                        assert!(b.selected.windows(2).all(|p| p[0].position < p[1].position));
                        for window in items.windows(w) {
                            let expected = window
                                .iter()
                                .min_by_key(|x| (x.hash, x.key, x.position))
                                .unwrap();
                            assert!(b.selected.contains(expected));
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn ties_collisions_and_invalid_windows() {
        let items: Vec<_> = [2, 1, 1, 3]
            .into_iter()
            .enumerate()
            .map(|(i, key)| Kmer {
                key,
                hash: 0,
                position: i as u32,
            })
            .collect();
        assert_eq!(
            deque(&items, 2)
                .unwrap()
                .selected
                .iter()
                .map(|x| x.position)
                .collect::<Vec<_>>(),
            [1, 2]
        );
        assert_eq!(
            naive(&items, 2).unwrap().selected,
            deque(&items, 2).unwrap().selected
        );
        assert!(deque(&items, 0).is_err());
        assert!(naive(&items, 0).is_err());
        assert_eq!(deque(&items, 1).unwrap().selected, items);
    }
}
