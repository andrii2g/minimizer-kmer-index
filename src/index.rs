//! Sorted exact keys and packed postings. No hash-table or per-key allocation.
use crate::dna::Kmer;

#[derive(Debug)]
pub struct Index {
    keys: Vec<u64>,
    offsets: Vec<u64>,
    positions: Vec<u32>,
}

impl Index {
    pub fn build(items: &[Kmer]) -> Self {
        let mut pairs: Vec<_> = items.iter().map(|x| (x.key, x.position)).collect();
        pairs.sort_unstable();
        pairs.dedup();
        let mut keys = Vec::new();
        let mut offsets = Vec::new();
        let mut positions = Vec::with_capacity(pairs.len());
        for (key, position) in pairs {
            if keys.last() != Some(&key) {
                keys.push(key);
                offsets.push(positions.len() as u64);
            }
            positions.push(position);
        }
        offsets.push(positions.len() as u64);
        Self {
            keys,
            offsets,
            positions,
        }
    }

    pub fn lookup(&self, key: u64) -> &[u32] {
        match self.keys.binary_search(&key) {
            Ok(i) => &self.positions[self.offsets[i] as usize..self.offsets[i + 1] as usize],
            Err(_) => &[],
        }
    }

    pub fn distinct_keys(&self) -> usize {
        self.keys.len()
    }
    pub fn postings(&self) -> usize {
        self.positions.len()
    }

    /// Array payload only: excludes Vec headers, allocator metadata and build scratch.
    pub fn payload_bytes(&self) -> usize {
        8 * self.keys.len() + 8 * self.offsets.len() + 4 * self.positions.len()
    }

    /// Array capacities only, with the same exclusions as payload_bytes.
    pub fn allocated_bytes(&self) -> usize {
        8 * self.keys.capacity() + 8 * self.offsets.capacity() + 4 * self.positions.capacity()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn groups_exact_keys_even_with_colliding_hashes() {
        let items = [
            Kmer {
                key: 7,
                hash: 0,
                position: 3,
            },
            Kmer {
                key: 2,
                hash: 0,
                position: 5,
            },
            Kmer {
                key: 7,
                hash: 0,
                position: 1,
            },
            Kmer {
                key: 7,
                hash: 0,
                position: 3,
            },
        ];
        let index = Index::build(&items);
        assert_eq!(index.lookup(7), [1, 3]);
        assert_eq!(index.lookup(2), [5]);
        assert!(index.lookup(0).is_empty());
        assert_eq!(index.postings(), 3);
        assert_eq!(index.distinct_keys(), 2);
        assert_eq!(index.payload_bytes(), 52);
        assert!(index.allocated_bytes() >= index.payload_bytes());
        assert!(Index::build(&[]).lookup(7).is_empty());
    }
}
