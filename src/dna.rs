//! Two-bit k-mers and reproducible synthetic data (not a cryptographic RNG).

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Kmer {
    pub key: u64,
    pub hash: u64,
    pub position: u32,
}

/// SplitMix64 finalizer. Its odd multipliers and xor-shifts form a permutation.
pub fn mix64(mut value: u64) -> u64 {
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d049bb133111eb);
    value ^ (value >> 31)
}

/// Encode A=0, C=1, G=2, T=3. All input is validated, even if shorter than k.
pub fn kmers(sequence: &[u8], k: usize) -> Result<Vec<Kmer>, String> {
    if !(1..=32).contains(&k) {
        return Err("k must be in 1..=32".into());
    }
    if sequence.len() > u32::MAX as usize {
        return Err("sequence length exceeds u32 position limit".into());
    }
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1_u64 << (2 * k)) - 1
    };
    let mut key = 0_u64;
    let mut result = Vec::with_capacity(sequence.len().saturating_sub(k - 1));
    for (i, base) in sequence.iter().enumerate() {
        let bits = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                return Err(format!(
                    "invalid DNA base at position {i}; expected A/C/G/T"
                ));
            }
        };
        key = ((key << 2) | bits) & mask;
        if i + 1 >= k {
            result.push(Kmer {
                key,
                hash: mix64(key),
                position: (i + 1 - k) as u32,
            });
        }
    }
    Ok(result)
}

pub struct Generator(u64);

impl Generator {
    pub fn new(seed: u64) -> Self {
        Self(seed)
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        mix64(self.0)
    }

    pub fn dna(&mut self, length: usize) -> Vec<u8> {
        (0..length)
            .map(|_| b"ACGT"[(self.next_u64() & 3) as usize])
            .collect()
    }

    /// Each chosen replacement differs from the source base.
    pub fn mutate(&mut self, sequence: &[u8], rate: f64) -> Result<Vec<u8>, String> {
        if !rate.is_finite() || !(0.0..=1.0).contains(&rate) {
            return Err("mutation rate must be finite and in 0..=1".into());
        }
        sequence
            .iter()
            .map(|base| {
                let original = b"ACGT"
                    .iter()
                    .position(|b| b == base)
                    .ok_or_else(|| "mutation input must contain only A/C/G/T".to_string())?;
                let probability = (self.next_u64() >> 11) as f64 / ((1_u64 << 53) as f64);
                Ok(if probability < rate {
                    b"ACGT"[(original + 1 + (self.next_u64() % 3) as usize) % 4]
                } else {
                    *base
                })
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rolling_encoding_matches_independent_slices() {
        let dna = Generator::new(42).dna(100);
        for k in 1..=32 {
            for item in kmers(&dna, k).unwrap() {
                let expected = dna[item.position as usize..item.position as usize + k]
                    .iter()
                    .fold(0_u64, |v, b| {
                        (v << 2) | b"ACGT".iter().position(|x| x == b).unwrap() as u64
                    });
                assert_eq!(item.key, expected);
            }
        }
        assert_eq!(kmers(&[b'T'; 32], 32).unwrap()[0].key, u64::MAX);
    }

    #[test]
    fn validates_boundaries_and_mutations() {
        assert!(kmers(b"AN", 3).is_err());
        assert!(kmers(b"A", 0).is_err());
        assert!(kmers(b"A", 33).is_err());
        assert!(kmers(b"AC", 3).unwrap().is_empty());
        let input = b"ACGTACGT";
        assert_eq!(Generator::new(7).mutate(input, 0.0).unwrap(), input);
        assert!(
            Generator::new(7)
                .mutate(input, 1.0)
                .unwrap()
                .iter()
                .zip(input)
                .all(|(a, b)| a != b)
        );
        assert!(Generator::new(7).mutate(input, f64::NAN).is_err());
        assert_eq!(Generator::new(7).dna(100), Generator::new(7).dna(100));
    }
}
