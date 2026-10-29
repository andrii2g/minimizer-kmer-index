//! Small-input validation independent of the experiment timing path.
use crate::{
    dna::{Generator, kmers},
    index::Index,
    minimizer, query,
};

pub fn run(seed: u64) -> Result<usize, String> {
    let mut rng = Generator::new(seed);
    let mut cases = 0;
    for _ in 0..30 {
        let reference = rng.dna(100);
        for k in [1, 3, 7, 15, 32] {
            let all = kmers(&reference, k)?;
            let full = Index::build(&all);
            for w in [1, 2, 5, 10, 40] {
                let selected = minimizer::deque(&all, w)?.selected;
                if selected != minimizer::naive(&all, w)?.selected {
                    return Err("selector disagreement".into());
                }
                let index = Index::build(&selected);
                for rate in [0.0, 0.05, 0.25, 1.0] {
                    let q = rng.mutate(&reference[20..80], rate)?;
                    let qk = kmers(&q, k)?;
                    let qs = minimizer::deque(&qk, w)?.selected;
                    let baseline = query::candidates(&full, &qk, &reference, &q);
                    let small = query::candidates(&index, &qs, &reference, &q);
                    for candidate in &small {
                        if !baseline.iter().any(|x| x.start == candidate.start) {
                            return Err(
                                "minimizer candidate missing from all-k-mer baseline".into()
                            );
                        }
                    }
                    for limit in [0, 1, 5] {
                        let exact = query::exhaustive_matches(&reference, &q, limit);
                        for c in &baseline {
                            if (c.mismatches <= limit) != exact.contains(&c.start) {
                                return Err("Hamming verification disagreement".into());
                            }
                        }
                        // Without an intact k-mer, errors separate at most (errors+1)
                        // runs of k-1 good bases: q <= (errors+1)*k - 1.
                        if q.len() >= (limit + 1) * k
                            && exact
                                .iter()
                                .any(|p| !baseline.iter().any(|x| x.start == *p))
                        {
                            return Err("all-k-mer oracle missed a guaranteed seed".into());
                        }
                    }
                    if rate == 0.0 && q.len() >= k + w - 1 && !small.iter().any(|x| x.start == 20) {
                        return Err("exact-window shared minimizer guarantee failed".into());
                    }
                    cases += 1;
                }
            }
        }
    }
    Ok(cases)
}

#[cfg(test)]
mod tests {
    #[test]
    fn runtime_verifier_passes() {
        assert_eq!(super::run(42).unwrap(), 3000);
    }
}
