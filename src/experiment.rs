//! Reproducible experiments; timings exclude data generation and oracle work.
use crate::{
    dna::{Generator, kmers},
    index::Index,
    minimizer::{self, Selection},
    query,
};
use std::{hint::black_box, time::Instant};

#[derive(Clone, Debug)]
pub struct Config {
    pub length: usize,
    pub query_length: usize,
    pub queries: usize,
    pub k: usize,
    pub w: usize,
    pub mutation_rate: f64,
    pub seed: u64,
    pub repeats: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            length: 100_000,
            query_length: 200,
            queries: 100,
            k: 15,
            w: 10,
            mutation_rate: 0.02,
            seed: 42,
            repeats: 3,
        }
    }
}

impl Config {
    pub fn validate(&self) -> Result<(), String> {
        if !(1..=32).contains(&self.k) {
            return Err("k must be in 1..=32".into());
        }
        if self.w == 0 {
            return Err("w must be positive".into());
        }
        if self.length > u32::MAX as usize {
            return Err("reference exceeds u32 position limit".into());
        }
        if self.query_length < self.k || self.query_length > self.length {
            return Err("require k <= query-length <= length".into());
        }
        if self.queries == 0 || self.repeats == 0 {
            return Err("queries and repeats must be positive".into());
        }
        if !self.mutation_rate.is_finite() || !(0.0..=1.0).contains(&self.mutation_rate) {
            return Err("mutation rate must be finite and in 0..=1".into());
        }
        Ok(())
    }
}

#[derive(Debug)]
pub struct Metrics {
    pub selected: usize,
    pub kmers: usize,
    pub distinct_keys: usize,
    pub all_distinct_keys: usize,
    pub bytes: usize,
    pub all_bytes: usize,
    pub allocated: usize,
    pub all_allocated: usize,
    pub naive_comparisons: u64,
    pub deque_comparisons: u64,
    pub pushes: u64,
    pub pops: u64,
    pub query_hashes: usize,
    pub query_selected: usize,
    pub candidates: usize,
    pub all_candidates: usize,
    pub mapped: usize,
    pub all_mapped: usize,
    pub shared: usize,
    pub surviving: usize,
    pub extract_ms: f64,
    pub naive_ms: f64,
    pub deque_ms: f64,
    pub build_ms: f64,
    pub all_build_ms: f64,
    pub query_extract_ms: f64,
    pub query_select_ms: f64,
    pub lookup_ms: f64,
    pub all_lookup_ms: f64,
    pub verify_ms: f64,
    pub all_verify_ms: f64,
}

fn elapsed(start: Instant) -> f64 {
    start.elapsed().as_secs_f64() * 1000.0
}

fn measure(
    mut f: impl FnMut() -> Result<Selection, String>,
    repeats: usize,
) -> Result<(Selection, f64), String> {
    black_box(f()?); // warm-up
    let mut times = Vec::with_capacity(repeats);
    let mut last = None;
    for _ in 0..repeats {
        let start = Instant::now();
        let result = black_box(f()?);
        times.push(elapsed(start));
        last = Some(result);
    }
    times.sort_by(f64::total_cmp);
    Ok((last.unwrap(), times[times.len() / 2]))
}

pub fn run(config: &Config) -> Result<Metrics, String> {
    config.validate()?;
    let reference = Generator::new(config.seed).dna(config.length);
    let start = Instant::now();
    let all = kmers(&reference, config.k)?;
    let extract_ms = elapsed(start);
    let (naive, naive_ms) = measure(|| minimizer::naive(&all, config.w), config.repeats)?;
    let (deque, deque_ms) = measure(|| minimizer::deque(&all, config.w), config.repeats)?;
    if naive.selected != deque.selected {
        return Err("selector disagreement".into());
    }
    let start = Instant::now();
    let index = Index::build(&deque.selected);
    let build_ms = elapsed(start);
    let start = Instant::now();
    let all_index = Index::build(&all);
    let all_build_ms = elapsed(start);
    let mut metrics = Metrics {
        selected: deque.selected.len(),
        kmers: all.len(),
        distinct_keys: index.distinct_keys(),
        all_distinct_keys: all_index.distinct_keys(),
        bytes: index.payload_bytes(),
        all_bytes: all_index.payload_bytes(),
        allocated: index.allocated_bytes(),
        all_allocated: all_index.allocated_bytes(),
        naive_comparisons: naive.work.comparisons,
        deque_comparisons: deque.work.comparisons,
        pushes: deque.work.pushes,
        pops: deque.work.pops,
        query_hashes: 0,
        query_selected: 0,
        candidates: 0,
        all_candidates: 0,
        mapped: 0,
        all_mapped: 0,
        shared: 0,
        surviving: 0,
        extract_ms,
        naive_ms,
        deque_ms,
        build_ms,
        all_build_ms,
        query_extract_ms: 0.0,
        query_select_ms: 0.0,
        lookup_ms: 0.0,
        all_lookup_ms: 0.0,
        verify_ms: 0.0,
        all_verify_ms: 0.0,
    };
    let mut origins = Generator::new(config.seed ^ 0xabcde);
    for i in 0..config.queries {
        let origin =
            (origins.next_u64() % (config.length - config.query_length + 1) as u64) as usize;
        // Separate mutation RNG keeps query origins identical across rates and w.
        let query = Generator::new(config.seed.wrapping_add(i as u64)).mutate(
            &reference[origin..origin + config.query_length],
            config.mutation_rate,
        )?;
        let start = Instant::now();
        let q = kmers(&query, config.k)?;
        metrics.query_extract_ms += elapsed(start);
        let start = Instant::now();
        let qs = minimizer::deque(&q, config.w)?.selected;
        metrics.query_select_ms += elapsed(start);
        metrics.query_hashes += q.len();
        metrics.query_selected += qs.len();
        let start = Instant::now();
        let candidates = query::candidate_starts(&index, &qs, reference.len(), query.len());
        metrics.lookup_ms += elapsed(start);
        let start = Instant::now();
        let all_candidates = query::candidate_starts(&all_index, &q, reference.len(), query.len());
        metrics.all_lookup_ms += elapsed(start);
        metrics.candidates += candidates.len();
        metrics.all_candidates += all_candidates.len();
        metrics.mapped += usize::from(candidates.iter().any(|x| x.0 as usize == origin));
        metrics.all_mapped += usize::from(all_candidates.iter().any(|x| x.0 as usize == origin));
        let start = Instant::now();
        for &(position, _) in &candidates {
            black_box(query::hamming(
                &reference[position as usize..position as usize + query.len()],
                &query,
            ));
        }
        metrics.verify_ms += elapsed(start);
        let start = Instant::now();
        for &(position, _) in &all_candidates {
            black_box(query::hamming(
                &reference[position as usize..position as usize + query.len()],
                &query,
            ));
        }
        metrics.all_verify_ms += elapsed(start);
        // Oracle metrics are outside timed paths; only the true reference slice is needed.
        let recall = query::seed_recall(
            &all[origin..origin + q.len()],
            &deque.selected,
            &q,
            &qs,
            origin as u32,
        );
        metrics.shared += recall.shared_selected;
        metrics.surviving += recall.surviving;
    }
    Ok(metrics)
}

impl Metrics {
    pub const CSV_HEADER: &str = "length,query_length,queries,k,w,mutation_rate,seed,repeats,reference_hashes,query_hashes,selected,query_selected,density,distinct_keys,all_distinct_keys,index_bytes,all_index_bytes,allocated_bytes,all_allocated_bytes,memory_reduction,candidates,all_candidates,mapped,all_mapped,mapping_recall,all_mapping_recall,shared_seeds,surviving_seeds,seed_recall,naive_comparisons,deque_comparisons,deque_pushes,deque_pops,extract_ms,naive_ms,deque_ms,build_ms,all_build_ms,query_extract_ms,query_select_ms,lookup_ms,all_lookup_ms,verify_ms,all_verify_ms";

    pub fn csv(&self, c: &Config) -> String {
        let recall = if self.surviving == 0 {
            "NA".into()
        } else {
            format!("{:.6}", self.shared as f64 / self.surviving as f64)
        };
        let fields = vec![
            c.length.to_string(),
            c.query_length.to_string(),
            c.queries.to_string(),
            c.k.to_string(),
            c.w.to_string(),
            c.mutation_rate.to_string(),
            c.seed.to_string(),
            c.repeats.to_string(),
            self.kmers.to_string(),
            self.query_hashes.to_string(),
            self.selected.to_string(),
            self.query_selected.to_string(),
            format!("{:.6}", self.selected as f64 / self.kmers as f64),
            self.distinct_keys.to_string(),
            self.all_distinct_keys.to_string(),
            self.bytes.to_string(),
            self.all_bytes.to_string(),
            self.allocated.to_string(),
            self.all_allocated.to_string(),
            format!("{:.6}", 1.0 - self.bytes as f64 / self.all_bytes as f64),
            self.candidates.to_string(),
            self.all_candidates.to_string(),
            self.mapped.to_string(),
            self.all_mapped.to_string(),
            format!("{:.6}", self.mapped as f64 / c.queries as f64),
            format!("{:.6}", self.all_mapped as f64 / c.queries as f64),
            self.shared.to_string(),
            self.surviving.to_string(),
            recall,
            self.naive_comparisons.to_string(),
            self.deque_comparisons.to_string(),
            self.pushes.to_string(),
            self.pops.to_string(),
        ];
        let times = [
            self.extract_ms,
            self.naive_ms,
            self.deque_ms,
            self.build_ms,
            self.all_build_ms,
            self.query_extract_ms,
            self.query_select_ms,
            self.lookup_ms,
            self.all_lookup_ms,
            self.verify_ms,
            self.all_verify_ms,
        ];
        format!(
            "{},{}",
            fields.join(","),
            times
                .iter()
                .map(|x| format!("{x:.6}"))
                .collect::<Vec<_>>()
                .join(",")
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn experiment_counts_and_csv_are_consistent() {
        let c = Config {
            length: 200,
            query_length: 60,
            queries: 5,
            k: 7,
            w: 5,
            mutation_rate: 0.0,
            repeats: 1,
            ..Config::default()
        };
        let m = run(&c).unwrap();
        assert_eq!(m.mapped, 5);
        assert_eq!(m.all_mapped, 5);
        assert_eq!(m.kmers, 194);
        assert_eq!(m.query_hashes, 270);
        assert_eq!(m.surviving, 270);
        assert!(m.bytes < m.all_bytes);
        assert_eq!(
            m.csv(&c).split(',').count(),
            Metrics::CSV_HEADER.split(',').count()
        );
        let empty = run(&Config {
            w: 1000,
            ..c.clone()
        })
        .unwrap();
        assert_eq!(empty.selected, 0);
        assert_eq!(empty.mapped, 0);
        assert!(
            Config {
                mutation_rate: f64::NAN,
                ..c
            }
            .validate()
            .is_err()
        );
    }
}
