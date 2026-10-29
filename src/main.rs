use minimizer_kmer_index::{
    experiment::{self, Config, Metrics},
    report, verify,
};
use std::{collections::BTreeMap, fs, path::Path, process::ExitCode};

const HELP: &str = "minimizer-kmer-index — deterministic DNA minimizer experiments

Commands:
  demo        Run one experiment and print its measurements
  bench       Sweep window sizes and mutation rates; write CSV
  visualize   Write a standalone interactive HTML report (query <= 160 bases)
  verify      Cross-check selectors, candidates and exhaustive Hamming matches

Options (demo/bench/visualize):
  --length N          Reference bases (default 100000; visualize 96)
  --query-length N    Query bases (default 200; visualize 96)
  --queries N         Query count (default 100; visualize 1)
  --k N               K-mer bases, 1..32 (default 15; visualize 5)
  --w N               K-mers per window (default 10; visualize 8)
  --mutation RATE     Substitution probability, 0..1 (default .02; visualize .08)
  --seed N            Data seed (default 42; also accepted by verify)
  --repeats N         Selector timing repeats after warm-up (default 3)
  --output PATH       bench: artifacts/bench.csv; visualize: artifacts/report.html
  --windows LIST      bench only: comma-separated w values (default 1,5,10,20,40)
  --rates LIST        bench only: comma-separated rates (default 0,.01,.05,.1)

--w and --mutation override their respective sweep lists when no list is given.
Reference and query extraction compute one hash per k-mer. Both selectors reuse
the reference fingerprints. Selection timings are medians; other timings are
single runs or totals over queries. No reverse complements or indels.
";

fn parse<T: std::str::FromStr>(value: &str, name: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("invalid value for {name}: {value}"))
}

fn write(path: &str, content: &str) -> Result<(), String> {
    if let Some(parent) = Path::new(path)
        .parent()
        .filter(|x| !x.as_os_str().is_empty())
    {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    fs::write(path, content).map_err(|e| format!("cannot write {path}: {e}"))
}

fn execute(args: Vec<String>) -> Result<(), String> {
    let Some(command) = args.first().map(String::as_str) else {
        print!("{HELP}");
        return Ok(());
    };
    if ["help", "--help", "-h"].contains(&command) {
        print!("{HELP}");
        return Ok(());
    }
    if !["demo", "bench", "visualize", "verify"].contains(&command) {
        return Err(format!("unknown command: {command}"));
    }
    let mut opts = BTreeMap::new();
    let mut rest = args[1..].iter();
    while let Some(key) = rest.next() {
        if key == "--help" {
            print!("{HELP}");
            return Ok(());
        }
        let common = [
            "--length",
            "--query-length",
            "--queries",
            "--k",
            "--w",
            "--mutation",
            "--seed",
            "--repeats",
        ];
        let allowed = if command == "verify" {
            key == "--seed"
        } else {
            common.contains(&key.as_str())
                || (command != "demo" && key == "--output")
                || (command == "bench" && ["--windows", "--rates"].contains(&key.as_str()))
        };
        if !allowed {
            return Err(format!("unknown option for {command}: {key}"));
        }
        let value = rest
            .next()
            .ok_or_else(|| format!("missing value for {key}"))?;
        if opts.insert(key.as_str(), value.as_str()).is_some() {
            return Err(format!("duplicate option: {key}"));
        }
    }
    let mut c = if command == "visualize" {
        Config {
            length: 96,
            query_length: 96,
            queries: 1,
            k: 5,
            w: 8,
            mutation_rate: 0.08,
            ..Config::default()
        }
    } else {
        Config::default()
    };
    macro_rules! option {
        ($name:literal,$field:ident) => {
            if let Some(value) = opts.get($name) {
                c.$field = parse(value, $name)?;
            }
        };
    }
    option!("--length", length);
    option!("--query-length", query_length);
    option!("--queries", queries);
    option!("--k", k);
    option!("--w", w);
    option!("--mutation", mutation_rate);
    option!("--seed", seed);
    option!("--repeats", repeats);
    if command == "verify" {
        println!(
            "Verified {} small-input cases (seed {}).",
            verify::run(c.seed)?,
            c.seed
        );
        return Ok(());
    }
    c.validate()?;
    if command == "bench" {
        let windows: Vec<usize> = if let Some(value) = opts.get("--windows") {
            value
                .split(',')
                .map(|x| parse(x, "--windows"))
                .collect::<Result<_, _>>()?
        } else if opts.contains_key("--w") {
            vec![c.w]
        } else {
            vec![1, 5, 10, 20, 40]
        };
        let rates: Vec<f64> = if let Some(value) = opts.get("--rates") {
            value
                .split(',')
                .map(|x| parse(x, "--rates"))
                .collect::<Result<_, _>>()?
        } else if opts.contains_key("--mutation") {
            vec![c.mutation_rate]
        } else {
            vec![0.0, 0.01, 0.05, 0.1]
        };
        // Validate the complete grid before starting expensive work.
        for &w in &windows {
            for &rate in &rates {
                Config {
                    w,
                    mutation_rate: rate,
                    ..c.clone()
                }
                .validate()?;
            }
        }
        let mut csv = format!("{}\n", Metrics::CSV_HEADER);
        for w in windows {
            for &rate in &rates {
                let config = Config {
                    w,
                    mutation_rate: rate,
                    ..c.clone()
                };
                let m = experiment::run(&config)?;
                csv.push_str(&m.csv(&config));
                csv.push('\n');
                eprintln!(
                    "w={w}, mutation={rate}: selected {}, mapped {}/{}",
                    m.selected, m.mapped, c.queries
                );
            }
        }
        let path = opts
            .get("--output")
            .copied()
            .unwrap_or("artifacts/bench.csv");
        write(path, &csv)?;
        println!("Wrote {path}");
    } else if command == "visualize" {
        let path = opts
            .get("--output")
            .copied()
            .unwrap_or("artifacts/report.html");
        write(path, &report::render(&c)?)?;
        println!("Wrote {path}");
    } else {
        let m = experiment::run(&c)?;
        println!(
            "Reference {} bases; {} queries of {} bases; k={}, w={}, mutation={}, seed={}",
            c.length, c.queries, c.query_length, c.k, c.w, c.mutation_rate, c.seed
        );
        println!(
            "Selected {}/{} reference k-mers ({:.2}%)",
            m.selected,
            m.kmers,
            100.0 * m.selected as f64 / m.kmers as f64
        );
        println!(
            "Index payload {} vs {} bytes; reduction {:.2}%",
            m.bytes,
            m.all_bytes,
            100.0 * (1.0 - m.bytes as f64 / m.all_bytes as f64)
        );
        println!(
            "Allocated arrays {} vs {} bytes",
            m.allocated, m.all_allocated
        );
        println!(
            "Candidates {} vs {}; true starts found {}/{} vs {}/{}",
            m.candidates, m.all_candidates, m.mapped, c.queries, m.all_mapped, c.queries
        );
        println!(
            "Seed recall: {} / {} surviving true-position anchors",
            m.shared, m.surviving
        );
        println!(
            "Hashes: {} reference + {} query (shared by selectors)",
            m.kmers, m.query_hashes
        );
        println!(
            "Comparisons naive/deque: {}/{}; deque pushes/pops: {}/{}",
            m.naive_comparisons, m.deque_comparisons, m.pushes, m.pops
        );
        println!(
            "Selection median ms naive/deque: {:.3}/{:.3}",
            m.naive_ms, m.deque_ms
        );
        println!("\n{}\n{}", Metrics::CSV_HEADER, m.csv(&c));
    }
    Ok(())
}

fn main() -> ExitCode {
    match execute(std::env::args().skip(1).collect()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("error: {error}\nUse --help for usage.");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rejects_invalid_cli_without_running_work() {
        for args in [
            vec!["nope"],
            vec!["demo", "--k", "0"],
            vec!["demo", "--w"],
            vec!["demo", "--seed", "1", "--seed", "2"],
            vec!["verify", "--output", "x"],
            vec!["bench", "--rates", "NaN"],
            vec!["visualize", "--output", "x", "--query-length", "999"],
        ] {
            assert!(execute(args.into_iter().map(String::from).collect()).is_err());
        }
    }
}
