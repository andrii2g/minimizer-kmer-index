# minimizer-kmer-index

A dependency-free Rust demonstration of the memory–sensitivity trade-off in DNA
minimizer indexing. Compare a window-scanning selector against a monotonic deque,
build a compact position index, and recover mutated queries through shared seeds.

## Run

Install current stable Rust (the project uses edition 2024), then:

```sh
cargo run --release -- demo
cargo run --release -- verify
cargo run --release -- bench
cargo run --release -- visualize
```

`bench` writes `artifacts/bench.csv`. Open `artifacts/report.html` in a browser for
the visualization. These generated outputs are ignored by Git. No network access,
external datasets, services or third-party crates are needed to run the project.

See the checked-in [example report](docs/example.html), [missed-match report](docs/missed.html)
and [example benchmark CSV](docs/example-bench.csv). Download the HTML reports and
open them locally; GitHub's source viewer does not execute their JavaScript.

```sh
# One experiment: larger windows and more substitutions
cargo run --release -- demo --length 1000000 --query-length 200 --queries 100 --k 15 --w 40 --mutation 0.1 --seed 42

# A custom parameter sweep
cargo run --release -- bench --windows 1,10,40 --rates 0,0.05,0.1 --output artifacts/custom.csv

# An aligned sequence with a visibly missed source match
cargo run --release -- visualize --k 15 --w 40 --mutation 0.08 --seed 2 --output artifacts/missed.html

# Complete CLI options
cargo run -- --help
```

## Algorithm contract

- Input is uppercase `A/C/G/T`. Invalid bases are errors, including in inputs shorter than k.
- `1 <= k <= 32`; two-bit encoding uses `A=0, C=1, G=2, T=3` in a `u64`.
- Positions are zero-based `u32`; references longer than `u32::MAX` bases are rejected.
- Each k-mer is encoded with a rolling shift/mask and ranked using the fixed SplitMix64 finalizer.
  Ranking is by `(hash, exact_key)`, with the leftmost occurrence winning ties.
  Index lookups use **exact encoded keys**, not hashes. Both sides must use the same k.
- A window contains **w k-mers**, spanning `k + w - 1` bases. `w=0` is invalid.
  Fewer than w k-mers yields no selected minimizers; no partial windows are emitted.
- Consecutive windows selecting the same position produce one indexed occurrence.
  Different occurrences of an identical k-mer remain separate positions.
- Sequence generation uses a seeded SplitMix64 generator. A substitution always changes
  the original base. Query origins stay the same across w and mutation-rate sweeps.

| Stage | Implementation | Cost |
|---|---|---|
| Extract | Rolling two-bit encoding, one fixed hash per k-mer | O(n) time, O(n) stored fingerprints |
| Naïve selection | Scan every complete window | O(nw) time |
| Deque selection | Expire old entries; pop strictly larger ranks | O(n) time, O(w) auxiliary queue |
| Index build | Sort and group exact-key/position pairs | O(m log m) time, O(m) scratch |
| Lookup | Binary search keys, enumerate postings, vote for starts | Depends on hits; repetitive keys can be expensive |
| Verification | Hamming distance at every seed-supported start | O(candidates × query length) |

Here n is the number of k-mers and m is the number of selected positions. These
are educational in-memory implementations: linear-time selection does not imply
constant-memory extraction or index construction.

The index uses three arrays: sorted `u64` keys, `u64` posting offsets, and `u32`
positions. The all-k-mer baseline uses exactly the same layout. Candidate starts
come from `reference_position - query_position`, must fit the reference, and rank
by distinct supporting query positions, then start position. Repeated input seeds
cannot inflate votes. Every candidate can be verified with Hamming distance;
the library also exposes exhaustive scanning as an independent oracle.

## Metrics

The demo prints a summary and one CSV row. The benchmark defaults to 20 combinations:
`w=1,5,10,20,40` and substitution rates `0,.01,.05,.1`. `--w` or `--mutation` limits
that dimension unless an explicit `--windows` or `--rates` list is supplied.

| Metric | Meaning |
|---|---|
| Selected / density | Unique selected positions divided by all k-mer positions |
| Index payload | `8*keys + 8*offsets + 4*positions`, including the final offset |
| Allocated bytes | Sum of array capacities times element sizes |
| Memory reduction | `1 - minimizer_payload / all_kmer_payload` |
| Candidates | Sum of distinct in-bounds voted starts over queries, before filtering |
| Mapping recall | Fraction of known query origins appearing among candidates |
| Seed recall | Shared selected anchors at the true positions / all surviving exact k-mer anchors |
| Hash operations | One per extracted reference/query k-mer; shared by selectors |
| Comparisons | Rank comparisons during selection; queue pushes/pops are also counted |

Seed recall is a pooled ratio over queries and is `NA` if no exact source anchors
survive. It is normally much lower than mapping recall: one shared anchor can
recover a source start. All-k-mer mapping recall is reported separately; even
that baseline can miss a heavily mutated query with no intact k-mer.

Payload and capacity figures exclude `Vec` headers, allocator metadata, stored
DNA, extraction/selection buffers, build scratch, and query voting structures.
They measure index storage, **not process peak memory**.

Selector timings use the same precomputed fingerprints, a warm-up and the median
of `--repeats` runs (default 3; upper middle sample for an even repeat count).
Extraction and index build are single-run timings. Query extraction, selection,
lookup/voting and Hamming verification are separate totals over all queries.
Synthetic generation and oracle recall accounting are outside measured phases.
Small-window scans can be faster despite doing more comparisons; use release
builds, larger windows and repeated runs before drawing runtime conclusions.

## Example results

Measured on Windows with Rust 1.97.1, seed 42, reference length 100,000, 100 queries
of length 200, k=15 and a 10% substitution rate. Full raw output, including runtime,
is in [example-bench.csv](docs/example-bench.csv). These are synthetic observations,
not a sensitivity guarantee for biological datasets.

| w | Reference positions retained | Index payload reduction | True starts recovered |
|---:|---:|---:|---:|
| 1 | 100.00% | 0.00% | 100/100 |
| 5 | 33.40% | 66.60% | 99/100 |
| 10 | 18.15% | 81.85% | 98/100 |
| 20 | 9.52% | 90.49% | 93/100 |
| 40 | 4.94% | 95.06% | 74/100 |

The default 2% substitution demo with w=10 recovered 100/100 source starts with
81.85% index payload reduction. The checked-in missed-match report uses seed 2,
k=15, w=40 and an 8% substitution rate: 18 exact source k-mers survive, but none
are selected on both sides. The all-k-mer baseline finds the true start while
the minimizer index misses it.

Random-order minimizer density is approximately
`2/(w+1)` under suitable assumptions about repeats; it is not a bound for arbitrary
sequences. Two strings sharing an exact substring of at least `k+w-1` bases have
a shared minimizer when selection and tie-breaking are consistent. Substitutions
can remove this condition.

## Visualization

The standalone HTML report shows every k-mer, selected minimizers, substitutions,
true-position anchor connections, and a histogram of candidate-start votes.
Move the window slider to see the two window winners. Tooltips reveal keys' ranks,
positions, vote counts and candidate Hamming distances.

The diagram shows the first generated query aligned to its source substring.
Reference minimizers are selected from the full reference, including flanking
windows. Cards summarize all configured queries; the diagram and histogram show
only the first. Visualizations require query length <=160, reference length <=10,000
and at least one complete query window. Experiments can use larger references or
windows without these display constraints.

## Validate and extend

```sh
cargo fmt --check
cargo test
cargo clippy --all-targets -- -D warnings
cargo run --release -- verify --seed 42
```

Tests cover exhaustive small DNA selector equivalence, independent slice encoding,
k=32, ties, forced hash collisions, repeated seeds, short/invalid inputs, mutation
boundaries, index postings, candidate votes and CLI output. The runtime verifier
checks 3,000 generated cases against the all-k-mer index and exhaustive Hamming
scanning. Exhaustive matches must appear in the all-k-mer candidates only when
the query length and mismatch limit guarantee an intact k-mer; minimizer results
are deliberately allowed to lose candidates.

Library entry points live in `dna`, `minimizer`, `index`, `query`, `experiment`,
`verify` and `report`. The MVP supports a single forward-strand reference and
substitutions. Reverse-complement canonicalization, FASTA input, indels, alignment
chaining, repetitive-seed filtering and parallelism are future extensions.

Background: [weighted minimizer mapping](https://pmc.ncbi.nlm.nih.gov/articles/PMC7355284/)
and [theory of local k-mer selection](https://academic.oup.com/bioinformatics/article/38/20/4659/6432031).

Licensed under [MIT](LICENSE).
