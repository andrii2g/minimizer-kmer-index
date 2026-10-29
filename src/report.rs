//! Standalone HTML/SVG report; no network, JavaScript dependencies or user HTML.
use crate::{
    dna::{Generator, Kmer, kmers},
    experiment::{self, Config},
    index::Index,
    minimizer, query,
};
use std::fmt::Write;

fn winners(items: &[Kmer], w: usize, offset: u32) -> String {
    if items.len() < w {
        return String::new();
    }
    items
        .windows(w)
        .map(|window| {
            let best = window
                .iter()
                .min_by_key(|x| (x.hash, x.key, x.position))
                .unwrap();
            (best.position - offset).to_string()
        })
        .collect::<Vec<_>>()
        .join(",")
}

pub fn render(c: &Config) -> Result<String, String> {
    c.validate()?;
    if c.query_length > 160 || c.length > 10_000 {
        return Err("visualize requires query-length <= 160 and reference length <= 10000; use bench for larger inputs".into());
    }
    if c.w > c.query_length - c.k + 1 {
        return Err(
            "visualize requires at least one complete query window (k+w-1 <= query-length)".into(),
        );
    }
    let reference = Generator::new(c.seed).dna(c.length);
    let mut origins = Generator::new(c.seed ^ 0xabcde);
    let origin = (origins.next_u64() % (c.length - c.query_length + 1) as u64) as usize;
    let source = &reference[origin..origin + c.query_length];
    let query = Generator::new(c.seed).mutate(source, c.mutation_rate)?;
    let all = kmers(&reference, c.k)?;
    let selected = minimizer::deque(&all, c.w)?.selected;
    let q = kmers(&query, c.k)?;
    let qs = minimizer::deque(&q, c.w)?.selected;
    let index = Index::build(&selected);
    let candidates = query::candidates(&index, &qs, &reference, &query);
    let local_all = &all[origin..origin + q.len()];
    let local_selected: Vec<_> = selected
        .iter()
        .filter(|x| (x.position as usize) >= origin && (x.position as usize) < origin + q.len())
        .copied()
        .collect();
    let recall = query::seed_recall(local_all, &selected, &q, &qs, origin as u32);
    let m = experiment::run(c)?;
    let pitch = 16;
    let left = 140;
    let x = |position: usize| left + position * pitch;
    let ref_dna_y = 65 + c.k * 5;
    let ref_seed_y = ref_dna_y + 35;
    let query_seed_y = ref_seed_y + 110;
    let query_dna_y = query_seed_y + 40;
    let query_all_y = query_dna_y + 35;
    let height = query_all_y + c.k * 5 + 35;
    let width = x(c.query_length) + 30;
    let mut svg = format!(
        "<svg id='alignment' viewBox='0 0 {width} {height}' role='img' aria-label='Aligned reference and mutated query k-mers'><rect id='window' x='{left}' y='28' width='{}' height='{}' fill='#c7d2fe' opacity='.3'/>",
        (c.k + c.w - 1) * pitch,
        height - 40
    );
    for (label, y) in [
        ("All reference k-mers", 45),
        ("Reference bases", ref_dna_y),
        ("Selected reference", ref_seed_y),
        ("Selected query", query_seed_y),
        ("Query bases", query_dna_y),
        ("All query k-mers", query_all_y),
    ] {
        write!(svg, "<text x='8' y='{y}' class='label'>{label}</text>").unwrap();
    }
    for (items, offset, base_y) in [(local_all, origin, 45), (&q, 0, query_all_y)] {
        for item in items {
            let position = item.position as usize - offset;
            let y = base_y + (position % c.k) * 5;
            let dna = if offset == origin && base_y == 45 {
                source
            } else {
                &query
            };
            write!(svg,"<line x1='{}' x2='{}' y1='{y}' y2='{y}' stroke='#94a3b8' stroke-width='2'><title>position {}; {} ; rank {}</title></line>",x(position),x(position+c.k)-2,item.position,String::from_utf8_lossy(&dna[position..position+c.k]),item.hash).unwrap();
        }
    }
    for (i, (&a, &b)) in source.iter().zip(&query).enumerate() {
        for (base, y) in [(a, ref_dna_y), (b, query_dna_y)] {
            if a != b {
                write!(
                    svg,
                    "<rect x='{}' y='{}' width='15' height='20' fill='#fee2e2'/>",
                    x(i),
                    y - 15
                )
                .unwrap();
            }
            write!(
                svg,
                "<text x='{}' y='{y}' class='base'>{}</text>",
                x(i) + 2,
                base as char
            )
            .unwrap();
        }
        if i % 10 == 0 {
            write!(
                svg,
                "<text x='{}' y='{}' class='tick'>{i}</text>",
                x(i),
                ref_dna_y + 16
            )
            .unwrap();
        }
    }
    for item in &local_selected {
        let position = item.position as usize - origin;
        let shared = qs
            .iter()
            .any(|q| q.position as usize == position && q.key == item.key);
        let color = if shared { "#0d9488" } else { "#d97706" };
        write!(svg,"<rect x='{}' y='{ref_seed_y}' width='{}' height='7' rx='2' fill='{color}'><title>reference position {}; {}</title></rect>",x(position),c.k*pitch-2,item.position,if shared {"shared anchor"} else {"no selected anchor at the true query position"}).unwrap();
        if shared {
            write!(svg,"<line x1='{}' x2='{}' y1='{}' y2='{query_seed_y}' stroke='#0d9488' opacity='.55'/>",x(position)+c.k*pitch/2,x(position)+c.k*pitch/2,ref_seed_y+7).unwrap();
        }
    }
    for item in &qs {
        let shared = local_selected
            .iter()
            .any(|r| r.position as usize == origin + item.position as usize && r.key == item.key);
        write!(svg,"<rect x='{}' y='{query_seed_y}' width='{}' height='7' rx='2' fill='{}'><title>query position {}; rank {}</title></rect>",x(item.position as usize),c.k*pitch-2,if shared {"#0d9488"} else {"#7c3aed"},item.position,item.hash).unwrap();
    }
    write!(svg,"<circle id='ref-winner' cy='{}' r='5' fill='#1d4ed8'/><circle id='query-winner' cy='{}' r='5' fill='#1d4ed8'/></svg>",ref_seed_y+15,query_seed_y-8).unwrap();
    let max_votes = candidates.iter().map(|c| c.votes).max().unwrap_or(1);
    let mut histogram = String::from(
        "<svg viewBox='0 0 1000 190' role='img' aria-label='Candidate start votes'><line x1='30' x2='980' y1='150' y2='150' stroke='#64748b'/>",
    );
    for candidate in &candidates {
        let xpos = 30.0 + candidate.start as f64 / (c.length - c.query_length + 1) as f64 * 940.0;
        let bar = candidate.votes as f64 / max_votes as f64 * 120.0;
        write!(histogram,"<rect x='{xpos:.2}' y='{:.2}' width='5' height='{bar:.2}' fill='{}'><title>start {}; {} votes; {} mismatches</title></rect>",150.0-bar,if candidate.start as usize == origin {"#0d9488"} else {"#94a3b8"},candidate.start,candidate.votes,candidate.mismatches).unwrap();
    }
    write!(histogram,"<text x='30' y='175' class='label'>0</text><text x='850' y='175' class='label'>reference start {}</text></svg>",c.length-c.query_length).unwrap();
    let mutation_count = source.iter().zip(&query).filter(|(a, b)| a != b).count();
    let status = if candidates.iter().any(|x| x.start as usize == origin) {
        "True start recovered"
    } else {
        "True start missed: no shared selected anchor"
    };
    let cards = format!(
        "<div class='card'><b>{}/{}</b><span>Reference minimizers / all k-mers</span></div><div class='card'><b>{:.1}%</b><span>Index payload reduction</span></div><div class='card'><b>{}/{}</b><span>First-query shared / surviving anchors</span></div><div class='card'><b>{}/{}</b><span>Queries with true-start candidates</span></div>",
        m.selected,
        m.kmers,
        100.0 * (1.0 - m.bytes as f64 / m.all_bytes as f64),
        recall.shared_selected,
        recall.surviving,
        m.mapped,
        c.queries
    );
    let reference_winners = winners(local_all, c.w, origin as u32);
    let query_winners = winners(&q, c.w, 0);
    let parameters = format!(
        "k={}, w={} k-mers ({} bases/window), seed={}, requested mutation rate={}, actual first-query substitutions={}/{}. Reference={} bases; query={} bases; displayed reference origin={}. All positions are zero-based; forward strand only.",
        c.k,
        c.w,
        c.k + c.w - 1,
        c.seed,
        c.mutation_rate,
        mutation_count,
        query.len(),
        reference.len(),
        query.len(),
        origin
    );
    let csv = format!("{}\n{}", crate::experiment::Metrics::CSV_HEADER, m.csv(c));
    Ok(include_str!("report.html")
        .replace("{{PARAMETERS}}", &parameters)
        .replace("{{CARDS}}", &cards)
        .replace("{{ALIGNMENT}}", &svg)
        .replace("{{HISTOGRAM}}", &histogram)
        .replace("{{STATUS}}", status)
        .replace("{{MAX_WINDOW}}", &(q.len() - c.w).to_string())
        .replace("{{REF_WINNERS}}", &reference_winners)
        .replace("{{QUERY_WINNERS}}", &query_winners)
        .replace("{{ORIGIN}}", &origin.to_string())
        .replace("{{K}}", &c.k.to_string())
        .replace("{{CSV}}", &csv))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn report_is_self_contained_and_checks_limits() {
        let c = Config {
            length: 60,
            query_length: 60,
            queries: 1,
            k: 5,
            w: 4,
            repeats: 1,
            ..Config::default()
        };
        let html = render(&c).unwrap();
        assert!(!html.contains("{{"));
        assert!(!html.contains("https://"));
        assert!(html.contains("id='alignment'"));
        assert!(html.contains("const refWinners = ["));
        assert!(
            render(&Config {
                w: 100,
                ..c.clone()
            })
            .is_err()
        );
        assert!(
            render(&Config {
                length: 200,
                query_length: 170,
                ..c
            })
            .is_err()
        );
    }
}
