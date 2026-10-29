use std::{fs, process::Command};

fn binary() -> Command {
    Command::new(env!("CARGO_BIN_EXE_minimizer-kmer-index"))
}

#[test]
fn commands_produce_consistent_files_and_fail_cleanly() {
    let directory = std::env::temp_dir().join(format!("minimizer-cli-{}", std::process::id()));
    fs::create_dir_all(&directory).unwrap();
    let csv = directory.join("bench.csv");
    let html = directory.join("report.html");
    let result = binary()
        .args([
            "bench",
            "--length",
            "300",
            "--query-length",
            "60",
            "--queries",
            "3",
            "--k",
            "5",
            "--windows",
            "1,5",
            "--rates",
            "0,0.1",
            "--repeats",
            "1",
            "--output",
        ])
        .arg(&csv)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let text = fs::read_to_string(&csv).unwrap();
    let rows: Vec<_> = text.lines().collect();
    assert_eq!(rows.len(), 5);
    assert!(
        rows.iter()
            .all(|row| row.split(',').count() == rows[0].split(',').count())
    );
    let result = binary()
        .args(["visualize", "--output"])
        .arg(&html)
        .output()
        .unwrap();
    assert!(result.status.success());
    assert!(
        fs::read_to_string(&html)
            .unwrap()
            .contains("id=\"window-slider\"")
    );
    let result = binary()
        .args(["demo", "--mutation", "NaN"])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("mutation rate"));
    fs::remove_file(csv).unwrap();
    fs::remove_file(html).unwrap();
    fs::remove_dir(directory).unwrap();
}
