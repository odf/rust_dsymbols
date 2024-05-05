use std::io::stdin;

use rust_dsymbols::euclidicity::is_euclidean;
use rust_dsymbols::dsyms::PartialDSym;
use rust_dsymbols::euclidicity::Euclidean::*;


fn main() {
    run();
}


#[cfg(not(feature = "pprof"))]
fn run() {
    run_test();
}


#[cfg(feature = "pprof")]
fn run() {
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(1000)
        .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        .build().unwrap();

    run_test();

    if let Ok(report) = guard.report().build() {
        let file = std::fs::File::create("flamegraph.svg").unwrap();
        report.flamegraph(file).unwrap();
    };
}


fn run_test() {
    let mut i = 0;
    let mut good = vec![];
    let mut ambiguous = vec![];

    for line in stdin().lines() {
        i += 1;

        if let Some(ds) = line.ok().and_then(|s|
            s.parse::<PartialDSym>().ok()
        ) {
            match is_euclidean(&ds) {
                Yes => {
                    let s = "simplified cover recognized".to_string();
                    println!("#Symbol {i} is good: {s}");
                    println!("{ds}");
                    good.push(i);
                },
                No(s) => println!("#Symbol {i} is bad: {s}"),
                Maybe(s, out) => {
                    println!("#Symbol {i} is undecided: {s}");
                    println!("#??? {out}");
                    ambiguous.push(i);
                },
            }
        }
    }

    println!(
        "### {} good symbols: {}",
        good.len(),
        good.iter().map(|i| i.to_string()).collect::<Vec<_>>().join(" ")
    );
    println!(
        "### {} ambiguous symbols: {}",
        ambiguous.len(),
        ambiguous.iter().map(|i| i.to_string()).collect::<Vec<_>>().join(" ")
    );
}
