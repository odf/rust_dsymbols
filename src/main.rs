use std::io::{Write, stdout};

use rust_dsymbols::dset_generators::DSets;
use rust_dsymbols::dsym_generators::{DSyms, Geometries};


fn main() {
    let args: Vec<_> = std::env::args().collect();

    // TODO add some error handling
    let max_size = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(16);
    let geometries = args.get(2).and_then(|s| s.parse().ok())
        .unwrap_or(Geometries::All);

    call_generate(max_size, geometries);
}


#[cfg(not(feature = "pprof"))]
fn call_generate(n: usize, geometries: Geometries) {
    generate_binary(n, geometries);
}


#[cfg(feature = "pprof")]
fn call_generate(n: usize, geometries: Geometries) {
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(1000)
        .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        .build().unwrap();

    generate_binary(n, geometries);

    if let Ok(report) = guard.report().build() {
        let file = std::fs::File::create("flamegraph.svg").unwrap();
        report.flamegraph(file).unwrap();
    };
}


fn generate_binary(n: usize, geometries: Geometries) {
    let mut count: u64 = 0;
    let mut previous = vec![];

    for dset in DSets::new(2, n) {
        for dsym in DSyms::new(&dset, geometries) {
            count += 1;

            let code = dsym.to_binary().unwrap();

            for i in 0..=code.len() {
                if i >= code.len()
                    || i >= previous.len()
                    || code[i] != previous[i]
                {
                    stdout().write_all(&[(i + 1) as u8]).unwrap();
                    stdout().write_all(&code[i..]).unwrap();

                    previous = code;
                    break;
                }
            }
        }
    }

    eprintln!("{} symbols generated.", count);
}
