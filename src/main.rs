use std::io::{Write, stdout};

use rust_dsymbols::dset_generators::DSets;
use rust_dsymbols::dsym_generators::DSyms;


fn main() {
    if let Some(arg) = std::env::args().nth(1) {
        let n: usize = arg.parse().unwrap();

        generate_binary(n);
        //generate_binary_profiled(n);
    } else {
        panic!("Expected an argument.");
    }
}


#[cfg(pprof)]
fn generate_binary_profiled(n: usize) {
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(1000)
        .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        .build().unwrap();

    generate_binary(n);

    if let Ok(report) = guard.report().build() {
        let file = std::fs::File::create("flamegraph.svg").unwrap();
        report.flamegraph(file).unwrap();
    };
}


fn generate_binary(n: usize) {
    let mut count: u64 = 0;
    let mut previous = vec![];

    for dset in DSets::new(2, n) {
        for dsym in DSyms::new(&dset) {
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
