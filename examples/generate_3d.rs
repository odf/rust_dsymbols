use rust_dsymbols::generators::dset_generators::DSets;


fn main() {
    let args: Vec<_> = std::env::args().collect();

    // TODO add some error handling
    let max_size = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(16);

    call_generate(max_size);
}


#[cfg(not(feature = "pprof"))]
fn call_generate(n: usize) {
    generate_binary(n);
}


#[cfg(feature = "pprof")]
fn call_generate(n: usize) {
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
    let mut count: usize = 0;

    for dset in DSets::new(3, n) {
        count += 1;
        // TODO call D-symbol generator here once it's implemented
        println!("{}", dset);
    }

    eprintln!("{} D-sets generated.", count);
}
