use rust_dsymbols::dset_generators::DSets;
use rust_dsymbols::dsym_generators::DSyms;


fn main() {
    let mut count = 0;

    if let Some(arg) = std::env::args().nth(1) {
        let n: usize = arg.parse().unwrap();

        for dset in DSets::new(2, n) {
            for dsym in DSyms::new(&dset) {
                println!("{}", dsym);
                count += 1;
            }
        }

        println!("Found {} in total.", count);
        println!();
    } else {
        panic!("Expected an argument.");
    }
    println!();
}
