use rust_dsymbols::dset_generators::DSets;
use rust_dsymbols::dsets::DSet;


fn main() {
    let mut count = 0;

    if let Some(arg) = std::env::args().nth(1) {
        let n: usize = arg.parse().unwrap();

        for ds in DSets::new(2, n) {
            if ds.size() == n {
                println!("{}", ds);
                count += 1;
            }
        }

        println!("Found {} in total.", count);
    } else {
        panic!("Expected an argument.");
    }
}
