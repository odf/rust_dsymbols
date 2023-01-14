use rust_dsymbols::dset_generators::DSets;
use rust_dsymbols::dsets::{DSet, OrientedCover};
use rust_dsymbols::dsyms::PartialDSym;


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
        println!();
    } else {
        panic!("Expected an argument.");
    }

    count = 0;
    for ds in DSets::new(2, 4) {
        count += 1;
        if count == 25 {
            let mut ds = PartialDSym::new(&ds);
            ds.set_v(0, 1, 2);
            ds.set_v(1, 1, 3);
            ds.set_v(1, 2, 4);
            println!("  Example symbol: {}", &ds);
            let cov = ds.oriented_cover().unwrap();
            println!("  Oriented cover: {}", &cov);
            println!("  Oriented cover is oriented: {}", cov.is_oriented());
            let map = cov.morphism(&ds, 1);
            println!("  Covering map: {:?}", &map);
        }
    }
    println!();
}
