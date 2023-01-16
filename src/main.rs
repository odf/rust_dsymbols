use std::str;

use rust_dsymbols::dset_generators::DSets;
use rust_dsymbols::dsets::DSet;
use rust_dsymbols::dsym_generators::DSyms;


fn main() {
    if let Some(arg) = std::env::args().nth(1) {
        let n: usize = arg.parse().unwrap();
        let mut previous = vec![];

        for dset in DSets::new(2, n) {
            for dsym in DSyms::new(&dset) {
                let code = dsym.to_binary().unwrap();

                for i in 0..=code.len() {
                    if i >= code.len()
                        || i >= previous.len()
                        || code[i] != previous[i]
                    {
                        print!(
                            "{}{}",
                            str::from_utf8(&[(i + 1) as u8]).unwrap(),
                            str::from_utf8(&code[i..]).unwrap()
                        );

                        previous = code;
                        break;
                    }
                }
            }
        }
    } else {
        panic!("Expected an argument.");
    }
}
