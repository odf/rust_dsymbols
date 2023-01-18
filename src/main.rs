use std::io::{Write, stdout};

use rust_dsymbols::dset_generators::DSets;
use rust_dsymbols::dsyms::DSym;
use rust_dsymbols::dsym_generators::DSyms;


fn main() {
    if let Some(arg) = std::env::args().nth(1) {
        let n: usize = arg.parse().unwrap();

        generate_binary(n);
    } else {
        panic!("Expected an argument.");
    }
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
