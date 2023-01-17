use std::str;

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
    let mut count = 0;
    let mut previous = vec![];

    for dset in DSets::new(2, n) {
        for dsym in DSyms::new(&dset) {
            let code = dsym.to_binary().unwrap();
            count += 1;

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

    eprintln!("{} symbols generated.", count);
}
