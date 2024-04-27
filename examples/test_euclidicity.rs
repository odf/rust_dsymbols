use std::io::stdin;

use rust_dsymbols::euclidicity::is_euclidean;
use rust_dsymbols::dsyms::PartialDSym;
use rust_dsymbols::euclidicity::Euclidean::*;


fn main() {
    let mut i = 0;
    for line in stdin().lines() {
        i += 1;

        if let Some(ds) = line.ok().and_then(|s|
            s.parse::<PartialDSym>().ok()
        ) {
            match is_euclidean(&ds) {
                Yes => {
                    println!(
                        "#Symbol {i} is good: simplified cover recognized"
                    );
                    println!("{ds}");
                },
                No(s) => println!("#Symbol {i} is bad: {s}"),
                Maybe(s, _) => println!("#Symbol {i} is undecided: {s}"),
            }
        }
    }
}
