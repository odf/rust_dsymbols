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
                Yes => println!("{i} good"),
                No(s) => println!("{i} bad: {s}"),
                Maybe(s, _) => println!("{i} undecided: {s}"),
            }
        }
    }
}
