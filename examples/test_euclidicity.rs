use std::io::stdin;

use rust_dsymbols::euclidicity::is_euclidean;
use rust_dsymbols::dsyms::PartialDSym;
use rust_dsymbols::euclidicity::Euclidean::*;


fn main() {
    for line in stdin().lines() {
        if let Some(ds) = line.ok().and_then(|s|
            s.parse::<PartialDSym>().ok()
        ) {
            println!("{ds}");
            match is_euclidean(&ds) {
                Yes => println!("good"),
                No(s) => println!("bad: {s}"),
                Maybe(s, _) => println!("undecided: {s}"),
            }
        }
    }
}
