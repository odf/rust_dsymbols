use rust_dsymbols::dset_generators::DSets;


fn main() {
    let mut count = 0;

    for ds in DSets::new(2, 4) {
        println!("{}", ds);
        count += 1;
    }

    println!("Found {} in total.", count);
}
