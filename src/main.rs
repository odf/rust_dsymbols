use rust_dsymbols::dsets::{DSet, PartialDSet, SimpleDSet};


fn report<T>(ds: &T, title: &str)
    where T: DSet + std::fmt::Display
{
    println!("{}: {}", title, ds);
    println!("  DSet of size {} and dim {}", ds.size(), ds.dim());
    println!("  partial orientation: {:?}", ds.partial_orientation());
    println!("  complete: {}", ds.is_complete());
    println!("  loopless: {}", ds.is_loopless());
    println!("  weakly oriented: {}", ds.is_weakly_oriented());
    println!("  oriented: {}", ds.is_oriented());
    println!("  0,1-orbits: {:?}", &ds.orbits(0, 1));
    println!("  1,2-orbits: {:?}", &ds.orbits(1, 2));
    println!("  0,2-orbits: {:?}", &ds.orbits(0, 2));
    println!("  automorphisms: {:?}", &ds.automorphisms());
    println!();
}


fn main() {
    let mut tmp = PartialDSet::new(4, 2);
    tmp.set(0, 1, 2);
    tmp.set(0, 3, 4);

    report(&tmp, "Unfinished DSet:");

    tmp.set(1, 1, 4);
    tmp.set(1, 2, 3);
    tmp.set(2, 1, 1);
    tmp.set(2, 2, 2);
    tmp.set(2, 3, 4);

    let ds = SimpleDSet::from(tmp);
    report(&ds, "Finished DSet:");

    let cov = SimpleDSet::from(ds.oriented_cover());
    report(&cov, "Oriented cover of finished set:");
}
