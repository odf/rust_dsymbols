use crate::derived::{as_partial_dsym, collapse, dual};
use crate::dsets::DSet;
use crate::dsyms::{DSym, PartialDSym};
use crate::fundamental_group::inner_edges;


fn merge_tiles(ds: &PartialDSym) -> PartialDSym {
    let inner = inner_edges(ds);
    let junk = inner.iter().cloned()
        .filter(|&(_, i)| i == 3)
        .flat_map(|(d, _)| ds.orbit([3], d));
    collapse(ds, junk, 3)
}


fn merge_facets(ds: &PartialDSym) -> PartialDSym {
    let reps = ds.orbit_reps([2, 3], 1..ds.size());
    let junk = reps.iter().cloned()
        .filter(|&d| ds.r(2, 3, d) == Some(2))
        .flat_map(|d| ds.orbit([2, 3], d));
    collapse(ds, junk, 2)
}


pub fn merge_all<T: DSym>(ds: &T) -> PartialDSym {
    let mut ds = as_partial_dsym(ds);

    for op in [
        merge_tiles, merge_facets, dual, merge_tiles, merge_facets, dual
    ] {
        ds = op(&ds);
    }

    ds
}


#[cfg(test)]
mod Test {
    use crate::delaney3d::pseudo_toroidal_cover;
    use crate::derived::minimal_image;

    use super::*;

    #[test]
    fn test_merge_all() {
        let test = |s: &str| {
            let ds = s.parse::<PartialDSym>().unwrap();
            eprintln!("{ds}");
            let cov = pseudo_toroidal_cover(&ds).unwrap();
            let out = minimal_image(&merge_all(&cov));
            assert_eq!(out.to_string(), "<1.1:1 3:1,1,1,1:4,3,4>");
        };

        test("<1.1:1 3:1,1,1,1:4,3,4>");
        test("<1.1:3 3:1 2 3,1 3,2 3,1 2 3:3 4,3,4 6>");
        test("<1.1:4 3:1 2 3 4,1 2 4,1 3 4,2 3 4:3 3 8,4 3,3 4>");
    }
}
