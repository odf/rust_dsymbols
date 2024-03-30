use crate::derived::{as_partial_dsym, build_set, build_sym_using_vs, collapse, dual};
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


fn merge_all(ds: &PartialDSym) -> PartialDSym {
    let mut ds = as_partial_dsym(ds);

    for op in [
        merge_tiles, merge_facets, dual, merge_tiles, merge_facets, dual
    ] {
        ds = op(&ds);
    }

    ds
}


fn fix_local_1_vertex(ds: &PartialDSym) -> PartialDSym {
    let mut ds = as_partial_dsym(ds);

    for c in ds.orbit_reps([1, 2], 1..ds.size()) {
        if ds.m(1, 2, c) == Some(1) {
            let d = ds.op(0, ds.op(1, c).unwrap()).unwrap();
            let e = ds.op(1, ds.op(0, c).unwrap()).unwrap();
            let f = ds.op(3, d).unwrap();
            let g = ds.op(3, e).unwrap();

            let d1 = ds.op(1, d).unwrap();
            let e1 = ds.op(1, e).unwrap();
            let f1 = ds.op(1, f).unwrap();
            let g1 = ds.op(1, g).unwrap();

            let op = |i, c| {
                if i == 1 {
                    if c == d { Some(e1) }
                    else if c == e { Some(d1) }
                    else if c == f { Some(g1) }
                    else if c == g { Some(f1) }
                    else if c == d1 { Some(e) }
                    else if c == e1 { Some(d) }
                    else if c == f1 { Some(g) }
                    else if c == g1 { Some(f) }
                    else { ds.op(i, c) }
                } else {
                    ds.op(i, c)
                }
            };

            ds = build_sym_using_vs(
                build_set(ds.size(), ds.dim(), op),
                |i, d| ds.v(i, i + 1, d)
            );

            return collapse(&ds, ds.orbit([0, 1, 3], c), 3);
        }
    }

    ds
}


pub fn simplify<T: DSym>(ds: &T) -> PartialDSym {
    // TODO very basic first version
    let mut ds = as_partial_dsym(ds);

    for op in [merge_all, fix_local_1_vertex] {
        let old = ds;
        ds = op(&old);
        if ds.size() == old.size() {
            break;
        }
    }

    ds
}


#[cfg(test)]
mod test {
    use crate::delaney3d::pseudo_toroidal_cover;
    use crate::fpgroups::invariants::abelian_invariants;
    use crate::fundamental_group::fundamental_group;

    use super::*;

    #[test]
    fn test_merge_all() {
        let test = |s: &str| {
            let ds = s.parse::<PartialDSym>().unwrap();
            eprintln!("{ds}");
            let cov = pseudo_toroidal_cover(&ds).unwrap();
            let out = simplify(&cov);

            assert_eq!(out.orbit_reps([0, 1, 2], 1..out.size()).len(), 1);
            assert_eq!(out.orbit_reps([1, 2, 3], 1..out.size()).len(), 1);

            let reps = out.orbit_reps([2, 3], 1..out.size());
            assert!(reps.iter().all(|&d| out.r(2, 3, d) != Some(2)));
            let reps = out.orbit_reps([0, 1], 1..out.size());
            assert!(reps.iter().all(|&d| out.r(0, 1, d) != Some(2)));
            let reps = out.orbit_reps([1, 2], 1..out.size());
            assert!(reps.iter().all(|&d| out.r(1, 2, d) != Some(1)));

            let fg = fundamental_group(&out);
            let nr_gens = fg.gen_to_edge.len();
            let rels: Vec<_> = fg.relators.iter().cloned().collect();
            let inv = abelian_invariants(nr_gens, &rels);
            assert_eq!(inv, vec![0, 0, 0]);
        };

        test("<1.4:1 3:1,1,1,1:4,3,4>");
        test("<2.1:2 3:1 2,1 2,1 2,2:3 3,3 4,4>");
        test("<513.5:2 3:2,1 2,1 2,2:4,2 4,6>");
        test("<513.8:2 3:2,1 2,1 2,2:6,2 3,6>");
        test("<3.3:3 3:1 2 3,1 2 3,1 3,2 3:3 3 4,4 4,3>");
        test("<167.3:3 3:1 2 3,1 3,2 3,1 2 3:3 4,3,4 6>");
        test("<184.4:3 3:1 2 3,1 3,2 3,1 3:4 6,3,3>");
        test("<23.14:4 3:1 2 3 4,1 2 4,1 3 4,2 3 4:3 3 8,4 3,3 4>");
        test("<71.3:4 3:1 2 3 4,1 2 4,1 3 4,2 4:3 3 6,3 3,4>");
        test("<514.7:4 3:2 4,1 2 3 4,1 2 3 4,3 4:4 4,2 4 4 3,4 4>");
        test("<553.3:4 3:2 4,1 2 3 4,3 4,2 4:4 6,2 6,4>");
        test("<45.2:5 3:1 2 3 5,1 2 4 5,1 3 4 5,2 3 4 5:3 3 3,3 3 3,6 4 4>");
        test("<45.7:5 3:1 2 3 5,1 2 4 5,1 3 4 5,2 3 4 5:3 3 3,4 3 3,6 3 3>");
        test("<45.12:5 3:1 2 3 5,1 2 4 5,1 3 4 5,2 3 4 5:3 3 6,4 3 3,3 4 4>");
        test("<54.2:5 3:1 2 3 5,1 2 4 5,1 3 5,2 3 4 5:3 3 3,3 4,3 6>");
        test("<54.4:5 3:1 2 3 5,1 2 4 5,1 3 5,2 3 4 5:3 3 3,4 4,3 4>");
        test("<222.77:5 3:1 2 4 5,1 3 5,2 3 4 5,1 5 4:4 12,3 2,3 4>");
    }
}
