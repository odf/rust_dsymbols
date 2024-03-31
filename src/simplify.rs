use std::collections::{HashMap, HashSet};

use crate::derived::{as_partial_dsym, build_set, build_sym_using_vs, dual};
use crate::dsets::DSet;
use crate::dsyms::{DSym, PartialDSym};
use crate::fundamental_group::inner_edges;


fn collapse<I>(ds: &PartialDSym, remove: I, connector: usize)
    -> Option<PartialDSym>
    where I: IntoIterator<Item=usize>
{
    let remove: HashSet<_> = remove.into_iter().collect();

    if remove.is_empty() {
        None
    } else {
        let mut src2img = vec![0; ds.size() + 1];
        let mut img2src = vec![0; ds.size() + 1];
        let mut next = 1;
        for d in 1..=ds.size() {
            if !remove.contains(&d) {
                src2img[d] = next;
                img2src[next] = d;
                next += 1;
            }
        }

        let op = |i, d| {
            let mut e = ds.op(i, img2src[d]).unwrap();
            if i != connector {
                while src2img[e] == 0 {
                    e = ds.op(i, ds.op(connector, e).unwrap()).unwrap();
                }
            }
            Some(src2img[e])
        };

        Some(build_sym_using_vs(
            build_set(ds.size() - remove.len(), ds.dim(), op),
            |i, d| ds.v(i, i + 1, img2src[d])
        ))
    }
}


fn reglue<I>(ds: &PartialDSym, pairs: I, index: usize) -> Option<PartialDSym>
    where I: IntoIterator<Item=(usize, usize)>
{
    let paired: HashMap<_, _> = pairs.into_iter()
        .flat_map(|(d, e)| [(d, e), (e, d)])
        .collect();

    if paired.is_empty() {
        None
    } else {
        let op = |i, d| {
            if i == index && paired.contains_key(&d) {
                Some(*paired.get(&d).unwrap())
            } else {
                ds.op(i, d)
            }
        };

        Some(build_sym_using_vs(
            build_set(ds.size(), ds.dim(), op),
            |i, d| ds.v(i, i + 1, d)
        ))
    }
}


fn merge_tiles(ds: &PartialDSym) -> Option<PartialDSym> {
    let inner = inner_edges(ds);
    let junk = inner.iter().cloned()
        .filter(|&(_, i)| i == 3)
        .flat_map(|(d, _)| ds.orbit([3], d));
    collapse(ds, junk, 3)
}


fn merge_facets(ds: &PartialDSym) -> Option<PartialDSym> {
    let reps = ds.orbit_reps([2, 3], 1..ds.size());
    let junk = reps.iter().cloned()
        .filter(|&d| ds.r(2, 3, d) == Some(2))
        .flat_map(|d| ds.orbit([2, 3], d));
    collapse(ds, junk, 2)
}


fn wrap_dual(ds: &PartialDSym) -> Option<PartialDSym> {
    Some(dual(ds))
}


fn merge_all(ds: &PartialDSym) -> Option<PartialDSym> {
    let mut ds = as_partial_dsym(ds);

    for op in [
        merge_tiles, merge_facets, wrap_dual,
        merge_tiles, merge_facets, wrap_dual
    ] {
        if let Some(out) = op(&ds) {
            ds = out;
        }
    }

    Some(ds)
}


fn fix_local_1_vertex(ds: &PartialDSym) -> Option<PartialDSym> {
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

            let tmp = reglue(
                &ds, [(d, e1), (e, d1), (f, g1), (g, f1)], 1
            ).unwrap();

            return collapse(&tmp, tmp.orbit([0, 1, 3], c), 3);
        }
    }

    None
}


pub fn simplify<T: DSym>(ds: &T) -> PartialDSym {
    // TODO very basic first version
    let ds = as_partial_dsym(ds);
    let ds = merge_all(&ds).or(Some(ds)).unwrap();

    for op in [fix_local_1_vertex] {
        if let Some(out) = op(&ds) {
            return out;
        }
    }

    ds
}


#[cfg(test)]
mod test {
    use crate::delaney2d::toroidal_cover;
    use crate::delaney3d::pseudo_toroidal_cover;
    use crate::derived::minimal_image;
    use crate::fpgroups::invariants::abelian_invariants;
    use crate::fundamental_group::fundamental_group;

    use super::*;


    #[test]
    fn test_collapse_2d() {
        let dsym = |s: &str| s.parse::<PartialDSym>().unwrap();

        let cov = toroidal_cover(&dsym("<1.1:1:1,1,1:3,6>"));
        let out = minimal_image(
            &collapse(&cov, cov.orbit([0, 2], 1), 2).unwrap()
        );

        assert_eq!(out, dsym("<1.1:1:1,1,1:4,4>"));
    }


    #[test]
    fn test_collapse_3d() {
        let dsym = |s: &str| s.parse::<PartialDSym>().unwrap();

        let ds = dsym("<1.1:4 3:1 2 3 4,1 2 4,1 3 4,2 3 4:3 3 8,4 3,3 4>");
        let cov = dual(&pseudo_toroidal_cover(&ds).unwrap());
        let remove = (1..=cov.size()).filter(|&d| cov.m(0, 1, d) == Some(3));
        let out = minimal_image(
            &collapse(&cov, remove, 3).unwrap()
        );

        assert_eq!(out, dsym("<1.1:1 3:1,1,1,1:4,3,4>"));
    }


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
