use std::collections::HashSet;

use crate::delaney3d::pseudo_toroidal_cover;
use crate::derived::{canonical, minimal_image, subsymbol};
use crate::dsets::DSet;
use crate::dsyms::{DSym, PartialDSym, SimpleDSym};
use crate::fpgroups::cosets::coset_tables;
use crate::fpgroups::invariants::abelian_invariants;
use crate::fpgroups::stabilizer::stabilizer;
use crate::fundamental_group::{fundamental_group, FundamentalGroup};
use crate::simplify::{find_small_tile_cut, simplify};


fn bad_subgroup_invariants(
    fg: &FundamentalGroup, index: usize, expected: Vec<usize>
) -> bool
{
    let nr_gens = fg.gen_to_edge.len();

    for table in coset_tables(nr_gens, &fg.relators, index) {
        let (sgens, srels) = stabilizer(0, fg.relators.clone(), &table);
        if abelian_invariants(sgens.len(), srels) != expected {
            return true;
        }
    }

    false
}


fn bad_connected_components(ds: &PartialDSym) -> bool {
    let mut seen_z3 = false;

    for d in ds.orbit_reps(0..=ds.dim(), 1..=ds.size()) {
        let component = subsymbol(ds, 0..ds.dim(), d);
        let fg = fundamental_group(&component);
        let nr_gens = fg.gen_to_edge.len();
        let invars = abelian_invariants(nr_gens, fg.relators.clone());

        if invars == [0, 0, 0] {
            if seen_z3 || bad_subgroup_invariants(&fg, 2, vec![0, 0, 0]) {
                return true;
            }
            seen_z3 = true;
        } else if invars == [] {
            if bad_subgroup_invariants(&fg, 5, vec![]) {
                return true;
            }
        } else {
            return true;
        }
    }

    false
}


fn bad_subgroup_count(fg: &FundamentalGroup, index: usize, expected: usize)
    -> bool
{
    let n = coset_tables(fg.gen_to_edge.len(), &fg.relators, index)
        .take(expected + 1)
        .count();

    n != expected
}


pub enum Euclidean {
    Yes,
    No(String),
    Maybe(String, SimpleDSym),
}


fn give_up(s: &str, ds: PartialDSym) -> Euclidean {
    Euclidean::Maybe(s.to_string(), SimpleDSym::from(ds))
}


fn fail(s: &str) -> Euclidean {
    Euclidean::No(s.to_string())
}


pub fn is_euclidean<T: DSym>(ds: &T) -> Euclidean {
    let good = HashSet::from([
        "<1.1:1 3:1,1,1,1:4,3,4>",
        "<1.1:8 3:2 4 6 8,6 3 5 7 8,2 7 8 5 6,4 3 6 8:3 4,5 3,3 4>",
    ]);

    if let Some(cov) = pseudo_toroidal_cover(ds) {
        if let Some(simp) = simplify(&cov) {
            let simp = canonical(&simp);
            let key = canonical(&minimal_image(&simp));

            if good.contains(&key.to_string()[..]) {
                Euclidean::Yes
            } else if !simp.is_connected() {
                if bad_connected_components(&simp) {
                    fail("cover is a non-trivial connected sum")
                } else {
                    give_up("cover is a (potentially trivial) connected sum", simp)
                }
            } else {
                let fg = fundamental_group(&simp);
                let invars = abelian_invariants(
                    fg.gen_to_edge.len(), fg.relators.clone()
                );

                if invars != [0, 0, 0] {
                    fail("cover has at least one handle")
                } else if fg.relators.is_empty() {
                    fail("cover has free fundamental group")
                } else if bad_subgroup_invariants(&fg, 2, vec![0, 0, 0]) {
                    fail("bad subgroups for cover")
                } else if bad_subgroup_count(&fg, 3, 21) {
                    fail("bad subgroup count for cover")
                } else if bad_subgroup_count(&fg, 4, 56) {
                    fail("bad subgroup count for cover")
                } else if let Some((d, cut)) = find_small_tile_cut(&simp) {
                    let msg = format!(
                        "unimplemented: split along {}-cut {:?}, glue at {}-face",
                        cut.len(), cut, simp.r(0, 1, d).unwrap()
                    );
                    give_up(&msg[..], simp)
                } else {
                    give_up("no decision found", key)
                }
            }
        } else {
            fail("cover is a lens space")
        }
    } else {
        fail("no pseudo-toroidal cover")
    }
}
