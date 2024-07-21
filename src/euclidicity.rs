use std::collections::BTreeSet;
use once_cell::sync::Lazy;

use crate::delaney3d::{orbifold_graph, pseudo_toroidal_cover};
use crate::derived::{canonical, minimal_image, subsymbol};
use crate::dsets::DSet;
use crate::dsyms::{DSym, PartialDSym, SimpleDSym};
use crate::fpgroups::cosets::coset_tables;
use crate::fpgroups::invariants::abelian_invariants;
use crate::fpgroups::stabilizer::stabilizer;
use crate::fundamental_group::{fundamental_group, FundamentalGroup};
use crate::simplify::simplify;


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


static INVARIANTS: Lazy<BTreeSet<String>> = Lazy::new(|| {
    include_str!(
        concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/data/euclideanInvariants.data"
        )
    )
        .split_whitespace()
        .into_iter()
        .filter(|s| s.len() > 0 && !s.starts_with("#"))
        .map(|s| s.to_string())
        .collect()
});


fn orbifold_invariant<T: DSym>(ds: &T) -> String {
    let (labels, edges) = orbifold_graph(ds);
    let fg = fundamental_group(ds);
    let nr_gens = fg.gen_to_edge.len();
    let invars = abelian_invariants(nr_gens, fg.relators.clone());

    let mut parts = vec![labels.len().to_string()];
    parts.extend(labels);
    parts.push(if ds.is_oriented() {
        "2".to_string()
    } else if ds.is_weakly_oriented() {
        "1".to_string()
    } else {
        "0".to_string()
    });
    parts.push(edges.len().to_string());
    parts.push(invars.len().to_string());
    parts.extend(invars.iter().map(|n| n.to_string()));
    parts.push("".to_string());

    parts.join("/")
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
    if !INVARIANTS.contains(&orbifold_invariant(ds)) {
        fail("orbifold invariants do not match")
    } else if let Some(cov) = pseudo_toroidal_cover(ds) {
        if let Some(simp) = simplify(&cov) {
            let simp = canonical(&simp);
            let key = canonical(&minimal_image(&simp));

            if key.to_string() == "<1.1:1 3:1,1,1,1:4,3,4>" {
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
                } else if bad_subgroup_count(&fg, 2, 8) {
                    fail("bad subgroup count for cover")
                } else if bad_subgroup_invariants(&fg, 2, vec![0, 0, 0]) {
                    fail("bad subgroups for cover")
                //} else if bad_subgroup_count(&fg, 3, 21) {
                //    fail("bad subgroup count for cover")
                //} else if bad_subgroup_count(&fg, 4, 56) {
                //    fail("bad subgroup count for cover")
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
