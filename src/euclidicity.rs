use crate::derived::subsymbol;
use crate::dsets::DSet;
use crate::dsyms::PartialDSym;
use crate::fpgroups::cosets::coset_tables;
use crate::fpgroups::invariants::abelian_invariants;
use crate::fpgroups::stabilizer::stabilizer;
use crate::fundamental_group::{fundamental_group, FundamentalGroup};


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
