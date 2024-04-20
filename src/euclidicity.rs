use crate::derived::subsymbol;
use crate::dsets::DSet;
use crate::dsyms::PartialDSym;
use crate::fpgroups::cosets::coset_tables;
use crate::fpgroups::invariants::abelian_invariants;
use crate::fpgroups::stabilizer::stabilizer;
use crate::fundamental_group::{fundamental_group, FundamentalGroup};


fn good_subgroup_invariants(
    fg: &FundamentalGroup, index: usize, expected: Vec<usize>
) -> bool
{
    let nr_gens = fg.gen_to_edge.len();
    let rels: Vec<_> = fg.relators.iter().cloned().collect();

    for table in coset_tables(nr_gens, &rels, index) {
        let (sgens, srels) = stabilizer(0, fg.relators.clone(), &table);
        if abelian_invariants(sgens.len(), srels) != expected {
            return false;
        }
    }

    true
}


fn good_connected_components(ds: &PartialDSym) -> bool {
    let mut seen_z3 = false;

    for d in ds.orbit_reps(0..=ds.dim(), 1..=ds.size()) {
        let component = subsymbol(ds, 0..ds.dim(), d);
        let fg = fundamental_group(&component);
        let nr_gens = fg.gen_to_edge.len();
        let invars = abelian_invariants(nr_gens, fg.relators.clone());

        if invars == [0, 0, 0] {
            if seen_z3 || !good_subgroup_invariants(&fg, 2, vec![0, 0, 0]) {
                return false;
            }
            seen_z3 = true;
        } else if invars == [] {
            if !good_subgroup_invariants(&fg, 5, vec![]) {
                return false;
            }
        } else {
            return false;
        }
    }

    true
}
