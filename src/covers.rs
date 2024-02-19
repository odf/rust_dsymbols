use std::collections::HashMap;

use crate::dsyms::*;

use crate::derived::cover;
use crate::fpgroups::cosets::{coset_table, CosetTable};
use crate::fpgroups::free_words::{FreeWord, Relator};
use crate::fundamental_group::fundamental_group;


fn trace_word(table: &CosetTable, start: usize, word: &FreeWord) -> usize {
    word.iter().fold(start, |row, g| table[row][g])
}


pub fn cover_for_table<T: DSym>(
    ds: &T,
    table: &CosetTable,
    edge_to_word: HashMap<(usize, usize), FreeWord>
)
    -> PartialDSym
{
    cover(
        ds,
        table.len(),
        |sheet, i, d| trace_word(table, sheet, &edge_to_word[&(d, i)])
    )
}


pub fn subgroup_cover<T: DSym>(ds: &T, subgens: &Vec<FreeWord>)
    -> PartialDSym
{
    let g = fundamental_group(ds);
    let table = coset_table(
        g.gen_to_edge.len(),
        &g.relators.iter().map(|w| Relator::from(w.clone())).collect(),
        subgens
    );
    cover_for_table(ds, &table, g.edge_to_word)
}


pub fn finite_universal_cover<T: DSym>(ds: &T) -> PartialDSym {
    subgroup_cover(ds, &vec![])
}


#[test]
fn test_finite_universal_cover_size() {
    use crate::dsets::DSet;

    let cov = |s: &str|
        finite_universal_cover(&s.parse::<PartialDSym>().unwrap());

    assert_eq!(cov("<1.1:1:1,1,1:3,3>").size(), 24);
    assert_eq!(cov("<1.1:1:1,1,1:4,3>").size(), 48);
    assert_eq!(cov("<1.1:1:1,1,1:3,5>").size(), 120);
    assert_eq!(cov("<1.1:1 3:1,1,1,1:3,3,3>").size(), 120);
    assert_eq!(cov("<1.1:1 3:1,1,1,1:4,3,3>").size(), 384);
}


#[test]
fn test_finite_universal_cover_round_trip() {
    use crate::derived::minimal_image;

    let check = |s: &str| {
        let ds = s.parse::<PartialDSym>().unwrap();
        assert_eq!(
            ds,
            minimal_image(&finite_universal_cover(&ds)).unwrap_or(ds.clone())
        );
    };

    check("<1.1:1:1,1,1:3,3>");
    check("<1.1:1:1,1,1:4,3>");
    check("<1.1:1:1,1,1:3,5>");
    check("<1.1:1 3:1,1,1,1:3,3,3>");
    check("<1.1:1 3:1,1,1,1:4,3,3>");
}


#[test]
fn test_finite_universal_cover_group() {
    let check = |s: &str| {
        let ds = s.parse::<PartialDSym>().unwrap();
        let g = fundamental_group(&finite_universal_cover(&ds));

        assert_eq!(g.gen_to_edge.len(), 0);
        assert_eq!(g.edge_to_word.len(), 0);
        assert_eq!(g.relators.len(), 0);
        assert_eq!(g.cones.len(), 0);
    };

    check("<1.1:1:1,1,1:3,3>");
    check("<1.1:1:1,1,1:4,3>");
    check("<1.1:1:1,1,1:3,5>");
    check("<1.1:1 3:1,1,1,1:3,3,3>");
    check("<1.1:1 3:1,1,1,1:4,3,3>");
}
