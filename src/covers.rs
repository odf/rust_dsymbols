use std::collections::BTreeMap;

use crate::dsyms::*;

use crate::derived::cover;
use crate::fpgroups::cosets::{coset_table, coset_tables, CosetTable};
use crate::fpgroups::free_words::FreeWord;
use crate::fundamental_group::fundamental_group;


fn trace_word(table: &CosetTable, start: usize, word: &FreeWord) -> usize {
    word.iter().fold(start, |row, g| table.get(row, *g).unwrap())
}


pub fn cover_for_table<T: DSym>(
    ds: &T,
    table: &CosetTable,
    edge_to_word: &BTreeMap<(usize, usize), FreeWord>
)
    -> PartialDSym
{
    cover(
        ds,
        table.len(),
        |sheet, i, d|
            trace_word(
                table,
                sheet,
                &edge_to_word.get(&(d, i)).unwrap_or(&FreeWord::empty())
            )
    )
}


pub fn subgroup_cover<T: DSym>(ds: &T, subgens: &Vec<FreeWord>)
    -> PartialDSym
{
    let g = fundamental_group(ds);
    let table = coset_table(g.nr_generators(), &g.relators, subgens);
    cover_for_table(ds, &table, &g.edge_to_word)
}


pub fn finite_universal_cover<T: DSym>(ds: &T) -> PartialDSym {
    subgroup_cover(ds, &vec![])
}


pub fn covers<T: DSym>(ds: &T, max_deg: usize) -> Vec<PartialDSym> {
    let mut result = Vec::new();
    let g = fundamental_group(ds);

    for table in coset_tables(g.nr_generators(), &g.relators, max_deg) {
        let cov = cover_for_table(ds, &table, &g.edge_to_word);
        result.push(cov);
    }

    result
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
    assert_eq!(cov("<1.1:2 3:2,2,2,2:4,3,3>").size(), 384);
}


#[test]
fn test_finite_universal_cover_round_trip() {
    use crate::derived::minimal_image;

    let check = |s: &str| {
        let ds = s.parse::<PartialDSym>().unwrap();
        assert_eq!(ds, minimal_image(&finite_universal_cover(&ds)));
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


#[test]
fn test_number_of_covers() {
    let covers = |s: &str, n: usize|
        covers(&s.parse::<PartialDSym>().unwrap(), n);

    assert_eq!(covers("<1.1:1 1:1,1:3>", 3).len(), 3);
    assert_eq!(covers("<1.1:2 1:2,2:3>", 3).len(), 2);

    assert_eq!(covers("<1.1:1:1,1,1:3,3>", 24).len(), 11);
}
