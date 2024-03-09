use std::collections::HashMap;
use std::iter::successors;

use crate::covers::cover_for_table;
use crate::derived::oriented_cover;
use crate::dsyms::{DSym, PartialDSym};
use crate::fpgroups::cosets::{
    core_table, coset_tables, intersection_table, CosetTable
};
use crate::fpgroups::free_words::FreeWord;
use crate::fpgroups::invariants::abelian_invariants;
use crate::fpgroups::stabilizer::stabilizer;
use crate::fundamental_group::{fundamental_group, FundamentalGroup};


fn point_groups() -> Vec<String> {
    ["z1", "z2", "z3", "z4", "v4", "s3", "z6", "d4", "d6", "a4", "s4"]
        .into_iter()
        .map(Into::into)
        .collect()
}


fn core_type_by_size(n: usize) -> String {
    match n {
         1 => "z1".into(),
         2 => "z2".into(),
         3 => "z3".into(),
         6 => "s3".into(),
         8 => "d4".into(),
        12 => "a4".into(),
        24 => "s4".into(),
         _ => panic!()
    }
}


fn is_fully_involutive(ct: &CosetTable) -> bool {
    ct.iter().all(|row| row.keys().all(|g| row[&g] == row[&-g]))
}


fn core_type(ct: &CosetTable) -> String {
    if ct.len() == 4 {
        if is_fully_involutive(ct) {
            "v4".into()
        } else {
            "z4".into()
        }
    } else {
        core_type_by_size(ct.len())
    }
}


fn degree(ct: &CosetTable, w: &FreeWord) -> usize {
    successors(
        Some((0, 0)),
        |&(i, row)| Some((i + 1, w.iter().fold(row, |a, g| ct[a][g])))
    )
        .skip(1)
        .skip_while(|&(_, row)| row != 0)
        .map(|(i, _)| i)
        .next()
        .unwrap()
}


fn flattens_all(ct: &CosetTable, cones: &Vec<(FreeWord, usize)>) -> bool {
    cones.iter().all(|(wd, deg)| degree(ct, wd) == *deg)
}


fn construct_candidates(fg: &FundamentalGroup)
    -> HashMap<String, Vec<CosetTable>>
{
    let nr_gens = fg.gen_to_edge.len();
    let rels: Vec<_> = fg.relators.iter().cloned().collect();
    let cones: Vec<_> = fg.cones.iter()
        .map(|(fw, &d)| (fw.clone(), d))
        .collect();

    let core_tables: Vec<_> = coset_tables(nr_gens, &rels, 4)
        .map(|ct| core_table(&ct))
        .collect();
    let cones2 = cones.iter().filter(|(_, d)| *d == 2).cloned().collect();
    let cones3 = cones.iter().filter(|(_, d)| *d == 3).cloned().collect();

    let mut result = HashMap::new();
    for p in point_groups() {
        result.insert(p, vec![]);
    }

    for table in core_tables.iter() {
        if flattens_all(table, &cones) {
            result.get_mut(&core_type(&table)).unwrap().push(table.clone());
        }
    }

    for ta in core_tables.iter().filter(|ct| flattens_all(ct, &cones3)) {
        for tb in core_tables.iter().filter(|ct| ct.len() == 2) {
            let tx = intersection_table(ta, tb);

            if flattens_all(&tx, &cones) {
                if ta.len() == 3 && tx.len() == 6 {
                    if flattens_all(&tb, &cones2) {
                        result.get_mut("z6").unwrap().push(tx);
                    }
                } else if ta.len() == 6 && tx.len() == 12 {
                    if !flattens_all(&tb, &cones2) {
                        result.get_mut("d6").unwrap().push(tx);
                    }
                }
            }
        }
    }

    result
}


pub fn pseudo_toroidal_cover<T: DSym>(ds: &T) -> Option<PartialDSym> {
    assert!(ds.dim() == 3, "must be three-dimensional");
    assert!(ds.is_complete(), "must be complete");

    for i in 0..ds.dim() {
        for d in ds.orbit_reps_2d(i, i + 1) {
            let v = ds.v(i, i + 1, d).unwrap();
            assert!(
                v <= 6 && v != 5,
                "violates the crystallographic restriction"
            );
        }
    }

    let ds = oriented_cover(ds);
    let fg = fundamental_group(&ds);
    let candidates = construct_candidates(&fg);

    for tp in point_groups() {
        for table in candidates[&tp].iter() {
            let rels: Vec<_> = fg.relators.iter().cloned().collect();
            let (sgens, srels) = stabilizer(0, &rels, table);
            let inv = abelian_invariants(sgens.len(), &srels);

            if inv == vec![0, 0, 0] {
                return Some(cover_for_table(&ds, table, &fg.edge_to_word));
            }
        }
    }

    None
}


#[test]
fn test_canonical() {
    let check_pseudo_toroidal = |src: &str, out: &str| {
        assert_eq!(
            src.parse::<PartialDSym>().map(|ds| pseudo_toroidal_cover(&ds)),
            out.parse::<PartialDSym>().map(|ds| Some(ds))
        );
    };

    check_pseudo_toroidal(
        "<1.1:1 3:1,1,1,1:4,3,4>",
        "<1.1:48 3:
        2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48,
        4 9 12 18 21 24 15 11 32 33 29 26 42 39 25 30 48 37 36 35 43 45 46 47,
        6 5 14 13 26 25 28 27 20 19 36 35 38 37 32 31 44 43 46 45 42 41 48 47,
        8 7 16 15 20 19 24 23 30 29 28 27 40 39 42 41 46 45 38 37 48 47 44 43:
        4 4 4 4 4 4,3 3 3 3 3 3 3 3,4 4 4 4 4 4>",
    );
    check_pseudo_toroidal(
        "<1.1:2 3:2,1 2,1 2,2:6,3 2,6>",
        "<1.1:96 3:
        2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48
        50 52 54 56 58 60 62 64 66 68 70 72 74 76 78 80 82 84 86 88 90 92 94 96,
        3 16 22 7 32 34 11 44 50 15 24 19 64 66 23 27 80 82 31 36 35 39 92 74
        43 52 47 96 58 51 55 88 70 59 94 63 68 67 71 86 75 90 79 84 83 87 91 95,
        7 8 17 18 19 20 39 40 45 46 55 56 33 34 63 64 69 70 75 76 41 42 87 88 65
        66 53 54 47 48 73 74 95 96 81 82 91 92 77 78 71 72 85 86 89 90 93 94,
        10 9 44 43 26 25 80 79 16 15 50 49 58 57 48 47 42 41 52 51 32 31 82 81
        78 77 84 83 54 53 88 87 66 65 92 91 64 63 94 93 96 95 90 89 76 75 86 85:
        6 6 6 6 6 6 6 6,3 2 2 2 3 2 2 3 2 3 3 2 2 3 2 2 3 3 2 2,6 6 6 6 6 6 6 6>",
    );
    check_pseudo_toroidal(
        "<1.1:2 3:1 2,1 2,1 2,2:3 3,3 4,4>",
        "<1.1:96 3:
        3 4 7 8 11 12 15 16 19 20 23 24 27 28 31 32 35 36 39 40 43 44 47 48 51
        52 55 56 59 60 63 64 67 68 71 72 75 76 79 80 83 84 87 88 91 92 95 96,
        7 8 17 18 19 20 31 32 37 38 43 44 45 46 59 60 65 66 71 72 73 74 39 40 87
        88 61 62 47 48 91 92 77 78 83 84 93 94 67 68 85 86 75 76 89 90 95 96,
        11 16 9 14 23 28 21 26 36 34 35 33 51 56 49 54 64 62 63 61 79 84 77 82
        59 76 57 74 75 60 73 58 83 80 81 78 88 86 87 85 91 96 89 94 95 92 93 90,
        2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50
        52 54 56 58 60 62 64 66 68 70 72 74 76 78 80 82 84 86 88 90 92 94 96:
        3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3,
        3 4 3 4 3 4 4 4 3 3 4 3 3 3,
        4 4 4 4 4 4 4 4 4 4 4 4>",
    );
}
