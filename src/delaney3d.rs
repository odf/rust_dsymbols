use std::iter::successors;

use crate::fpgroups::cosets::CosetTable;
use crate::fpgroups::free_words::FreeWord;


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


fn flattens_all(ct: &CosetTable, cones: Vec<(&FreeWord, usize)>) -> bool {
    cones.iter().all(|&(wd, deg)| degree(ct, wd) == deg)
}
