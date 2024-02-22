use num_rational::Rational64;

use crate::dsyms::*;


fn orbit_types_2d<T: DSym>(ds: &T) -> Vec<(Option<usize>, bool)> {
    let mut result = vec![];

    for i in 0..ds.dim() {
        for j in (i + 1)..=ds.dim() {
            for d in ds.orbit_reps_2d(i, j) {
                let loopless = ds.orbit([i, j].into_iter(), d).iter()
                    .all(|&e|
                        ds.op(i, e) != Some(e) && ds.op(j, e) != Some(e)
                    );
                result.push((ds.v(i, j, d), loopless))
            }
        }
    }

    result
}


pub fn curvature<T: DSym>(ds: &T) -> Rational64 {
    assert!(ds.dim() == 2, "must be two-dimensional");
    assert!(ds.is_complete(), "must be complete");

    let s: Rational64 = orbit_types_2d(ds).iter()
        .map(|&(v, loopless)|
            Rational64::new(if loopless { 2 } else { 1 }, v.unwrap() as i64)
        )
        .sum();

    s - Rational64::from(ds.size() as i64)
}


#[test]
fn test_curvature() {
    let curv = |s: &str| curvature(&s.parse::<PartialDSym>().unwrap());

    assert_eq!(curv("<1.1:1:1,1,1:3,3>"), Rational64::new(1, 6));
    assert_eq!(curv("<1.1:1:1,1,1:3,4>"), Rational64::new(1, 12));
    assert_eq!(curv("<1.1:1:1,1,1:3,5>"), Rational64::new(1, 30));
    assert_eq!(curv("<1.1:1:1,1,1:3,6>"), Rational64::new(0, 1));
    assert_eq!(curv("<1.1:1:1,1,1:3,7>"), Rational64::new(-1, 42));

    assert_eq!(
        curv("
            <1.1:12:
            1 3 5 8 10 11 12,2 3 6 7 10 12 11,1 4 5 9 7 11 10 12:
            3 3 3,6 3 3>
        "),
        Rational64::new(1, 1));

        assert_eq!(
            curv("
                <1.1:12:
                1 3 4 5 7 8 9 11 12,2 4 6 8 10 12,12 2 3 5 6 7 9 10 11:
                8 12 16,8 12 16>
            "),
            Rational64::new(-23, 6));
}
