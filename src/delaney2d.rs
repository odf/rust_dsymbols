use num_rational::Rational64;
use num_traits::{Signed, Zero};

use crate::derived::oriented_cover;
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


pub fn is_euclidean<T: DSym>(ds: &T) -> bool {
    curvature(ds).is_zero()
}


pub fn is_hyperbolic<T: DSym>(ds: &T) -> bool {
    curvature(ds).is_negative()
}


pub fn is_spherical<T: DSym>(ds: &T) -> bool {
    if !curvature(ds).is_positive() {
        false
    } else {
        let dso = oriented_cover(ds);
        let cones: Vec<_> = orbit_types_2d(&dso).iter()
            .map(|&(v, _)| v)
            .filter(|&v| v.is_some_and(|v| v > 1))
            .collect();

        match cones.len() {
            1 => false,
            2 => cones[0] == cones[1],
            _ => true,
        }
    }
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


#[test]
fn test_is_spherical() {
    let passes = |s: &str| is_spherical(&s.parse::<PartialDSym>().unwrap());

    assert!(passes("<1.1:1:1,1,1:3,3>"));
    assert!(passes("<1.1:1:1,1,1:3,4>"));
    assert!(passes("<1.1:1:1,1,1:3,5>"));
    assert!(!passes("<1.1:1:1,1,1:3,6>"));
    assert!(!passes("<1.1:1:1,1,1:3,7>"));

    assert!(passes("<1.1:2:2,1 2,1 2:2,1 1>"));
    assert!(!passes("<1.1:2:2,1 2,1 2:2,1 3>"));
    assert!(passes("<1.1:2:2,1 2,1 2:2,3 3>"));
    assert!(!passes("<1.1:2:2,1 2,1 2:2,3 4>"));
}
