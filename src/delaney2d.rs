use std::collections::HashSet;

use num_rational::Rational64;
use num_traits::{Signed, Zero};

use crate::dsets::*;
use crate::dsyms::*;
use crate::covers::covers;
use crate::derived::oriented_cover;


fn orbit_types_2d<T: DSym>(ds: &T) -> Vec<(usize, bool)> {
    let mut result = vec![];

    for i in 0..ds.dim() {
        for j in (i + 1)..=ds.dim() {
            for d in ds.orbit_reps_2d(i, j) {
                let loopless = ds.orbit([i, j].into_iter(), d).iter()
                    .all(|&e|
                        ds.op(i, e) != Some(e) && ds.op(j, e) != Some(e)
                    );
                result.push((ds.v(i, j, d).unwrap(), loopless))
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
            Rational64::new(if loopless { 2 } else { 1 }, v as i64)
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
            .filter(|&v| v > 1)
            .collect();

        match cones.len() {
            1 => false,
            2 => cones[0] == cones[1],
            _ => true,
        }
    }
}


pub fn toroidal_cover<T: DSym>(ds: &T) -> PartialDSym {
    assert!(ds.dim() == 2, "must be two-dimensional");
    assert!(is_euclidean(ds), "must be euclidean");

    let ds = &oriented_cover(ds);
    let degree = orbit_types_2d(ds).iter()
        .map(|&(v, _)| v)
        .max()
        .unwrap_or(1);

    for cov in covers(ds, degree) {
        if orbit_types_2d(&cov).iter().all(|&(v, _)| v == 1) {
            return cov;
        }
    }

    panic!("symbol is 2d euclidean, should have found a toroidal cover");
}


fn opposite<T: DSym>(ds: &T, i: usize, j: usize, d: usize) -> (usize, usize) {
    let mut k = i;
    let mut e = d;

    while ds.op(k, e).unwrap_or(e) != e {
        e = ds.op(k, e).unwrap();
        k = i + j - k;
    }

    (k, e)
}


fn trace_boundary<T: DSym>(ds: &T) -> Vec<Vec<usize>> {
    let ori = ds.partial_orientation();
    let mut result = vec![];
    let mut seen: HashSet<(usize, usize)> = HashSet::new();

    for i in 0..=ds.dim() {
        for d in 1..=ds.size() {
            if ds.op(i, d) != Some(d) || seen.contains(&(i, d)) {
                continue;
            }

            let mut corners = vec![];
            let mut j = i;
            let mut k = match ori[d] {
                Sign::PLUS => (i + 1) % 3,
                Sign::MINUS => (i + 2) % 3,
                Sign::ZERO => panic!("orientation should be total"),
            };
            let mut e = d;
            let mut nu = k;

            loop {
                let v = ds.v(j, k, e).unwrap();
                if v > 1 {
                    corners.push(v);
                }

                seen.insert((j, e));
                let (nu, e) = opposite(ds, k, j, e);
                k = 3 - j - k;
                j = nu;

                if seen.contains(&(j, e)) {
                    break;
                }
            }

            result.push(corners);
        }
    }
    result
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


#[test]
fn test_toroidal_cover() {
    let check = |s: &str| {
        let ds = s.parse::<PartialDSym>().unwrap();
        let cov = toroidal_cover(&ds);

        assert!(is_euclidean(&cov));
        assert!(cov.is_complete());
        assert!(cov.is_connected());
        assert!(cov.is_oriented());
        assert!(orbit_types_2d(&cov).iter().all(|&(v, _)| v == 1));
        assert!(cov.morphism(&ds, 1).is_some());
    };

    check("<1.1:1:1,1,1:3,6>");
    check("<1.1:3:1 2 3,1 3,2 3:4 8,3>");
    check("<1.1:8:2 4 6 8,8 3 5 7,6 5 8 7:4,4>");
    check("<1.1:8:2 4 6 8,8 3 5 7,5 6 8 7:4,4>");
    check("<1.1:5:2 4 5,1 2 3 5,3 4 5:8 3,8 3>");
    check("<1.1:4:2 4,1 3 4,3 4:4,4>");
    check("<1.1:4:1 3 4,2 4,4 2 3:4,4>")
}
