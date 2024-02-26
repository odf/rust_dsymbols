use std::collections::HashSet;
use std::iter::successors;

use num_rational::Rational64;
use num_traits::{Signed, Zero};

use crate::covers::covers;
use crate::derived::*;
use crate::dsets::*;
use crate::dsyms::*;


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


fn best_cyclic(corners: &Vec<usize>) -> Vec<usize> {
    (0..corners.len())
        .map(|i|
            corners[i..].iter().chain(corners[..i].iter()).cloned().collect()
        )
        .max()
        .unwrap_or(vec![])
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
            let mut nu;

            loop {
                let v = ds.v(j, k, e).unwrap();
                if v > 1 {
                    corners.push(v);
                }

                seen.insert((j, e));
                (nu, e) = opposite(ds, k, j, e);
                k = 3 - j - k;
                j = nu;

                if seen.contains(&(j, e)) {
                    break;
                }
            }

            result.push(best_cyclic(&corners));
        }
    }

    result.sort();
    result.reverse();
    result
}


fn euler_characteristic<T: DSym>(ds: &T) -> isize {
    let nr_orbits = |i, j| ds.orbit_reps_2d(i, j).len();
    let nr_loops = |i|
        (1..=ds.size()).filter(|&d| ds.op(i, d) == Some(d)).count();

    let nf = ds.size();
    let ne = (3 * nf + nr_loops(0) + nr_loops(1) + nr_loops(2)) / 2;
    let nv = nr_orbits(0, 1) + nr_orbits(0, 2) + nr_orbits(1, 2);

    (nf + nv) as isize - ne as isize
}


fn degree_list_as_string(vs: Vec<usize>) -> String {
    vs.iter()
        .map(|v| if *v < 10 { v.to_string() } else { format!("({})", v) })
        .collect::<Vec<_>>()
        .join("")
}


fn cone_degrees<T: DSym>(ds: &T) -> Vec<usize> {
    orbit_types_2d(ds).iter()
        .filter(|(v, c)| *c && *v > 1)
        .map(|&(v, _)| v)
        .collect()
}


pub fn orbifold_symbol<T: DSym>(ds: &T) -> String {
    assert!(ds.dim() == 2, "must be two-dimensional");
    assert!(ds.is_complete(), "must be complete");

    let bnd_comps = trace_boundary(ds);
    let chi = euler_characteristic(ds) + bnd_comps.len() as isize;
    let x = 2 - chi;

    let mut cones = cone_degrees(ds);
    cones.sort();
    cones.reverse();

    let mut parts = vec![degree_list_as_string(cones)];

    for corners in bnd_comps {
        parts.push("*".to_string());
        parts.push(degree_list_as_string(corners));
    }

    if ds.is_weakly_oriented() {
        parts.push(vec!["o"; (x as usize) / 2].join(""));
    } else {
        parts.push(vec!["x"; x as usize].join(""));
    }

    let s = parts.join("");

    (match &s[..] {
        "x" => "1x",
        "*" => "1*",
        "" => "1",
        s => s
    }).to_string()
}


fn split_along<T: DSym>(ds: &T, cut_elms: &Vec<usize>, cut_idx: usize)
    -> (PartialDSym, PartialDSym)
{
    let in_cut: HashSet<_> = cut_elms.iter()
        .flat_map(|&d| [(cut_idx, d), (cut_idx, ds.op(cut_idx, d).unwrap())])
        .collect();

    let op = |i, d| {
        if in_cut.contains(&(i, d)) { Some(d) } else { ds.op(i, d) }
    };
    let v = |i, d| ds.v(i, i + 1, d);
    let tmp = build_sym_using_vs(build_set(ds.size(), ds.dim(), op), v);

    (
        subsymbol(&tmp, 0..=ds.dim(), cut_elms[0]),
        subsymbol(&tmp, 0..=ds.dim(), ds.op(1, cut_elms[0]).unwrap()),
    )
}


fn check_cones(cones: Vec<usize>, allow_2_cone: bool) -> bool {
    cones == vec![] || (allow_2_cone && cones == vec![2])
}


fn bad_spherical_cut<T: DSym>(
    ds: &T,
    cut_elms: &Vec<usize>,
    patch: &PartialDSym,
    rest: &PartialDSym,
    allow_2_cone: bool
)
    -> bool
{
    if euler_characteristic(ds) <= 0 {
        false
    } else if patch.size() == ds.size() {
        let mut vs = vec![
            ds.v(0, 1, cut_elms[0]),
            ds.v(1, 2, cut_elms[0])
        ];
        if cut_elms.len() > 2 {
            vs.push(ds.v(1, 2, cut_elms[1]));
        }
        let vs = vs.into_iter().flatten().filter(|&v| v > 1).collect();

        check_cones(vs, allow_2_cone)
    } else if patch.size() == ds.size() - cut_elms.len() {
        cut_elms.iter().all(|&d|
            ds.v(0, 1, d) == Some(1) && ds.v(1, 2, d) == Some(1)
        )
        &&
        check_cones(cone_degrees(rest), allow_2_cone) 
    } else {
        false
    }
}


fn cuts_off_disk<T: DSym>(
    ds: &T,
    cut_elms: &Vec<usize>,
    cut_idx: usize,
    allow_2_cone: bool
)
    -> bool
{
    let (patch, rest) = split_along(ds, cut_elms, cut_idx);

    patch.size() != cut_elms.len() &&
    !bad_spherical_cut(ds, cut_elms, &patch, &rest, allow_2_cone) &&
    patch.is_weakly_oriented() &&
    euler_characteristic(&patch) == 1 &&
    check_cones(cone_degrees(&patch), allow_2_cone)
}


pub fn is_pseudo_convex<T: DSym>(ds: &T) -> bool {
    assert!(ds.dim() == 2, "must be two-dimensional");
    assert!(ds.is_complete(), "must be complete");

    let ds = &oriented_cover(ds);
    let step = |i, j, d| ds.op(j, d).and_then(|e| ds.op(i, e));
    let ori = ds.partial_orientation();

    for a1 in (1..=ds.size()).filter(|&d| ori[d] == Sign::PLUS) {
        let mut seen1 = HashSet::from([a1]);

        for a2 in successors(ds.op(0, a1), |&d| step(0, 1, d)) {
            if seen1.contains(&a2) {
                break;
            }
            let mut seen2 = seen1.clone();
            seen2.insert(a2);
            seen1.extend([a2, ds.op(1, a2).unwrap()]);

            for b2 in successors(ds.op(2, a2), |&d| step(2, 1, d)) {
                if seen2.contains(&b2) {
                    if b2 == a1 && cuts_off_disk(ds, &vec![a1, a2], 1, true) {
                        return false;
                    }
                    break;
                }
                let mut seen3 = seen2.clone();
                seen3.insert(b2);
                seen2.extend([b2, ds.op(1, b2).unwrap()]);

                for b1 in successors(ds.op(0, b2), |&d| step(0, 1, d)) {
                    if seen3.contains(&b1) {
                        break;
                    }
                    seen3.extend([b1, ds.op(1, b1).unwrap()]);
                    let mut seen4 = seen3.clone();

                    for t in successors(ds.op(2, b1), |&d| step(2, 1, d)) {
                        if seen4.contains(&t) {
                            if t == a1 && cuts_off_disk(
                                ds, &vec![a1, a2, b2, b1], 1, false
                            ) {
                                return false;
                            }
                            break;
                        }
                        seen4.extend([t, ds.op(1, t).unwrap()]);
                    }
                }
            }
        }
    }

    true
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


#[test]
fn test_orbifold_symbol() {
    let orb = |s: &str| orbifold_symbol(&s.parse::<PartialDSym>().unwrap());

    assert_eq!(orb("<1.1:3:1 2 3,1 3,2 3:4 8,3>"), "*442");
    assert_eq!(orb("<1.1:1:1,1,1:5,3>"), "*532");
    assert_eq!(orb("<1.1:1:1,1,1:6,3>"), "*632");
    assert_eq!(orb("<1.1:1:1,1,1:7,3>"), "*732");
    assert_eq!(orb("<1.1:1:1,1,1:15,3>"), "*(15)32");
    assert_eq!(orb("<1.1:2:2,1 2,1 2:2,4 4>"), "*44");
    assert_eq!(orb("<1.1:2:2,1 2,1 2:2,4 5>"), "*54");
    assert_eq!(orb("<1.1:8:2 4 6 8,8 3 5 7,6 5 8 7:4,4>"), "o");
    assert_eq!(orb("<1.1:8:2 4 6 8,8 3 5 7,5 6 8 7:4,4>"), "xx");
    assert_eq!(orb("<1.1:5:2 4 5,1 2 3 5,3 4 5:8 3,8 3>"), "*442");
    assert_eq!(orb("<1.1:4:2 4,1 3 4,3 4:4,4>"), "*x");
    assert_eq!(orb("<1.1:4:1 3 4,2 4,4 2 3:4,4>"), "**");
    assert_eq!(orb("<1.1:4:1 3 4,2 4,4 2 3:8,12>"), "*3*2");
    assert_eq!(
        orb("
            <1.1:12:
            1 3 4 5 7 8 9 11 12,2 4 6 8 10 12,12 2 3 5 6 7 9 10 11:
            8 12 16,8 12 16>
        "),
        "*432*423"
    );
    assert_eq!(
        orb("
            <1.1:12:
            1 3 5 8 10 11 12,2 3 6 7 10 12 11,1 4 5 9 7 11 10 12:
            3 3 3,6 3 3>
        "),
        "*22"
    );
    assert_eq!(
        orb("
            <1.1:16:
            2 4 6 8 10 12 14 16,16 3 5 7 9 11 13 15,6 5 8 7 14 13 16 15:
            8,8>
        "),
        "oo"
    );
    assert_eq!(
        orb("
            <1.1:16:
            2 4 6 8 10 12 14 16,16 3 5 7 9 11 13 15,5 6 8 7 14 13 16 15:
            8,8>
        "),
        "xxxx"
    );
}


#[test]
fn test_is_pseudo_convex() {
    let passes = |s: &str| is_pseudo_convex(&s.parse::<PartialDSym>().unwrap());

    // has handles, cross-caps or multiple boundaries
    assert!(passes("<1.1:8:2 4 6 8,8 3 5 7,6 5 8 7:4,4>"));
    assert!(passes("<1.1:8:2 4 6 8,8 3 5 7,5 6 8 7:4,4>"));
    assert!(passes("<1.1:4:2 4,1 3 4,3 4:4,4>"));
    assert!(passes("
        <1.1:12:
        1 3 4 5 7 8 9 11 12,2 4 6 8 10 12,12 2 3 5 6 7 9 10 11:
        8 12 16,8 12 16>
    "));

    // disk
    assert!(passes("<1.1:3:1 2 3,1 3,2 3:4 8,3>"));
    assert!(passes("<1.1:1:1,1,1:5,3>"));
    assert!(passes("<1.1:2:2,1 2,1 2:2,4 4>"));
    assert!(!passes("<1.1:5:2 4 5,1 2 3 5,3 4 5:8 3,8 3>"));
    assert!(!passes("
        <1.1:12:
        1 3 5 8 10 11 12,2 3 6 7 10 12 11,1 4 5 9 7 11 10 12:
        3 3 3,6 3 3>
    "));
}
