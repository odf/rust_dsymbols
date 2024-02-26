use crate::dsets::*;
use crate::dsyms::*;
use crate::util::partitions::Partition;


pub fn build_set<F>(size: usize, dim: usize, op: F) -> PartialDSet
    where F: Fn(usize, usize) -> Option<usize>
{
    let mut dset = PartialDSet::new(size, dim);
    for i in 0..=dset.dim() {
        for d in 1..=dset.size() {
            if let Some(di) = op(i, d) {
                dset.set(i, d, di);
            }
        }
    }
    dset
}


pub fn build_sym_using_vs<F>(dset: PartialDSet, v: F) -> PartialDSym
    where F: Fn(usize, usize) -> Option<usize>
{
    let mut dsym: PartialDSym = dset.into();
    for i in 0..dsym.dim() {
        for d in dsym.orbit_reps_2d(i, i + 1) {
            if let Some(v) = v(i, d) {
                dsym.set_v(i, d, v);
            }
        }
    }
    dsym
}


pub fn build_sym_using_ms<F>(dset: PartialDSet, m: F) -> PartialDSym
    where F: Fn(usize, usize) -> Option<usize>
{
    let mut dsym: PartialDSym = dset.into();
    for i in 0..dsym.dim() {
        for d in dsym.orbit_reps_2d(i, i + 1) {
            if let Some(r) = dsym.r(i, i + 1, d) {
                if let Some(m) = m(i, d) {
                    dsym.set_v(i, d, m / r);
                }
            }
        }
    }
    dsym
}


pub fn canonical<T: DSym>(ds: &T) -> PartialDSym {
    let src2img = minimal_traversal_code(ds).get_map();

    let mut img2src = vec![0; ds.size() + 1];
    for d in 1..=ds.size() {
        img2src[src2img[d]] = d;
    }

    let op = |i, d| ds.op(i, img2src[d]).map(|e| src2img[e]);
    let v = |i, d| ds.v(i, i + 1, img2src[d]);

    build_sym_using_vs(build_set(ds.size(), ds.dim(), op), v)
}


pub fn dual<T: DSym>(ds: &T) -> PartialDSym {
    let n = ds.dim();

    build_sym_using_vs(
        build_set(ds.size(), n, |i, d| ds.op(n - i, d)),
        |i, d| ds.v(n - i - 1, n - i, d)
    )
}


pub fn subsymbol<T, I>(ds: &T, indices: I, seed: usize) -> PartialDSym
    where T: DSym, I: IntoIterator<Item=usize>
{
    let indices: Vec<_> = indices.into_iter().collect();
    let elements = ds.orbit(indices.iter().cloned(), seed);

    let mut src2img = vec![0; ds.size() + 1];
    let mut img2src = vec![0; ds.size() + 1];
    let mut next = 1;
    for d in 1..=ds.size() {
        if elements.contains(&d) {
            src2img[d] = next;
            img2src[next] = d;
            next += 1;
        }
    }

    build_sym_using_vs(
        build_set(
            elements.len(),
            indices.len() - 1,
            |i, d| ds.op(indices[i], img2src[d]).map(|e| src2img[e])
        ),
        |i, d| ds.v(indices[i], indices[i + 1], img2src[d])
    )
}


pub fn cover<T, F>(ds: &T, nr_sheets: usize, sheet_map: F) -> PartialDSym
    where
        T: DSym,
        F: Fn(usize, usize, usize) -> usize
{
    let sz = ds.size();
    let src = |d: usize| (d - 1) % sz + 1;
    let op = |i, d| ds.op(i, src(d))
        .map(|di| sz * sheet_map((d - src(d)) / sz, i, src(d)) + di);

    build_sym_using_ms(
        build_set(nr_sheets * sz, ds.dim(), op),
        |i, d| ds.m(i, i + 1, (d - 1) % ds.size() + 1)
    )
}


pub fn as_partial_dsym<T: DSym>(ds: &T) -> PartialDSym {
    let op = |i, d| ds.op(i, d);
    let v = |i, d| ds.v(i, i + 1, d);

    build_sym_using_vs(build_set(ds.size(), ds.dim(), op), v)
}


pub fn oriented_cover<T: DSym>(ds: &T) -> PartialDSym {
    if ds.is_oriented() {
        as_partial_dsym(ds)
    } else {
        let ori = ds.partial_orientation();
        let sheet_map = |k, i, d| {
            if ori[d] == ori[ds.op(i, d).unwrap()] { k ^ 1 } else { k }
        };
        cover(ds, 2, sheet_map)
    }
}


pub fn minimal_image<T: DSym>(ds: &T) -> PartialDSym {
    if ds.is_minimal() {
        as_partial_dsym(ds)
    } else {
        let p = (2..=ds.size())
            .fold(Partition::new(), |p, d| ds.fold(&p, 1, d).unwrap_or(p));

        let mut src2img = vec![0; ds.size() + 1];
        let mut img2src = vec![0; ds.size() + 1];
        let mut next = 1;
        for d in 1..=ds.size() {
            let e = p.find(&d);
            if src2img[e] == 0 {
                src2img[e] = next;
                img2src[next] = e;
                next += 1;
            }
            src2img[d] = src2img[e];
        }


        build_sym_using_ms(
            build_set(
                next - 1,
                ds.dim(),
                |i, d| ds.op(i, img2src[d]).map(|e| src2img[e])
            ),
            |i, d| ds.m(i, i + 1, img2src[d])
        )
    }
}


#[test]
fn test_oriented_cover() {
    let check_cover = |src: &str, out: &str| {
        assert_eq!(
            src.parse::<PartialDSym>().map(|ds| oriented_cover(&ds)),
            out.parse::<PartialDSym>(),
        );
    };

    check_cover(
        "<1.1:6:2 4 6,6 3 5,2 4 6:3,6>",
        "<1.1:6:2 4 6,6 3 5,2 4 6:3,6>",
    );
    check_cover(
        "<1.1:6:2 4 6,6 3 5,2 5 6:3,6>",
        "<1.1:12:2 4 6 8 10 12,6 3 5 12 9 11,2 11 12 9 10 8:3 3,6 6>",
    );
    check_cover(
        "<1.1:3:1 2 3,3 2,2 3:6 4,3>",
        "<1.1:6:4 5 6,3 5 6,2 6 5:6 4,3>",
    );
    check_cover(
        "<1.1:2 3:2,1 2,1 2,2:6,3 2,6>",
        "<1.1:4 3:2 4,3 4,3 4,2 4:6,3 2,6>"
    );
}


#[test]
fn test_canonical() {
    let check_canonical = |src: &str, out: &str| {
        assert_eq!(
            src.parse::<PartialDSym>().map(|ds| canonical(&ds)),
            out.parse::<PartialDSym>()
        );
    };

    check_canonical(
        "<1.1:3:1 2 3,3 2,2 3:6 4,3>",
        "<1.1:3:1 2 3,2 3,1 3:6 4,3>",
    );
    check_canonical(
        "<1.1:2 3:2,1 2,1 2,2:6,3 2,6>",
        "<1.1:2 3:2,1 2,1 2,2:6,2 3,6>",
    );
    check_canonical(
        "<1.1:24:
        2 4 6 8 10 12 14 16 18 20 22 24,
        16 3 5 7 9 11 13 15 24 19 21 23,
        10 9 20 19 14 13 22 21 24 23 18 17:
        8 4,3 3 3 3>",
        "<1.1:24:
        2 4 6 8 10 12 14 16 18 20 22 24,
        8 3 5 7 24 11 13 15 17 19 21 23,
        9 10 21 22 17 18 13 14 20 19 24 23:
        4 8,3 3 3 3>",
    )
}


#[test]
fn test_dual() {
    let check_dual = |src: &str, out: &str| {
        assert_eq!(
            src.parse::<PartialDSym>().map(|ds| dual(&ds)),
            out.parse::<PartialDSym>()
        );
    };

    check_dual(
        "<1.1:3:1 2 3,3 2,2 3:6 4,3>",
        "<1.1:3:2 3,3 2,1 2 3:3,6 4>",
    );
    check_dual(
        "<1.1:2 3:2,1 2,1 2,2:6,3 2,6>",
        "<1.1:2 3:2,1 2,1 2,2:6,3 2,6>",
    );
    check_dual(
        "<1.1:24:
        2 4 6 8 10 12 14 16 18 20 22 24,
        16 3 5 7 9 11 13 15 24 19 21 23,
        10 9 20 19 14 13 22 21 24 23 18 17:
        8 4,3 3 3 3>",
        "<1.1:24:
        10 9 20 19 14 13 22 21 24 23 18 17,
        16 3 5 7 9 11 13 15 24 19 21 23,
        2 4 6 8 10 12 14 16 18 20 22 24:
        3 3 3 3,8 4>",
    );
}


#[test]
fn test_subsymbol() {
    fn check_subsymbol<I>(src: &str, idcs: I, seed: usize, out: &str)
        where I: IntoIterator<Item=usize>
    {
        assert_eq!(
            src.parse::<PartialDSym>().map(|ds| subsymbol(&ds, idcs, seed)),
            out.parse::<PartialDSym>()
        );
    }

    check_subsymbol(
        "<1.1:3:1 2 3,3 2,2 3:6 4,3>", [0, 1], 1,
        "<1.1:2 1:1 2,2:6>"
    );
    check_subsymbol(
        "<1.1:3:1 2 3,3 2,2 3:6 4,3>", [0, 1], 2,
        "<1.1:1 1:1,1:4>"
    );
    check_subsymbol(
        "<1.1:3:1 2 3,3 2,2 3:6 4,3>", [1, 2], 1,
        "<1.1:3 1:3 2,2 3:3>"
    );
    check_subsymbol(
        "<1.1:2 3:2,1 2,1 2,2:6,3 2,6>", [0, 1, 2], 1,
        "<1.1:2:2,1 2,1 2:6,3 2>"
    )
}


#[test]
fn test_minimal() {
    let check_minimal = |src: &str, out: &str| {
        assert_eq!(
            src.parse::<PartialDSym>().map(|ds| minimal_image(&ds)),
            out.parse::<PartialDSym>()
        );
    };

    check_minimal(
        "<1.1:24:
        2 4 6 8 10 12 14 16 18 20 22 24,
        16 3 5 7 9 11 13 15 24 19 21 23,
        10 9 20 19 14 13 22 21 24 23 18 17:
        8 4,3 3 3 3>",
        "<1.1:3:1 2 3,2 3,1 3:8 4,3>",
    );
    check_minimal(
        "<1.1:2 3:2,1 2,1 2,2:6,3 2,6>",
        "<1.1:2 3:2,1 2,1 2,2:6,3 2,6>",
    );
    check_minimal(
        "<1.1:3:1 2 3,3 2,2 3:6 6,3>",
        "<1.1:1:1,1,1:6,3>",
    );
    check_minimal(
        "<1.1:3:1 2 3,3 2,2 3:6 4,3>",
        "<1.1:3:1 2 3,3 2,2 3:6 4,3>",
    );
}
