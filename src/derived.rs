use crate::dsets::*;
use crate::dsyms::*;


fn build_set<F>(size: usize, dim: usize, op: F) -> PartialDSet
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


fn build_sym_from_set<F>(dset: PartialDSet, v: F) -> PartialDSym
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


pub fn build_sym<F1, F2>(size: usize, dim: usize, op: F1, v: F2) -> PartialDSym
    where
        F1: Fn(usize, usize) -> Option<usize>,
        F2: Fn(usize, usize) -> Option<usize>
{
    build_sym_from_set(build_set(size, dim, op), v)
}


pub fn canonical<T: DSym>(ds: &T) -> PartialDSym {
    let src2img = minimal_traversal_code(ds).get_map();

    let mut img2src = vec![0; ds.size() + 1];
    for d in 1..=ds.size() {
        img2src[src2img[d]] = d;
    }
    let img2src = img2src; // shadow with immutable version

    build_sym(
        ds.size(),
        ds.dim(),
        |i, d| ds.op(i, img2src[d]).map(|e| src2img[e]),
        |i, d| ds.v(i, i + 1, img2src[d]),
    )
}


fn set_cover<T, F>(ds: &T, nr_sheets: usize, sheet_map: F) -> PartialDSet
    where
        T: DSet,
        F: Fn(usize, usize, usize) -> usize
{
    let sz = ds.size();
    let src = |d: usize| (d - 1) % sz + 1;
    let op = |i, d| ds.op(i, src(d))
        .map(|di| sz * sheet_map((d - src(d)) / sz, i, src(d)) + di);

    build_set(nr_sheets * sz, ds.dim(), op)
}


fn oriented_set_cover<T>(ds: &T) -> Option<PartialDSet>
    where T: DSet
{
    if ds.is_oriented() {
        None
    } else {
        let ori = ds.partial_orientation();
        let sheet_map = |k, i, d| {
            if ori[d] == ori[ds.op(i, d).unwrap()] { k ^ 1 } else { k }
        };
        Some(set_cover(ds, 2, sheet_map))
    }
}


fn cover<T, F>(ds: &T, nr_sheets: usize, sheet_map: F) -> PartialDSym
    where
        T: DSym,
        F: Fn(usize, usize, usize) -> usize
{
    let mut cov: PartialDSym = set_cover(ds, nr_sheets, sheet_map).into();

    for i in 0..=cov.dim() {
        for d in cov.orbit_reps_2d(i, i + 1) {
            if let Some(r) = cov.r(i, i + 1, d) {
                if let Some(m) = ds.m(i, i + 1, (d - 1) % ds.size() + 1) {
                    cov.set_v(i, d, m / r);
                }
            }
        }
    }

    cov
}


fn oriented_cover<T>(ds: &T) -> Option<PartialDSym>
    where T: DSym
{
    if ds.is_oriented() {
        None
    } else {
        let ori = ds.partial_orientation();
        let sheet_map = |k, i, d| {
            if ori[d] == ori[ds.op(i, d).unwrap()] { k ^ 1 } else { k }
        };
        Some(cover(ds, 2, sheet_map))
    }
}


pub trait OrientedCover<T> {
    fn oriented_cover(&self) -> Option<T>;
}


impl OrientedCover<PartialDSet> for PartialDSet {
    fn oriented_cover(&self) -> Option<PartialDSet> {
        oriented_set_cover(self)
    }
}

impl OrientedCover<SimpleDSet> for SimpleDSet {
    fn oriented_cover(&self) -> Option<SimpleDSet> {
        oriented_set_cover(self).map(Into::into)
    }
}


impl OrientedCover<PartialDSym> for PartialDSym {
    fn oriented_cover(&self) -> Option<PartialDSym> {
        oriented_cover(self)
    }
}

impl OrientedCover<SimpleDSym> for SimpleDSym {
    fn oriented_cover(&self) -> Option<SimpleDSym> {
        oriented_cover(self).map(Into::into)
    }
}


#[test]
fn test_oriented_cover() {
    let check_cover = |src: &str, cov: &str| {
        assert_eq!(
            src.parse::<PartialDSym>()
                .map(|ds| ds.oriented_cover().unwrap_or(ds)),
            cov.parse::<PartialDSym>(),
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
    let check_canonical = |src: &str, canon: &str| {
        assert_eq!(
            src.parse::<PartialDSym>().map(|ds| canonical(&ds)),
            canon.parse::<PartialDSym>()
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
