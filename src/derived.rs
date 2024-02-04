use crate::dsets::*;
use crate::dsyms::*;


fn oriented_set_cover<T>(ds: &T) -> Option<PartialDSet>
    where T: DSet
{
    if ds.is_oriented() {
        None
    } else {
        let sz = ds.size();
        let ori = ds.partial_orientation();
        let mut cov = PartialDSet::new(2 * ds.size(), ds.dim());

        for i in 0..=ds.dim() {
            for d in 1..=ds.size() {
                if let Some(di) = ds.op(i, d) {
                    if ori[di] != ori[d] {
                        cov.set(i, d, di);
                        cov.set(i, d + sz, di + sz);
                    } else {
                        cov.set(i, d, di + sz);
                        cov.set(i, d + sz, di);
                    }
                }
            }
        }

        Some(cov)
    }
}


fn oriented_cover<T>(ds: &T) -> Option<PartialDSym>
    where T: DSym
{
    if ds.is_oriented() {
        None
    } else {
        let mut cov: PartialDSym = oriented_set_cover(ds)?.into();

        for i in 0..cov.dim() {
            for d in cov.orbit_reps_2d(i, i + 1) {
                let e = (d - 1) % ds.size() + 1;
                let vd = ds.m(i, i + 1, e)? / cov.r(i, i + 1, d)?;
                cov.set_v(i, d, vd);
            }
        }

        Some(cov)
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
            src.parse::<PartialDSym>().ok().and_then(|ds|
                Some(ds.oriented_cover().unwrap_or(ds).to_string())
            ).unwrap(),
            cov
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
