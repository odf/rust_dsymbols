use std::fmt;

use crate::dsets::*;

pub trait DSym : DSet {
    fn get_r(&self, i: usize, d: usize) -> usize;
    fn get_v(&self, i: usize, d: usize) -> usize;

    fn r(&self, i: usize, j: usize, d: usize) -> usize {
        assert!(i <= self.dim());
        assert!(j <= self.dim());
        assert!(1 <= d && d <= self.size());

        if j == i + 1 {
            self.get_r(i, d)
        } else if i == j + 1 {
            self.get_r(j, d)
        } else if j != i && self.get(i, d) == self.get(j, d) {
            1
        } else {
            2
        }
    }

    fn v(&self, i: usize, j: usize, d: usize) -> usize {
        assert!(i <= self.dim());
        assert!(j <= self.dim());
        assert!(1 <= d && d <= self.size());

        if j == i + 1 {
            self.get_v(i, d)
        } else if i == j + 1 {
            self.get_v(j, d)
        } else if j != i && self.get(i, d) == self.get(j, d) {
            2
        } else {
            1
        }
    }
}


pub fn collect_orbits(ds: &SimpleDSet)
    -> (Vec<usize>, Vec<bool>, Vec<Vec<usize>>)
{
    let mut orbit_rs = vec![];
    let mut orbit_is_chain = vec![];
    let mut orbit_index = vec![vec![0; ds.size() + 1]; ds.dim()];
    let mut seen = vec![false; ds.size() + 1];

    for i in 0..ds.dim() {
        seen.fill(false);

        for d in 1..=ds.size() {
            if !seen[d] {
                let orbit_nr = orbit_rs.len();
                let mut e = d;
                let mut steps = 0;
                let mut is_chain = false;

                loop {
                    let ei = ds.get_unchecked(i, e);
                    is_chain |= ei == e;
                    orbit_index[i][ei] = orbit_nr;
                    seen[ei] = true;

                    e = ds.get_unchecked(i + 1, ei);
                    is_chain |= e == ei;
                    orbit_index[i][e] = orbit_nr;
                    seen[e] = true;

                    steps += 1;

                    if e == d {
                        break;
                    }
                }

                orbit_rs.push(steps);
                orbit_is_chain.push(is_chain);
            }
        }
    }

    (orbit_rs, orbit_is_chain, orbit_index)
}


#[derive(Clone)]
pub struct PartialDSym {
    dset: SimpleDSet,
    orbit_index: Vec<Vec<usize>>,
    orbit_rs: Vec<usize>,
    orbit_vs: Vec<usize>,
}

impl PartialDSym {
    pub fn new(dset: &SimpleDSet) -> PartialDSym {
        let (orbit_rs, _, orbit_index) = collect_orbits(&dset);
        let orbit_vs = vec![0; orbit_rs.len()];

        PartialDSym { dset: dset.clone(), orbit_index, orbit_rs, orbit_vs }
    }

    pub fn from(
        dset: SimpleDSet,
        orbit_index: Vec<Vec<usize>>,
        orbit_rs: Vec<usize>,
        orbit_vs: Vec<usize>
    ) -> PartialDSym
    {
        PartialDSym { dset, orbit_index, orbit_rs, orbit_vs }
    }

    pub fn set_v(&mut self, i: usize, d: usize, v: usize) {
        assert!(1 <= d);
        self.orbit_vs[self.orbit_index[i][d]] = v;
    }
}

impl DSet for PartialDSym {
    fn set_count(&self) -> usize {
        self.dset.set_count()
    }

    fn size(&self) -> usize {
        self.dset.size()
    }

    fn dim(&self) -> usize {
        self.dset.dim()
    }

    fn get(&self, i: usize, d: usize) -> Option<usize> {
        self.dset.get(i, d)
    }

    fn m(&self, i: usize, j: usize, d: usize) -> usize {
        self.r(i, j, d) * self.v(i, j, d)
    }

    fn is_complete(&self) -> bool {
        self.dset.is_complete() && self.orbit_vs.iter().all(|&v| v > 0)
    }
}

impl DSym for PartialDSym {
    fn get_r(&self, i: usize, d: usize) -> usize {
        self.orbit_rs[self.orbit_index[i][d]]
    }

    fn get_v(&self, i: usize, d: usize) -> usize {
        self.orbit_vs[self.orbit_index[i][d]]
    }
}

impl fmt::Display for PartialDSym {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        DSet::fmt(self, f)
    }
}


#[derive(Clone)]
pub struct SimpleDSym {
    dset: SimpleDSet,
    orbit_index: Vec<Vec<usize>>,
    orbit_rs: Vec<usize>,
    orbit_vs: Vec<usize>,
    counter: usize,
}

impl SimpleDSym {
    pub fn from_partial(ds: PartialDSym, counter: usize) -> Self {
        assert!(ds.is_complete());
        // TODO add more consistency checks here

        let PartialDSym { dset, orbit_index, orbit_rs, orbit_vs } = ds;
        SimpleDSym { dset, orbit_index, orbit_rs, orbit_vs, counter }
    }

    pub fn to_binary(&self) -> Option<Vec<u8>> {
        const END: usize = 0;

        if self.size() > 254 {
            None
        } else {
            let mut result = Vec::with_capacity(self.dim() * self.size());

            result.push(self.dim() as u8);

            for d in 1..=self.size() {
                for i in 0..=self.dim() {
                    let e = self.dset.get_unchecked(i, d);
                    if e >= d {
                        result.push(e as u8);
                    }
                }
            }

            for &v in &self.orbit_vs {
                result.push(v as u8);
            }

            result.push(END as u8);

            Some(result)
        }
    }
}

impl DSet for SimpleDSym {
    fn set_count(&self) -> usize {
        self.dset.set_count()
    }

    fn symbol_count(&self) -> usize {
        self.counter
    }

    fn size(&self) -> usize {
        self.dset.size()
    }

    fn dim(&self) -> usize {
        self.dset.dim()
    }

    fn get(&self, i: usize, d: usize) -> Option<usize> {
        self.dset.get(i, d)
    }

    fn m(&self, i: usize, j: usize, d: usize) -> usize {
        self.r(i, j, d) * self.v(i, j, d)
    }

    fn is_complete(&self) -> bool {
        // checked on creation
        true
    }
}

impl DSym for SimpleDSym {
    fn get_r(&self, i: usize, d: usize) -> usize {
        self.orbit_rs[self.orbit_index[i][d]]
    }

    fn get_v(&self, i: usize, d: usize) -> usize {
        self.orbit_vs[self.orbit_index[i][d]]
    }
}

impl From<PartialDSym> for SimpleDSym {
    fn from(value: PartialDSym) -> Self {
        SimpleDSym::from_partial(value, 1)
    }
}

impl fmt::Display for SimpleDSym {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        DSet::fmt(self, f)
    }
}


fn oriented_cover<T>(ds: &T, setcov: &SimpleDSet) -> Option<PartialDSym>
    where T: DSym
{
    if ds.is_oriented() {
        None
    } else {
        let mut cov = PartialDSym::new(setcov);

        for i in 0..cov.dim() {
            for d in cov.orbit_reps_2d(i, i + 1) {
                let e = (d - 1) % ds.size() + 1;
                let vd = ds.m(i, i + 1, e) / cov.r(i, i + 1, d);
                cov.set_v(i, d, vd);
            }
        }

        Some(cov)
    }
}


impl OrientedCover<PartialDSym> for PartialDSym {
    fn oriented_cover(&self) -> Option<PartialDSym> {
        self.dset.oriented_cover()
            .and_then(|setcov| oriented_cover(self, &setcov))
    }
}

impl OrientedCover<SimpleDSym> for SimpleDSym {
    fn oriented_cover(&self) -> Option<SimpleDSym> {
        self.dset.oriented_cover()
            .and_then(|setcov| oriented_cover(self, &setcov))
            .map(Into::into)
    }
}
