use std::fmt;

use crate::dsets::*;

pub trait DSym : DSet {
    fn r(&self, i: usize, j: usize, d: usize) -> usize;
    fn v(&self, i: usize, j: usize, d: usize) -> usize;

    fn m(&self, i: usize, j: usize, d: usize) -> usize {
        self.r(i, j, d) * self.v(i, j, d)
    }
}


pub fn collect_orbits(ds: &SimpleDSet) -> (Vec<usize>, Vec<Vec<usize>>) {
    let mut orbit_rs = vec![];
    let mut orbit_index = vec![vec![0; ds.size() + 1]; ds.dim()];
    let mut seen = vec![false; ds.size() + 1];

    for i in 0..ds.dim() {
        seen.fill(false);

        for d in 1..=ds.size() {
            if !seen[d] {
                let orbit_nr = orbit_rs.len();

                orbit_index[i][d] = orbit_nr;
                seen[d] = true;

                let mut e = d;
                let mut k = i;
                let mut steps = 0;

                loop {
                    e = ds.get(k, e).unwrap();
                    k = i + (i + 1) - k;
                    steps += 1;

                    orbit_index[i][e] = orbit_nr;
                    seen[e] = true;

                    if e == d && k == i {
                        break;
                    }
                }

                orbit_rs.push(steps / 2);
            }
        }
    }

    (orbit_rs, orbit_index)
}


#[derive(Clone)]
struct PartialDSym {
    dset: SimpleDSet,
    orbit_index: Vec<Vec<usize>>,
    orbit_rs: Vec<usize>,
    orbit_vs: Vec<usize>,
}

impl PartialDSym {
    fn new(dset: &SimpleDSet) -> PartialDSym {
        let (orbit_rs, orbit_index) = collect_orbits(&dset);
        let orbit_vs = vec![0; orbit_rs.len()];

        PartialDSym { dset: dset.clone(), orbit_index, orbit_rs, orbit_vs }
    }

    fn set_v(&mut self, i: usize, j: usize, d: usize, v: usize) {
        assert!(1 <= d);

        if j == i + 1 {
            self.orbit_vs[self.orbit_index[i][d]] = v;
        } else if i == j + 1 {
            self.orbit_vs[self.orbit_index[j][d]] = v;
        } else {
            panic!("Illegal index pair: {},{}", i, j);
        }
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
}

impl DSym for PartialDSym {
    fn r(&self, i: usize, j: usize, d: usize) -> usize {
        assert!(i <= self.dim());
        assert!(j <= self.dim());
        assert!(1 <= d && d <= self.size());

        if j == i + 1 {
            self.orbit_rs[self.orbit_index[i][d]]
        } else if i == j + 1 {
            self.orbit_rs[self.orbit_index[j][d]]
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
            self.orbit_vs[self.orbit_index[i][d]]
        } else if i == j + 1 {
            self.orbit_vs[self.orbit_index[j][d]]
        } else if j != i && self.get(i, d) == self.get(j, d) {
            2
        } else {
            1
        }
    }
}

impl fmt::Display for PartialDSym {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        (self as &dyn DSet).fmt(f)
    }
}


#[derive(Clone)]
struct SimpleDSym {
    dset: SimpleDSet,
    orbit_index: Vec<Vec<usize>>,
    orbit_rs: Vec<usize>,
    orbit_vs: Vec<usize>,
    counter: usize,
}

impl SimpleDSym {
    fn from_partial(ds: PartialDSym, counter: usize) -> Self {
        assert!(ds.is_complete());
        // TODO add more consistency checks here

        let PartialDSym { dset, orbit_index, orbit_rs, orbit_vs } = ds;
        SimpleDSym { dset, orbit_index, orbit_rs, orbit_vs, counter }
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
}

impl DSym for SimpleDSym {
    fn r(&self, i: usize, j: usize, d: usize) -> usize {
        assert!(i <= self.dim());
        assert!(j <= self.dim());
        assert!(1 <= d && d <= self.size());

        if j == i + 1 {
            self.orbit_rs[self.orbit_index[i][d]]
        } else if i == j + 1 {
            self.orbit_rs[self.orbit_index[j][d]]
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
            self.orbit_vs[self.orbit_index[i][d]]
        } else if i == j + 1 {
            self.orbit_vs[self.orbit_index[j][d]]
        } else if j != i && self.get(i, d) == self.get(j, d) {
            2
        } else {
            1
        }
    }
}

impl fmt::Display for SimpleDSym {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        (self as &dyn DSet).fmt(f)
    }
}
