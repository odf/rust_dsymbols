use std::{fmt, str::FromStr};

use crate::dsets::*;
use crate::parse_dsym;


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
                    let ei = ds.op_unchecked(i, e);
                    is_chain |= ei == e;
                    orbit_index[i][ei] = orbit_nr;
                    seen[ei] = true;

                    e = ds.op_unchecked(i + 1, ei);
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


pub trait DSym : DSet {
    fn r(&self, i: usize, j: usize, d: usize) -> Option<usize>;
    fn v(&self, i: usize, j: usize, d: usize) -> Option<usize>;
}


pub struct TraversalCode<'a, T: DSym> {
    ds: &'a T,
    traversal: Traversal<'a, T, std::array::IntoIter<usize, 1>>,
    next_element: usize,
    element_map: Vec<usize>,
    buffer: Vec<isize>,
}


impl<'a, T: DSym> TraversalCode<'a, T> {
    pub fn new(ds: &'a T, seed: usize) -> Self {
        Self {
            ds,
            traversal: ds.traversal(0..=ds.dim(), [seed].into_iter()),
            next_element: 1,
            element_map: vec![0; ds.size() + 1],
            buffer: vec![],
        }
    }

    fn advance(&mut self) -> bool {
        if let Some((maybe_i, di, d)) = self.traversal.next() {
            if self.element_map[d] == 0 {
                self.element_map[d] = self.next_element;
            }

            if let Some(i) = maybe_i {
                self.buffer.push(i as isize);
                self.buffer.push(self.element_map[di] as isize);
                self.buffer.push(self.element_map[d] as isize);
            } else {
                self.buffer.push(-1);
                self.buffer.push(self.element_map[di] as isize);
            }

            if self.element_map[d] == self.next_element {
                self.next_element += 1;

                for i in 0..self.ds.dim() {
                    self.buffer.push(self.ds.v(i, i + 1, d).unwrap() as isize);
                }
            }

            true
        } else {
            false
        }
    }

    pub fn get(&mut self, i: usize) -> Option<isize> {
        while self.buffer.len() < i && self.advance() {
        }
        self.buffer.get(i).cloned()
    }

    pub fn get_code(&mut self) -> Vec<isize> {
        while self.advance() {
        }
        self.buffer.clone()
    }

    pub fn get_map(&mut self) -> Vec<usize> {
        while self.advance() {
        }
        self.element_map.clone()
    }
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

    pub fn build(
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

    fn op(&self, i: usize, d: usize) -> Option<usize> {
        self.dset.op(i, d)
    }

    fn is_complete(&self) -> bool {
        self.dset.is_complete() && self.orbit_vs.iter().all(|&v| v > 0)
    }

    fn m(&self, i: usize, j: usize, d: usize) -> Option<usize> {
        Some(self.r(i, j, d)? * self.v(i, j, d)?)
    }
}

impl DSym for PartialDSym {
    fn r(&self, i: usize, j: usize, d: usize) -> Option<usize> {
        if i > self.dim() || j > self.dim() || d < 1 || d > self.size() {
            None
        } else if j == i {
            Some(1)
        } else if j == i + 1 {
            Some(self.orbit_rs[self.orbit_index[i][d]])
        } else if j == i - 1 {
            Some(self.orbit_rs[self.orbit_index[j][d]])
        } else if self.op(i, d) == self.op(j, d) {
            Some(1)
        } else {
            Some(2)
        }
    }

    fn v(&self, i: usize, j: usize, d: usize) -> Option<usize> {
        if i > self.dim() || j > self.dim() || d < 1 || d > self.size() {
            None
        } else if j == i {
            Some(1)
        } else if j == i + 1 {
            Some(self.orbit_vs[self.orbit_index[i][d]])
        } else if j == i - 1 {
            Some(self.orbit_vs[self.orbit_index[j][d]])
        } else if self.op(i, d) == self.op(j, d) {
            Some(2)
        } else {
            Some(1)
        }
    }
}

impl From<SimpleDSet> for PartialDSym {
    fn from(dset: SimpleDSet) -> Self {
        let (orbit_rs, _, orbit_index) = collect_orbits(&dset);
        let orbit_vs = vec![0; orbit_rs.len()];

        PartialDSym { dset, orbit_index, orbit_rs, orbit_vs }
    }
}

impl From<PartialDSet> for PartialDSym {
    fn from(dset: PartialDSet) -> Self {
        SimpleDSet::from(dset).into()
    }
}

impl fmt::Display for PartialDSym {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        DSet::fmt(self, f)
    }
}

impl FromStr for PartialDSym {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (_, spec) = parse_dsym::parse_dsymbol(s)?;

        if spec.size < 1 {
            Err("size must be at least 1".into())
        } else if spec.dim < 1 {
            Err("dimension must be at least 1".into())
        } else if spec.op_spec.len() != spec.dim as usize + 1 {
            Err("incorrect dimension for op specifications".into())
        } else if spec.m_spec.len() != spec.dim as usize {
            Err("incorrect dimension for degree specifications".into())
        } else {
            let mut dset = PartialDSet::new(spec.size, spec.dim);

            for i in 0..=spec.dim {
                let op_i = spec.op_spec.get(i).unwrap();
                let mut k = 0;

                for d in 1..=spec.size {
                    if dset.op_unchecked(i, d) == 0 {
                        let &di = op_i.get(k)
                            .ok_or("incomplete op spec".to_string())?;
                        dset.set(i, d, di);
                        k += 1;
                    }
                }

                if k < op_i.len() {
                    return Err("unused data in op spec".into());
                }
            }

            let mut dsym = PartialDSym::from(dset);

            for i in 0..spec.dim {
                let ms_i = spec.m_spec.get(i).unwrap();
                let mut k = 0;

                for d in 1..=spec.size {
                    if dsym.v(i, i + 1, d) == Some(0) {
                        let &m = ms_i.get(k)
                            .ok_or("incomplete degree spec".to_string())?;
                        let r = dsym.r(i, i + 1, d).unwrap(); 
                        if m % r != 0 {
                            return Err("illegal degree value".into());
                        }
                        dsym.set_v(i, d, m / r);
                        k += 1;
                    }
                }

                if k < ms_i.len() {
                    return Err("unused data in degree spec".into());
                }
            }

            Ok(dsym.into())
        }
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
                    let e = self.dset.op_unchecked(i, d);
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

    fn op(&self, i: usize, d: usize) -> Option<usize> {
        self.dset.op(i, d)
    }

    fn is_complete(&self) -> bool {
        // checked on creation
        true
    }

    fn m(&self, i: usize, j: usize, d: usize) -> Option<usize> {
        Some(self.r(i, j, d)? * self.v(i, j, d)?)
    }
}

impl DSym for SimpleDSym {
    fn r(&self, i: usize, j: usize, d: usize) -> Option<usize> {
        if i > self.dim() || j > self.dim() || d < 1 || d > self.size() {
            None
        } else if j == i {
            Some(1)
        } else if j == i + 1 {
            Some(self.orbit_rs[self.orbit_index[i][d]])
        } else if j == i - 1 {
            Some(self.orbit_rs[self.orbit_index[j][d]])
        } else if self.op(i, d) == self.op(j, d) {
            Some(1)
        } else {
            Some(2)
        }
    }

    fn v(&self, i: usize, j: usize, d: usize) -> Option<usize> {
        if i > self.dim() || j > self.dim() || d < 1 || d > self.size() {
            None
        } else if j == i {
            Some(1)
        } else if j == i + 1 {
            Some(self.orbit_vs[self.orbit_index[i][d]])
        } else if j == i - 1 {
            Some(self.orbit_vs[self.orbit_index[j][d]])
        } else if self.op(i, d) == self.op(j, d) {
            Some(2)
        } else {
            Some(1)
        }
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
    where T: DSet
{
    if ds.is_oriented() {
        None
    } else {
        let mut cov: PartialDSym = setcov.clone().into();

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


#[test]
fn test_parse_from_string() {
    let s = "<1.1:2 3:2,1 2,1 2,2:6,3 2,6>";
    let dsym: PartialDSym = s.parse().unwrap();

    assert_eq!(dsym.dim(), 3);
    assert_eq!(dsym.size(), 2);
    assert_eq!(dsym.op(0, 0), None);
    assert_eq!(dsym.op(0, 3), None);
    assert_eq!(dsym.op(4, 1), None);
    assert!(dsym.is_complete());
    assert!(!dsym.is_loopless());
    assert!(dsym.is_weakly_oriented());
    assert!(!dsym.is_oriented());
    assert!(dsym.is_minimal());
    assert_eq!(dsym.to_string(), s);
    assert_eq!(dsym.automorphisms().len(), 1);
    assert_eq!(dsym.dset.automorphisms().len(), 2);
    assert_eq!(
        dsym.oriented_cover().and_then(|dso| Some(dso.to_string())),
        Some("<1.1:4 3:2 4,3 4,3 4,2 4:6,3 2,6>".to_string())
    )
}


#[test]
fn test_is_minimal() {
    let is_minimal = |s: &str| s.parse::<PartialDSym>().unwrap().is_minimal();

    assert!(is_minimal("<1.1:3:2 3,1 3,1 2 3:3,4 6>"));
    assert!(!is_minimal("<1.1:3:2 3,1 3,1 2 3:3,6 6>"));
    assert!(is_minimal("<1.1:8:2 4 6 8,8 3 5 7,1 2 3 4 5 6 7 8:4,4 6 8 4>"));
    assert!(!is_minimal("<1.1:8:2 4 6 8,8 3 5 7,1 2 3 4 5 6 7 8:4,4 6 8 6>"));
}


#[test]
fn test_traversal_code() {
    let s = "<1.1:3:1 2 3,3 2,2 3:6 4,3>";
    let dsym : PartialDSym = s.parse().unwrap();

    assert_eq!(
        TraversalCode::new(&dsym, 1).get_code(),
        &[
            -1, 1, 3, 1,
            0, 1, 1,
            1, 1, 2, 3, 1,
            0, 2, 2,
            2, 1, 3, 4, 1,
            0, 3, 3,
            1, 3, 3,
            2, 2, 2,
        ]
    );
    assert_eq!(
        TraversalCode::new(&dsym, 2).get_code(),
        &[
            -1, 1, 4, 1,
            0, 1, 1,
            1, 1, 1,
            2, 1, 2, 3, 1,
            0, 2, 2,
            1, 2, 3, 3, 1,
            0, 3, 3,
            2, 3, 3,
        ]
    );
    assert_eq!(
        TraversalCode::new(&dsym, 3).get_code(),
        &[
            -1, 1, 3, 1,
            0, 1, 1,
            1, 1, 2, 3, 1,
            0, 2, 2,
            2, 1, 1,
            2, 2, 3, 4, 1,
            0, 3, 3,
            1, 3, 3,
        ]
    );
}


#[test]
fn test_traversal_map() {
    let s = "<1.1:3:1 2 3,3 2,2 3:6 4,3>";
    let dsym : PartialDSym = s.parse().unwrap();

    assert_eq!(TraversalCode::new(&dsym, 1).get_map(), &[0, 1, 3, 2]);
    assert_eq!(TraversalCode::new(&dsym, 2).get_map(), &[0, 2, 1, 3]);
    assert_eq!(TraversalCode::new(&dsym, 3).get_map(), &[0, 2, 3, 1]);
}
