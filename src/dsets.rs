use std::{collections::{VecDeque, HashSet}, fmt};


#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum Sign {
    PLUS,
    MINUS,
    ZERO,
}

use Sign::*;


#[derive(Debug)]
pub struct Orbit2d {
    elements: Vec<usize>,
    i: usize,
    j: usize,
    is_chain: bool,
}

impl Orbit2d {
    fn len(&self) -> usize {
        self.elements.len()
    }

    fn r(&self) -> usize {
        if self.is_chain { self.len() } else { (self.len() + 1) / 2 }
    }

    fn v_min(&self) -> usize {
        match self.r() {
            1 => 3,
            2 => 2,
            _ => 1,
        }
    }
}


pub trait DSet {
    fn size(&self) -> usize;
    fn dim(&self) -> usize;
    fn get(&self, _i: usize, _d: usize) -> Option<usize>;

    fn set_count(&self) -> usize { 1 }
    fn symbol_count(&self) -> usize { 1 }
    fn m(&self, _i: usize, _j: usize, _d: usize) -> usize { 0 }


    fn partial_orientation(&self) -> Vec<Sign> {
        let mut sgn = vec![ZERO; self.size() + 1];
        let mut queue = VecDeque::new();

        sgn[1] = Sign::PLUS;
        queue.push_back(1);

        while let Some(d) = queue.pop_front() {
            for i in 0..=self.dim() {
                if let Some(di) = self.get(i, d) {
                    if sgn[di] == ZERO {
                        sgn[di] = if sgn[d] == PLUS { MINUS } else { PLUS };
                        queue.push_back(di);
                    }
                }
            }
        }

        sgn
    }


    fn is_complete(&self) -> bool {
        for i in 0..=self.dim() {
            for d in 1..=self.size() {
                if self.get(i, d) == None {
                    return false;
                }
            }
        }
        true
    }


    fn is_loopless(&self) -> bool {
        for i in 0..=self.dim() {
            for d in 1..=self.size() {
                if self.get(i, d) == Some(d) {
                    return false;
                }
            }
        }
        true
    }


    fn is_weakly_oriented(&self) -> bool {
        let ori = self.partial_orientation();

        for i in 0..=self.dim() {
            for d in 1..=self.size() {
                if let Some(di) = self.get(i, d) {
                    if di != d && ori[d] != ZERO && ori[di] == ori[d] {
                        return false;
                    }
                }
            }
        }
        true
    }


    fn is_oriented(&self) -> bool {
        self.is_loopless() && self.is_weakly_oriented()
    }


    fn oriented_cover(&self) -> PartialDSet {
        if self.is_oriented() {
            let mut ds = PartialDSet::new(self.size(), self.dim());

            for i in 0..=self.dim() {
                for d in 1..=self.size() {
                    if let Some(di) = self.get(i, d) {
                        ds.set(i, d, di);
                    }
                }
            }

            ds
        } else {
            let sz = self.size();
            let ori = self.partial_orientation();
            let mut cov = PartialDSet::new(2 * self.size(), self.dim());

            for i in 0..=self.dim() {
                for d in 1..=self.size() {
                    if let Some(di) = self.get(i, d) {
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

            cov
        }
    }


    fn orbits(&self, i: usize, j: usize) -> Vec<Orbit2d> {
        let mut orbits = vec![];
        let mut seen = HashSet::new();

        for d in 1..=self.size() {
            if !seen.contains(&d) {
                seen.insert(d);

                let mut elements = vec![d];
                let mut is_chain = false;
                let mut e = d;
                let mut k = i;

                loop {
                    if let Some(ek) = self.get(k, e) {
                        is_chain = is_chain || ek == e;
                        e = ek;
                    }
                    k = i + j - k;

                    if !seen.contains(&e) {
                        seen.insert(e);
                        elements.push(e);
                    }

                    if e == d && k == i {
                        break;
                    }
                }

                orbits.push(Orbit2d { elements, i, j, is_chain });
            }
        }

        orbits
    }


    fn morphism(&self, other: &dyn DSet, img0: usize)
        -> Option<Vec<usize>>
    {
        let mut m = vec![0; self.size() + 1];
        let mut queue = VecDeque::new();

        m[1] = img0;
        queue.push_back((1, img0));

        while let Some((d, e)) = queue.pop_front() {
            for i in 0..=self.dim() {
                let di = self.get(i, d).unwrap_or(0);
                let ei = other.get(i, e).unwrap_or(0);

                if di > 0 || ei > 0 {
                    if m[di] == 0 {
                        m[di] = ei;
                        queue.push_back((di, ei));
                    } else if m[di] != ei {
                        return None;
                    }
                }
            }
        }

        Some(m)
    }


    fn automorphisms(&self) -> Vec<Vec<usize>>
        where Self: Sized
    {
        let mut result = vec![];

        for d in 1..=self.size() {
            if let Some(map) = self.morphism(self, d) {
                result.push(map);
            }
        }

        result
    }


    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut chunks = vec!();

        chunks.push(format!("<{}.{}:", self.set_count(), self.symbol_count()));

        if self.dim() != 2 {
            chunks.push(format!("{} {}:", self.size(), self.dim()));
        } else {
            chunks.push(format!("{}:", self.size()));
        }

        for i in 0..=self.dim() {
            if i > 0 {
                chunks.push(",".to_string());
            }
            for d in 1..=self.size() {
                let e = self.get(i, d).unwrap_or(0);
                if e == 0 || e >= d {
                    if d > 1 {
                        chunks.push(" ".to_string());
                    }
                    chunks.push(e.to_string());
                }
            }
        }
        chunks.push(":".to_string());

        for i in 0..self.dim() {
            if i > 0 {
                chunks.push(",".to_string());
            }

            for orb in self.orbits(i, i + 1) {
                let &d = orb.elements.first().unwrap();
                if d > 1 {
                    chunks.push(" ".to_string());
                }
                chunks.push(self.m(i, i + 1, d).to_string());
            }
        }
        chunks.push(">".to_string());

        write!(f, "{}", chunks.join(""))
    }
}


#[derive(Clone)]
pub struct PartialDSet {
    _size: usize,
    _dim: usize,
    _op: Vec<usize>,
}

impl PartialDSet {
    pub fn new(size: usize, dim: usize) -> PartialDSet {
        assert!(size >= 1);
        assert!(dim >= 1);

        let op = vec![0; size * (dim + 1)];
        PartialDSet { _size: size, _dim: dim, _op: op }
    }

    fn idx(&self, i: usize, d: usize) -> usize {
        i * self._size + d - 1
    }

    pub fn set(&mut self, i: usize, d: usize, e: usize) {
        assert!(i <= self._dim);
        assert!(1 <= d && d <= self._size);
        assert!(1 <= e && e <= self._size);

        let kd = self.idx(i, d);
        let ke = self.idx(i, e);
        self._op[kd] = e;
        self._op[ke] = d;
    }

    pub fn grow(&mut self, count: usize) {
        self._size += count;
        self._op.append(&mut vec![0 as usize; count * (self.dim() + 1)]);
    }
}

impl DSet for PartialDSet {
    fn size(&self) -> usize {
        self._size
    }

    fn dim(&self) -> usize {
        self._dim
    }

    fn get(&self, i: usize, d: usize) -> Option<usize> {
        if i > self._dim || d < 1 || d > self._size {
            None
        } else {
            match self._op[self.idx(i, d)] {
                0 => None,
                di => Some(di)
            }
        }
    }
}

impl fmt::Display for PartialDSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        (self as &dyn DSet).fmt(f)
    }
}


pub struct SimpleDSet {
    _size: usize,
    _dim: usize,
    _op: Vec<usize>,
}

impl SimpleDSet {
    fn idx(&self, i: usize, d: usize) -> usize {
        i * self._size + d - 1
    }
}

impl DSet for SimpleDSet {
    fn size(&self) -> usize {
        self._size
    }

    fn dim(&self) -> usize {
        self._dim
    }

    fn get(&self, i: usize, d: usize) -> Option<usize> {
        if i > self._dim || d < 1 || d > self._size {
            None
        } else {
            Some(self._op[self.idx(i, d)])
        }
    }
}

impl From<PartialDSet> for SimpleDSet {
    fn from(ds: PartialDSet) -> Self {
        assert!(ds.is_complete());
        // TODO add more consistency checks here

        let PartialDSet { _size, _dim, _op } = ds;
        SimpleDSet { _size, _dim, _op }
    }
}

impl From<&PartialDSet> for SimpleDSet {
    fn from(ds: &PartialDSet) -> Self {
        assert!(ds.is_complete());
        // TODO add more consistency checks here

        let PartialDSet { _size, _dim, _op } = ds;
        SimpleDSet { _size: *_size, _dim: *_dim, _op: _op.clone() }
    }
}

impl fmt::Display for SimpleDSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        (self as &dyn DSet).fmt(f)
    }
}
