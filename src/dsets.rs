use std::collections::VecDeque;
use std::fmt;


#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum Sign {
    PLUS,
    MINUS,
    ZERO,
}

use Sign::*;


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


    fn orbit_reps_2d(&self, i: usize, j: usize) -> Vec<usize> {
        let mut result = vec![];
        let mut seen = vec![false; self.size() + 1];

        for d in 1..=self.size() {
            if !seen[d] {
                result.push(d);
                seen[d] = true;

                let mut e = d;
                let mut k = i;

                loop {
                    e = self.get(k, e).unwrap_or(e);
                    k = i + j - k;
                    seen[e] = true;

                    if e == d && k == i {
                        break;
                    }
                }
            }
        }

        result
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

            for d in self.orbit_reps_2d(i, i + 1) {
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
    size: usize,
    dim: usize,
    op: Vec<usize>,
}

impl PartialDSet {
    pub fn new(size: usize, dim: usize) -> PartialDSet {
        assert!(size >= 1);
        assert!(dim >= 1);

        let op = vec![0; size * (dim + 1)];
        PartialDSet { size, dim, op }
    }

    fn idx(&self, i: usize, d: usize) -> usize {
        (d - 1) * (self.dim + 1) + i
    }

    pub fn set(&mut self, i: usize, d: usize, e: usize) {
        assert!(i <= self.dim);
        assert!(1 <= d && d <= self.size);
        assert!(1 <= e && e <= self.size);

        let kd = self.idx(i, d);
        let ke = self.idx(i, e);
        self.op[kd] = e;
        self.op[ke] = d;
    }

    pub fn grow(&mut self, count: usize) {
        self.size += count;
        self.op.append(&mut vec![0 as usize; count * (self.dim() + 1)]);
    }
}

impl DSet for PartialDSet {
    fn size(&self) -> usize {
        self.size
    }

    fn dim(&self) -> usize {
        self.dim
    }

    fn get(&self, i: usize, d: usize) -> Option<usize> {
        if i > self.dim || d < 1 || d > self.size {
            None
        } else {
            match self.op[self.idx(i, d)] {
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


#[derive(Clone)]
pub struct SimpleDSet {
    size: usize,
    dim: usize,
    op: Vec<usize>,
    counter: usize,
}

impl SimpleDSet {
    fn idx(&self, i: usize, d: usize) -> usize {
        (d - 1) * (self.dim + 1) + i
    }

    pub fn from_partial(ds: PartialDSet, counter: usize) -> Self {
        assert!(ds.is_complete());
        // TODO add more consistency checks here

        let PartialDSet { size, dim, op } = ds;
        SimpleDSet { size, dim, op, counter }
    }
}

impl DSet for SimpleDSet {
    fn set_count(&self) -> usize {
        self.counter
    }

    fn size(&self) -> usize {
        self.size
    }

    fn dim(&self) -> usize {
        self.dim
    }

    fn get(&self, i: usize, d: usize) -> Option<usize> {
        if i > self.dim || d < 1 || d > self.size {
            None
        } else {
            Some(self.op[self.idx(i, d)])
        }
    }
}

impl fmt::Display for SimpleDSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        (self as &dyn DSet).fmt(f)
    }
}
