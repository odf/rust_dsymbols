use core::cmp::Ordering;
use std::cell::UnsafeCell;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fmt;
use std::ops::Neg;

use num_rational::BigRational;

use crate::geometry::traits::Array2d;
use crate::geometry::vec_matrix::VecMatrix;
use crate::geometry::modular_solver_new::solve;


#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VectorLabelledEdge {
    head: usize,
    tail: usize,
    shift: VecMatrix<i64>
}


impl VectorLabelledEdge {
    pub fn make(head: usize, tail: usize, shift: VecMatrix<i64>) -> Self {
        assert_eq!(shift.nr_columns(), 1);
        Self { head, tail, shift }
    }

    pub fn dim(&self) -> usize {
        self.shift.nr_rows()
    }

    pub fn canonical(&self) -> Self {
        if self.tail < self.head {
            return -self;
        } else if self.tail == self.head {
            for i in 0..self.dim() {
                if self.shift[i][0] < 0 {
                    return -self;
                }
            }
        }

        self.clone()
    }
}


impl fmt::Display for VectorLabelledEdge {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {}, [", self.head, self.tail)?;

        for i in 0..self.shift.nr_rows() {
            if i > 0 {
                write!(f, ", ")?;
            }
            write!(f, "{}", self.shift[i][0])?;
        }

        write!(f, "])")?;

        Ok(())
    }
}


impl PartialOrd for VectorLabelledEdge {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        assert_eq!(self.shift.nr_rows(), other.shift.nr_rows());

        match self.head.partial_cmp(&other.head) {
            Some(Ordering::Equal) => {}
            ord => return ord,
        }

        match self.tail.partial_cmp(&other.tail) {
            Some(Ordering::Equal) => {}
            ord => return ord,
        }

        for i in 0..self.shift.nr_rows() {
            match self.shift[i][0].partial_cmp(&other.shift[i][0]) {
                Some(Ordering::Equal) => {}
                ord => return ord,
            }
        }

        Some(Ordering::Equal)
    }
}


impl Ord for VectorLabelledEdge {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}


impl Neg for VectorLabelledEdge {
    type Output = VectorLabelledEdge;

    fn neg(self) -> Self::Output {
        VectorLabelledEdge::make(self.tail, self.head, -&self.shift)
    }
}


impl Neg for &VectorLabelledEdge {
    type Output = VectorLabelledEdge;

    fn neg(self) -> Self::Output {
        VectorLabelledEdge::make(self.tail, self.head, -&self.shift)
    }
}


pub struct PeriodicGraph {
    edges: Vec<VectorLabelledEdge>,
    vertices: Vec<usize>,
    incidences: HashMap<usize, Vec<VectorLabelledEdge>>,
    positions: UnsafeCell<HashMap<usize, VecMatrix<BigRational>>>
}


impl<I> From<I> for PeriodicGraph
where I: IntoIterator<Item=VectorLabelledEdge>
{
    fn from(edges: I) -> PeriodicGraph {
        let edges: Vec<_> = edges.into_iter()
            .map(|e| e.canonical())
            .collect::<BTreeSet<_>>().into_iter() // sorts and deduplicates
            .collect();

        assert!(edges.len() > 0);
        let d = edges[0].dim();
        assert!(edges.iter().all(|e| e.dim() == d));

        let mut vertices = BTreeSet::new();
        for e in edges.iter() {
            vertices.insert(e.head);
            vertices.insert(e.tail);
        }
        let vertices: Vec<_> = vertices.into_iter().collect();

        let mut incidences: HashMap<_, Vec<_>> = HashMap::new();
        for e in edges.iter() {
            incidences.entry(e.head).or_default().push(e.clone());
            incidences.entry(e.tail).or_default().push(-e);
        }

        let positions = UnsafeCell::new(HashMap::new());

        PeriodicGraph { edges, vertices, incidences, positions }
    }
}


impl PeriodicGraph {
    pub fn dim(&self) -> usize {
        self.edges[0].dim()
    }

    pub fn edges(&self) -> &Vec<VectorLabelledEdge> {
        &self.edges
    }

    pub fn vertices(&self) -> &Vec<usize> {
        &self.vertices
    }

    pub fn incidences(&self, v: usize) -> Option<&Vec<VectorLabelledEdge>> {
        self.incidences.get(&v)
    }

    pub fn position(&self, v: usize) -> VecMatrix<BigRational> {
        assert!(self.incidences.contains_key(&v));

        if let Some(output) = unsafe { (*self.positions.get()).get(&v) } {
            output.clone()
        } else {
            let positions = barycentric_placement(self);
            let output = positions[&v].clone();
            unsafe { *self.positions.get() = positions };
            output
        }
    }
}


fn barycentric_placement(g: &PeriodicGraph) -> HashMap<usize, VecMatrix<BigRational>> {
    let verts = g.vertices();
    let vidcs: BTreeMap<_, _> = verts.iter().enumerate()
        .map(|(i, &e)| (e, i)).collect();

    let n = verts.len();
    let d = g.dim();

    let mut a = VecMatrix::<i64>::new(n, n);
    let mut t = VecMatrix::<i64>::new(n, d);

    a[0][0] = 1;

    for i in 1..n {
        for ngb in g.incidences(verts[i]).unwrap() {
            let j = vidcs[&ngb.tail];
            a[i][j] -= 1;
            a[i][i] += 1;

            let s = &ngb.shift;
            for k in 0..d {
                t[i][k] = t[i][k] + s[k][0];
            }
        }
    }

    let p = solve(&a, &t).unwrap();

    let mut result = HashMap::new();
    for i in 0..n {
        result.insert(verts[i], p.submatrix([i], 0..d));
    }

    result
}
