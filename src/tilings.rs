use std::collections::BTreeMap;
use std::ops::Neg;

use crate::dsyms::DSym;
use crate::fpgroups::invariants::relator_as_vector;
use crate::fundamental_group::fundamental_group;
use crate::geometry::traits::Array2d;
use crate::geometry::vec_matrix::VecMatrix;


type EdgeVectors = BTreeMap<(usize, usize), VecMatrix<i64>>;


#[derive(Clone)]
struct VectorLabelledEdge {
    head: usize,
    tail: usize,
    shift: VecMatrix<i64>
}


impl VectorLabelledEdge {
    fn from(head: usize, tail: usize, shift: VecMatrix<i64>) -> Self {
        assert_eq!(shift.nr_columns(), 1);
        Self { head, tail, shift }
    }

    fn dim(&self) -> usize {
        self.shift.nr_rows()
    }

    fn canonical(&self) -> Self {
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


impl Neg for VectorLabelledEdge {
    type Output = VectorLabelledEdge;

    fn neg(self) -> Self::Output {
        VectorLabelledEdge::from(self.tail, self.head, -&self.shift)
    }
}


impl Neg for &VectorLabelledEdge {
    type Output = VectorLabelledEdge;

    fn neg(self) -> Self::Output {
        VectorLabelledEdge::from(self.tail, self.head, -&self.shift)
    }
}


struct Skeleton {
    chamber_to_node: Vec<usize>,
    edge_translations: EdgeVectors,
    corner_shifts: EdgeVectors,
    edges: Vec<VectorLabelledEdge>
}


impl Skeleton {
    fn of<T: DSym>(cov: &T) -> Skeleton {
        todo!()
    }
}


fn skeleton_edge<T: DSym>(
    d: usize,
    cov: &T,
    e2t: &EdgeVectors,
    c2s: &EdgeVectors,
    c2v: &BTreeMap<usize, usize>
)
    -> VectorLabelledEdge
{
    let e = cov.op(0, d).unwrap();
    let sd = &c2s[&(d, 0)];
    let se = &c2s[&(e, 0)];

    let shift = if let Some(t) = e2t.get(&(d, 0)) {
        se + t - sd
    } else {
        se - sd
    };

    VectorLabelledEdge::from(c2v[&d], c2v[&e], shift)
}


fn chamber_to_node<T: DSym>(cov: &T) -> BTreeMap<usize, usize> {
    let mut result = BTreeMap::new();
    let reps = cov.orbit_reps(1..=cov.dim(), cov.elements());

    for (i, &d) in reps.iter().enumerate() {
        for e in cov.orbit(1..=cov.dim(), d) {
            result.insert(e, i + 1);
        }
    }

    result
}


fn edge_translations<T: DSym>(cov: &T) -> EdgeVectors
{
    let fg = fundamental_group(cov);
    let nr_gens = fg.nr_generators();

    let mut mat = VecMatrix::new(fg.relators.len(), nr_gens);
    for (i, w) in fg.relators.iter().enumerate() {
        mat[i].copy_from_slice(&relator_as_vector(nr_gens, w));
    }

    let nul = mat.null_space_matrix();

    BTreeMap::from_iter(fg.edge_to_word.iter().map(|(&(d, i), w)| (
        (d, i),
        (VecMatrix::from(relator_as_vector(nr_gens, w)) * &nul).transpose()
    )))
}


fn corner_shifts<T: DSym>(cov: &T, e2t: &EdgeVectors) -> EdgeVectors {
    let zero = VecMatrix::new(3, 1);

    let mut result = BTreeMap::new();

    for i in cov.indices() {
        let idcs = cov.indices().filter(|&k| k != i);
        for (maybe_k, d, dk) in cov.traversal(idcs, cov.elements()) {
            let shift = if let Some(k) = maybe_k {
                result.get(&(d, i)).unwrap() - e2t.get(&(d, k)).unwrap_or(&zero)
            } else {
                zero.clone()
            };
            result.insert((dk, i), shift);
        }
    }

    result
}


mod test {
    use crate::dsets::DSet;
    use crate::delaney3d::pseudo_toroidal_cover;
    use crate::dsyms::PartialDSym;

    use super::*;

    #[test]
    fn test_edge_translations() {
        fn test(spec: &str) {
            let ds = spec.parse::<PartialDSym>().unwrap();
            let cov = pseudo_toroidal_cover(&ds).unwrap();
            let e2t = edge_translations(&cov);

            for i in 0..=cov.dim() {
                for j in i..=cov.dim() {
                    for d in cov.orbit_reps([i, j], 1..=cov.size()) {
                        let mut s = VecMatrix::new(3, 1);
                        let mut e = d;
                        let mut k = i;
                        loop {
                            if let Some(v) = e2t.get(&(e, k)) {
                                s = s + v;
                            }
                            e = cov.op(k, e).unwrap();
                            k = i + j - k;

                            if e == d && k == i {
                                break;
                            }
                        }

                        assert!(s.is_zero());
                    }
                }
            }
        }

        test("<1.1:1 3:1,1,1,1:4,3,4>");
        test("<1.1:2 3:2,1 2,1 2,2:6,3 2,6>");
        test("<1.1:2 3:1 2,1 2,1 2,2:3 3,3 4,4>");
        test("<1.1:6 3:2 4 6,1 2 3 5 6,3 4 5 6,2 3 4 5 6:6 4,2 3 3,8 4 4>");
    }


    #[test]
    fn test_corner_shifts() {
        fn test(spec: &str) {
            let ds = spec.parse::<PartialDSym>().unwrap();
            let cov = pseudo_toroidal_cover(&ds).unwrap();
            let e2t = edge_translations(&cov);
            let c2s = corner_shifts(&cov, &e2t);

            for i in cov.indices() {
                for j in cov.indices() {
                    if j != i {
                        for d in cov.elements() {
                            let sd = &c2s[&(d, i)];
                            let sdj = &c2s[&(cov.op(j, d).unwrap(), i)];
                            if let Some(t) = e2t.get(&(d, j)) {
                                assert_eq!(sd, &(sdj + t));
                            } else {
                                assert_eq!(sd, sdj);
                            }
                        }
                    }
                }
            }
        }

        test("<1.1:1 3:1,1,1,1:4,3,4>");
        test("<1.1:2 3:2,1 2,1 2,2:6,3 2,6>");
        test("<1.1:2 3:1 2,1 2,1 2,2:3 3,3 4,4>");
        test("<1.1:6 3:2 4 6,1 2 3 5 6,3 4 5 6,2 3 4 5 6:6 4,2 3 3,8 4 4>");
    }


    #[test]
    fn test_chamber_to_node() {
        fn test(spec: &str) {
            let ds = spec.parse::<PartialDSym>().unwrap();
            let cov = pseudo_toroidal_cover(&ds).unwrap();
            let c2n = chamber_to_node(&cov);

            for d in cov.elements() {
                for i in 1..=cov.dim() {
                    assert_eq!(c2n[&d], c2n[&cov.op(i, d).unwrap()]);
                }
            }

            let reps = cov.orbit_reps(1..=cov.dim(), cov.elements());
            for d in &reps {
                assert!(c2n[&d] > 0);
                for e in &reps {
                    if d != e {
                        assert_ne!(c2n[d], c2n[e]);
                    }
                }
            }
        }

        test("<1.1:1 3:1,1,1,1:4,3,4>");
        test("<1.1:2 3:2,1 2,1 2,2:6,3 2,6>");
        test("<1.1:2 3:1 2,1 2,1 2,2:3 3,3 4,4>");
        test("<1.1:6 3:2 4 6,1 2 3 5 6,3 4 5 6,2 3 4 5 6:6 4,2 3 3,8 4 4>");
    }
}
