use std::collections::BTreeMap;

use num_rational::BigRational;
use num_traits::{FromPrimitive, Signed};

use crate::dsets::Sign;
use crate::geometry::traits::{Array2d, Entry, ScalarPtr};
use crate::pgraphs::{PeriodicGraph, VectorLabelledEdge};

use crate::dsyms::DSym;
use crate::fpgroups::invariants::relator_as_vector;
use crate::fundamental_group::fundamental_group;
use crate::geometry::vec_matrix::VecMatrix;


type EdgeVectors<T> = BTreeMap<(usize, usize), VecMatrix<T>>;


pub struct Skeleton {
    pub chamber_to_node: Vec<usize>,
    pub edge_translations: EdgeVectors<i64>,
    pub corner_shifts: EdgeVectors<i64>,
    pub graph: PeriodicGraph
}


impl Skeleton {
    pub fn of<T: DSym>(cov: &T) -> Skeleton {
        let chamber_to_node = chamber_to_node(cov);
        let edge_translations = edge_translations(cov);
        let corner_shifts = corner_shifts(cov, &edge_translations);

        let indices = cov.indices().filter(|&i| i != 1);

        let graph = cov.orbit_reps(indices, cov.elements()).iter()
            .map(|&d| skeleton_edge(
                d, cov, &edge_translations, &corner_shifts, &chamber_to_node
            ))
            .into();

        Skeleton { chamber_to_node, edge_translations, corner_shifts, graph }
    }
}


fn skeleton_edge<T: DSym>(
    d: usize,
    cov: &T,
    e2t: &EdgeVectors<i64>,
    c2s: &EdgeVectors<i64>,
    c2n: &Vec<usize>
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

    VectorLabelledEdge::make(c2n[d], c2n[e], shift)
}


fn chamber_to_node<T: DSym>(cov: &T) -> Vec<usize> {
    let mut result = vec![0; cov.size() + 1];
    let reps = cov.orbit_reps(1..=cov.dim(), cov.elements());

    for (i, &d) in reps.iter().enumerate() {
        for e in cov.orbit(1..=cov.dim(), d) {
            result[e] = i + 1;
        }
    }

    result
}


fn edge_translations<T: DSym>(cov: &T) -> EdgeVectors<i64>
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


fn corner_shifts<T: DSym>(cov: &T, e2t: &EdgeVectors<i64>) -> EdgeVectors<i64> {
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


pub fn chamber_positions<T: DSym>(cov: &T, skel: &Skeleton)
    -> EdgeVectors<BigRational>
{
    let graph = &skel.graph;
    let nodes = &skel.chamber_to_node;
    let shifts = &skel.corner_shifts;

    let shift = |d, i| shifts[&(d, i)].upcast().unwrap();

    let corners: BTreeMap<_, _> = cov.elements().map(|d|
        (d, graph.position(nodes[d]) + shift(d, 0))
    ).collect();

    let mut result = BTreeMap::new();
    for (d, p) in &corners {
        result.insert((*d, 0), p.clone());
    }

    for i in 1..=cov.dim() {
        for d0 in cov.orbit_reps(0..i, cov.elements()) {
            let orb = cov.orbit(0..i, d0);

            let mut sum = VecMatrix::zero(cov.dim(), 1);
            for d in &orb {
                sum = sum + &corners[d] - shift(*d, i);
            }
            let center = sum / BigRational::from_usize(orb.len()).unwrap();

            for d in &orb {
                result.insert((*d, i), &center + shift(*d, i));
            }
        }
    }

    result
}


fn chamber_basis<T>(pos: &EdgeVectors<T>, d: usize)
    -> VecMatrix<T>
    where T: Entry + Clone, for <'a> &'a T: ScalarPtr<T>
{
    let c0 = &pos[&(d, 0)];
    let dim = c0.nr_rows();

    let mut result = VecMatrix::zero(dim, dim);
    for i in 0..dim {
        let row = (&pos[&(d, i + 1)] - c0).transpose();
        result[i].clone_from_slice(&row[0]);
    }

    result
}


fn normalized_orientation<D>(cov: &D, skel: &Skeleton) -> Vec<Sign>
    where D: DSym,
{
    let pos = chamber_positions(cov, skel);

    let mut vol = BigRational::from_i64(0).unwrap();
    for d in cov.elements() {
        vol = vol + chamber_basis(&pos, d).determinant();
    }

    let ori = cov.partial_orientation();
    if vol.is_negative() {
        ori.iter().map(|s|
            match s {
                Sign::PLUS => Sign::MINUS,
                Sign::MINUS => Sign::PLUS,
                Sign::ZERO => Sign::ZERO,
            }
        ).collect()
    } else {
        ori
    }
}


fn non_degenerate_chamber<D, T>(cov: &D, pos: &EdgeVectors<T>)
    -> Option<usize>
    where
        D: DSym,
        T: Entry + Clone, for <'a> &'a T: ScalarPtr<T>
{
    cov.elements().find(|&d| !chamber_basis(pos, d).determinant().is_zero())
}


fn orbit_indices<T, I1, I2>(ds: &T, indices: I1, seeds: I2) -> Vec<usize>
    where T: DSym, I1: IntoIterator<Item=usize>, I2: IntoIterator<Item=usize>
{
    let mut result = vec![0; ds.size() + 1];
    let mut next_index = 0;

    for (i, _, d) in ds.traversal(indices, seeds) {
        if i.is_none() {
            next_index += 1;
        }
        result[d] = next_index - 1;
    }

    result
}


fn tile_surfaces<S, D, I>(
    cov: &D,
    skel: &Skeleton,
    vertex_pos: &BTreeMap<usize, VecMatrix<S>>,
    seeds: I
)
    -> Vec<(Vec<VecMatrix<S>>, Vec<Vec<usize>>)>
    where
        S: Entry + Clone + FromPrimitive, for <'a> &'a S: ScalarPtr<S>,
        D: DSym,
        I: IntoIterator<Item=usize>
{
    let dim = cov.dim();
    let ori = normalized_orientation(cov, skel);
    let mut result = vec![];

    for d0 in seeds {
        let tile = cov.orbit(0..dim, d0);
        let corner_indices = orbit_indices(cov, 1..dim, tile.clone());

        let pos = cov.orbit_reps(1..dim, tile.clone()).iter().map(|&d|
            &vertex_pos[&skel.chamber_to_node[d]] +
            skel.corner_shifts[&(d, 0)].upcast().unwrap()
        ).collect();

        let mut faces = vec![];
        for d in cov.orbit_reps([0, 1], tile) {
            let e = match ori[d] {
                Sign::MINUS => cov.op(0, d).unwrap(),
                _ => d
            };
            faces.push(
                cov.orbit_2d(0, 1, e).into_iter()
                    .step_by(2)
                    .map(|d| corner_indices[d])
                    .collect()
            );
        }

        result.push((pos, faces))
    }

    result
}


#[cfg(test)]
mod test {
    use crate::derived::canonical;
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
                    assert_eq!(c2n[d], c2n[cov.op(i, d).unwrap()]);
                }
            }

            let reps = cov.orbit_reps(1..=cov.dim(), cov.elements());
            for &d in &reps {
                assert!(c2n[d] > 0);
                for &e in &reps {
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

    #[test]
    fn test_skeleton() {
        fn test(spec: &str, nv: usize, ne: usize) {
            let ds = spec.parse::<PartialDSym>().unwrap();
            let cov = pseudo_toroidal_cover(&ds).unwrap();
            let skel = Skeleton::of(&cov);

            assert_eq!(skel.graph.vertices().len(), nv);
            assert_eq!(skel.graph.edges().len(), ne);
        }

        test("<1.1:1 3:1,1,1,1:4,3,4>", 1, 3);
        test("<1.1:2 3:2,1 2,1 2,2:6,3 2,6>", 2, 4);
        test("<1.1:2 3:1 2,1 2,1 2,2:3 3,3 4,4>", 1, 6);
        test("<1.1:6 3:2 4 6,1 2 3 5 6,3 4 5 6,2 3 4 5 6:6 4,2 3 3,8 4 4>", 4, 12);
    }

    #[test]
    fn test_chamber_positions() {
        fn test(spec: &str) {
            let ds = spec.parse::<PartialDSym>().unwrap();
            let cov = canonical(&pseudo_toroidal_cover(&ds).unwrap());
            let skel = Skeleton::of(&cov);
            let pos = chamber_positions(&cov, &skel);

            for i in 1..=cov.dim() {
                let mut sum: VecMatrix<BigRational> = VecMatrix::zero(cov.dim(), 1);
                let mut always_zero = true;

                for d in cov.elements() {
                    let n = BigRational::from_usize(
                        cov.orbit(1..i, d).iter().count()
                    ).unwrap();
                    sum = sum + (&pos[&(d, i)] - &pos[&(d, 0)]) / n;
                    always_zero &= sum.is_zero();

                    for k in 0..i {
                        assert_eq!(
                            pos[&(d, i)],
                            pos[&(cov.op(k, d).unwrap_or(d), i)]
                        )
                    }
                }

                assert!(sum.is_zero());
                assert!(!always_zero);
            }
        }

        test("<1.1:1 3:1,1,1,1:4,3,4>");
        test("<1.1:2 3:2,1 2,1 2,2:6,3 2,6>");
        test("<1.1:3 3:1 3,2 3,1 2 3,1 2 3:3,4 4,4 3 3>");
        test("<1.1:3 3:1 2 3,1 2 3,2 3,1 3:4 3 3,4 4,3>");
    }

    #[test]
    fn test_tile_surfaces() {
        fn test(spec: &str) {
            let ds = spec.parse::<PartialDSym>().unwrap();
            let cov = canonical(&pseudo_toroidal_cover(&ds).unwrap());
            let skel = Skeleton::of(&cov);
            let pos = skel.graph.vertices().iter()
                .map(|&v| (v, skel.graph.position(v)))
                .collect();
            let surf = tile_surfaces(&cov, &skel, &pos, [1]);

            for (pos, faces) in surf {
                for p in pos {
                    for i in 0..p.nr_rows() {
                        print!("  {}", p[(i, 0)]);
                    }
                    println!();
                }
                for f in faces {
                    println!("{f:?}");
                }
                println!();
            }
        }

        test("<1.1:1 3:1,1,1,1:4,3,4>");
        test("<1.1:2 3:2,1 2,1 2,2:6,3 2,6>");
    }
}
