use std::collections::BTreeMap;

use crate::dsyms::DSym;
use crate::fpgroups::invariants::relator_as_vector;
use crate::fundamental_group::fundamental_group;
use crate::geometry::vec_matrix::VecMatrix;


type EdgeVectors = BTreeMap<(usize, usize), VecMatrix<i64>>;


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


fn corner_shifts<T: DSym>(cov: &T, e2t: EdgeVectors) -> EdgeVectors {
    let zero = VecMatrix::new(3, 1);

    let mut result = BTreeMap::new();

    for i in 0..=cov.dim() {
        let idcs = (0..=cov.dim()).filter(|&k| k != i);
        for (maybe_k, d, dk) in cov.traversal(idcs, 0..=cov.size()) {
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


#[test]
fn test_edge_translations() {
    use crate::dsets::DSet;
    use crate::delaney3d::pseudo_toroidal_cover;
    use crate::dsyms::PartialDSym;

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
