use std::collections::BTreeMap;

use crate::dsyms::DSym;
use crate::fpgroups::invariants::relator_as_vector;
use crate::fundamental_group::fundamental_group;
use crate::geometry::vec_matrix::VecMatrix;


fn edge_translations<T: DSym>(cov: &T)
    -> BTreeMap<(usize, usize), VecMatrix<i64>>
{
    let fg = fundamental_group(cov);
    let nr_gens = fg.nr_generators();

    let mut mat = VecMatrix::new(fg.relators.len(), nr_gens);
    for (i, w) in fg.relators.iter().enumerate() {
        let row = relator_as_vector(nr_gens, w);
        for j in 0..row.len() {
            mat[i][j] = row[j] as i64;
        }
    }

    let nul = mat.null_space_matrix();

    let mut result = BTreeMap::new();
    for ((d, i), w) in fg.edge_to_word {
        let vec: Vec<_> = relator_as_vector(nr_gens, &w).iter()
            .map(|&x| x as i64)
            .collect();
        result.insert((d, i), (VecMatrix::from(vec) * &nul).transpose());
    }

    result
}


#[test]
fn test_edge_translations() {
    use crate::delaney3d::pseudo_toroidal_cover;
    use crate::dsyms::PartialDSym;
    use crate::geometry::traits::Array2d;

    let ds = "<1.1:1 3:1,1,1,1:4,3,4>".parse::<PartialDSym>().unwrap();
    let cov = pseudo_toroidal_cover(&ds).unwrap();
    let e2t = edge_translations(&cov);

    for (_, t) in e2t {
        assert_eq!(t.nr_rows(), 3);
        assert_eq!(t.nr_columns(), 1);
    }
}
