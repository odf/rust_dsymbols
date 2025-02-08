use num_traits::Zero;

use super::prime_residue_classes::PrimeResidueClass;
use super::traits::{Array2d, Entry, Scalar, ScalarPtr};


impl<const P: i64> Scalar for PrimeResidueClass<P> {}
impl<const P: i64> ScalarPtr<PrimeResidueClass<P>> for &PrimeResidueClass<P> {}


impl<const P: i64> Entry for PrimeResidueClass<P> {
    fn can_divide(a: &Self, b: &Self) -> bool {
        !b.is_zero()
    }

    fn pivot_row<M: Array2d<Self>>(col: usize, row0: usize, a: &M)
        -> Option<usize>
    {
        for row in row0..a.nr_rows() {
            if !a[(row, col)].is_zero() {
                return Some(row);
            }
        }

        None
    }

    fn clear_col<A: Array2d<Self>, B: Array2d<Self>>(
        col: usize, row1: usize, row2: usize, a: &mut A, x: Option<&mut B>
    ) {
        let f = a[(row1, col)] / a[(row2, col)];
        a[(row1, col)] = Self::zero();

        for k in (col + 1)..a.nr_columns() {
            a[(row1, k)] = a[(row1, k)] - a[(row2, k)] * f;
        }

        if let Some(x) = x {
            for k in 0..x.nr_columns() {
                x[(row1, k)] = x[(row1, k)] - x[(row2, k)] * f;
            }
        }
    }
}


mod property_based_tests {
    use super::*;
    use proptest::prelude::*;
    use proptest::collection::vec;
    use crate::geometry::vec_matrix::VecMatrix;

    fn matrix_from_values<T>(v: &[T], m: usize) -> VecMatrix<T>
        where T: Scalar + Clone
    {
        let n = v.len() / m;
        let mut result = VecMatrix::new(n, m);

        let mut k = 0;
        for i in 0..n {
            for j in 0..m{
                result[i][j] = v[k].clone();
                k += 1;
            }
        }

        result
    }

    fn entry<T>(size: i32)
        -> impl Strategy<Value=T>
        where T: From<i32> + std::fmt::Debug
    {
        (0..size).prop_map(T::from)
    }

    fn sized_matrix<T>(size: i32, n: usize, m: usize)
        -> impl Strategy<Value=VecMatrix<T>>
        where T: Scalar + Clone + From<i32> + std::fmt::Debug + 'static
    {
        vec(entry(size), n * m).prop_map(move |v| matrix_from_values(&v, m))
    }

    fn matrix<T>(entry_size: i32, dmin: usize, dmax: usize)
        -> impl Strategy<Value=VecMatrix<T>>
        where T: Scalar + Clone + From<i32> + std::fmt::Debug + 'static
    {
        (dmin..=dmax).prop_flat_map(move |n|
            sized_matrix(entry_size, n, n)
        )
    }

    fn equations<T>(entry_size: i32, dmin: usize, dmax: usize)
        -> impl Strategy<Value=(VecMatrix<T>, VecMatrix<T>)>
        where T: Scalar + Clone + From<i32> + std::fmt::Debug + 'static
    {
        (dmin..=dmax).prop_flat_map(move |n|
            (sized_matrix(entry_size, n, n), sized_matrix(entry_size, n, 1))
        )
    }

    fn test_matrix<T>(m: &VecMatrix<T>)
        where
            T: Entry + Clone + PartialEq + std::fmt::Debug,
            for <'a> &'a T: ScalarPtr<T>
    {
        let n = m.nr_rows();
        let zero = VecMatrix::new(n, 1);
        let one = VecMatrix::identity(n);

        assert_eq!(m.determinant().is_zero(), m.rank() < n);
        assert_eq!(m.null_space().len(), n - m.rank());

        for v in m.null_space() {
            assert_eq!(m * v, zero);
        }

        let nul = m.null_space_matrix();
        assert_eq!((m * &nul), VecMatrix::new(n, nul.nr_columns()));

        if let Some(inv) = m.inverse() {
            assert_eq!(m * inv, one);
        }

        assert_eq!(m.inverse().is_some(), m.rank() == m.nr_rows());
    }

    fn test_solver<T>(m: &VecMatrix<T>, v: &VecMatrix<T>)
        where
            T: Entry + Clone + PartialEq + std::fmt::Debug,
            for <'a> &'a T: ScalarPtr<T>
    {
        let b = m * v;
        if let Some(sol) = m.solve(&b) {
            assert_eq!(m * sol, b);
        }
    }

    const P_SMALL: i64 = 61;
    const P_LARGE: i64 = 3_037_000_493;

    proptest! {
        #[test]
        fn test_matrix_small(
            m in matrix::<PrimeResidueClass<P_SMALL>>(3, 2, 8)
        ) {
            test_matrix(&m);
        }

        #[test]
        fn test_matrix_large(
            m in matrix::<PrimeResidueClass<P_LARGE>>(1_000_000_000, 2, 16)
        ) {
            test_matrix(&m);
        }

        #[test]
        fn test_solver_small(
            (m, v) in equations::<PrimeResidueClass<P_SMALL>>(3, 2, 8)
        ) {
            test_solver(&m, &v);
        }

        #[test]
        fn test_solver_large(
            (m, v) in equations::<PrimeResidueClass<P_LARGE>>(1_000_000_000, 2, 16)
        ) {
            test_solver(&m, &v);
        }
    }
}
