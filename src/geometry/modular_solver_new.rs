use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Zero;

use super::prime_residue_classes::PrimeResidueClass;
use super::traits::{Array2d, Entry, Scalar, ScalarPtr};
use super::vec_matrix::VecMatrix;


impl<const P: i64> Scalar for PrimeResidueClass<P> {}
impl<const P: i64> ScalarPtr<PrimeResidueClass<P>> for &PrimeResidueClass<P> {}


impl<const P: i64> Entry for PrimeResidueClass<P> {
    fn can_divide(_: &Self, b: &Self) -> bool {
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


fn number_of_p_adic_steps_needed(
    a: &VecMatrix<i64>, b: &VecMatrix<i64>, prime: i64
) -> u64
{
    fn column_norm(a: &VecMatrix<i64>, j: usize) -> f64 {
        (0..a.nr_rows())
            .map(|i| (a[i][j] as f64).powf(2.0))
            .sum::<f64>()
            .sqrt()
    }

    let mut log_norms: Vec<_> = (0..a.nr_columns())
        .map(|j| column_norm(a, j).ln())
        .collect();
    log_norms.push((0..b.nr_columns())
        .map(|j| column_norm(b, j).ln())
        .max_by(|a, b| a.total_cmp(b))
        .unwrap());

    log_norms.sort_by(|a, b| a.total_cmp(b));

    let log_delta: f64 = log_norms.iter().skip(1).sum();
    let golden_ratio = (1.0 + (5 as f64).sqrt()) / 2.0;

    (2.0 * (log_delta + golden_ratio.ln()) / (prime as f64).ln()).ceil() as u64
}


fn rational_reconstruction(s: &BigInt, h: &BigInt) -> BigRational {
    let (mut u, mut u1) = (h.clone(), s.clone());
    let (mut v, mut v1) = (BigInt::from(0), BigInt::from(1));
    let mut sign = BigInt::from(1);

    while &u1.pow(2) > &h {
        let (q, r) = (&u / &u1, &u % &u1);

        (u, u1) = (u1, r);

        let v1_next = v + q * &v1;
        (v, v1) = (v1, v1_next);

        sign = -sign;
    }

    BigRational::new(sign * u1, v1)
}


const PRIME: i64 = 9999991;


pub fn solve(a: &VecMatrix<i64>, b: &VecMatrix<i64>)
    -> Option<VecMatrix<BigRational>>
{
    if let Some(c) = a.to::<PrimeResidueClass<PRIME>>().inverse() {
        let nr_steps = number_of_p_adic_steps_needed(a, b, PRIME);
        let nrows = b.nr_rows();
        let ncols = b.nr_columns();

        let mut p = BigInt::from(1);
        let mut b = b.clone();
        let mut s = VecMatrix::<BigInt>::new(nrows, ncols);

        for step in 0..nr_steps {
            let x = (&c * b.to()).to();
            for i in 0..nrows {
                for j in 0..ncols {
                    s[i][j] += &p * x[i][j];
                }
            }

            p *= PRIME;

            if step + 1 < nr_steps {
                let ax = a * x;
                for i in 0..nrows {
                    for j in 0..ncols {
                        b[i][j] = (b[i][j] - ax[i][j]) / PRIME;
                    }
                }
            }
        }

        let mut result = VecMatrix::new(nrows, ncols);
        for i in 0..nrows {
            for j in 0..ncols {
                result[i][j] = rational_reconstruction(&s[i][j], &p);
            }
        }
        Some(result)
    } else {
        None
    }
}


mod property_based_tests {
    use super::*;
    use proptest::prelude::*;
    use proptest::collection::vec;

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

    fn test_modular_solver(m: &VecMatrix<i64>, v: &VecMatrix<i64>) {
        let convert = |m: &VecMatrix<i64>| m.to::<BigInt>().to();

        let b = m * v;
        if let Some(sol) = solve(m, &b) {
            assert_eq!(convert(m) * sol, convert(&b));
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

        #[test]
        fn test_modular_solver_small(
            (m, v) in equations::<i64>(3, 2, 6)
        ) {
            test_modular_solver(&m, &v);
        }

        #[test]
        fn test_modular_solver_large(
            (m, v) in equations::<i64>(1_000_000_000, 2, 10)
        ) {
            test_modular_solver(&m, &v);
        }
    }
}
