use std::ops::{Add, Div, Index, IndexMut, Mul};
use num_traits::Zero;

use crate::geometry::traits::{Entry, Scalar, ScalarPtr, Array2d};


#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Matrix<T, const N: usize, const M: usize> {
    data: [[T; M]; N]
}


impl<T, const N: usize, const M: usize> Array2d<T> for Matrix<T, N, M> {
    fn nr_rows(&self) -> usize {
        N
    }

    fn nr_columns(&self) -> usize {
        M
    }
}


impl<T, const N: usize, const M: usize>
    Index<usize> for Matrix<T, N, M>
{
    type Output = [T; M];

    fn index(&self, index: usize) -> &Self::Output {
        assert!(index < N);
        &self.data[index]
    }
}


impl<T, const N: usize, const M: usize>
    Index<(usize, usize)> for Matrix<T, N, M>
{
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (i, j) = index;
        assert!(i < N);
        assert!(j < M);

        &self.data[i][j]
    }
}


impl<T, const N: usize, const M: usize>
    IndexMut<usize> for Matrix<T, N, M>
{
    fn index_mut(&mut self, index: usize) -> &mut [T; M] {
        assert!(index < N);
        &mut self.data[index]
    }
}


impl<T, const N: usize, const M: usize>
    IndexMut<(usize, usize)> for Matrix<T, N, M>
{
    fn index_mut(&mut self, index: (usize, usize)) -> &mut T {
        let (i, j) = index;
        assert!(i < N);
        assert!(j < M);

        &mut self.data[i][j]
    }
}


impl<T: Scalar, const N: usize, const M: usize>
    From<[[T; M]; N]> for Matrix<T, N, M>
{
    fn from(data: [[T; M]; N]) -> Self {
        Self { data }
    }
}


impl<T: Scalar, const M: usize>
    From<[T; M]> for Matrix<T, 1, M>
{
    fn from(data: [T; M]) -> Self {
        Self { data: [data] }
    }
}


impl<T: Scalar + Clone, const M: usize>
    From<&[T; M]> for Matrix<T, 1, M>
{
    fn from(data: &[T; M]) -> Self {
        Self { data: [data.clone()] }
    }
}


impl<T: Scalar>
    From<T> for Matrix<T, 1, 1>
{
    fn from(data: T) -> Self {
        Self { data: [[data]] }
    }
}


impl<T: Scalar + Clone , const N: usize, const M: usize> Matrix<T, N, M> {
    pub fn new() -> Self {
        Matrix::from(
            core::array::from_fn(|_|
                core::array::from_fn(|_|
                    T::zero()
                )
            )
        )
    }

    pub fn transpose(&self) -> Matrix<T, M, N> {
        let mut result = Matrix::new();
        for i in 0..M {
            for j in 0..N {
                result[i][j] = self[j][i].clone();
            }
        }
        result
    }

    pub fn swap_rows(&mut self, i: usize, j: usize) {
        assert!(i < N);
        assert!(j < N);
        assert_ne!(i, j);

        self.data.swap(i, j)
    }

    pub fn hstack<const L: usize, const S: usize>(self, rhs: &Matrix<T, N, L>)
        -> Matrix<T, N, S>
    {
        assert_eq!(S, N + M);

        let mut result = Matrix::new();

        for i in 0..N {
            for j in 0..M {
                result[i][j] = self[i][j].clone();
            }
            for j in 0..L {
                result[i][M + j] = rhs[i][j].clone();
            }
        }
        result
    }

    pub fn vstack<const L: usize, const S: usize>(self, rhs: &Matrix<T, L, M>)
        -> Matrix<T, S, M>
    {
        assert_eq!(S, N + L);

        let mut result = Matrix::new();

        for j in 0..M {
            for i in 0..N {
                result[i][j] = self[i][j].clone();
            }
            for i in 0..L {
                result[N + i][j] = rhs[i][j].clone();
            }
        }
        result
    }

    pub fn submatrix<I, J, const K: usize, const L: usize>(
        &self, rows: I, columns: J
    )
        -> Matrix<T, K, L>
        where
            I: IntoIterator<Item=usize>,
            J: IntoIterator<Item=usize>
    {
        let rows: Vec<_> = rows.into_iter().collect();
        let columns: Vec<_> = columns.into_iter().collect();

        assert!(rows.iter().all(|&i| i < N));
        assert!(columns.iter().all(|&j| j < M));
        assert_eq!(rows.len(), K);
        assert_eq!(columns.len(), L);

        let mut result = Matrix::new();

        for i in 0..K {
            for j in 0..L {
                result[i][j] = self[rows[i]][columns[j]].clone();
            }
        }

        result
    }
}


impl<T: Scalar + Clone, const N: usize> Matrix<T, N, N> {
    pub fn identity() -> Self {
        let mut result = Matrix::new();
        for i in 0..N {
            result[i][i] = T::one();
        }
        result
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize>
    Add<&Matrix<T, N, M>> for &Matrix<T, N, M>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = Matrix<T, N, M>;

    fn add(self, rhs: &Matrix<T, N, M>) -> Self::Output {
        let mut result = Matrix::new();

        for i in 0..N {
            for j in 0..M {
                result[i][j] = &self[i][j] + &rhs[i][j];
            }
        }

        result
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize>
    Add<Matrix<T, N, M>> for Matrix<T, N, M>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = Matrix<T, N, M>;

    fn add(self, rhs: Matrix<T, N, M>) -> Self::Output {
        &self + &rhs
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize>
    Add<[[T; M]; N]> for Matrix<T, N, M>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = Matrix<T, N, M>;

    fn add(self, rhs: [[T; M]; N]) -> Self::Output {
        self + Matrix::from(rhs)
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize> Zero for Matrix<T, N, M>
    where for <'a> &'a T: ScalarPtr<T>
{
    fn zero() -> Self {
        Matrix::new()
    }

    fn is_zero(&self) -> bool {
        for i in 0..N {
            for j in 0..M {
                if !self[i][j].is_zero() {
                    return false
                }
            }
        }

        true
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize, const L: usize>
    Mul<&Matrix<T, M, L>> for &Matrix<T, N, M>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = Matrix<T, N, L>;

    fn mul(self, rhs: &Matrix<T, M, L>) -> Self::Output {
        let mut result = Matrix::new();

        for i in 0..N {
            for j in 0..L {
                let mut x = T::zero();
                for k in 0..M {
                    x = x + &self[i][k] * &rhs[k][j];
                }
                result[i][j] = x;
            }
        }

        result
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize, const L: usize>
    Mul<&Matrix<T, M, L>> for Matrix<T, N, M>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = Matrix<T, N, L>;

    fn mul(self, rhs: &Matrix<T, M, L>) -> Self::Output {
        &self * rhs
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize, const L: usize>
    Mul<Matrix<T, M, L>> for &Matrix<T, N, M>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = Matrix<T, N, L>;

    fn mul(self, rhs: Matrix<T, M, L>) -> Self::Output {
        self * &rhs
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize, const L: usize>
    Mul<Matrix<T, M, L>> for Matrix<T, N, M>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = Matrix<T, N, L>;

    fn mul(self, rhs: Matrix<T, M, L>) -> Self::Output {
        &self * &rhs
    }
}


impl<T: Scalar + Clone, const M: usize, const L: usize>
    Mul<&Matrix<T, M, L>> for &[T; M]
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = [T; L];

    fn mul(self, rhs: &Matrix<T, M, L>) -> Self::Output {
        (Matrix::from(self) * rhs)[0].clone()
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize, const L: usize>
    Mul<Matrix<T, M, L>> for [[T; M]; N]
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = Matrix<T, N, L>;

    fn mul(self, rhs: Matrix<T, M, L>) -> Self::Output {
        Matrix::from(self) * rhs
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize, const L: usize>
    Mul<[[T; L]; M]> for Matrix<T, N, M>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = Matrix<T, N, L>;

    fn mul(self, rhs: [[T; L]; M]) -> Self::Output {
        self * Matrix::from(rhs)
    }
}


impl<const N: usize, const M: usize>
    Mul<Matrix<i64, N, M>> for i64
{
    type Output = Matrix<i64, N, M>;

    fn mul(self, rhs: Matrix<i64, N, M>) -> Self::Output {
        let mut result = Matrix::new();

        for i in 0..N {
            for j in 0..M {
                result[i][j] = self * rhs[i][j];
            }
        }
        result
    }
}


impl<const N: usize, const M: usize>
    Mul<Matrix<f64, N, M>> for f64
{
    type Output = Matrix<f64, N, M>;

    fn mul(self, rhs: Matrix<f64, N, M>) -> Self::Output {
        let mut result = Matrix::new();

        for i in 0..N {
            for j in 0..M {
                result[i][j] = self * rhs[i][j];
            }
        }
        result
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize>
    Mul<&T> for &Matrix<T, N, M>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = Matrix<T, N, M>;

    fn mul(self, rhs: &T) -> Self::Output {
        let mut result = Matrix::new();

        for i in 0..N {
            for j in 0..M {
                result[i][j] = &self[i][j] * rhs;
            }
        }
        result
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize>
    Mul<&T> for Matrix<T, N, M>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = Self;

    fn mul(self, rhs: &T) -> Self::Output {
        &self * rhs
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize>
    Mul<T> for Matrix<T, N, M>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        &self * &rhs
    }
}


#[derive(Debug)]
pub struct RowEchelonMatrix<T: Entry, const N: usize, const M: usize> {
    multiplier: Matrix<T, N, N>,
    result: Matrix<T, N, M>,
    columns: [usize; N],
    rank: usize,
    nr_swaps: usize
}


impl<T: Entry + Clone, const N: usize, const M: usize>
    RowEchelonMatrix<T, N, M>
{
    pub fn new(m: &Matrix<T, N, M>) -> Self {
        let mut u = m.clone();
        let mut s = Matrix::identity();
        let mut row = 0;
        let mut nr_swaps = 0;
        let mut cols = [N; N];

        for col in 0..M {
            if let Some(pr) = Entry::pivot_row(col, row, &u) {
                if pr != row {
                    u.swap_rows(pr, row);
                    s.swap_rows(pr, row);
                    nr_swaps += 1;
                }

                for r in (row + 1)..N {
                    Entry::clear_col(col, r, row, &mut u, Some(&mut s));
                }

                cols[row] = col;
                row += 1;
            }
        }

        RowEchelonMatrix {
            multiplier: s,
            result: u,
            columns: cols,
            rank: row,
            nr_swaps
        }
    }
}


impl<T: Entry + Clone, const N: usize, const M: usize> Matrix<T, N, M>
    where for <'a> &'a T: ScalarPtr<T>
{
    fn rank(&self) -> usize {
        RowEchelonMatrix::new(self).rank
    }

    fn null_space(&self) -> Vec<Matrix<T, M, 1>> {
        let re = RowEchelonMatrix::new(&self.transpose());
        let s = re.multiplier.transpose();

        (re.rank..M).map(|i| s.submatrix(0..M, [i])).collect()
    }

    fn solve<const K: usize>(&self, rhs: &Matrix<T, N, K>)
        -> Option<Matrix<T, M, K>>
        where T: Div<T, Output=T>
    {
        let re = RowEchelonMatrix::new(self);
        let y = re.multiplier * rhs;

        if !(re.rank..N).all(|i| (0..K).all(|j| y[i][j].is_zero())) {
            return None;
        }

        let mut result = Matrix::zero();

        for row in (0..re.rank).rev() {
            let a = &re.result[row] * &result;
            let b = &y[row];
            let x = &re.result[row][re.columns[row]];
            for k in 0..K {
                let t = &b[k] - &a[k];
                if Entry::can_divide(&t, x) {
                    result[re.columns[row]][k] = &t / x;
                } else {
                    return None;
                }
            }
        }

        Some(result)
    }
}


impl<T: Entry + Clone, const N: usize> Matrix<T, N, N>
    where for <'a> &'a T: ScalarPtr<T>
{
    fn determinant(&self) -> T {
        match self.nr_rows() {
            0 => T::one(),
            1 => self[0][0].clone(),
            2 => {
                &self[0][0] * &self[1][1] - &self[0][1] * &self[1][0]
            },
            3 => {
                &self[0][0] * (&self[1][1] * &self[2][2]) +
                &self[0][1] * (&self[1][2] * &self[2][0]) +
                &self[0][2] * (&self[1][0] * &self[2][1]) -
                &self[0][2] * (&self[1][1] * &self[2][0]) -
                &self[0][0] * (&self[1][2] * &self[2][1]) -
                &self[0][1] * (&self[1][0] * &self[2][2])
            },
            _ => {
                let re = RowEchelonMatrix::new(self);
                let mut result = T::one();
                for i in 0..N {
                    result = &result * &re.result[i][i];
                }
                if re.nr_swaps % 2 == 0 { result } else { -result }
            }
        }
    }

    fn inverse(&self) -> Option<Self>
        where T: Div<T, Output=T>
    {
        self.solve(&Matrix::identity())
    }
}


#[test]
fn test_matrix_copy() {
    let m1 = Matrix::from([[1, 2], [3, 4]]);
    let mut m2 = m1;
    m2[0][0] = 0;
    assert_eq!(m1 + m2, Matrix::from([[1, 4], [6, 8]]));
    assert_eq!(m1 + (-1 * m2), Matrix::from([[1, 0], [0, 0]]));
}


#[test]
fn test_matrix_indexing() {
    let mut m = Matrix::from([[1.0, 1.0], [0.0, 1.0]]);

    m[1][0] = m[1][0] + 4.0;

    assert_eq!(m, [[1.0, 1.0], [4.0, 1.0]].into());
    assert_eq!(m[1], [4.0, 1.0]);
}


#[test]
fn test_matrix_identity() {
    let m = Matrix::<i64, 3, 3>::identity();
    assert_eq!(m, [[1, 0, 0], [0, 1, 0], [0, 0, 1]].into());
    assert_eq!(Matrix::identity(), [[1, 0], [0, 1]].into());
}


#[test]
fn test_matrix_add() {
    assert_eq!(
        Matrix::from([[1, 2], [3, 4]]) + Matrix::from([[1, 2], [3, 4]]),
        Matrix::from([[2, 4], [6, 8]])
    );

    assert_eq!(
        Matrix::from([[1, 2], [3, 4]]) + [[4, 3], [2, 1]],
        Matrix::from([[5, 5], [5, 5]])
    );
}


#[test]
fn test_matrix_zero() {
    assert_eq!(Matrix::zero(), Matrix::from([[0; 3]; 2]));
    assert!(Matrix::from([[0; 4]; 5]).is_zero());
}


#[test]
fn test_matrix_mul() {
    assert_eq!(
        (Matrix::from([[1, 2, 3]]) * [[3], [2], [1]])[0][0],
        10
    );
    assert_eq!(
        Matrix::from([1, 2, 3]) * Matrix::from([3, 2, 1]).transpose(),
        10.into()
    );
    assert_eq!(
        [[1, 2, 3]] * Matrix::from([3, 2, 1]).transpose(),
        [[10]].into()
    );
    assert_eq!(
        Matrix::from([[1.0, 1.0], [0.0, 1.0]]) * [[1.0, 2.0], [0.0, 1.0]],
        [[1.0, 3.0], [0.0, 1.0]].into()
    );
    assert_eq!(
        (
            Matrix::from([[1.0, 1.0], [0.0, 1.0]]) *
            [[1.0, 2.0], [0.0, 1.0]]
        )[0][1],
        3.0
    );
    assert_eq!(
        (
            Matrix::from([[1.0, 1.0], [0.0, 1.0]]) *
            [[1.0, 2.0], [0.0, 1.0]]
        )[0],
        [1.0, 3.0]
    );
    assert_eq!(
        (Matrix::from([[1.0, 1.0]]) * [[1.0, 2.0], [0.0, 1.0]])[0],
        [1.0, 3.0]
    );
    assert_eq!(
        Matrix::from([[1, 2], [3, 4]]) * 2,
        Matrix::from([[2, 4], [6, 8]])
    );
    assert_eq!(
        Matrix::from([[1.0, 2.0], [3.0, 4.0]]) * 2.0,
        Matrix::from([[2.0, 4.0], [6.0, 8.0]])
    );
    assert_eq!(
        2.5 * Matrix::from([[1.0, 2.0], [3.0, 4.0]]),
        Matrix::from([[2.5, 5.0], [7.5, 10.0]])
    );
}


#[test]
fn test_matrix_stack() {
    assert_eq!(
        Matrix::from([[1, 2, 3]]).vstack(&Matrix::from([[4, 5, 6], [7, 8, 9]])),
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]].into()
    );
    assert_eq!(
        Matrix::from([[1, 2], [3, 4]])
            .hstack(&Matrix::from([[5, 6], [7, 8]])),
        [[1, 2, 5, 6], [3, 4, 7, 8]].into()
    );
}


#[test]
fn test_matrix_submatrix() {
    assert_eq!(
        Matrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).submatrix(0..2, [0, 2]),
        Matrix::from([[1, 3], [4, 6]])
    )
}


#[test]
fn test_matrix_determinant() {
    assert_eq!(
        Matrix::from([[1.0, 0.3, 0.7], [0.0, 2.0, 1.2], [0.0, 0.0, 0.25]])
            .determinant(),
        0.5
    );

    assert_eq!(
        Matrix::from([
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 0, 1],
            [0, 0, 1, 0]
        ]).determinant(),
        -1
    );

    assert_eq!(
        Matrix::from([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 0.0]
        ]).determinant(),
        -1.0
    );

    assert_eq!(Matrix::<f64, 1, 1>::identity().determinant(), 1.0);
    assert_eq!(Matrix::<f64, 2, 2>::identity().determinant(), 1.0);
    assert_eq!(Matrix::<f64, 3, 3>::identity().determinant(), 1.0);
    assert_eq!(Matrix::<f64, 4, 4>::identity().determinant(), 1.0);

    assert_eq!((Matrix::<f64, 1, 1>::identity() * 2.0).determinant(), 2.0);
    assert_eq!((Matrix::<f64, 2, 2>::identity() * 2.0).determinant(), 4.0);
    assert_eq!((Matrix::<f64, 3, 3>::identity() * 2.0).determinant(), 8.0);
    assert_eq!((Matrix::<f64, 4, 4>::identity() * 2.0).determinant(), 16.0);
}

#[test]
fn test_matrix_nullspace() {
    let a = Matrix::from([[0.0, 0.0], [0.0, 0.0]]);
    let n = a.null_space();
    assert_eq!(n.len(), 2);
    for v in a.null_space() {
        assert_eq!(&a * v, Matrix::from([[0.0], [0.0]]));
    }

    let a = Matrix::from([[-25, 48], [50, -96]]);
    let n = a.null_space();
    assert_eq!(n.len(), 2 - a.rank());
    for v in a.null_space() {
        assert_eq!(&a * v, Matrix::from([[0], [0]]));
    }

    let a = Matrix::from([[1, 2], [3, 4]]);
    let n = a.null_space();
    assert_eq!(n, vec![]);

    let a = Matrix::from([[1.0, 2.0], [3.0, 6.0]]);
    let n = a.null_space();
    assert_eq!(n.len(), 1);
    for v in a.null_space() {
        assert_eq!(&a * v, Matrix::from([[0.0], [0.0]]));
    }

    let a = Matrix::from([[1.0, 2.0], [3.0, 4.0]]);
    let n = a.null_space();
    assert_eq!(n, vec![]);

    let a = Matrix::from([[0.0, 1.0, 0.0]]);
    let n = a.null_space();
    assert_eq!(n.len(), 2);
    for v in n {
        assert_eq!(&a * v, Matrix::from([[0.0]]));
    }
}


#[test]
fn test_matrix_solve_i64() {
    fn test<const N: usize, const M: usize, const L: usize>(
        a: [[i64; M]; N], b: [[i64; L]; M]
    )
    {
        let a = Matrix::from(a);
        let b = Matrix::from(b);
        let ab = &a * &b;

        assert_eq!(&a * a.solve(&ab).unwrap(), ab);

        if N == M && a.rank() == N {
            assert_eq!(a.solve(&ab), Some(b));
        }
    }

    test([[1, 2], [3, 4]], [[1, 0], [1, -3]]);
    test([[1, 2], [3, 6]], [[1], [1]]);
}


#[test]
fn test_matrix_solve_f64() {
    use crate::geometry::traits::allclose;

    fn test<const N: usize, const M: usize, const L: usize>(
        a: [[f64; M]; N], b: [[f64; L]; M]
    )
    {
        let a = Matrix::from(a);
        let b = Matrix::from(b);
        let ab = &a * &b;

        assert!(allclose(&(&a * a.solve(&ab).unwrap()), &ab, 1.0e-9, 1.0e-9));

        if N == M && a.rank() == N {
            assert!(allclose(&a.solve(&ab).unwrap(), &b, 1.0e-9, 1.0e-9));
        }
    }

    test([[1.0, 2.0], [3.0, 4.0]], [[1.0, 0.0], [1.0, -3.0]]);
    test([[1.0, 2.0], [3.0, 6.0]], [[1.0], [1.0]]);
}


#[test]
fn test_matrix_inverse_i64() {
    fn test<const N: usize>(a: [[i64; N]; N])
    {
        let a = Matrix::from(a);

        if a.rank() == N {
            assert_eq!(&a * a.inverse().unwrap(), Matrix::identity());
        } else {
            assert_eq!(a.inverse(), None);
        }
    }

    test([[1, 2], [3, 5]]);

    // Inverse exists, but is not integral:
    assert_eq!(Matrix::from([[1, 2], [3, 4]]).inverse(), None);
}


#[test]
fn test_matrix_inverse_f64() {
    use crate::geometry::traits::allclose;

    fn test<const N: usize>(a: [[f64; N]; N])
    {
        let a = Matrix::from(a);

        if a.rank() == N {
            assert!(allclose(
                &(&a * a.inverse().unwrap()),
                &Matrix::identity(),
                1.0e-9,
                1.0e-9
            ));
        } else {
            assert_eq!(a.inverse(), None);
        }
    }

    test([[1.0, 2.0], [3.0, 4.0]]);
    test([[1.0, 2.0], [3.0, 6.0]]);
}


mod test_big_rational {
    use super::*;
    use num_bigint::BigInt;
    use num_traits::{One, FromPrimitive};
    use num_rational::BigRational;

    fn matrix<const N: usize, const M: usize>(m: [[i64; M]; N])
        -> Matrix<BigRational, N, M>
    {
        let mut result = Matrix::new();

        for i in 0..N {
            for j in 0..M {
                result[(i, j)] = BigRational::from_i64(m[i][j]).unwrap();
            }
        }

        result
    }

    #[test]
    fn test_determinant() {
        assert_eq!(
            matrix([
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 0, 1],
                [0, 0, 1, 0]
            ]).determinant(),
            -BigRational::one()
        );

        assert_eq!(
            matrix([[7, 8, 8], [4, 5, 6], [1, 2, 3]]).determinant(),
            BigRational::from_i64(-3).unwrap()
        );

        let one = &BigRational::one();
        let two = &BigRational::from_i64(2).unwrap();
        let half = &(one / two);

        let i1 = &Matrix::<BigRational, 1, 1>::identity();
        let i2 = &Matrix::<BigRational, 2, 2>::identity();
        let i3 = &Matrix::<BigRational, 3, 3>::identity();
        let i4 = &Matrix::<BigRational, 4, 4>::identity();

        assert_eq!(i1.determinant(), *one);
        assert_eq!(i2.determinant(), *one);
        assert_eq!(i3.determinant(), *one);
        assert_eq!(i4.determinant(), *one);

        assert_eq!((i1 * two).determinant(), two.pow(1));
        assert_eq!((i2 * two).determinant(), two.pow(2));
        assert_eq!((i3 * two).determinant(), two.pow(3));
        assert_eq!((i4 * two).determinant(), two.pow(4));

        assert_eq!((i1 * half).determinant(), half.pow(1));
        assert_eq!((i2 * half).determinant(), half.pow(2));
        assert_eq!((i3 * half).determinant(), half.pow(3));
        assert_eq!((i4 * half).determinant(), half.pow(4));
    }

    #[test]
    fn test_rank() {
        assert_eq!(matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]).rank(), 3);
        assert_eq!(matrix([[1, 2, 3], [2, 4, 6], [3, 6, 9]]).rank(), 1);
        assert_eq!(matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).rank(), 2);
        assert_eq!(matrix([[7, 8, 8], [4, 5, 6], [1, 2, 3]]).rank(), 3);
    }

    #[test]
    fn test_solve() {
        let half = &BigRational::new(BigInt::from(1), BigInt::from(2));

        let a = matrix([[1, 2], [3, 4]]);
        let x = matrix([[1, 0], [1, -3]]);
        let b = &a * &x;
        assert_eq!(b, matrix([[3, -6], [7, -12]]));
        let s = a.solve(&b).unwrap();
        assert_eq!(s, x);

        let a = matrix([[1, 2], [3, 6]]);
        let x = matrix([[1], [1]]);
        let b = &a * &x;
        assert_eq!(b, matrix([[3], [9]]));
        let s = a.solve(&b).unwrap();
        assert_eq!(&a * s, b);

        let a = matrix([[1, 2], [3, 5]]);
        let x = matrix([[1], [2]]) * half;
        let b = &a * &x;
        assert_eq!(b, matrix([[5], [13]]) * half);
        let s = a.solve(&b).unwrap();
        assert_eq!(s, x);
    }

    #[test]
    fn test_nullspace() {
        let a = matrix([[1, 2], [3, 4]]);
        let n = a.null_space();
        for v in &n {
            assert_eq!(&a * v, Matrix::zero());
        }
        assert_eq!(n.len(), 0);

        let a = matrix([[1, 2], [3, 6]]);
        let n = a.null_space();
        for v in &n {
            assert_eq!(&a * v, Matrix::zero());
        }
        assert_eq!(n.len(), 1);
    }

    #[test]
    fn test_inverse() {
        let a = matrix([[1, 2], [3, 4]]);
        assert_eq!(&a * &a.inverse().unwrap(), Matrix::identity());

        let a = matrix([[1, 2], [3, 5]]);
        assert_eq!(&a * &a.inverse().unwrap(), Matrix::identity());

        let a = matrix([[1, 2], [3, 6]]);
        assert!(a.inverse().is_none());
    }
}


mod property_based_tests {
    use super::*;
    use num_rational::BigRational;
    use num_traits::FromPrimitive;
    use proptest::prelude::*;
    use proptest::collection::vec;
    use crate::geometry::traits::allclose;

    trait FromI32 {
        fn from(i: i32) -> Self;
    }

    impl FromI32 for i64 {
        fn from(i: i32) -> Self {
            i as Self
        }
    }

    impl FromI32 for f64 {
        fn from(i: i32) -> Self {
            i as Self
        }
    }

    impl FromI32 for BigRational {
        fn from(i: i32) -> Self {
            BigRational::from_i32(i).unwrap()
        }
    }

    fn matrix_from_values<T, const N: usize, const M: usize>(v: &[T])
        -> Matrix<T, N, M>
        where T: Scalar + Clone
    {
        let mut result = Matrix::new();

        let mut k = 0;
        for i in 0..N {
            for j in 0..M {
                result[i][j] = v[k].clone();
                k += 1;
            }
        }

        result
    }

    fn entry<T>(size: i32)
        -> impl Strategy<Value=T>
        where T: FromI32 + std::fmt::Debug
    {
        (0..size).prop_map(|i: i32| T::from(i))
    }

    fn matrix<T, const N: usize, const M: usize>(size: i32)
        -> impl Strategy<Value=Matrix<T, N, M>>
        where T: Scalar + Clone + FromI32 + std::fmt::Debug + 'static
    {
        vec(entry(size), N * M).prop_map(|v| matrix_from_values(&v))
    }

    fn test_numerical_matrix<const N: usize>(m: &Matrix<f64, N, N>) {
        let zero = Matrix::new();
        let one = Matrix::identity();

        assert_eq!(m.determinant().is_zero(), m.rank() < N);
        assert_eq!(m.null_space().len(), N - m.rank());

        for v in m.null_space() {
            assert!(allclose(&(m * v), &zero, 1e-9, 1e-9));
        }

        if let Some(inv) = m.inverse() {
            assert!(allclose(&(m * inv), &one, 1e-9, 1e-9));
        }

        assert_eq!(m.inverse().is_some(), m.rank() == N);
    }

    fn test_exact_matrix<T, const N: usize>(m: &Matrix<T, N, N>)
        where
            T: Entry + Clone + PartialEq + std::fmt::Debug,
            for <'a> &'a T: ScalarPtr<T>
    {
        let zero = Matrix::new();
        let one = Matrix::identity();

        assert_eq!(m.determinant().is_zero(), m.rank() < N);
        assert_eq!(m.null_space().len(), N - m.rank());

        for v in m.null_space() {
            assert_eq!(m * v, zero);
        }

        if let Some(inv) = m.inverse() {
            assert_eq!(m * inv, one);
        }
    }

    fn test_rational_matrix<T, const N: usize>(m: &Matrix<T, N, N>)
        where
            T: Entry + Clone + PartialEq + std::fmt::Debug,
            for <'a> &'a T: ScalarPtr<T>
    {
        test_exact_matrix(m);

        assert_eq!(m.inverse().is_some(), m.rank() == N);
    }

    fn test_solver<T, const N: usize, const M: usize>(
        m: &Matrix<T, N, M>, v: &Matrix<T, M, 1>
    )
        where
            T: Entry + Clone + PartialEq + std::fmt::Debug,
            for <'a> &'a T: ScalarPtr<T>
    {
        let b = m * v;
        if let Some(sol) = m.solve(&b) {
            assert_eq!(m * sol, b);
        }
    }

    proptest! {
        #[test]
        fn test_matrix_2i(m in matrix::<i64, 2, 2>(10)) {
            test_exact_matrix(&m);
        }

        #[test]
        fn test_matrix_2i_small(m in matrix::<i64, 2, 2>(2)) {
            test_exact_matrix(&m);
        }

        #[test]
        fn test_solver_2i(
            m in matrix::<i64, 2, 2>(10),
            v in matrix::<i64, 2, 1>(10)
        ) {
            test_solver(&m, &v);
        }

        #[test]
        fn test_solver_2i_small(
            m in matrix::<i64, 2, 2>(2),
            v in matrix::<i64, 2, 1>(2)
        ) {
            test_solver(&m, &v);
        }

        #[test]
        fn test_matrix_3i(m in matrix::<i64, 3, 3>(10)) {
            test_exact_matrix(&m);
        }

        #[test]
        fn test_matrix_3i_small(m in matrix::<i64, 3, 3>(2)) {
            test_exact_matrix(&m);
        }

        #[test]
        fn test_solver_3i(
            m in matrix::<i64, 3, 3>(10),
            v in matrix::<i64, 3, 1>(10)
        ) {
            test_solver(&m, &v);
        }

        #[test]
        fn test_solver_3i_small(
            m in matrix::<i64, 3, 3>(2),
            v in matrix::<i64, 3, 1>(2)
        ) {
            test_solver(&m, &v);
        }

        #[test]
        fn test_matrix_4i(m in matrix::<i64, 4, 4>(10)) {
            test_exact_matrix(&m);
        }

        #[test]
        fn test_matrix_4i_small(m in matrix::<i64, 4, 4>(2)) {
            test_exact_matrix(&m);
        }

        #[test]
        fn test_solver_4i(
            m in matrix::<i64, 4, 4>(10),
            v in matrix::<i64, 4, 1>(10)
        ) {
            test_solver(&m, &v);
        }

        #[test]
        fn test_solver_4i_small(
            m in matrix::<i64, 4, 4>(2),
            v in matrix::<i64, 4, 1>(2)
        ) {
            test_solver(&m, &v);
        }
    }

    proptest! {
        #[test]
        fn test_matrix_2q(m in matrix::<BigRational, 2, 2>(10)) {
            test_rational_matrix(&m);
        }

        #[test]
        fn test_matrix_2q_small(m in matrix::<BigRational, 2, 2>(2)) {
            test_rational_matrix(&m);
        }

        #[test]
        fn test_solver_2q(
            m in matrix::<BigRational, 2, 2>(10),
            v in matrix::<BigRational, 2, 1>(10)
        ) {
            test_solver(&m, &v);
        }

        #[test]
        fn test_solver_2q_small(
            m in matrix::<BigRational, 2, 2>(2),
            v in matrix::<BigRational, 2, 1>(2)
        ) {
            test_solver(&m, &v);
        }

        #[test]
        fn test_matrix_3q(m in matrix::<BigRational, 3, 3>(10)) {
            test_rational_matrix(&m);
        }

        #[test]
        fn test_matrix_3q_small(m in matrix::<BigRational, 3, 3>(2)) {
            test_rational_matrix(&m);
        }

        #[test]
        fn test_solver_3q(
            m in matrix::<BigRational, 3, 3>(10),
            v in matrix::<BigRational, 3, 1>(10)
        ) {
            test_solver(&m, &v);
        }

        #[test]
        fn test_solver_3q_small(
            m in matrix::<BigRational, 3, 3>(2),
            v in matrix::<BigRational, 3, 1>(2)
        ) {
            test_solver(&m, &v);
        }

        #[test]
        fn test_matrix_4q(m in matrix::<BigRational, 4, 4>(10)) {
            test_rational_matrix(&m);
        }

        #[test]
        fn test_matrix_4q_small(m in matrix::<BigRational, 4, 4>(2)) {
            test_rational_matrix(&m);
        }
    }

    proptest! {
        #[test]
        fn test_matrix_2f(m in matrix::<f64, 2, 2>(10)) {
            test_numerical_matrix(&m);
        }

        #[test]
        fn test_matrix_2f_small(m in matrix::<f64, 2, 2>(2)) {
            test_numerical_matrix(&m);
        }

        #[test]
        fn test_matrix_3f(m in matrix::<f64, 3, 3>(10)) {
            test_numerical_matrix(&m);
        }

        #[test]
        fn test_matrix_3f_small(m in matrix::<f64, 3, 3>(2)) {
            test_numerical_matrix(&m);
        }

        #[test]
        fn test_matrix_4f(m in matrix::<f64, 4, 4>(10)) {
            test_numerical_matrix(&m);
        }

        #[test]
        fn test_matrix_4f_small(m in matrix::<f64, 4, 4>(2)) {
            test_numerical_matrix(&m);
        }
    }
}
