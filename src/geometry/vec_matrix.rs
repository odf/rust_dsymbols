use std::ops::{Add, Div, Index, IndexMut, Mul};

use crate::geometry::traits::{Entry, Scalar, ScalarPtr, Array2d};
use crate::geometry::matrix::Matrix;


#[derive(Clone, Debug, PartialEq)]
pub struct VecMatrix<T> {
    data: Vec<T>,
    nr_rows: usize,
    nr_cols: usize,
}


impl<T> Array2d<T> for VecMatrix<T> {
    fn nr_rows(&self) -> usize {
        self.nr_rows
    }

    fn nr_columns(&self) -> usize {
        self.nr_cols
    }
}


impl<T> Index<usize> for VecMatrix<T>
{
    type Output = [T];

    fn index(&self, index: usize) -> &Self::Output {
        assert!(index < self.nr_rows);
        &self.data[index * self.nr_cols .. (index + 1) * self.nr_cols]
    }
}


impl<T> Index<(usize, usize)> for VecMatrix<T>
{
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (i, j) = index;
        assert!(i < self.nr_rows);
        assert!(j < self.nr_cols);

        &self.data[i * self.nr_cols + j]
    }
}


impl<T> IndexMut<usize> for VecMatrix<T>
{
    fn index_mut(&mut self, index: usize) -> &mut [T] {
        assert!(index < self.nr_rows);
        &mut self.data[index * self.nr_cols .. (index + 1) * self.nr_cols]
    }
}


impl<T> IndexMut<(usize, usize)> for VecMatrix<T>
{
    fn index_mut(&mut self, index: (usize, usize)) -> &mut T {
        let (i, j) = index;
        assert!(i < self.nr_rows);
        assert!(j < self.nr_cols);

        &mut self.data[i * self.nr_cols + j]
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize>
    From<Matrix<T, N, M>> for VecMatrix<T>
{
    fn from(data: Matrix<T, N, M>) -> Self {
        let mut result = VecMatrix::new(N, M);

        for i in 0..N {
            result[i].clone_from_slice(&data[i]);
        }

        result
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize>
    From<[[T; M]; N]> for VecMatrix<T>
{
    fn from(data: [[T; M]; N]) -> Self {
        Matrix::from(data).into()
    }
}


impl<T: Scalar + Clone, const M: usize> From<[T; M]> for VecMatrix<T> {
    fn from(data: [T; M]) -> Self {
        Matrix::from(data).into()
    }
}


impl<T: Scalar + Clone> From<T> for VecMatrix<T> {
    fn from(data: T) -> Self {
        Matrix::from(data).into()
    }
}


impl<T: Scalar + Clone> From<&[T]> for VecMatrix<T> {
    fn from(data: &[T]) -> Self {
        let mut result = VecMatrix::new(1, data.len());

        result[0].clone_from_slice(&data);

        result
    }
}


impl<T: Scalar + Clone> From<Vec<T>> for VecMatrix<T> {
    fn from(data: Vec<T>) -> Self {
        let mut result = VecMatrix::new(1, data.len());

        result[0].clone_from_slice(&data);

        result
    }
}


impl<T: Scalar + Clone> VecMatrix<T> {
    pub fn new(nr_rows: usize, nr_cols: usize) -> Self {
        VecMatrix {
            data: vec![T::zero(); nr_rows * nr_cols],
            nr_rows,
            nr_cols
        }
    }

    pub fn transpose(&self) -> VecMatrix<T> {
        let mut result = VecMatrix::new(self.nr_cols, self.nr_rows);

        for i in 0..self.nr_cols {
            for j in 0..self.nr_rows {
                result[i][j] = self[j][i].clone();
            }
        }

        result
    }

    pub fn swap_rows(&mut self, i: usize, j: usize) {
        assert!(i < self.nr_rows);
        assert!(j < self.nr_rows);
        assert_ne!(i, j);

        let ri = i * self.nr_cols;
        let rj = j * self.nr_cols;

        for k in 0..self.nr_cols {
            self.data.swap(ri + k, rj + k);
        }
    }

    pub fn hstack(self, rhs: &VecMatrix<T>) -> VecMatrix<T> {
        assert_eq!(self.nr_rows, rhs.nr_rows);

        let m = self.nr_cols;
        let mut result = VecMatrix::new(self.nr_rows, self.nr_cols + rhs.nr_cols);

        for i in 0..self.nr_rows {
            result[i][0..m].clone_from_slice(&self[i]);
            result[i][m..(m + rhs.nr_cols)].clone_from_slice(&rhs[i]);
        }
        result
    }

    pub fn vstack(self, rhs: &VecMatrix<T>) -> VecMatrix<T> {
        assert_eq!(self.nr_cols, rhs.nr_cols);

        let mut result = VecMatrix::new(self.nr_rows + rhs.nr_rows, self.nr_cols);

        for i in 0..self.nr_rows {
            result[i].clone_from_slice(&self[i]);
        }
        for i in 0..rhs.nr_rows {
            result[self.nr_rows + i].clone_from_slice(&rhs[i]);
        }

        result
    }

    pub fn submatrix<I, J>(&self, rows: I, columns: J) -> VecMatrix<T>
        where
            I: IntoIterator<Item=usize>,
            J: IntoIterator<Item=usize>
    {
        let rows: Vec<_> = rows.into_iter().collect();
        let columns: Vec<_> = columns.into_iter().collect();

        assert!(rows.iter().all(|&i| i < self.nr_rows));
        assert!(columns.iter().all(|&j| j < self.nr_cols));

        let mut result = VecMatrix::new(rows.len(), columns.len());

        for i in 0..rows.len() {
            for j in 0..columns.len() {
                result[i][j] = self[rows[i]][columns[j]].clone();
            }
        }

        result
    }
}


impl<T: Scalar + Clone> VecMatrix<T> {
    pub fn identity(dim: usize) -> Self {
        let mut result = VecMatrix::new(dim, dim);
        for i in 0..dim {
            result[i][i] = T::one();
        }
        result
    }
}


impl<T: Scalar + Clone> Add<&VecMatrix<T>> for &VecMatrix<T>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = VecMatrix<T>;

    fn add(self, rhs: &VecMatrix<T>) -> Self::Output {
        assert_eq!(self.nr_rows, rhs.nr_rows);
        assert_eq!(self.nr_cols, rhs.nr_cols);

        let mut result = VecMatrix::new(self.nr_rows, self.nr_cols);

        for i in 0..self.nr_rows {
            for j in 0..self.nr_cols {
                result[i][j] = &self[i][j] + &rhs[i][j];
            }
        }

        result
    }
}


impl<T: Scalar + Clone> Add<&VecMatrix<T>> for VecMatrix<T>
where for <'a> &'a T: ScalarPtr<T>
{
    type Output = VecMatrix<T>;

    fn add(self, rhs: &VecMatrix<T>) -> Self::Output {
        &self + rhs
    }
}


impl<T: Scalar + Clone> Add<VecMatrix<T>> for &VecMatrix<T>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = VecMatrix<T>;

    fn add(self, rhs: VecMatrix<T>) -> Self::Output {
        self + &rhs
    }
}


impl<T: Scalar + Clone> Add<VecMatrix<T>> for VecMatrix<T>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = VecMatrix<T>;

    fn add(self, rhs: VecMatrix<T>) -> Self::Output {
        self + &rhs
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize>
    Add<[[T; M]; N]> for &VecMatrix<T>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = VecMatrix<T>;

    fn add(self, rhs: [[T; M]; N]) -> Self::Output {
        self + VecMatrix::from(rhs)
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize>
    Add<[[T; M]; N]> for VecMatrix<T>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = VecMatrix<T>;

    fn add(self, rhs: [[T; M]; N]) -> Self::Output {
        self + VecMatrix::from(rhs)
    }
}


impl<T: Scalar + Clone> VecMatrix<T> {
    pub fn zero(nr_rows: usize, nr_cols: usize) -> Self {
        VecMatrix::new(nr_rows, nr_cols)
    }

    pub fn is_zero(&self) -> bool {
        for i in 0..self.nr_rows {
            for j in 0..self.nr_cols {
                if !self[i][j].is_zero() {
                    return false
                }
            }
        }

        true
    }
}


impl<T: Scalar + Clone> Mul<&VecMatrix<T>> for &VecMatrix<T>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = VecMatrix<T>;

    fn mul(self, rhs: &VecMatrix<T>) -> Self::Output {
        let mut result = VecMatrix::new(self.nr_rows, rhs.nr_cols);

        for i in 0..self.nr_rows {
            for j in 0..rhs.nr_cols {
                let mut x = T::zero();
                for k in 0..self.nr_cols {
                    x = x + &self[i][k] * &rhs[k][j];
                }
                result[i][j] = x;
            }
        }

        result
    }
}


impl<T: Scalar + Clone> Mul<&VecMatrix<T>> for VecMatrix<T>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = VecMatrix<T>;

    fn mul(self, rhs: &VecMatrix<T>) -> Self::Output {
        &self * rhs
    }
}


impl<T: Scalar + Clone> Mul<VecMatrix<T>> for &VecMatrix<T>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = VecMatrix<T>;

    fn mul(self, rhs: VecMatrix<T>) -> Self::Output {
        self * &rhs
    }
}


impl<T: Scalar + Clone> Mul<VecMatrix<T>> for VecMatrix<T>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = VecMatrix<T>;

    fn mul(self, rhs: VecMatrix<T>) -> Self::Output {
        self * &rhs
    }
}


impl<T: Scalar + Clone> Mul<&VecMatrix<T>> for &[T]
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = Vec<T>;

    fn mul(self, rhs: &VecMatrix<T>) -> Self::Output {
        assert_eq!(self.len(), rhs.nr_rows());
        (VecMatrix::from(self) * rhs)[0].into()
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize>
    Mul<&VecMatrix<T>> for [[T; M]; N]
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = VecMatrix<T>;

    fn mul(self, rhs: &VecMatrix<T>) -> Self::Output {
        VecMatrix::from(self) * rhs
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize>
    Mul<VecMatrix<T>> for [[T; M]; N]
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = VecMatrix<T>;

    fn mul(self, rhs: VecMatrix<T>) -> Self::Output {
        VecMatrix::from(self) * rhs
    }
}


impl<T: Scalar + Clone, const M: usize, const L: usize>
    Mul<[[T; L]; M]> for VecMatrix<T>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = VecMatrix<T>;

    fn mul(self, rhs: [[T; L]; M]) -> Self::Output {
        self * VecMatrix::from(rhs)
    }
}


impl Mul<VecMatrix<i64>> for i64 {
    type Output = VecMatrix<i64>;

    fn mul(self, rhs: VecMatrix<i64>) -> Self::Output {
        let mut result = VecMatrix::new(rhs.nr_rows, rhs.nr_cols);

        for i in 0..rhs.nr_rows {
            for j in 0..rhs.nr_cols {
                result[i][j] = self * rhs[i][j];
            }
        }

        result
    }
}


impl Mul<VecMatrix<f64>> for f64 {
    type Output = VecMatrix<f64>;

    fn mul(self, rhs: VecMatrix<f64>) -> Self::Output {
        let mut result = VecMatrix::new(rhs.nr_rows, rhs.nr_cols);

        for i in 0..rhs.nr_rows {
            for j in 0..rhs.nr_cols {
                result[i][j] = self * rhs[i][j];
            }
        }

        result
    }
}


impl<T: Scalar + Clone> Mul<&T> for &VecMatrix<T>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = VecMatrix<T>;

    fn mul(self, rhs: &T) -> Self::Output {
        let mut result = VecMatrix::new(self.nr_rows, self.nr_cols);

        for i in 0..self.nr_rows {
            for j in 0..self.nr_cols {
                result[i][j] = &self[i][j] * rhs;
            }
        }

        result
    }
}


impl<T: Scalar + Clone> Mul<&T> for VecMatrix<T>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = Self;

    fn mul(self, rhs: &T) -> Self::Output {
        &self * rhs
    }
}


impl<T: Scalar + Clone> Mul<T> for VecMatrix<T>
    where for <'a> &'a T: ScalarPtr<T>
{
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        &self * &rhs
    }
}


#[derive(Debug)]
pub struct RowEchelonVecMatrix<T: Entry> {
    multiplier: VecMatrix<T>,
    result: VecMatrix<T>,
    columns: Vec<usize>,
    rank: usize,
    nr_swaps: usize
}


impl<T: Entry + Clone> RowEchelonVecMatrix<T> {
    pub fn new(m: &VecMatrix<T>) -> Self {
        let mut u = m.clone();
        let mut s = VecMatrix::identity(m.nr_rows());
        let mut row = 0;
        let mut nr_swaps = 0;
        let mut cols = vec![m.nr_rows(); m.nr_rows()];

        for col in 0..m.nr_columns() {
            if let Some(pr) = Entry::pivot_row(col, row, &u) {
                if pr != row {
                    u.swap_rows(pr, row);
                    s.swap_rows(pr, row);
                    nr_swaps += 1;
                }

                for r in (row + 1)..m.nr_rows() {
                    Entry::clear_col(col, r, row, &mut u, Some(&mut s));
                }

                cols[row] = col;
                row += 1;
            }
        }

        RowEchelonVecMatrix {
            multiplier: s,
            result: u,
            columns: cols,
            rank: row,
            nr_swaps
        }
    }
}


impl<T: Entry + Clone> VecMatrix<T>
    where for <'a> &'a T: ScalarPtr<T>
{
    pub fn rank(&self) -> usize {
        RowEchelonVecMatrix::new(self).rank
    }

    pub fn null_space(&self) -> Vec<VecMatrix<T>> {
        let re = RowEchelonVecMatrix::new(&self.transpose());
        let s = re.multiplier.transpose();

        (re.rank..self.nr_columns())
            .map(|i| s.submatrix(0..s.nr_rows(), [i]))
            .collect()
    }

    pub fn null_space_matrix(&self) -> VecMatrix<T> {
        let re = RowEchelonVecMatrix::new(&self.transpose());
        let s = re.multiplier.transpose();

        s.submatrix(0..s.nr_rows(), re.rank..self.nr_columns())
    }

    pub fn solve(&self, rhs: &VecMatrix<T>) -> Option<VecMatrix<T>>
        where T: Div<T, Output=T>
    {
        let re = RowEchelonVecMatrix::new(self);
        let y = re.multiplier * rhs;

        if !(re.rank..self.nr_rows()).all(|i|
            (0..rhs.nr_columns()).all(|j| y[i][j].is_zero())
        ) {
            return None;
        }

        let mut result = VecMatrix::zero(self.nr_columns(), rhs.nr_columns());

        for row in (0..re.rank).rev() {
            let a = &re.result[row] * &result;
            let b = &y[row];
            let x = &re.result[row][re.columns[row]];
            for k in 0..rhs.nr_columns() {
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


impl<T: Entry + Clone> VecMatrix<T>
    where for <'a> &'a T: ScalarPtr<T>
{
    pub fn determinant(&self) -> T {
        assert_eq!(self.nr_rows(), self.nr_columns());

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
                let re = RowEchelonVecMatrix::new(self);
                let mut result = T::one();
                for i in 0..self.nr_rows {
                    result = &result * &re.result[i][i];
                }
                if re.nr_swaps % 2 == 0 { result } else { -result }
            }
        }
    }

    pub fn inverse(&self) -> Option<Self>
        where T: Div<T, Output=T>
    {
        assert_eq!(self.nr_rows(), self.nr_columns());

        self.solve(&VecMatrix::identity(self.nr_rows()))
    }
}


#[test]
fn test_matrix_clone() {
    let m1 = VecMatrix::from([[1, 2], [3, 4]]);
    let mut m2 = m1.clone();
    m2[0][0] = 0;
    assert_eq!(&m1 + &m2, VecMatrix::from([[1, 4], [6, 8]]));
    assert_eq!(m1 + (-1 * m2), VecMatrix::from([[1, 0], [0, 0]]));
}


#[test]
fn test_matrix_indexing() {
    let mut m = VecMatrix::from([[1.0, 1.0], [0.0, 1.0]]);

    m[1][0] = m[1][0] + 4.0;

    assert_eq!(m, [[1.0, 1.0], [4.0, 1.0]].into());
    assert_eq!(m[1], [4.0, 1.0]);
}


#[test]
fn test_matrix_identity() {
    let m = VecMatrix::<i64>::identity(3);
    assert_eq!(m, [[1, 0, 0], [0, 1, 0], [0, 0, 1]].into());
    assert_eq!(VecMatrix::identity(2), [[1, 0], [0, 1]].into());
}


#[test]
fn test_matrix_add() {
    assert_eq!(
        &VecMatrix::from([[1, 2], [3, 4]]) + &VecMatrix::from([[1, 2], [3, 4]]),
        VecMatrix::from([[2, 4], [6, 8]])
    );
    assert_eq!(
        &VecMatrix::from([[1, 2], [3, 4]]) + VecMatrix::from([[1, 2], [3, 4]]),
        VecMatrix::from([[2, 4], [6, 8]])
    );
    assert_eq!(
        VecMatrix::from([[1, 2], [3, 4]]) + &VecMatrix::from([[1, 2], [3, 4]]),
        VecMatrix::from([[2, 4], [6, 8]])
    );

    assert_eq!(
        VecMatrix::from([[1, 2], [3, 4]]) + VecMatrix::from([[1, 2], [3, 4]]),
        VecMatrix::from([[2, 4], [6, 8]])
    );

    assert_eq!(
        &VecMatrix::from([[1, 2], [3, 4]]) + [[4, 3], [2, 1]],
        VecMatrix::from([[5, 5], [5, 5]])
    );
    assert_eq!(
        VecMatrix::from([[1, 2], [3, 4]]) + [[4, 3], [2, 1]],
        VecMatrix::from([[5, 5], [5, 5]])
    );
}


#[test]
fn test_matrix_stack() {
    assert_eq!(
        VecMatrix::from([[1, 2, 3]])
            .vstack(&VecMatrix::from([[4, 5, 6], [7, 8, 9]])),
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]].into()
    );
    assert_eq!(
        VecMatrix::from([[1, 2], [3, 4]])
            .hstack(&VecMatrix::from([[5, 6], [7, 8]])),
        [[1, 2, 5, 6], [3, 4, 7, 8]].into()
    );
}


#[test]
fn test_matrix_submatrix() {
    assert_eq!(
        VecMatrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).submatrix(0..2, [0, 2]),
        VecMatrix::from([[1, 3], [4, 6]])
    )
}


#[test]
fn test_matrix_zero() {
    assert_eq!(VecMatrix::zero(2, 3), VecMatrix::from([[0; 3]; 2]));
    assert!(VecMatrix::from([[0; 4]; 5]).is_zero());
}


#[test]
fn test_matrix_mul() {
    assert_eq!(
        (VecMatrix::from([[1, 2, 3]]) * [[3], [2], [1]])[0][0],
        10
    );
    assert_eq!(
        VecMatrix::from([1, 2, 3]) * VecMatrix::from([3, 2, 1]).transpose(),
        10.into()
    );
    assert_eq!(
        [[1, 2, 3]] * VecMatrix::from([3, 2, 1]).transpose(),
        [[10]].into()
    );
    assert_eq!(
        VecMatrix::from([[1.0, 1.0], [0.0, 1.0]]) * [[1.0, 2.0], [0.0, 1.0]],
        [[1.0, 3.0], [0.0, 1.0]].into()
    );
    assert_eq!(
        (
            VecMatrix::from([[1.0, 1.0], [0.0, 1.0]]) *
            [[1.0, 2.0], [0.0, 1.0]]
        )[0][1],
        3.0
    );
    assert_eq!(
        (
            VecMatrix::from([[1.0, 1.0], [0.0, 1.0]]) *
            [[1.0, 2.0], [0.0, 1.0]]
        )[0],
        [1.0, 3.0]
    );
    assert_eq!(
        (VecMatrix::from([[1.0, 1.0]]) * [[1.0, 2.0], [0.0, 1.0]])[0],
        [1.0, 3.0]
    );
    assert_eq!(
        VecMatrix::from([[1, 2], [3, 4]]) * 2,
        VecMatrix::from([[2, 4], [6, 8]])
    );
    assert_eq!(
        VecMatrix::from([[1.0, 2.0], [3.0, 4.0]]) * 2.0,
        VecMatrix::from([[2.0, 4.0], [6.0, 8.0]])
    );
    assert_eq!(
        2.5 * VecMatrix::from([[1.0, 2.0], [3.0, 4.0]]),
        VecMatrix::from([[2.5, 5.0], [7.5, 10.0]])
    );
}


#[test]
fn test_matrix_determinant() {
    assert_eq!(
        VecMatrix::from([[1.0, 0.3, 0.7], [0.0, 2.0, 1.2], [0.0, 0.0, 0.25]])
            .determinant(),
        0.5
    );

    assert_eq!(
        VecMatrix::from([
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 0, 1],
            [0, 0, 1, 0]
        ]).determinant(),
        -1
    );

    assert_eq!(
        VecMatrix::from([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 0.0]
        ]).determinant(),
        -1.0
    );

    assert_eq!(VecMatrix::<f64>::identity(1).determinant(), 1.0);
    assert_eq!(VecMatrix::<f64>::identity(2).determinant(), 1.0);
    assert_eq!(VecMatrix::<f64>::identity(3).determinant(), 1.0);
    assert_eq!(VecMatrix::<f64>::identity(4).determinant(), 1.0);

    assert_eq!((VecMatrix::<f64>::identity(1) * 2.0).determinant(), 2.0);
    assert_eq!((VecMatrix::<f64>::identity(2) * 2.0).determinant(), 4.0);
    assert_eq!((VecMatrix::<f64>::identity(3) * 2.0).determinant(), 8.0);
    assert_eq!((VecMatrix::<f64>::identity(4) * 2.0).determinant(), 16.0);
}

#[test]
fn test_matrix_nullspace() {
    let a = VecMatrix::from([[1, 2], [3, 4]]);
    let n = a.null_space();
    assert_eq!(n, vec![]);

    let a = VecMatrix::from([[1.0, 2.0], [3.0, 6.0]]);
    let n = a.null_space();
    assert_eq!(n.len(), 1);
    for v in a.null_space() {
        assert_eq!(&a * v, VecMatrix::from([[0.0], [0.0]]));
    }
    assert!((&a * a.null_space_matrix()).is_zero());

    let a = VecMatrix::from([[1.0, 2.0], [3.0, 4.0]]);
    let n = a.null_space();
    assert_eq!(n, vec![]);

    let a = VecMatrix::from([[0.0, 1.0, 0.0]]);
    let n = a.null_space();
    assert_eq!(n.len(), 2);
    for v in n {
        assert_eq!(&a * v, VecMatrix::from([[0.0]]));
    }
    assert!((&a * a.null_space_matrix()).is_zero());
}


#[test]
fn test_matrix_solve_i64() {
    fn test<const N: usize, const M: usize, const L: usize>(
        a: [[i64; M]; N], b: [[i64; L]; M]
    )
    {
        let a = VecMatrix::from(a);
        let b = VecMatrix::from(b);
        let ab = &a * &b;

        assert_eq!(&a * a.solve(&ab).unwrap(), ab);

        if N == M && a.rank() == N {
            assert_eq!(a.solve(&ab), Some(b));
        }
    }

    test([[1, 2], [3, 4]], [[1, 0], [1, -3]]);
    test([[1, 2], [3, 6]], [[1], [1]]);
    test([[0, 0], [0, 1]], [[0], [-1]]);
}


#[test]
fn test_matrix_solve_f64() {
    use crate::geometry::traits::allclose;

    fn test<const N: usize, const M: usize, const L: usize>(
        a: [[f64; M]; N], b: [[f64; L]; M]
    )
    {
        let a = VecMatrix::from(a);
        let b = VecMatrix::from(b);
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
        let a = VecMatrix::from(a);

        if a.rank() == N {
            assert_eq!(&a * a.inverse().unwrap(), VecMatrix::identity(N));
        } else {
            assert_eq!(a.inverse(), None);
        }
    }

    test([[1, 2], [3, 5]]);

    // Inverse exists, but is not integral:
    assert_eq!(VecMatrix::from([[1, 2], [3, 4]]).inverse(), None);
}


#[test]
fn test_matrix_inverse_f64() {
    use crate::geometry::traits::allclose;

    fn test<const N: usize>(a: [[f64; N]; N])
    {
        let a = VecMatrix::from(a);

        if a.rank() == N {
            assert!(allclose(
                &(&a * a.inverse().unwrap()),
                &VecMatrix::identity(N),
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
        -> VecMatrix<BigRational>
    {
        let mut result = VecMatrix::new(N, M);

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

        let i1 = &VecMatrix::<BigRational>::identity(1);
        let i2 = &VecMatrix::<BigRational>::identity(2);
        let i3 = &VecMatrix::<BigRational>::identity(3);
        let i4 = &VecMatrix::<BigRational>::identity(4);

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
        assert_eq!(n.len(), 0);

        let a = matrix([[1, 2], [3, 6]]);
        let n = a.null_space();
        assert_eq!(n.len(), 1);
        for v in &n {
            assert_eq!(&a * v, VecMatrix::zero(2, 1));
        }
        assert!((&a * a.null_space_matrix()).is_zero());
    }

    #[test]
    fn test_inverse() {
        let a = matrix([[1, 2], [3, 4]]);
        assert_eq!(&a * &a.inverse().unwrap(), VecMatrix::identity(2));

        let a = matrix([[1, 2], [3, 5]]);
        assert_eq!(&a * &a.inverse().unwrap(), VecMatrix::identity(2));

        let a = matrix([[1, 2], [3, 6]]);
        assert!(a.inverse().is_none());
    }
}


mod property_based_tests {
    use super::*;
    use num_rational::BigRational;
    use num_traits::{FromPrimitive, Zero};
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
        where T: FromI32 + std::fmt::Debug
    {
        (0..size).prop_map(T::from)
    }

    fn sized_matrix<T>(size: i32, n: usize, m: usize)
        -> impl Strategy<Value=VecMatrix<T>>
        where T: Scalar + Clone + FromI32 + std::fmt::Debug + 'static
    {
        vec(entry(size), n * m).prop_map(move |v| matrix_from_values(&v, m))
    }

    fn matrix<T>(entry_size: i32, dmin: usize, dmax: usize)
        -> impl Strategy<Value=VecMatrix<T>>
        where T: Scalar + Clone + FromI32 + std::fmt::Debug + 'static
    {
        (dmin..=dmax).prop_flat_map(move |n|
            sized_matrix(entry_size, n, n)
        )
    }

    fn equations<T>(entry_size: i32, dmin: usize, dmax: usize)
        -> impl Strategy<Value=(VecMatrix<T>, VecMatrix<T>)>
        where T: Scalar + Clone + FromI32 + std::fmt::Debug + 'static
    {
        (dmin..=dmax).prop_flat_map(move |n|
            (sized_matrix(entry_size, n, n), sized_matrix(entry_size, n, 1))
        )
    }

    fn test_numerical_matrix(m: &VecMatrix<f64>) {
        let n = m.nr_rows();
        let zero = VecMatrix::new(n, 1);
        let one = VecMatrix::identity(n);

        assert_eq!(m.determinant().is_zero(), m.rank() < n);
        assert_eq!(m.null_space().len(), n - m.rank());

        for v in m.null_space() {
            assert!(allclose(&(m * v), &zero, 1e-9, 1e-9));
        }

        let nul = m.null_space_matrix();
        assert!(allclose(
            &(m * &nul),
            &VecMatrix::new(n, nul.nr_columns()),
            1e-9, 1e-9
        ));

        if let Some(inv) = m.inverse() {
            assert!(allclose(&(m * inv), &one, 1e-9, 1e-9));
        }

        assert_eq!(m.inverse().is_some(), m.rank() == n);
    }

    fn test_exact_matrix<T>(m: &VecMatrix<T>)
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
    }

    fn test_rational_matrix<T>(m: &VecMatrix<T>)
        where
            T: Entry + Clone + PartialEq + std::fmt::Debug,
            for <'a> &'a T: ScalarPtr<T>
    {
        test_exact_matrix(m);

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

    fn test_solver_numerical(m: &VecMatrix<f64>, v: &VecMatrix<f64>) {
        let b = m * v;
        if let Some(sol) = m.solve(&b) {
            assert!(allclose(&(m * sol), &b, 1e-9, 1e-9));
        }
    }

    proptest! {
        #[test]
        fn test_matrix_int(m in matrix::<i64>(10, 2, 4)) {
            test_exact_matrix(&m);
        }

        #[test]
        fn test_matrix_int_small(m in matrix::<i64>(3, 2, 4)) {
            test_exact_matrix(&m);
        }

        #[test]
        fn test_solver_int((m, v) in equations::<i64>(10, 2, 4)) {
            test_solver(&m, &v);
        }

        #[test]
        fn test_solver_int_small((m, v) in equations::<i64>(3, 2, 4)) {
            test_solver(&m, &v);
        }
    }

    proptest! {
        #[test]
        fn test_matrix_rational(m in matrix::<BigRational>(10, 2, 6)) {
            test_rational_matrix(&m);
        }

        #[test]
        fn test_matrix_rational_small(m in matrix::<BigRational>(3, 2, 6)) {
            test_rational_matrix(&m);
        }

        #[test]
        fn test_solver_rational((m, v) in equations::<BigRational>(10, 2, 6)) {
            test_solver(&m, &v);
        }

        #[test]
        fn test_solver_rational_small((m, v) in equations::<BigRational>(3, 2, 6)) {
            test_solver(&m, &v);
        }
    }

    proptest! {
        #[test]
        fn test_matrix_float(m in matrix::<f64>(10, 2, 6)) {
            test_numerical_matrix(&m);
        }

        #[test]
        fn test_matrix_float_small(m in matrix::<f64>(3, 2, 6)) {
            test_numerical_matrix(&m);
        }

        #[test]
        fn test_solver_float((m, v) in equations::<f64>(10, 2, 6)) {
            test_solver_numerical(&m, &v);
        }

        #[test]
        fn test_solver_float_small((m, v) in equations::<f64>(3, 2, 6)) {
            test_solver_numerical(&m, &v);
        }
    }
}
