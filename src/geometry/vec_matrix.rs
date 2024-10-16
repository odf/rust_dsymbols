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
            for j in 0..M {
                result[i][j] = data[i][j].clone();
            }
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

        for j in 0..data.len() {
            result[0][j] = data[j].clone();
        }

        result
    }
}


impl<T: Scalar + Clone> From<Vec<T>> for VecMatrix<T> {
    fn from(data: Vec<T>) -> Self {
        let mut result = VecMatrix::new(1, data.len());

        for j in 0..data.len() {
            result[0][j] = data[j].clone();
        }

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

    pub fn get_row(&self, i: usize) -> VecMatrix<T> {
        assert!(i < self.nr_rows);
        VecMatrix::from(&self[i])
    }

    pub fn set_row(&mut self, i: usize, row: &VecMatrix<T>) {
        assert_eq!(row.nr_rows, 1);
        assert_eq!(row.nr_cols, self.nr_cols);
        assert!(i < self.nr_rows);

        for j in 0..self.nr_cols {
            self[i][j] = row[0][j].clone();
        }
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

    pub fn get_column(&self, j: usize) -> VecMatrix<T> {
        assert!(j < self.nr_cols);

        let mut result = VecMatrix::new(self.nr_rows, 1);

        for i in 0..self.nr_rows {
            result[i][0] = self[i][j].clone();
        }

        result
    }

    pub fn set_column(&mut self, j: usize, col: &VecMatrix<T>) {
        assert_eq!(col.nr_cols, 1);
        assert_eq!(col.nr_rows, self.nr_rows);
        assert!(j < self.nr_cols);

        for i in 0..self.nr_rows {
            self[i][j] = col[i][0].clone();
        }
    }

    pub fn hstack(self, rhs: &VecMatrix<T>) -> VecMatrix<T> {
        assert_eq!(self.nr_rows, rhs.nr_rows);

        let mut result = VecMatrix::new(self.nr_rows, self.nr_cols + rhs.nr_cols);

        for i in 0..self.nr_rows {
            for j in 0..self.nr_cols {
                result[i][j] = self[i][j].clone();
            }
            for j in 0..rhs.nr_cols {
                result[i][self.nr_cols + j] = rhs[i][j].clone();
            }
        }
        result
    }

    pub fn vstack(self, rhs: &VecMatrix<T>) -> VecMatrix<T> {
        assert_eq!(self.nr_cols, rhs.nr_cols);

        let mut result = VecMatrix::new(self.nr_rows + rhs.nr_rows, self.nr_cols);

        for j in 0..self.nr_cols {
            for i in 0..self.nr_rows {
                result[i][j] = self[i][j].clone();
            }
            for i in 0..rhs.nr_rows {
                result[self.nr_rows + i][j] = rhs[i][j].clone();
            }
        }
        result
    }

    pub fn submatrix<I, J>(self, rows: I, columns: J) -> VecMatrix<T>
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

    fn is_zero(&self) -> bool {
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
    fn rank(&self) -> usize {
        RowEchelonVecMatrix::new(self).rank
    }

    fn null_space(&self) -> Vec<VecMatrix<T>> {
        let re = RowEchelonVecMatrix::new(&self.transpose());
        let s = re.multiplier;

        (re.rank..self.nr_columns())
            .map(|i| s.get_row(i).transpose())
            .collect()
    }

    fn solve(&self, rhs: &VecMatrix<T>) -> Option<VecMatrix<T>>
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
            let a = re.result.get_row(row) * &result;
            let b = y.get_row(row);
            let x = &re.result[row][re.columns[row]];
            for k in 0..rhs.nr_columns() {
                let t = &b[0][k] - &a[0][k];
                if Entry::can_divide(&t, &x) {
                    result[row][k] = &t / x;
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
                let re = RowEchelonVecMatrix::new(self);
                let mut result = T::one();
                for i in 0..self.nr_rows {
                    result = &result * &re.result[i][i];
                }
                if re.nr_swaps % 2 == 0 { result } else { -result }
            }
        }
    }

    fn inverse(&self) -> Option<Self>
        where T: Div<T, Output=T>
    {
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

    m.set_row(0, &(m.get_row(0) * 3.0));
    m[1][0] = m[1][0] + 4.0;

    assert_eq!(m, [[3.0, 3.0], [4.0, 1.0]].into());
}


#[test]
fn test_matrix_row_column_manipulation() {
    let mut m = VecMatrix::from([[1.0, 1.0], [0.0, 1.0]]);

    m.set_row(0, &(m.get_row(0) * 2.0));
    m.set_column(1, &(m.get_column(1) * 3.0));

    assert_eq!(m, [[2.0, 6.0], [0.0, 3.0]].into());
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

    let a = VecMatrix::from([[1.0, 2.0], [3.0, 4.0]]);
    let n = a.null_space();
    assert_eq!(n, vec![]);

    let a = VecMatrix::from([[0.0, 1.0, 0.0]]);
    let n = a.null_space();
    assert_eq!(n.len(), 2);
    for v in n {
        assert_eq!(&a * v, VecMatrix::from([[0.0]]));
    }
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
}
