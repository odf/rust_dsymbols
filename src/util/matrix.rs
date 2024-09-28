use std::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};
use num_traits::{One, Zero};


pub trait Scalar:
    Zero + One + Mul<Output=Self> + Add<Output=Self> + Neg<Output=Self>
{
}

impl Scalar for f64 {}
impl Scalar for i64 {}


#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Matrix<T, const N: usize, const M: usize> {
    data: [[T; M]; N]
}


impl<T: Scalar, const N: usize, const M: usize> Matrix<T, N, M> {
    pub fn nr_rows(&self) -> usize {
        N
    }

    pub fn nr_columns(&self) -> usize {
        M
    }
}


impl<T: Scalar, const N: usize, const M: usize>
    Index<usize> for Matrix<T, N, M>
{
    type Output = [T; M];

    fn index(&self, index: usize) -> &Self::Output {
        assert!(index < N);
        &self.data[index]
    }
}


impl<T: Scalar, const N: usize, const M: usize>
    IndexMut<usize> for Matrix<T, N, M>
{
    fn index_mut(&mut self, index: usize) -> &mut [T; M] {
        assert!(index < N);
        &mut self.data[index]
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


impl<T: Scalar>
    From<T> for Matrix<T, 1, 1>
{
    fn from(data: T) -> Self {
        Self { data: [[data]] }
    }
}


impl<T: Scalar + Copy , const N: usize, const M: usize> Matrix<T, N, M> {
    pub fn new() -> Self {
        Matrix::from([[T::zero(); M]; N])
    }

    pub fn transpose(&self) -> Matrix<T, M, N> {
        let mut result = Matrix::new();
        for i in 0..M {
            for j in 0..N {
                result[i][j] = self[j][i];
            }
        }
        result
    }

    pub fn get_row(&self, i: usize) -> Matrix<T, 1, M> {
        assert!(i < N);
        Matrix::from([self[i]])
    }

    pub fn set_row(&mut self, i: usize, row: Matrix<T, 1, M>) {
        assert!(i < N);
        self[i] = row[0];
    }

    pub fn get_column(&self, j: usize) -> Matrix<T, N, 1> {
        assert!(j < M);
        let mut result = Matrix::new();
        for i in 0..N {
            result[i][0] = self[i][j];
        }
        result
    }

    pub fn set_column(&mut self, j: usize, column: Matrix<T, N, 1>) {
        assert!(j < M);
        for i in 0..N {
            self[i][j] = column[i][0];
        }
    }

    pub fn hstack<const L: usize, const S: usize>(self, rhs: Matrix<T, N, L>)
        -> Matrix<T, N, S>
    {
        assert_eq!(S, N + M);

        let mut result = Matrix::new();

        for i in 0..N {
            for j in 0..M {
                result[i][j] = self[i][j];
            }
            for j in 0..L {
                result[i][M + j] = rhs[i][j];
            }
        }
        result
    }

    pub fn vstack<const L: usize, const S: usize>(self, rhs: Matrix<T, L, M>)
        -> Matrix<T, S, M>
    {
        assert_eq!(S, N + L);

        let mut result = Matrix::new();

        for j in 0..M {
            for i in 0..N {
                result[i][j] = self[i][j];
            }
            for i in 0..L {
                result[N + i][j] = rhs[i][j];
            }
        }
        result
    }

    pub fn submatrix<I, J, const K: usize, const L: usize>(
        self, rows: I, columns: J
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
                result[i][j] = self[rows[i]][columns[j]];
            }
        }

        result
    }
}


impl<T: Scalar + Copy, const N: usize> Matrix<T, N, N> {
    pub fn identity() -> Self {
        let mut result = Matrix::new();
        for i in 0..N {
            result[i][i] = T::one();
        }
        result
    }
}


impl<T: Scalar + Copy, const N: usize, const M: usize>
    Add<[[T; M]; N]> for Matrix<T, N, M>
{
    type Output = Matrix<T, N, M>;

    fn add(self, rhs: [[T; M]; N]) -> Self::Output {
        let mut result = Matrix::new();

        for i in 0..N {
            for j in 0..M {
                result[i][j] = self[i][j] + rhs[i][j];
            }
        }

        result
    }
}


impl<T: Scalar + Copy, const N: usize, const M: usize>
    Add<Matrix<T, N, M>> for Matrix<T, N, M>
{
    type Output = Matrix<T, N, M>;

    fn add(self, rhs: Matrix<T, N, M>) -> Self::Output {
        self + rhs.data
    }
}


impl<T: Scalar + Copy, const N: usize, const M: usize> Zero for Matrix<T, N, M> {
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


impl<T: Scalar + Copy, const N: usize, const M: usize, const L: usize>
    Mul<Matrix<T, M, L>> for [[T; M]; N]
{
    type Output = [[T; L]; N];

    fn mul(self, rhs: Matrix<T, M, L>) -> Self::Output {
        let mut result = Matrix::new();

        for i in 0..N {
            for j in 0..L {
                let mut x = T::zero();
                for k in 0..M {
                    x = x + self[i][k] * rhs[k][j];
                }
                result[i][j] = x;
            }
        }

        result.data
    }
}


impl<T: Scalar + Copy, const N: usize, const M: usize, const L: usize>
    Mul<[[T; L]; M]> for Matrix<T, N, M>
{
    type Output = Matrix<T, N, L>;

    fn mul(self, rhs: [[T; L]; M]) -> Self::Output {
        let mut result = Matrix::new();

        for i in 0..N {
            for j in 0..L {
                let mut x = T::zero();
                for k in 0..M {
                    x = x + self[i][k] * rhs[k][j];
                }
                result[i][j] = x;
            }
        }
        result
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


impl<T: Scalar + Copy, const N: usize, const M: usize>
    Mul<T> for Matrix<T, N, M>
{
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        let mut result = Matrix::new();

        for i in 0..N {
            for j in 0..M {
                result[i][j] = self[i][j] * rhs;
            }
        }
        result
    }
}


impl<T: Scalar + Copy, const N: usize, const M: usize, const L: usize>
    Mul<Matrix<T, M, L>> for Matrix<T, N, M>
{
    type Output = Matrix<T, N, L>;

    fn mul(self, rhs: Matrix<T, M, L>) -> Self::Output {
        self * rhs.data
    }
}


pub fn gcdx<T>(a: T, b: T) -> (T, T, T, T, T) // TODO return a struct?
    where T:
        Copy + Zero + One +
        Div<Output=T> + Sub<Output=T> + Mul<Output=T>
{
    let (mut a, mut a_next) = (a, b);
    let (mut r, mut r_next) = (T::one(), T::zero());
    let (mut s, mut s_next) = (T::zero(), T::one());

    while !a_next.is_zero() {
        let q = a / a_next;
        (a, a_next) = (a_next, a - q * a_next);
        (r, r_next) = (r_next, r - q * r_next);
        (s, s_next) = (s_next, s - q * s_next);
    }

    (a, r, s, r_next, s_next)
}


pub trait Entry: Scalar + Sub<Output=Self> + std::fmt::Display {
    fn clear_column<const N: usize, const M: usize>(
        col: usize, v: &mut [Self; N], b: &mut [Self; N],
        vx: Option<&mut [Self; M]>, bx: Option<&mut [Self; M]>
    );
    fn normalize_column<const N: usize>(col: usize, v: &mut [Self; N]);
    fn reduce_column<const N: usize>(
        col: usize, v: &mut [Self; N], b: &[Self; N]
    );
    fn can_divide(a: Self, b: Self) -> bool;
}


impl Entry for f64 {
    fn clear_column<const N: usize, const M: usize>(
        col: usize, v: &mut [Self; N], b: &mut [Self; N],
        vx: Option<&mut [Self; M]>, bx: Option<&mut [Self; M]>
    ) {
        let f = v[col] / b[col];
        v[col] = 0.0;

        for k in (col + 1)..v.len() {
            v[k] -= b[k] * f;
        }

        if let Some(vx) = vx {
            if let Some(bx) = bx {
                for k in 0..vx.len() {
                    vx[k] -= bx[k] * f;
                }
            }
        }
    }

    fn normalize_column<const N: usize>(col: usize, v: &mut [Self; N]) {
        let f = v[col];
        v[col] = 1.0;

        for k in (col + 1)..v.len() {
            v[k] /= f;
        }
    }

    fn reduce_column<const N: usize>(
        col: usize, v: &mut [Self; N], b: &[Self; N]
    ) {
        let f = v[col];
        v[col] = 0.0;

        for k in (col + 1)..v.len() {
            v[k] -= b[k] * f;
        }
    }

    fn can_divide(a: Self, b: Self) -> bool {
        b != 0.0 && (a == 0.0 || (a / b * b / a - 1.0).abs() <= f64::EPSILON)
    }
}


impl Entry for i64 {
    fn clear_column<const N: usize, const M: usize>(
        col: usize, v: &mut [Self; N], b: &mut [Self; N],
        vx: Option<&mut [Self; M]>, bx: Option<&mut [Self; M]>
    ) {
        let (_, r, s, t, u) = gcdx(b[col], v[col]);
        let det = r * u - s * t;

        for k in col..v.len() {
            let tmp = det * (b[k] * r + v[k] * s);
            v[k] = b[k] * t + v[k] * u;
            b[k] = tmp;
        }

        if let Some(vx) = vx {
            if let Some(bx) = bx {
                for k in 0..vx.len() {
                    let tmp = det * (bx[k] * r + vx[k] * s);
                    vx[k] = bx[k] * t + vx[k] * u;
                    bx[k] = tmp;
                }
            }
        }
    }

    fn normalize_column<const N: usize>(col: usize, v: &mut [Self; N]) {
        if v[col] < 0 {
            for k in col..v.len() {
                v[k] = -v[k];
            }
        }
    }

    fn reduce_column<const N: usize>(
        col: usize, v: &mut [Self; N], b: &[Self; N]
    ) {
        let f = v[col] / b[col] - (
            if v[col] < 0 { 1 } else { 0 }
        );

        if f != 0 {
            for k in col..v.len() {
                v[k] -= b[k] * f;
            }
        }
    }

    fn can_divide(a: Self, b: Self) -> bool {
        b != 0 && a / b * b == a
    }
}


fn pivot_column<T: Zero>(v: &[T]) -> Option<usize> {
    v.iter().position(|x| !x.is_zero())
}


#[derive(Debug, PartialEq)]
struct Basis<T: Entry, const N: usize> {
    vectors: Matrix<T, N, N>,
    rank: usize
}


impl<T: Copy + Entry, const N: usize> Basis<T, N> {
    fn new() -> Self {
        Basis {
            vectors: Matrix::from([[T::zero(); N]; N]),
            rank: 0
        }
    }

    fn extend(&mut self, v: &[T; N]) {
        let mut v = v.clone();

        for i in 0..self.rank {
            let mut b = &mut self.vectors[i];

            if let Some(col) = pivot_column(&v) {
                let col_b = pivot_column(b).unwrap();

                if col < col_b {
                    if (N - i) % 2 > 0 {
                        for j in 0..N {
                            v[j] = -v[j];
                        }
                    }
                    for k in (i..self.rank).rev() {
                        self.vectors[k + 1] = self.vectors[k];
                    }
                    self.vectors[i] = v;
                    self.rank += 1;
                    return;
                } else if col == col_b {
                    Entry::clear_column::<N, N>(col, &mut v, &mut b, None, None);
                }
            } else {
                break;
            }
        }

        if pivot_column(&v).is_some() {
            self.vectors[self.rank] = v;
            self.rank += 1;
        }
    }

    fn reduce(&mut self) {
        let mut col = 0;
        for row in 0..self.rank {
            while self.vectors[row][col].is_zero() {
                col += 1;
            }

            Entry::normalize_column(col, &mut self.vectors[row]);

            let b = self.vectors[row];
            for i in 0..row {
                Entry::reduce_column(col, &mut self.vectors[i], &b);
            }
        }
    }

    fn rank(&self) -> usize {
        self.rank
    }

    fn vectors(&self) -> Vec<[T; N]> {
        (0..self.rank).map(|i| self.vectors[i]).collect()
    }
}


pub struct RowEchelonMatrix<T: Entry, const N: usize, const M: usize> {
    original: Matrix<T, N, M>,
    multiplier: Matrix<T, N, N>,
    result: Matrix<T, N, M>,
    columns: [usize; N],
    rank: usize
}


impl<T: Entry + Copy, const N: usize, const M: usize>
    RowEchelonMatrix<T, N, M>
{
    pub fn from(m: Matrix<T, N, M>) -> Self {
        let mut u = m.clone();
        let mut s = Matrix::identity();
        let mut row = 0;
        let mut cols = [N; N];

        for col in 0..M {
            let pivot_row = (row..N).find(|&r| !u[r][col].is_zero());

            if let Some(pr) = pivot_row {
                if pr != row {
                    (u[pr], u[row]) = (u[row], u[pr]);
                    (s[pr], s[row]) = (s[row], s[pr]);
                }

                let mut vu = u[row].clone();
                let mut vs = s[row].clone();

                for r in (row + 1)..N {
                    Entry::clear_column(
                        col,
                        &mut u[r], &mut vu,
                        Some(&mut s[r]), Some(&mut vs)
                    );
                }

                u[row] = vu;
                s[row] = vs;

                cols[row] = col;
                row += 1;
            }
        }

        RowEchelonMatrix {
            original: m,
            multiplier: s,
            result: u,
            columns: cols,
            rank: row
        }
    }
}


impl<T: Entry + Copy, const N: usize, const M: usize> Matrix<T, N, M> {
    fn rank(&self) -> usize {
        RowEchelonMatrix::from(self.clone()).rank
    }

    fn null_space(&self) -> Vec<Matrix<T, M, 1>> {
        let re = RowEchelonMatrix::from(self.transpose());
        let s = re.multiplier;

        (re.rank..M).map(|i| Matrix::from(s[i]).transpose()).collect()
    }

    fn solve<const K: usize>(&self, rhs: Matrix<T, N, K>)
        -> Option<Matrix<T, M, K>>
        where T: Div<T, Output=T>
    {
        let re = RowEchelonMatrix::from(self.clone());
        let y = re.multiplier * rhs;

        if !(re.rank..N).all(|i| (0..K).all(|j| y[i][j].is_zero())) {
            return None;
        }

        let mut result = Matrix::zero();

        for row in (0..re.rank).rev() {
            let a = (Matrix::from(re.result[row]) * result)[0];
            let b = y[row];
            let x = re.result[row][re.columns[row]];
            for k in 0..K {
                let t = b[k] - a[k];
                if Entry::can_divide(t, x) {
                    result[row][k] = t / x;
                } else {
                    return None;
                }
            }
        }

        Some(result)
    }
}


impl<T: Entry + Copy, const N: usize> Matrix<T, N, N> {
    fn determinant(&self) -> T {
        match self.nr_rows() {
            0 => T::one(),
            1 => self[0][0],
            2 => {
                self[0][0] * self[1][1] - self[0][1] * self[1][0]
            },
            3 => {
                self[0][0] * self[1][1] * self[2][2] +
                self[0][1] * self[1][2] * self[2][0] +
                self[0][2] * self[1][0] * self[2][1] -
                self[0][2] * self[1][1] * self[2][0] -
                self[0][0] * self[1][2] * self[2][1] -
                self[0][1] * self[1][0] * self[2][2]
            },
            _ => {
                let re = RowEchelonMatrix::from(self.clone());
                (0..N).map(|i| re.result[i][i])
                    .reduce(|a, b| a * b)
                    .unwrap_or(T::zero())
            }
        }
    }

    fn inverse(&self) -> Option<Self>
        where T: Div<T, Output=T>
    {
        self.solve(Matrix::identity())
    }
}


#[test]
fn test_matrix_indexing() {
    let mut m = Matrix::from([[1.0, 1.0], [0.0, 1.0]]);

    m[0] = (Matrix::from(m[0]) * 3.0)[0];
    m[1][0] = m[1][0] + 4.0;

    assert_eq!(m, [[3.0, 3.0], [4.0, 1.0]].into());
}


#[test]
fn test_matrix_row_column_manipulation() {
    let mut m = Matrix::from([[1.0, 1.0], [0.0, 1.0]]);

    m.set_row(0, m.get_row(0) * 2.0);
    m.set_column(1, m.get_column(1) * 3.0);

    assert_eq!(m, [[2.0, 6.0], [0.0, 3.0]].into());
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
        [[10]]
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
        Matrix::from([[1, 2, 3]]).vstack(Matrix::from([[4, 5, 6], [7, 8, 9]])),
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]].into()
    );
    assert_eq!(
        Matrix::from([[1, 2], [3, 4]])
            .hstack(Matrix::from([[5, 6], [7, 8]])),
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
    let a = Matrix::from([[1, 2], [3, 4]]);
    let n = a.null_space();
    assert_eq!(n, vec![]);

    let a = Matrix::from([[1.0, 2.0], [3.0, 6.0]]);
    let n = a.null_space();
    assert_eq!(n.len(), 1);
    for v in a.null_space() {
        assert_eq!(a * v, Matrix::from([[0.0], [0.0]]));
    }

    let a = Matrix::from([[1.0, 2.0], [3.0, 4.0]]);
    let n = a.null_space();
    assert_eq!(n, vec![]);

    let a = Matrix::from([[0.0, 1.0, 0.0]]);
    let n = a.null_space();
    assert_eq!(n.len(), 2);
    for v in n {
        assert_eq!(a * v, Matrix::from([[0.0]]));
    }
}

#[test]
fn test_matrix_solve() {
    fn test<T, const N: usize, const M: usize, const L: usize>(
        a: [[T; M]; N], b: [[T; L]; M]
    )
        where T: Entry + Copy + std::fmt::Debug + PartialEq + Div<T, Output=T>
    {
        let a = Matrix::from(a);
        let b = Matrix::from(b);
        let ab = a * b;
        assert_eq!(a * a.solve(ab).unwrap(), ab);

        if N == M && a.rank() == N {
            assert_eq!(a.solve(ab), Some(b));
        }
    }

    test([[1.0, 2.0], [3.0, 4.0]], [[1.0, 0.0], [1.0, -3.0]]);
    test([[1, 2], [3, 4]], [[1, 0], [1, -3]]);
    test([[1.0, 2.0], [3.0, 6.0]], [[1.0], [1.0]]);
    test([[1, 2], [3, 6]], [[1], [1]]);
}

#[test]
fn test_matrix_inverse() {
    fn test<T, const N: usize>(a: [[T; N]; N])
        where T: Entry + Copy + std::fmt::Debug + PartialEq + Div<T, Output=T>
    {
        let a = Matrix::from(a);

        if a.rank() == N {
            assert_eq!(a * a.inverse().unwrap(), Matrix::identity());
        } else {
            assert_eq!(a.inverse(), None);
        }
    }

    test([[1.0, 2.0], [3.0, 4.0]]);
    test([[1.0, 2.0], [3.0, 6.0]]);

    test([[1, 2], [3, 5]]);

    // Inverse exists, but is not integral:
    assert_eq!(Matrix::from([[1, 2], [3, 4]]).inverse(), None);
}
