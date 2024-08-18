use std::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};

use num_traits::{One, Zero};


pub trait Scalar:
    Zero + One + Mul<Output=Self> + Add<Output=Self> + Neg<Output=Self>
{
}


#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Matrix<T, const N: usize, const M: usize> {
    data: [[T; M]; N]
}


impl<T: Scalar + Copy , const N: usize, const M: usize> Matrix<T, N, M> {
    pub fn nr_rows(&self) -> usize {
        N
    }

    pub fn nr_columns(&self) -> usize {
        M
    }

    pub fn transpose(&self) -> Matrix<T, M, N> {
        let mut result = [[T::zero(); N]; M];
        for i in 0..M {
            for j in 0..N {
                result[i][j] = self.data[j][i];
            }
        }
        Matrix::from(result)
    }

    pub fn get_row(&self, i: usize) -> Matrix<T, 1, M> {
        assert!(i < N);
        Matrix::from([self.data[i]])
    }

    pub fn set_row(&mut self, i: usize, row: Matrix<T, 1, M>) {
        assert!(i < N);
        self.data[i] = row.data[0];
    }

    pub fn get_column(&self, j: usize) -> Matrix<T, N, 1> {
        assert!(j < M);
        let mut result = [[T::zero(); 1]; N];
        for i in 0..N {
            result[i][0] = self.data[i][j];
        }
        Matrix::from(result)
    }

    pub fn set_column(&mut self, j: usize, column: Matrix<T, N, 1>) {
        assert!(j < M);
        for i in 0..N {
            self.data[i][j] = column.data[i][0];
        }
    }

    pub fn hstack<const L: usize, const S: usize>(self, rhs: Matrix<T, N, L>)
        -> Matrix<T, N, S>
    {
        assert_eq!(S, N + M);

        let mut result = [[T::zero(); S]; N];

        for i in 0..N {
            for j in 0..M {
                result[i][j] = self.data[i][j];
            }
            for j in 0..L {
                result[i][M + j] = rhs[i][j];
            }
        }
        Matrix::from(result)
    }

    pub fn vstack<const L: usize, const S: usize>(self, rhs: Matrix<T, L, M>)
        -> Matrix<T, S, M>
    {
        assert_eq!(S, N + L);

        let mut result = [[T::zero(); M]; S];

        for j in 0..M {
            for i in 0..N {
                result[i][j] = self.data[i][j];
            }
            for i in 0..L {
                result[N + i][j] = rhs[i][j];
            }
        }
        Matrix::from(result)
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

        let mut result = [[T::zero(); L]; K];

        for i in 0..K {
            for j in 0..L {
                result[i][j] = self.data[rows[i]][columns[j]];
            }
        }

        Matrix::from(result)
    }
}


impl<T: Scalar + Copy, const N: usize> Matrix<T, N, N> {
    pub fn identity() -> Self {
        let mut data = [[T::zero(); N]; N];
        for i in 0..N {
            data[i][i] = T::one();
        }
        Self { data }
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


impl<T: Scalar, const N: usize, const M: usize>
    Index<(usize, usize)> for Matrix<T, N, M>
{
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        assert!(index.0 < N);
        assert!(index.0 < M);
        &self.data[index.0][index.1]
    }
}


impl<T: Scalar, const N: usize, const M: usize>
    IndexMut<(usize, usize)> for Matrix<T, N, M>
{
    fn index_mut(&mut self, index: (usize, usize)) -> &mut T {
        assert!(index.0 < N);
        assert!(index.0 < M);
        &mut self.data[index.0][index.1]
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


impl<T: Scalar + Copy, const N: usize, const M: usize, const L: usize>
    Mul<[[T; L]; M]> for Matrix<T, N, M>
{
    type Output = Matrix<T, N, L>;

    fn mul(self, rhs: [[T; L]; M]) -> Self::Output {
        let mut result = [[T::zero(); L]; N];

        for i in 0..N {
            for j in 0..L {
                let mut x = T::zero();
                for k in 0..M {
                    x = x + self.data[i][k] * rhs[k][j];
                }
                result[i][j] = x;
            }
        }
        Matrix::from(result)
    }
}


impl<T: Scalar + Copy, const N: usize, const M: usize>
    Mul<T> for Matrix<T, N, M>
{
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        let mut result = [[T::zero(); M]; N];

        for i in 0..N {
            for j in 0..M {
                result[i][j] = self.data[i][j] * rhs;
            }
        }
        Matrix::from(result)
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


pub trait Entry: Scalar {
    fn clear_column<const N: usize>(
        col: usize, v: &mut [Self; N], b: &mut [Self; N]
    );
    fn normalize_column<const N: usize>(
        col: usize, v: &mut [Self; N]
    );
    fn reduce_column<const N: usize>(
        col: usize, v: &mut [Self; N], b: &[Self; N]
    );
    fn solve_row<const N: usize>(
        a: &[Self; N], x: &Vec<&[Self; N]>, b: &[Self; N]
    )
        -> Option<[Self; N]>;
}


impl Scalar for f64 {}

impl Entry for f64 {
    fn clear_column<const N: usize>(
        col: usize, v: &mut [Self; N], b: &mut [Self; N]
    ) {
        let f = v[col] / b[col];
        v[col] = 0.0;

        for k in (col + 1)..v.len() {
            v[k] -= b[k] * f;
        }
    }

    fn normalize_column<const N: usize>(
        col: usize, v: &mut [Self; N]
    ) {
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

    fn solve_row<const N: usize>(
        _a: &[Self; N], _x: &Vec<&[Self; N]>, b: &[Self; N]
    )
        -> Option<[Self; N]>
    {
        Some(b.clone())
    }
}


impl Scalar for i64 {}

impl Entry for i64 {
    fn clear_column<const N: usize>(
        col: usize, v: &mut [Self; N], b: &mut [Self; N]
    ) {
        let (_, r, s, t, u) = gcdx(b[col], v[col]);
        let det = r * u - s * t;

        for k in col..v.len() {
            let tmp = det * (b[k] * r + v[k] * s);
            v[k] = b[k] * t + v[k] * u;
            b[k] = tmp;
        }
    }

    fn normalize_column<const N: usize>(
        col: usize, v: &mut [Self; N]
    ) {
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

    fn solve_row<const N: usize>(
        a: &[Self; N], x: &Vec<&[Self; N]>, b: &[Self; N]
    )
        -> Option<[Self; N]>
    {
        let k = x.len();
        let mut result = [0; N];

        for col in 0..b.len() {
            let mut t = b[col];
            for i in 0..k {
                t -= a[i] * x[i][col];
            }
            if t % a[k] == 0 {
                result[col] = t / a[k];
            } else {
                return None;
            }
        }
        Some(result)
    }
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
        let pivot_column = |v: &[T]| { v.iter().position(|x| !x.is_zero()) };

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
                    self.vectors[i] = v;
                    return;
                } else if col == col_b {
                    Entry::clear_column(col, &mut v, &mut b);
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
            while self.vectors[(row, col)].is_zero() {
                col += 1;
            }

            Entry::normalize_column(col, &mut self.vectors[row]);

            let b = self.vectors[row];
            for i in 0..row {
                Entry::reduce_column(col, &mut self.vectors[i], &b);
            }
        }
    }
}


impl<T: Entry + Copy , const N: usize, const M: usize> Matrix<T, N, M> {
    fn reduced_basis(&self) -> Basis<T, M> {
        let mut b = Basis::new();
        for i in 0..self.nr_rows() {
            b.extend(&self[i]);
        }
        b.reduce();
        b
    }
}


#[test]
fn test_matrix_indexing() {
    let mut m = Matrix::from([[1.0, 1.0], [0.0, 1.0]]);

    m[0] = (Matrix::from(m[0]) * 3.0)[0];
    m[(1, 0)] = m[(1, 0)] + 4.0;

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
fn test_matrix_mul() {
    assert_eq!(
        (Matrix::from([[1, 2, 3]]) * [[3], [2], [1]])[(0, 0)],
        10
    );
    assert_eq!(
        Matrix::from([1, 2, 3]) * Matrix::from([3, 2, 1]).transpose(),
        10.into()
    );
    assert_eq!(
        Matrix::from([[1.0, 1.0], [0.0, 1.0]]) * [[1.0, 2.0], [0.0, 1.0]],
        [[1.0, 3.0], [0.0, 1.0]].into()
    );
    assert_eq!(
        (
            Matrix::from([[1.0, 1.0], [0.0, 1.0]]) *
            [[1.0, 2.0], [0.0, 1.0]]
        )[(0, 1)],
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
fn test_matrix_reduced_basis() {
    assert_eq!(
        Matrix::from([[1, 4, 7], [2, 5, 8], [3, 6, 8]]).reduced_basis(),
        Basis {
            vectors: Matrix::from([[1, 1, 0], [0, 3, 0], [0, 0, 1]]),
            rank: 3
        }
    );
    assert_eq!(
        Matrix::from([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 8.0]])
            .reduced_basis(),
        Basis { vectors: Matrix::identity(), rank: 3 }
    );
    assert_eq!(
        Matrix::from([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
            .reduced_basis(),
        Basis {
            vectors: Matrix::from(
                [[1.0, 0.0, -1.0], [0.0, 1.0, 2.0], [0.0, 0.0, 0.0]]
            ),
            rank: 2
        }
    );
}
