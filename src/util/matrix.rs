use std::ops::{Add, Index, IndexMut, Mul};

use num_traits::{One, Zero};


pub trait Scalar: Zero + One + Mul<Output=Self> + Add<Output=Self> {
}


#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Matrix<T, const N: usize, const M: usize> {
    data: [[T; M]; N]
}


impl<T: Scalar + Copy , const N: usize, const M: usize> Matrix<T, N, M> {
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
    fn identity() -> Self {
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


impl Scalar for i64 {}
impl Scalar for f64 {}


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
