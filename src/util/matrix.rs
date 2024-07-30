use std::ops::{Add, Index, Mul};

use num_traits::Zero;


trait Scalar: Zero + Copy + Mul<Output=Self> + Add<Output=Self> {
}


#[derive(Copy, Clone, Debug, PartialEq)]
struct Matrix<T, const N: usize, const M: usize> {
    data: [[T; M]; N]
}


impl<T: Scalar, const N: usize, const M: usize> Matrix<T, N, M> {
    fn new(data: [[T; M]; N]) -> Self {
        Self { data }
    }

    fn transpose(self) -> Matrix<T, M, N> {
        let mut result = [[T::zero(); N]; M];
        for i in 0..M {
            for j in 0..N {
                result[i][j] = self.data[j][i];
            }
        }
        Matrix::new(result)
    }
}


impl<T: Scalar, const N: usize, const M: usize>
    Index<(usize, usize)> for Matrix<T, N, M>
{
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.data[index.0][index.1]
    }
}


impl<T: Scalar, const N: usize, const M: usize>
    Index<usize> for Matrix<T, N, M>
{
    type Output = [T; M];

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}


trait AsMatrix<T, const N: usize, const M: usize> {
    fn matrix(self) -> Matrix<T, N, M>;
}


impl<T: Scalar, const N: usize, const M: usize>
    AsMatrix<T, N, M> for [[T; M]; N]
{
    fn matrix(self) -> Matrix<T, N, M> {
        Matrix::new(self)
    }
}


impl<T: Scalar, const M: usize> AsMatrix<T, 1, M> for [T; M] {
    fn matrix(self) -> Matrix<T, 1, M> {
        Matrix::new([self])
    }
}


impl<T: Scalar> AsMatrix<T, 1, 1> for T {
    fn matrix(self) -> Matrix<T, 1, 1> {
        Matrix::new([[self]])
    }
}


impl<T: Scalar> Matrix<T, 1, 1> {
    fn scalar(self) -> T {
        self.data[0][0]
    }
}


impl<T: Scalar, const M: usize> Matrix<T, 1, M> {
    fn array(self) -> [T; M] {
        self.data[0]
    }
}


impl<T: Scalar, const N: usize, const M: usize, const L: usize>
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
        Matrix::new(result)
    }
}


impl<T: Scalar, const N: usize, const M: usize, const L: usize>
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
fn test_matrix_mul() {
    assert_eq!(
        (Matrix::new([[1, 2, 3]]) * [[3], [2], [1]]).scalar(),
        10
    );
    assert_eq!(
        [1, 2, 3].matrix() * [3, 2, 1].matrix().transpose(),
        10.matrix()
    );
    assert_eq!(
        Matrix::new([[1.0, 1.0], [0.0, 1.0]]) * [[1.0, 2.0], [0.0, 1.0]],
        Matrix::new([[1.0, 3.0], [0.0, 1.0]])
    );
    assert_eq!(
        (
            Matrix::new([[1.0, 1.0], [0.0, 1.0]]) *
            [[1.0, 2.0], [0.0, 1.0]]
        )[(0, 1)],
        3.0
    );
    assert_eq!(
        (
            Matrix::new([[1.0, 1.0], [0.0, 1.0]]) *
            [[1.0, 2.0], [0.0, 1.0]]
        )[0],
        [1.0, 3.0]
    );
    assert_eq!(
        (Matrix::new([[1.0, 1.0]]) * [[1.0, 2.0], [0.0, 1.0]]).array(),
        [1.0, 3.0]
    );
}
