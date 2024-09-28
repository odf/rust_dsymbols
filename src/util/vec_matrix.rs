use std::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};
use num_traits::{One, Zero};


pub trait Scalar:
    Zero + One + Mul<Output=Self> + Add<Output=Self> + Neg<Output=Self>
{
}

impl Scalar for f64 {}
impl Scalar for i64 {}


#[derive(Clone, Debug, PartialEq)]
pub struct VecMatrix<T> {
    data: Vec<T>,
    nr_rows: usize,
    nr_cols: usize,
}


impl<T: Scalar> VecMatrix<T> {
    pub fn nr_rows(&self) -> usize {
        self.nr_rows
    }

    pub fn nr_columns(&self) -> usize {
        self.nr_cols
    }
}


impl<T: Scalar> Index<usize> for VecMatrix<T>
{
    type Output = [T];

    fn index(&self, index: usize) -> &Self::Output {
        assert!(index < self.nr_rows);
        &self.data[index * self.nr_cols .. (index + 1) * self.nr_cols]
    }
}


impl<T: Scalar> IndexMut<usize> for VecMatrix<T>
{
    fn index_mut(&mut self, index: usize) -> &mut [T] {
        assert!(index < self.nr_rows);
        &mut self.data[index * self.nr_cols .. (index + 1) * self.nr_cols]
    }
}


impl<T: Scalar + Clone, const N: usize, const M: usize>
    From<[[T; M]; N]> for VecMatrix<T>
{
    fn from(data: [[T; M]; N]) -> Self {
        let mut result = VecMatrix::new(N, M);

        for i in 0..N {
            for j in 0..M {
                result[i][j] = data[i][j].clone();
            }
        }

        result
    }
}


impl<T: Scalar + Clone, const M: usize> From<[T; M]> for VecMatrix<T> {
    fn from(data: [T; M]) -> Self {
        let mut result = VecMatrix::new(1, M);

        for j in 0..M {
            result[0][j] = data[j].clone();
        }

        result
    }
}


impl<T: Scalar + Clone> From<T> for VecMatrix<T> {
    fn from(data: T) -> Self {
        let mut result = VecMatrix::new(1, 1);
        result[0][0] = data.clone();
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

        let mut result = VecMatrix::new(1, self.nr_cols);

        for j in 0..self.nr_cols {
            result[0][j] = self[i][j].clone();
        }

        result
    }

    pub fn set_row(&mut self, i: usize, row: VecMatrix<T>) {
        assert_eq!(row.nr_rows, 1);
        assert_eq!(row.nr_cols, self.nr_cols);
        assert!(i < self.nr_rows);

        for j in 0..self.nr_cols {
            self[i][j] = row[0][j].clone();
        }
    }

    pub fn get_col(&self, j: usize) -> VecMatrix<T> {
        assert!(j < self.nr_cols);

        let mut result = VecMatrix::new(self.nr_rows, 1);

        for i in 0..self.nr_rows {
            result[i][0] = self[i][j].clone();
        }

        result
    }

    pub fn set_col(&mut self, j: usize, col: VecMatrix<T>) {
        assert_eq!(col.nr_cols, 1);
        assert_eq!(col.nr_rows, self.nr_rows);
        assert!(j < self.nr_cols);

        for i in 0..self.nr_rows {
            self[i][j] = col[i][0].clone();
        }
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


impl<T: Scalar + Copy> Mul<VecMatrix<T>> for VecMatrix<T> {
    type Output = VecMatrix<T>;

    fn mul(self, rhs: VecMatrix<T>) -> Self::Output {
        let mut result = VecMatrix::new(self.nr_rows, rhs.nr_cols);

        for i in 0..self.nr_rows {
            for j in 0..rhs.nr_cols {
                let mut x = T::zero();
                for k in 0..self.nr_cols {
                    x = x + self[i][k] * rhs[k][j];
                }
                result[i][j] = x;
            }
        }

        result
    }
}


impl<T: Scalar + Copy, const N: usize, const M: usize>
    Mul<VecMatrix<T>> for [[T; M]; N]
{
    type Output = VecMatrix<T>;

    fn mul(self, rhs: VecMatrix<T>) -> Self::Output {
        VecMatrix::from(self) * rhs
    }
}


impl<T: Scalar + Copy, const M: usize, const L: usize>
    Mul<[[T; L]; M]> for VecMatrix<T>
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


impl<T: Scalar + Copy> Mul<T> for VecMatrix<T> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        let mut result = VecMatrix::new(self.nr_rows, self.nr_cols);

        for i in 0..self.nr_rows {
            for j in 0..self.nr_cols {
                result[i][j] = self[i][j] * rhs;
            }
        }

        result
    }
}


#[test]
fn test_matrix_indexing() {
    let mut m = VecMatrix::from([[1.0, 1.0], [0.0, 1.0]]);

    m.set_row(0, m.get_row(0) * 3.0);
    m[1][0] = m[1][0] + 4.0;

    assert_eq!(m, [[3.0, 3.0], [4.0, 1.0]].into());
}


#[test]
fn test_matrix_submatrix() {
    assert_eq!(
        VecMatrix::from([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).submatrix(0..2, [0, 2]),
        VecMatrix::from([[1, 3], [4, 6]])
    )
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
