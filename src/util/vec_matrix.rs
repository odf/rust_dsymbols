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


impl<T: Scalar + Clone> VecMatrix<T> {
    pub fn nr_rows(&self) -> usize {
        self.nr_rows
    }

    pub fn nr_columns(&self) -> usize {
        self.nr_cols
    }
}


impl<T: Scalar> Index<(usize, usize)> for VecMatrix<T>
{
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        assert!(index.0 < self.nr_rows);
        assert!(index.1 < self.nr_cols);
        &self.data[index.0 * self.nr_cols + index.1]
    }
}


impl<T: Scalar> IndexMut<(usize, usize)> for VecMatrix<T>
{
    fn index_mut(&mut self, index: (usize, usize)) -> &mut T {
        assert!(index.0 < self.nr_rows);
        assert!(index.1 < self.nr_cols);
        &mut self.data[index.0 * self.nr_cols + index.1]
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
