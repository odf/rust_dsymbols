use num_traits::Zero;

use crate::geometry::matrix::Matrix;

use super::traits::{Scalar, ScalarPtr, gcdx};


pub trait Entry: Scalar {
    fn clear_column(col: usize, v: &mut [Self], b: &mut [Self]);
    fn normalize_column(col: usize, v: &mut [Self]);
    fn reduce_column(col: usize, v: &mut [Self], b: &[Self]);
}


impl Entry for f64 {
    fn clear_column(col: usize, v: &mut [Self], b: &mut [Self]) {
        let f = v[col] / b[col];
        v[col] = 0.0;

        for k in (col + 1)..v.len() {
            v[k] -= b[k] * f;
        }
    }

    fn normalize_column(col: usize, v: &mut [Self]) {
        let f = v[col];
        v[col] = 1.0;

        for k in (col + 1)..v.len() {
            v[k] /= f;
        }
    }

    fn reduce_column(col: usize, v: &mut [Self], b: &[Self]) {
        let f = v[col];
        v[col] = 0.0;

        for k in (col + 1)..v.len() {
            v[k] -= b[k] * f;
        }
    }
}


impl Entry for i64 {
    fn clear_column(col: usize, v: &mut [Self], b: &mut [Self]) {
        assert_eq!(v.len(), b.len());

        let (_, r, s, t, u) = gcdx(b[col], v[col]);
        let det = r * u - s * t;

        for k in col..v.len() {
            let tmp = det * (b[k] * r + v[k] * s);
            v[k] = b[k] * t + v[k] * u;
            b[k] = tmp;
        }
    }

    fn normalize_column(col: usize, v: &mut [Self]) {
        if v[col] < 0 {
            for k in col..v.len() {
                v[k] = -v[k];
            }
        }
    }

    fn reduce_column(col: usize, v: &mut [Self], b: &[Self]) {
        assert_eq!(v.len(), b.len());

        let f = v[col] / b[col] - (
            if v[col] < 0 { 1 } else { 0 }
        );

        if f != 0 {
            for k in col..v.len() {
                v[k] -= b[k] * f;
            }
        }
    }
}

fn pivot_column<T: Zero>(v: &[T]) -> Option<usize> {
    v.iter().position(|x| !x.is_zero())
}


#[derive(Debug, PartialEq)]
pub struct Basis<T: Entry, const N: usize> {
    vectors: Matrix<T, N, N>,
    rank: usize
}


impl<T: Copy + Entry, const N: usize> Basis<T, N>
    where for <'a> &'a T: ScalarPtr<T>
{
    pub fn new() -> Self {
        Basis {
            vectors: Matrix::new(),
            rank: 0
        }
    }

    pub fn extend(&mut self, v: &[T; N]) {
        let mut v = v.clone();

        for i in 0..self.rank {
            let b = &mut self.vectors[i];

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
                    Entry::clear_column(col, &mut v, b);
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

    pub fn reduce(&mut self) {
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

    pub fn rank(&self) -> usize {
        self.rank
    }

    pub fn vectors(&self) -> Vec<[T; N]> {
        (0..self.rank).map(|i| self.vectors[i]).collect()
    }
}
