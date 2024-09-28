use crate::util::matrix::Matrix;
use crate::util::entries::Entry;

use num_traits::Zero;


pub(crate) fn pivot_column<T: Zero>(v: &[T]) -> Option<usize> {
    v.iter().position(|x| !x.is_zero())
}


#[derive(Debug, PartialEq)]
pub(crate) struct Basis<T: Entry, const N: usize> {
    pub(crate) vectors: Matrix<T, N, N>,
    pub(crate) rank: usize
}


impl<T: Copy + Entry, const N: usize> Basis<T, N> {
    pub(crate) fn new() -> Self {
        Basis {
            vectors: Matrix::from([[T::zero(); N]; N]),
            rank: 0
        }
    }

    pub(crate) fn extend(&mut self, v: &[T; N]) {
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

    pub(crate) fn reduce(&mut self) {
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

    pub(crate) fn rank(&self) -> usize {
        self.rank
    }

    pub(crate) fn vectors(&self) -> Vec<[T; N]> {
        (0..self.rank).map(|i| self.vectors[i]).collect()
    }
}
