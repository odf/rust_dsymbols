use core::f64;
use std::ops::{Add, Div, Index, IndexMut, Mul, Neg, Sub};
use num_rational::BigRational;
use num_traits::{One, Zero, Signed};


pub trait Scalar:
    Zero + One +
    Neg<Output=Self> +
    Add<Output=Self> +
    Sub<Output=Self> +
    Mul<Output=Self> +
    Div<Output=Self>
{
}

impl Scalar for BigRational {}
impl Scalar for f64 {}
impl Scalar for i64 {}


pub trait ScalarPtr<T>:
    Sized +
    Neg<Output=T> +
    Add<Output=T> + Add<T, Output=T> +
    Sub<Output=T> + Sub<T, Output=T> +
    Mul<Output=T> + Mul<T, Output=T> +
    Div<Output=T> + Div<T, Output=T>
{
}

impl ScalarPtr<BigRational> for &BigRational {}
impl ScalarPtr<f64> for &f64 {}
impl ScalarPtr<i64> for &i64 {}

pub trait Array2d<T>:
    Index<(usize, usize), Output=T> +
    IndexMut<(usize, usize), Output=T>
{
    fn nr_rows(&self) -> usize;
    fn nr_columns(&self) -> usize;
}


pub fn allclose<M: Array2d<f64>>(a: &M, b: &M, rtol: f64, atol: f64) -> bool {
    assert_eq!(a.nr_rows(), b.nr_rows());
    assert_eq!(a.nr_columns(), b.nr_columns());

    for i in 0..a.nr_rows() {
        for j in 0..a.nr_columns() {
            if (a[(i, j)] - b[(i, j)]).abs() > (atol + rtol * b[(i, j)].abs()) {
                return false;
            }
        }
    }

    true
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
    fn can_divide(a: &Self, b: &Self) -> bool;
    fn pivot_row<M: Array2d<Self>>(col: usize, row0: usize, a: &M)
        -> Option<usize>;
    fn clear_col<A: Array2d<Self>, B: Array2d<Self>>(
        col: usize, row1: usize, row2: usize, a: &mut A, x: Option<&mut B>
    );
}


impl Entry for BigRational {
    fn can_divide(a: &Self, b: &Self) -> bool {
        !b.is_zero()
    }

    fn pivot_row<M: Array2d<Self>>(col: usize, row0: usize, a: &M)
        -> Option<usize>
    {
        let mut best_row = row0;

        for row in (row0 + 1)..a.nr_rows() {
            if a[(row, col)].abs() > a[(best_row, col)].abs() {
                best_row = row;
            }
        }

        if a[(best_row, col)].is_zero() { None } else { Some(best_row) }
    }

    fn clear_col<A: Array2d<Self>, B: Array2d<Self>>(
        col: usize, row1: usize, row2: usize, a: &mut A, x: Option<&mut B>
    ) {
        let f = &a[(row1, col)] / &a[(row2, col)];
        a[(row1, col)] = Self::zero();

        for k in (col + 1)..a.nr_columns() {
            a[(row1, k)] = &a[(row1, k)] - &a[(row2, k)] * &f;
        }

        if let Some(x) = x {
            for k in 0..x.nr_columns() {
                x[(row1, k)] = &x[(row1, k)] - &x[(row2, k)] * &f;
            }
        }
    }
}


impl Entry for f64 {
    fn can_divide(a: &Self, b: &Self) -> bool {
        *b != 0.0 && (*a == 0.0 || (a / b * b / a - 1.0).abs() <= f64::EPSILON)
    }

    fn pivot_row<M: Array2d<Self>>(col: usize, row0: usize, a: &M)
        -> Option<usize>
    {
        let mut best_row = row0;

        for row in (row0 + 1)..a.nr_rows() {
            if a[(row, col)].abs() > a[(best_row, col)].abs() {
                best_row = row;
            }
        }

        if a[(best_row, col)] != 0.0 { Some(best_row) } else { None }
    }

    fn clear_col<A: Array2d<Self>, B: Array2d<Self>>(
        col: usize, row1: usize, row2: usize, a: &mut A, x: Option<&mut B>
    ) {
        let eps = 1e-12;

        let f = a[(row1, col)] / a[(row2, col)];
        a[(row1, col)] = 0.0;

        for k in (col + 1)..a.nr_columns() {
            let s = a[(row1, k)];
            let t = s - a[(row2, k)] * f;
            a[(row1, k)] = if t.abs() < eps * s.abs() { 0.0 } else { t };
        }

        if let Some(x) = x {
            for k in 0..x.nr_columns() {
                let s = x[(row1, k)];
                let t = s - x[(row2, k)] * f;
                x[(row1, k)] = if t.abs() < eps * s.abs() { 0.0 } else { t };
            }
        }
    }
}


impl Entry for i64 {
    fn can_divide(a: &Self, b: &Self) -> bool {
        *b != 0 && a / b * b == *a
    }

    fn pivot_row<M: Array2d<Self>>(col: usize, row0: usize, a: &M)
        -> Option<usize>
    {
        let mut best_row = row0;

        for row in (row0 + 1)..a.nr_rows() {
            let x = a[(row, col)];
            let y = a[(best_row, col)];
            if x != 0 && (y == 0 || x.abs() < y.abs()) {
                best_row = row;
            }
        }

        if a[(best_row, col)] != 0 { Some(best_row) } else { None }
    }

    fn clear_col<A: Array2d<Self>, B: Array2d<Self>>(
        col: usize, row1: usize, row2: usize, a: &mut A, x: Option<&mut B>
    ) {
        let (_, r, s, t, u) = gcdx(a[(row2, col)], a[(row1, col)]);
        let det = r * u - s * t;

        for k in col..a.nr_columns() {
            let tmp = det * (a[(row2, k)] * r + a[(row1, k)] * s);
            a[(row1, k)] = a[(row2, k)] * t + a[(row1, k)] * u;
            a[(row2, k)] = tmp;
        }

        if let Some(x) = x {
            for k in 0..x.nr_columns() {
                let tmp = det * (x[(row2, k)] * r + x[(row1, k)] * s);
                x[(row1, k)] = x[(row2, k)] * t + x[(row1, k)] * u;
                x[(row2, k)] = tmp;
            }
        }
    }
}
