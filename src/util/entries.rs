use std::ops::{Add, Div, Mul, Neg, Sub};
use num_traits::{One, Zero};


pub trait Scalar:
    Zero + One + Mul<Output=Self> + Add<Output=Self> + Neg<Output=Self>
{
}

impl Scalar for f64 {}
impl Scalar for i64 {}


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


pub trait Entry: Scalar + Sub<Output=Self> {
    fn clear_column(
        col: usize, v: &mut [Self], b: &mut [Self],
        vx: Option<&mut [Self]>, bx: Option<&mut [Self]>
    );
    fn normalize_column(col: usize, v: &mut [Self]);
    fn reduce_column(col: usize, v: &mut [Self], b: &[Self]);
    fn can_divide(a: Self, b: Self) -> bool;
}


impl Entry for f64 {
    fn clear_column(
        col: usize, v: &mut [Self], b: &mut [Self],
        vx: Option<&mut [Self]>, bx: Option<&mut [Self]>
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

    fn can_divide(a: Self, b: Self) -> bool {
        b != 0.0 && (a == 0.0 || (a / b * b / a - 1.0).abs() <= f64::EPSILON)
    }
}


impl Entry for i64 {
    fn clear_column(
        col: usize, v: &mut [Self], b: &mut [Self],
        vx: Option<&mut [Self]>, bx: Option<&mut [Self]>
    ) {
        assert_eq!(v.len(), b.len());

        let (_, r, s, t, u) = gcdx(b[col], v[col]);
        let det = r * u - s * t;

        for k in col..v.len() {
            let tmp = det * (b[k] * r + v[k] * s);
            v[k] = b[k] * t + v[k] * u;
            b[k] = tmp;
        }

        if let Some(vx) = vx {
            if let Some(bx) = bx {
                assert_eq!(vx.len(), bx.len());

                for k in 0..vx.len() {
                    let tmp = det * (bx[k] * r + vx[k] * s);
                    vx[k] = bx[k] * t + vx[k] * u;
                    bx[k] = tmp;
                }
            }
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

    fn can_divide(a: Self, b: Self) -> bool {
        b != 0 && a / b * b == a
    }
}
