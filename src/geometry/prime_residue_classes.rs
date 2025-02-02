use std::ops::{Add, Div, Mul, Neg, Sub};

use num_traits::{One, Zero};


#[derive(Copy, Clone, Debug, PartialEq, Eq)]
struct PrimeResidueClass<const P: i64> {
    value: i64
}


impl<const P: i64> PrimeResidueClass<P> {
    fn valid() -> bool {
        if P < 2 || P as f64 > (i64::MAX as f64).sqrt() {
            false
        } else {
            for n in 2.. {
                if n * n > P {
                    break;
                } else if P % n == 0 {
                    return false;
                }
            }
            true
        }
    }

    fn inverse(self) -> Self {
        let (mut t, mut t1) = (0, 1);
        let (mut r, mut r1) = (P, self.value);

        while r1 != 0 {
            let q = r / r1;
            (t, t1) = (t1, t - q * t1);
            (r, r1) = (r1, r - q * r1);
        }

        assert_eq!(r, 1);

        t.into()
    }
}


impl<const P: i64> From<i64> for PrimeResidueClass<P> {
    fn from(n: i64) -> Self {
        PrimeResidueClass {
            value: if n >= 0 { n % P } else { n % P + P }
        }
    }
}


impl<const P: i64> Add<PrimeResidueClass<P>> for PrimeResidueClass<P> {
    type Output = Self;

    fn add(self, rhs: PrimeResidueClass<P>) -> Self::Output {
        (self.value + rhs.value).into()
    }
}


impl<const P: i64> Sub<PrimeResidueClass<P>> for PrimeResidueClass<P> {
    type Output = Self;

    fn sub(self, rhs: PrimeResidueClass<P>) -> Self::Output {
        (self.value - rhs.value).into()
    }
}


impl<const P: i64> Mul<PrimeResidueClass<P>> for PrimeResidueClass<P> {
    type Output = Self;

    fn mul(self, rhs: PrimeResidueClass<P>) -> Self::Output {
        (self.value * rhs.value).into()
    }
}


impl<const P: i64> Div<PrimeResidueClass<P>> for PrimeResidueClass<P> {
    type Output = Self;

    fn div(self, rhs: PrimeResidueClass<P>) -> Self::Output {
        self * rhs.inverse()
    }
}


impl<const P: i64> Neg for PrimeResidueClass<P> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        (-self.value).into()
    }
}


impl<const P: i64> Zero for PrimeResidueClass<P> {
    fn zero() -> Self {
        0.into()
    }

    fn is_zero(&self) -> bool {
        self.value == 0
    }
}


impl<const P: i64> One for PrimeResidueClass<P> {
    fn one() -> Self {
        1.into()
    }

    fn is_one(&self) -> bool {
        self.value == 1
    }
}


mod test {
    use proptest::prelude::*;
    use proptest::collection::vec;

    use super::*;

    #[test]
    fn test_valid() {
        assert!(PrimeResidueClass::<2>::valid());
        assert!(PrimeResidueClass::<61>::valid());
        assert!(PrimeResidueClass::<9999991>::valid());
        assert!(PrimeResidueClass::<3037000493>::valid()); // largest prime we can use
    }


    #[test]
    fn test_invalid() {
        assert!(!PrimeResidueClass::<1>::valid());
        assert!(!PrimeResidueClass::<77>::valid());
        assert!(!PrimeResidueClass::<9999992>::valid());
        assert!(!PrimeResidueClass::<3037000507>::valid()); // prime, but too large
    }


    #[test]
    fn test_add_sub() {
        let a = PrimeResidueClass::<61>::from(54);
        assert_eq!(a + 47.into() - 47.into(), a);
    }


    #[test]
    fn test_div_mul() {
        let a = PrimeResidueClass::<61>::from(54);
        assert_eq!(a * 47.into() / 47.into(), a);
    }

    proptest! {
        #[test]
        fn test_binomial(n in 0..61) {
            let n = PrimeResidueClass::<61>::from(n as i64);
            assert_eq!((n - 1.into()) * (n + 1.into()), n * n - 1.into());
        }

        #[test]
        fn test_binomial_large(n in 0..3037000493u32) {
            let n = PrimeResidueClass::<3037000493>::from(n as i64);
            assert_eq!((n - 1.into()) * (n + 1.into()), n * n - 1.into());
        }
    }
}
