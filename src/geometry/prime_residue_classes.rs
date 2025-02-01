use std::ops::Add;

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


#[test]
fn test_valid() {
    assert!(PrimeResidueClass::<2>::valid());
    assert!(PrimeResidueClass::<61>::valid());
    assert!(PrimeResidueClass::<9999991>::valid());
}


#[test]
fn test_invalid() {
    assert!(!PrimeResidueClass::<1>::valid());
    assert!(!PrimeResidueClass::<77>::valid());
    assert!(!PrimeResidueClass::<9999992>::valid());
}
