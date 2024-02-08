use std::cmp::Ordering;
use std::ops::Mul;


fn normalized(ws: &[&[isize]]) -> Vec<isize> {
    let mut buffer = vec![];

    for &w in ws {
        for &x in w {
            if buffer.last().is_some_and(|y| x == -y) {
                buffer.pop();
            } else if x != 0 {
                buffer.push(x);
            }
        }
    }

    buffer
}


#[derive(Debug, Eq, PartialEq)]
pub struct FreeWord {
    w: Vec<isize>
}


impl FreeWord {
    pub fn new(w: &[isize]) -> Self {
        Self { w: normalized(&[w]) }
    }

    pub fn empty() -> Self {
        Self::new(&[])
    }

    pub fn inverse(&self) -> Self {
        Self::new(&self.w.iter().rev().map(|x| -x).collect::<Vec<_>>())
    }

    pub fn raised_to(&self, m: isize) -> Self {
        if m < 0 {
            self.inverse().raised_to(-m)
        } else {
            (0..m).fold(FreeWord::empty(), |a, _| a * self)
        }
    }

    pub fn commutator(&self, other: &FreeWord) -> Self {
        self * other * self.inverse() * other.inverse()
    }

    pub fn rotated(&self, i: isize) -> Self {
        let n = self.w.len() as isize;
        let i = (n - i.rem_euclid(n)) as usize;
        let w_iter = self.w[i..].iter().chain(self.w[..i].iter()).cloned();
        Self::new(&w_iter.collect::<Vec<_>>())
    }
}


impl Mul<&FreeWord> for &FreeWord {
    type Output = FreeWord;

    fn mul(self, rhs: &FreeWord) -> Self::Output {
        FreeWord::new(&normalized(&[&self.w, &rhs.w]))
    }
}


impl Mul<FreeWord> for &FreeWord {
    type Output = FreeWord;

    fn mul(self, rhs: FreeWord) -> Self::Output {
        self * &rhs
    }
}


impl Mul<&FreeWord> for FreeWord {
    type Output = FreeWord;

    fn mul(self, rhs: &FreeWord) -> Self::Output {
        &self * rhs
    }
}


impl Mul<FreeWord> for FreeWord {
    type Output = FreeWord;

    fn mul(self, rhs: FreeWord) -> Self::Output {
        &self * &rhs
    }
}


impl PartialOrd for &FreeWord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(&other))
    }
}


impl Ord for &FreeWord {
    fn cmp(&self, other: &Self) -> Ordering {
        for i in 0..(self.w.len().min(other.w.len())) {
            let x = self.w[i];
            let y = other.w[i];

            if x != y {
                if x > 0 && y > 0 {
                    return x.cmp(&y);
                } else {
                    return y.cmp(&x);
                }
            }
        }

        self.w.len().cmp(&other.w.len())
    }
}


impl From<&[isize]> for FreeWord {
    fn from(value: &[isize]) -> Self {
        Self { w: normalized(&[value]) }
    }
}


impl From<Vec<isize>> for FreeWord {
    fn from(value: Vec<isize>) -> Self {
        (&value[..]).into()
    }
}


#[test]
fn test_freeword_normalized() {
    assert_eq!(normalized(&[&[]]), &[]);
    assert_eq!(normalized(&[&[1, 2]]), &[1, 2]);
    assert_eq!(normalized(&[&[1, 2], &[-2, -1]]), &[]);
    assert_eq!(normalized(&[&[1, 2, -2], &[1, 2]]), &[1, 1, 2]);
}


#[test]
fn test_freeword_inverse() {
    assert_eq!(FreeWord::empty().inverse(), FreeWord::empty());
    assert_eq!(
        FreeWord::new(&[1, 2, 1, -2]).inverse(),
        FreeWord::new(&[2, -1, -2, -1])
    );
}


#[test]
fn test_freeword_raise_to() {
    assert_eq!(FreeWord::empty().raised_to(0), FreeWord::empty());
    assert_eq!(FreeWord::empty().raised_to(-3), FreeWord::empty());
    assert_eq!(FreeWord::new(&[1]).raised_to(0), FreeWord::empty());
    assert_eq!(FreeWord::new(&[1]).raised_to(1), FreeWord::new(&[1]));

    assert_eq!(
        FreeWord::new(&[1]).raised_to(-3),
        FreeWord::new(&[-1, -1, -1])
    );
    assert_eq!(
        FreeWord::new(&[1, 2, 1, -2, -1]).raised_to(3),
        FreeWord::new(&[1, 2, 1, 1, 1, -2, -1])
    );
}


#[test]
fn test_freeword_mul() {
    assert_eq!(
        FreeWord::new(&[1, 2, 3]) * FreeWord::new(&[-3, -2, 1]),
        FreeWord::new(&[1, 1])
    );
}


#[test]
fn test_freeword_commutator() {
    assert_eq!(
        FreeWord::new(&[1, 2]).commutator(&FreeWord::new(&[3, 2])),
        FreeWord::new(&[1, 2, 3, -1, -2, -3])
    )
}


#[test]
fn test_freeword_cmp() {
    assert!(&FreeWord::empty() == &FreeWord::new(&[]));
    assert!(&FreeWord::empty() < &FreeWord::new(&[1]));
    assert!(&FreeWord::new(&[2]) > &FreeWord::new(&[1]));
    assert!(&FreeWord::new(&[1]) < &FreeWord::new(&[-1]));
    assert!(&FreeWord::new(&[-1]) > &FreeWord::new(&[2]));
    assert!(&FreeWord::new(&[-1]) < &FreeWord::new(&[-2]));
    assert!(&FreeWord::new(&[1, 2, 3, -1]) < &FreeWord::new(&[1, 2, 3, -2]));
    assert!(&FreeWord::new(&[1, 2, 3]) < &FreeWord::new(&[1, 2, 3, -2]));
}


#[test]
fn test_freeword_rotated() {
    assert_eq!(
        FreeWord::new(&[1, 2, 3]).rotated(0),
        FreeWord::new(&[1, 2, 3])
    );
    assert_eq!(
        FreeWord::new(&[1, 2, 3]).rotated(2),
        FreeWord::new(&[2, 3, 1])
    );
    assert_eq!(
        FreeWord::new(&[1, 2, 3]).rotated(-5),
        FreeWord::new(&[3, 1, 2])
    );
    assert_eq!(
        FreeWord::new(&[1, 2, -1]).rotated(1),
        FreeWord::new(&[2])
    );
}
