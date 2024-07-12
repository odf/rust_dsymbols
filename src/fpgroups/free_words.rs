use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::ops::{Index, Mul, MulAssign};


fn normalized<I>(w: I) -> Vec<isize> where I: IntoIterator<Item=isize> {
    let mut buffer = Vec::with_capacity(32);

    for x in w.into_iter() {
        if buffer.last().is_some_and(|y| x == -y) {
            buffer.pop();
        } else if x != 0 {
            buffer.push(x);
        }
    }

    buffer
}


#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct FreeWord {
    w: Vec<isize>
}


impl FreeWord {
    pub fn new<I>(w: I) -> Self where I: IntoIterator<Item=isize> {
        Self { w: normalized(w) }
    }

    pub fn empty() -> Self {
        Self::new([])
    }

    pub fn len(&self) -> usize {
        self.w.len()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, isize> {
        self.w.iter()
    }

    pub fn inverse(&self) -> Self {
        Self::new(self.w.iter().rev().map(|x| -x))
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
        let i = i.rem_euclid(n) as usize;

        let r = std::iter::empty()
            .chain(self.w.iter().skip(i))
            .chain(self.w.iter().take(i))
            .cloned();

        Self::new(r)
    }
}

fn mul(lhs: &[isize], rhs: &[isize]) -> Vec<isize> {
    lhs.iter().chain(rhs.iter()).cloned().collect()
}


impl Index<usize> for FreeWord {
    type Output = isize;

    fn index(&self, index: usize) -> &Self::Output {
        &self.w[index]
    }
}


impl Mul<&FreeWord> for &FreeWord {
    type Output = FreeWord;

    fn mul(self, rhs: &FreeWord) -> Self::Output {
        FreeWord::new(mul(&self.w, &rhs.w))
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


impl Mul<isize> for &FreeWord {
    type Output = FreeWord;

    fn mul(self, rhs: isize) -> Self::Output {
        FreeWord::new(mul(&self.w, &[rhs]))
    }
}


impl Mul<isize> for FreeWord {
    type Output = FreeWord;

    fn mul(self, rhs: isize) -> Self::Output {
        FreeWord::new(mul(&self.w, &[rhs]))
    }
}


impl MulAssign<&FreeWord> for FreeWord {
    fn mul_assign(&mut self, rhs: &FreeWord) {
        self.w = mul(&self.w, &rhs.w);
    }
}


impl PartialOrd for FreeWord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(&other))
    }
}


impl Ord for FreeWord {
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


impl<I> From<I> for FreeWord where I: IntoIterator<Item=isize> {
    fn from(value: I) -> Self {
        Self::new(value)
    }
}


pub fn relator_permutations(fw: &FreeWord) -> BTreeSet<FreeWord> {
    if fw.w.len() == 0 {
        BTreeSet::from([fw.clone()])
    } else {
        let mut result = BTreeSet::from([]);

        for i in 0..fw.w.len() {
            let w = fw.rotated(i as isize);
            result.insert(w.inverse());
            result.insert(w);
        }
        result
    }
}


pub fn relator_representative(fw: &FreeWord) -> FreeWord {
    if fw.w.len() == 0 {
        fw.clone()
    } else {
        let mut best = fw.clone();

        for i in 0..fw.w.len() {
            let w = fw.rotated(i as isize);
            let winv = w.inverse();
            if winv < best {
                best = winv;
            }
            if w < best {
                best = w;
            }
        }
        best
    }
}


#[test]
fn test_freeword_creation() {
    assert_eq!(FreeWord::new([]).w, &[]);
    assert_eq!(FreeWord::new([1, 2]).w, &[1, 2]);
    assert_eq!(FreeWord::from(vec![1, 2, -2, -1]).w, &[]);
    assert_eq!(FreeWord::from([1, 2, -2, 1, 2]).w, &[1, 1, 2]);
}


#[test]
fn test_freeword_len() {
    assert_eq!(FreeWord::from([1, 2, 3, -3, 1]).len(), 3);
}


#[test]
fn test_freeword_index() {
    assert_eq!(FreeWord::from([1, 2, 3])[1], 2);
}


#[test]
fn test_freeword_inverse() {
    assert_eq!(FreeWord::empty().inverse(), FreeWord::empty());
    assert_eq!(
        FreeWord::new([1, 2, 1, -2]).inverse(),
        FreeWord::new([2, -1, -2, -1])
    );
}


#[test]
fn test_freeword_raise_to() {
    assert_eq!(FreeWord::empty().raised_to(0), FreeWord::empty());
    assert_eq!(FreeWord::empty().raised_to(-3), FreeWord::empty());
    assert_eq!(FreeWord::new([1]).raised_to(0), FreeWord::empty());
    assert_eq!(FreeWord::new([1]).raised_to(1), FreeWord::new([1]));

    assert_eq!(
        FreeWord::new([1]).raised_to(-3),
        FreeWord::new([-1, -1, -1])
    );
    assert_eq!(
        FreeWord::new([1, 2, 1, -2, -1]).raised_to(3),
        FreeWord::new([1, 2, 1, 1, 1, -2, -1])
    );
}


#[test]
fn test_freeword_mul() {
    assert_eq!(
        FreeWord::new([1, 2, 3]) * FreeWord::new([-3, -2, 1]),
        FreeWord::new([1, 1])
    );
}


#[test]
fn test_freeword_commutator() {
    assert_eq!(
        FreeWord::new([1, 2]).commutator(&FreeWord::new([3, 2])),
        FreeWord::new([1, 2, 3, -1, -2, -3])
    )
}


#[test]
fn test_freeword_cmp() {
    assert!(&FreeWord::empty() == &FreeWord::new([]));
    assert!(&FreeWord::empty() < &FreeWord::new([1]));
    assert!(&FreeWord::new([2]) > &FreeWord::new([1]));
    assert!(&FreeWord::new([1]) < &FreeWord::new([-1]));
    assert!(&FreeWord::new([-1]) > &FreeWord::new([2]));
    assert!(&FreeWord::new([-1]) < &FreeWord::new([-2]));
    assert!(&FreeWord::new([1, 2, 3, -1]) < &FreeWord::new([1, 2, 3, -2]));
    assert!(&FreeWord::new([1, 2, 3]) < &FreeWord::new([1, 2, 3, -2]));
}


#[test]
fn test_freeword_rotated() {
    assert_eq!(
        FreeWord::new([1, 2, 3]).rotated(0),
        FreeWord::new([1, 2, 3])
    );
    assert_eq!(
        FreeWord::new([1, 2, 3]).rotated(-2),
        FreeWord::new([2, 3, 1])
    );
    assert_eq!(
        FreeWord::new([1, 2, 3]).rotated(5),
        FreeWord::new([3, 1, 2])
    );
    assert_eq!(
        FreeWord::new([1, 2, -1]).rotated(-1),
        FreeWord::new([2])
    );
}


#[test]
fn test_relator_representative() {
    assert_eq!(
        relator_representative(&FreeWord::new([])),
        FreeWord::new([])
    );
    assert_eq!(
        relator_representative(&FreeWord::new([3, 2, 1])),
        FreeWord::new([1, 3, 2])
    );
    assert_eq!(
        relator_representative(&FreeWord::new([3, 2, -1])),
        FreeWord::new([1, -2, -3])
    );
    assert_eq!(
        relator_representative(&FreeWord::new([3, -1, 2, -1])),
        FreeWord::new([1, -2, 1, -3])
    );
}


#[test]
fn test_relator_permutations() {
    let perms = |w: Vec<isize>| {
        relator_permutations(&FreeWord::new(w)).iter()
            .cloned()
            .collect::<Vec<_>>()
    };

    assert_eq!(perms(vec![]), vec![FreeWord::new([])]);
    assert_eq!(perms(vec![1]), vec![FreeWord::new([1]), FreeWord::new([-1])]);

    assert_eq!(
        perms(vec![1, 2]),
        vec![
            FreeWord::new([1, 2]),
            FreeWord::new([2, 1]),
            FreeWord::new([-1, -2]),
            FreeWord::new([-2, -1]),
        ]
    );
    assert_eq!(
        perms(vec![1, 2, 1, 2]),
        vec![
            FreeWord::new([1, 2, 1, 2]),
            FreeWord::new([2, 1, 2, 1]),
            FreeWord::new([-1, -2, -1, -2]),
            FreeWord::new([-2, -1, -2, -1]),
        ]
    );
}
