use std::iter as iter;


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
            let w: Vec<_> = iter::repeat(&self.w[..])
                .take(m as usize)
                .collect();
            Self::new(&normalized(&w))
        }
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
