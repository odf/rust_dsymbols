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


pub fn empty() -> Vec<isize> {
    vec![]
}


pub fn word(w: &[isize]) -> Vec<isize> {
    normalized(&[w])
}


pub fn inverse(w: &[isize]) -> Vec<isize> {
    word(w).iter().rev().map(|x| -x).collect()
}


pub fn raised_to(m: isize, w: &[isize]) -> Vec<isize> {
    if m < 0 {
        raised_to(-m, &inverse(w))
    } else {
        normalized(&iter::repeat(w).take(m as usize).collect::<Vec<_>>())
    }
}


#[test]
fn test_fp_words_normalized() {
    assert_eq!(normalized(&[&[]]), empty());
    assert_eq!(normalized(&[&[1, 2]]), word(&[1, 2]));
    assert_eq!(normalized(&[&[1, 2], &[-2, -1]]), empty());
    assert_eq!(normalized(&[&[1, 2, -2], &[1, 2]]), word(&[1, 1, 2]));
}


#[test]
fn test_fp_words_inverse() {
    assert_eq!(inverse(&empty()), empty());
    assert_eq!(inverse(&[1, 2, 1, -2]), &[2, -1, -2, -1]);
}


#[test]
fn test_fp_words_raise_to() {
    assert_eq!(raised_to(0, &empty()), empty());
    assert_eq!(raised_to(-3, &empty()), empty());
    assert_eq!(raised_to(0, &[1]), empty());
    assert_eq!(raised_to(1, &[1]), &[1]);
    assert_eq!(raised_to(-3, &[1]), &[-1, -1, -1]);
    assert_eq!(raised_to(3, &[1, 2, 1, -2, -1]), &[1, 2, 1, 1, 1, -2, -1]);
}
