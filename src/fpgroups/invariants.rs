use std::ops::{AddAssign, SubAssign};

use num_traits::{One, Zero};

use super::free_words::FreeWord;


fn gcdx(a: isize, b: isize) -> (isize, isize, isize, isize, isize)
{
    let (mut a, mut a_next) = (a, b);
    let (mut r, mut r_next) = (1, 0);
    let (mut s, mut s_next) = (0, 1);

    while a_next != 0 {
        let q = a / a_next;
        (a, a_next) = (a_next, a - q * a_next);
        (r, r_next) = (r_next, r - q * r_next);
        (s, s_next) = (s_next, s - q * s_next);
    }

    (a, r, s, r_next, s_next)
}


fn find_pivot(mat: &Vec<Vec<isize>>, start: usize) -> (usize, usize) {
    let (n, m) = (mat.len(), mat[0].len());
    let mut row = start;
    let mut col = start;
    let mut min = isize::max_value();

    for r in start..n {
        for c in start..m {
            let v = mat[r][c].abs();
            if v != 0 && v < min {
                (row, col, min) = (r, c, v)
            }
        }
    }
    (row, col)
}


fn move_pivot_in_place(
    mat: &mut Vec<Vec<isize>>,
    target: usize,
    (row, col): (usize, usize)
)
{
    let (n, m) = (mat.len(), mat[0].len());

    if row != target {
        for c in 0..m {
            (mat[row][c], mat[target][c]) = (mat[target][c], mat[row][c]);
        }
    }
    if col != target {
        for r in 0..n {
            (mat[r][col], mat[r][target]) = (mat[r][target], mat[r][col]);
        }
    }
}


fn clear_later_rows_in_place(mat: &mut Vec<Vec<isize>>, i: usize) -> usize {
    let (n, m) = (mat.len(), mat[0].len());
    let mut count = 0;

    for row in (i + 1)..n {
        let (e, f) = (mat[i][i], mat[row][i]);

        if e != 0 && f % e == 0 {
            let x = f / e;
            for col in i..m {
                mat[row][col] -= x * mat[i][col];
            }
        } else if f != 0 {
            let (_, a, b, c, d) = gcdx(e, f);
            for col in i..m {
                let (v, w) = (mat[i][col], mat[row][col]);
                mat[i][col] = v * a + w * b;
                mat[row][col] = v * c + w * d;
            }
            count += 1;
        }
    }

    count
}


fn clear_later_cols_in_place(mat: &mut Vec<Vec<isize>>, i: usize) -> usize {
    let (n, m) = (mat.len(), mat[0].len());
    let mut count = 0;

    for col in (i + 1)..m {
        let (e, f) = (mat[i][i], mat[i][col]);

        if e != 0 && f % e == 0 {
            let x = f / e;
            for row in i..n {
                mat[row][col] -= x * mat[row][i];
            }
        } else if f != 0 {
            let (_, a, b, c, d) = gcdx(e, f);
            for row in i..n {
                let (v, w) = (mat[row][i], mat[row][col]);
                mat[row][i] = v * a + w * b;
                mat[row][col] = v * c + w * d;
            }
            count += 1;
        }
    }

    count
}


fn diagonalize_in_place(mat: &mut Vec<Vec<isize>>) {
    let (n, m) = (mat.len(), mat[0].len());

    for i in 0.. n.min(m) {
        let (row, col) = find_pivot(mat, i);

        if mat[row][col] != 0 {
            move_pivot_in_place(mat, i, (row, col));
            loop {
                clear_later_rows_in_place(mat, i);
                if clear_later_cols_in_place(mat, i) == 0 {
                    break;
                }
            }
        }

        mat[i][i] = mat[i][i].abs();
    }
}


pub fn relator_as_vector<T>(nr_gens: usize, w: &FreeWord)
    -> Vec<T>
    where T: Zero + One + Clone + AddAssign + SubAssign
{
    let mut row = vec![T::zero(); nr_gens];

    for &g in w.iter() {
        if g < 0 {
            row[(-g - 1) as usize] -= T::one();
        } else {
            row[(g - 1) as usize] += T::one();
        }
    }

    row
}


pub fn abelian_invariants<'a, I>(nr_gens: usize, rels: I)
    -> Vec<usize>
    where I: IntoIterator<Item=&'a FreeWord>
{
    let mut mat: Vec<_> = rels.into_iter()
        .map(|w| relator_as_vector(nr_gens, w))
        .collect();

    if nr_gens == 0 {
        return vec![];
    } else if mat.len() == 0 {
        return vec![0; nr_gens];
    }

    diagonalize_in_place(&mut mat);

    let n = mat.len().min(nr_gens);
    let mut factors: Vec<_> = (0..n).map(|i| mat[i][i]).collect();

    for i in 0..n {
        for j in (i + 1)..n {
            let (a, b) = (factors[i], factors[j]);
            if a != 0 && b % a != 0 {
                let (g, _, _, _, _) = gcdx(a, b);
                factors[i] = g;
                factors[j] = a / g * b;
            }
        }
    }

    let mut result: Vec<_> = factors.iter().cloned()
        .filter(|x| *x != 1)
        .chain(std::iter::repeat(0).take(nr_gens - n))
        .map(|x| x.abs() as usize)
        .collect();
    result.sort();

    result
}


#[test]
fn test_diagonalize_in_place() {
    fn check<const N: usize, const M: usize>(
        mat: [[isize; M]; N],
        out: [[isize; M]; N],
    ) {
        let mut mat: Vec<Vec<_>> =
            mat.iter().map(|row| row.iter().cloned().collect()).collect();
        let out: Vec<Vec<_>> =
            out.iter().map(|row| row.iter().cloned().collect()).collect();

        diagonalize_in_place(&mut mat);
        assert_eq!(mat, out);
    }

    check(
        [[2, 2, 3], [3, 3, 4], [2, 1, 3]],
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    );
    check(
        [[1, 2], [3, 4]],
        [[1, 0], [0, 2]],
    );
    check(
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
        [[1, 0, 0], [0, 3, 0], [0, 0, 0]],
    );
    check(
        [[2, 0, 0], [0, 2, 0], [0, 0, 2], [2, 2, 0], [2, 0, 2], [0, 2, 2]],
        [[2, 0, 0], [0, 2, 0], [0, 0, 2], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
    )
}


#[test]
fn test_abelian_invariants() {
    fn check(nr_gens: usize, rels: &[&[isize]], invariants: &[usize])
    {
        let rels: Vec<_> = rels.iter().map(|&r| Vec::from(r).into()).collect();
        assert_eq!(abelian_invariants(nr_gens, &rels), invariants);
    }

    check(3,
        &[&[1, 2, -1, -2], &[1, 3, -1, -3], &[2, 3, -2, -3]],
        &[0, 0, 0]
    );
    check(3,
        &[&[1, 1], &[2, 2], &[3, 3], &[1, 2, 1, 2], &[1, 3, 1, 3], &[2, 3, 2, 3]],
        &[2, 2, 2]
    );
}
