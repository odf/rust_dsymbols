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


fn find_pivot<const N: usize, const M: usize>(
    mat: &[[isize; M]; N],
    start: usize
)
    -> (usize, usize)
{
    let mut row = start;
    let mut col = start;
    let mut min = mat[row][col].abs();

    for r in start..N {
        for c in start..M {
            let v = mat[r][c].abs();
            if v != 0 && v < min {
                (row, col, min) = (r, c, v)
            }
        }
    }
    (row, col)
}


fn move_pivot_in_place<const N: usize, const M: usize>(
    mat: &mut [[isize; M]; N],
    target: usize,
    (row, col): (usize, usize)
)
{
    if row != target {
        (mat[row], mat[target]) = (mat[target], mat[row]);
    }
    if col != target {
        for r in 0..N {
            (mat[r][col], mat[r][target]) = (mat[r][target], mat[r][col]);
        }
    }
}


fn clear_later_rows_in_place<const N: usize, const M: usize>(
    mat: &mut [[isize; M]; N],
    i: usize
)
    -> usize
{
    let mut count = 0;

    for row in (i + 1)..N {
        let (e, f) = (mat[i][i], mat[row][i]);

        if e != 0 && f % e == 0 {
            let x = f / e;
            for col in i..M {
                mat[row][col] -= x * mat[i][col];
            }
        } else if f != 0 {
            let (_, a, b, c, d) = gcdx(e, f);
            for col in i..M {
                let (v, w) = (mat[i][col], mat[row][col]);
                mat[i][col] = v * a + w * b;
                mat[row][col] = v * c + w * d;
            }
            count += 1;
        }
    }

    count
}


fn clear_later_cols_in_place<const N: usize, const M: usize>(
    mat: &mut [[isize; M]; N],
    i: usize
)
    -> usize
{
    let mut count = 0;

    for col in (i + 1)..M {
        let (e, f) = (mat[i][i], mat[i][col]);

        if e != 0 && f % e == 0 {
            let x = f / e;
            for row in i..N {
                mat[row][col] -= x * mat[row][i];
            }
        } else if f != 0 {
            let (_, a, b, c, d) = gcdx(e, f);
            for row in i..N {
                let (v, w) = (mat[row][i], mat[row][col]);
                mat[row][i] = v * a + w * b;
                mat[row][col] = v * c + w * d;
            }
            count += 1;
        }
    }

    count
}


fn diagonalize_in_place<const N: usize, const M: usize>(
    mat: &mut [[isize; M]; N]
) {
    for i in 0.. N.min(M) {
        let (row, col) = find_pivot(mat, i);
        let val = mat[row][col].abs();

        if val != 0 {
            loop {
                move_pivot_in_place(mat, i, (row, col));
                clear_later_rows_in_place(mat, i);
                if clear_later_cols_in_place(mat, i) == 0 {
                    break;
                }
            }
        }

        mat[i][i] = mat[i][i].abs();
    }
}


#[test]
fn test_diagonalize_in_place() {
    fn check<const N: usize, const M: usize>(
        mat: [[isize; M]; N],
        out: [[isize; M]; N],
    ) {
        let mut tmp = mat.clone();
        diagonalize_in_place(&mut tmp);
        assert_eq!(tmp, out);
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
