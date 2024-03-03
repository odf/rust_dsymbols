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
