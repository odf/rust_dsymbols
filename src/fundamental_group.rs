use std::collections::{HashMap, HashSet, VecDeque};

use crate::dsyms::*;


type Ridge = (usize, usize, usize);

struct Boundary<'a, T: DSym> {
    ds: &'a T,
    opposite: HashMap<Ridge, (Ridge, usize)>,
}


impl<'a, T: DSym> Boundary<'a, T> {
    fn new(ds: &'a T) -> Self {
        let mut opposite = HashMap::new();
        for d in 1..=ds.size() {
            for i in 0..=ds.dim() {
                for j in 0..=ds.dim() {
                    opposite.insert((d, i, j), ((d, j, i), 1));
                }
            }
        }
        Self { ds, opposite }
    }

    fn opposite(&self, d: usize, i: usize, j: usize)
        -> Option<(Ridge, usize)>
    {
        self.opposite.get(&(d, i, j)).cloned()
    }

    fn glue(&mut self, d: usize, i: usize) -> Vec<Ridge> {
        let di = self.ds.op(i, d).unwrap();
        let mut result = Vec::new();

        for j in (0..=self.ds.dim()).filter(|&j| j != i) {
            if let Some((d_opp, d_cnt)) = self.opposite(d, i, j) {
                let (e, k, _) = d_opp;

                if d == di {
                    self.opposite.insert(d_opp, ((0, 0, 0), d_cnt));
                    self.opposite.remove(&(d, i, j));

                    if self.ds.op(k, e) == Some(e) {
                        result.push(d_opp)
                    }
                } else {
                    let (di_opp, di_cnt) = self.opposite(di, i, j).unwrap();
                    let count = d_cnt + di_cnt;

                    self.opposite.insert(d_opp, (di_opp, count));
                    self.opposite.insert(di_opp, (d_opp, count));
                    self.opposite.remove(&(d, i, j));
                    self.opposite.remove(&(di, i, j));

                    if self.ds.op(k, e) != Some(e) {
                        result.push(d_opp)
                    }
                }
            }
        }

        result
    }

    fn glue_recursively(&mut self, todo: Vec<Ridge>) -> Vec<Ridge> {
        let mut todo = VecDeque::from(todo);
        let mut result = Vec::new();

        while let Some(next) = todo.pop_front() {
            let (d, i, j) = next;
            let t = if self.ds.op(i, d) == Some(d) { 1 } else { 2 };
            let m = self.ds.m(i, j, d).unwrap() * t;

            if i == j || self.opposite(d, i, j).is_some_and(|(_, n)| n == m) {
                todo.extend(self.glue(d, i));
                result.push(next);
            }
        }

        result
    }
}


fn spanning_tree<T: DSym>(ds: T) -> Vec<(usize, usize)> {
    let mut seen = HashSet::new();
    let mut result = Vec::new();

    for (maybe_i, d, di) in ds.full_traversal() {
        if !seen.contains(&di) {
            if let Some(i) = maybe_i {
                result.push((d, i));
            }
            seen.insert(di);
        }
    }

    result
}


#[test]
fn test_spanning_tree() {
    let tree = |s: &str| spanning_tree(s.parse::<PartialDSym>().unwrap());

    assert_eq!(
        tree("
            <1.1:24:
            2 4 6 8 10 12 14 16 18 20 22 24,
            16 3 5 7 9 11 13 15 24 19 21 23,
            10 9 20 19 14 13 22 21 24 23 18 17:
            8 4,3 3 3 3
            >
        ").len(),
        23
    );
    assert_eq!(tree("<1.1:3:1 2 3,1 3,2 3:4 8,3>"), vec![(1, 2), (2, 1)]);
    assert_eq!(tree("<1.1:2 3:2,1 2,1 2,2:6,3 2,6>"), vec![(1, 0)]);
    assert_eq!(
        tree("<1.1:8:2 4 6 8,8 3 5 7,1 2 3 4 5 6 7 8:4,4 6 8 4>"),
        vec![(1, 0), (2, 1), (3, 0), (4, 1), (5, 0), (6, 1), (7, 0)]
    );
}
