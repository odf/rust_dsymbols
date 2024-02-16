use std::collections::{HashMap, VecDeque};

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
