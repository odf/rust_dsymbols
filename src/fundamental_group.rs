use std::collections::HashMap;

use crate::dsyms::*;


struct Boundary<'a, T: DSym> {
    ds: &'a T,
    opposite: HashMap<(usize, usize, usize), (usize, usize)>,
}


impl<'a, T: DSym> Boundary<'a, T> {
    fn new(ds: &'a T) -> Self {
        let mut opposite = HashMap::new();
        for d in 1..=ds.size() {
            for i in 0..=ds.dim() {
                for j in 0..=ds.dim() {
                    opposite.insert((d, i, j), (d, 1));
                }
            }
        }
        Self { ds, opposite }
    }

    fn opposite(&self, d: usize, i: usize, j: usize)
        -> Option<(usize, usize)>
    {
        self.opposite.get(&(d, i, j)).cloned()
    }

    fn glue(&mut self, d: usize, i: usize) {
        let e = self.ds.op(i, d).unwrap();

        for j in (0..=self.ds.dim()).filter(|&j| j != i) {
            if let Some((d_opp, d_cnt)) = self.opposite(d, i, j) {
                if d == e {
                    if d_cnt % 2 == 0 {
                        self.opposite.insert((d_opp, i, j), (0, d_cnt));
                    } else {
                        self.opposite.insert((d_opp, j, i), (0, d_cnt));
                    }
                    self.opposite.remove(&(d, i, j));
                } else {
                    let (e_opp, e_cnt) = self.opposite(e, i, j).unwrap();
                    let count = d_cnt + e_cnt;

                    if d_cnt % 2 == 0 {
                        self.opposite.insert((d_opp, i, j), (e_opp, count));
                    } else {
                        self.opposite.insert((d_opp, j, i), (e_opp, count));
                    }
                    if e_cnt % 2 == 0 {
                        self.opposite.insert((e_opp, i, j), (d_opp, count));
                    } else {
                        self.opposite.insert((e_opp, j, i), (d_opp, count));
                    }
                    self.opposite.remove(&(d, i, j));
                    self.opposite.remove(&(e, i, j));
                }
            }
        }
    }
}
