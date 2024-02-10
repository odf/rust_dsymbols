use std::collections::{HashMap, VecDeque};

use crate::util::partitions::Partition;


#[derive(Clone)]
struct CosetTable {
    nr_gens: usize,
    table: Vec<Vec<usize>>,
    part: Partition<usize>
}


impl CosetTable {
    fn new(nr_gens: usize) -> Self {
        Self {
            nr_gens,
            table: vec![vec![], vec![0; 2 * nr_gens + 1]],
            part: Partition::new(),
        }
    }

    fn all_gens(&self) -> Vec<isize> {
        let tmp: Vec<_> = (1..=self.nr_gens).map(|g| g as isize).collect();
        tmp.iter().cloned().chain(tmp.iter().map(|g| -g)).collect()
    }

    fn len(&self) -> usize {
        // don't count the initial dummy row
        self.table.len() - 1
    }

    fn get(&self, c: usize, g: isize) -> usize {
        let n = self.nr_gens as isize;

        if c == 0 || c > self.len() || g < -n || g > n {
            // using 0 instead of None since this is not user-facing code
            0
        } else {
            // generator inverses are negative numbers, so offset second index
            self.table[c][(g + n) as usize]
        }
    }

    fn set(&mut self, c: usize, g: isize, d: usize) {
        let n = self.nr_gens as isize;
        assert!(c > 0);
        assert!(d > 0);
        assert!(g >= -n && g <= n);

        while c < self.len() {
            self.table.push(vec![0; 2 * self.nr_gens + 1]);
        }
        self.table[c][(g + n) as usize] = d;
    }

    fn join(&mut self, c: usize, d: usize, g: isize) {
        self.set(c, g, d);
        self.set(d, -g, c);
    }

    fn canon(&self, c: usize) -> usize {
        self.part.find(&c)
    }

    fn identify(&mut self, a: usize, b: usize) {
        let mut queue: VecDeque<(usize, usize)> = VecDeque::from([(a, b)]);

        while let Some((a, b)) = queue.pop_back() {
            let a = self.canon(a);
            let b = self.canon(b);

            if a != b {
                self.part.unite(&a, &b);

                for g in self.all_gens() {
                    let ag = self.get(a, g);
                    let bg = self.get(b, g);
                    if ag == 0 {
                        self.set(a, g, bg);
                    } else {
                        if bg != 0 && self.canon(ag) != self.canon(bg) {
                            queue.push_front((ag, bg));
                        }
                        self.set(b, g, ag);
                    }
                }
            }
        }
    }

    fn compact(&self) -> Vec<HashMap<isize, isize>> {
        let mut to_idx = vec![0; self.len() + 1];
        let mut i = 0;
        for k in 1..=self.len() {
            if self.canon(k) == k {
                to_idx[k] = i;
                i += 1;
            }
        }
        let to_idx = to_idx;

        let mut result = vec![];
        for k in 1..=self.len() {
            if to_idx[k] != 0 {
                let mut row = HashMap::new();
                for g in self.all_gens() {
                    row.insert(g, to_idx[self.canon(self.get(k, g))]);
                }
                result.push(row);
            }
        }

        result
    }
}
