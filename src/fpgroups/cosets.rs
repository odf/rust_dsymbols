use std::collections::{BTreeSet, HashMap, VecDeque};

use crate::util::partitions::Partition;

use super::free_words::{FreeWord, Relator};


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

        while c > self.len() {
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

    fn merge(&mut self, a: usize, b: usize) {
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


fn scan(table: &CosetTable, w: &FreeWord, start: usize, limit: usize)
    -> (usize, usize)
{
    let mut row = start;

    for index in 0..limit {
        let next = table.get(row, w[index]);
        if next == 0 {
            return (row, index);
        } else {
            row = next
        }
    }
    (row, limit)
}


fn scan_both_ways(table: &CosetTable, w: &FreeWord, start: usize)
    -> (usize, usize, usize, isize)
{
    let n = w.len();
    let (head, i) = scan(table, w, start, n);
    let (tail, j) = scan(table, &w.inverse(), start, n - i);
    (head, tail, n - i - j, w[i])
}


fn scan_and_connect(table: &mut CosetTable, w: &FreeWord, start: usize) {
    let (head, tail, gap, c) = scan_both_ways(table, w, start);

    if gap == 1 {
        table.join(head, tail, c);
    } else if gap == 0 && head != tail {
        table.merge(head, tail);
    }
}


pub fn coset_table(
    nr_gens: usize, relators: &Vec<Relator>, subgroup_gens: &Vec<FreeWord>
) -> Vec<HashMap<isize, isize>>
{
    let mut rels = BTreeSet::new();
    for rel in relators {
        rels.extend(rel.permutations());
    }
    let rels = rels;

    let mut table = CosetTable::new(nr_gens);
    let mut i = 0;

    while i < table.len() {
        i += 1;
        if i != table.canon(i) {
            continue;
        }

        for g in table.all_gens() {
            if table.get(i, g) == 0 {
                let n = table.len();
                assert!(n < 100_000, "Reached coset table limit of 100_000");

                table.join(i, n, g);
                for w in &rels {
                    scan_and_connect(&mut table, w, n);
                }

                for w in subgroup_gens {
                    let c = table.canon(0);
                    scan_and_connect(&mut table, w, c);
                }
            }
        }
    }

    table.compact()
}
