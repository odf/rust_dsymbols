use core::fmt;
use std::collections::{BTreeSet, HashMap, VecDeque};

use crate::util::partitions::Partition;

use super::free_words::{FreeWord, Relator};


type CosetTable = Vec<HashMap<isize, usize>>;


#[derive(Clone)]
struct DynamicCosetTable {
    nr_gens: usize,
    nr_rows: usize,
    table: HashMap<(usize, isize), usize>,
    part: Partition<usize>
}


impl DynamicCosetTable {
    fn new(nr_gens: usize) -> Self {
        Self {
            nr_gens,
            nr_rows: 1,
            table: HashMap::new(),
            part: Partition::new(),
        }
    }

    fn all_gens(&self) -> Vec<isize> {
        let tmp: Vec<_> = (1..=self.nr_gens).map(|g| g as isize).collect();
        tmp.iter().cloned().chain(tmp.iter().map(|g| -g)).collect()
    }

    fn nr_rows(&self) -> usize {
        self.nr_rows
    }

    fn get(&self, c: usize, g: isize) -> Option<&usize> {
        self.table.get(&(c, g))
    }

    fn set(&mut self, c: usize, g: isize, d: usize) {
        let n = self.nr_gens as isize;
        assert!(c > 0);
        assert!(d > 0);
        assert!(g >= -n && g <= n);

        self.nr_rows = self.nr_rows.max(c).max(d);
        self.table.insert((c, g), d);
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

        while let Some((a, b)) = queue.pop_front() {
            let a = self.canon(a);
            let b = self.canon(b);

            if a != b {
                self.part.unite(&a, &b);

                for g in self.all_gens() {
                    if let Some(&ag) = self.get(a, g) {
                        self.set(b, g, ag);
                        if let Some(&bg) = self.get(b, g) {
                            queue.push_back((ag, bg));
                        }
                    } else if let Some(&bg) = self.get(b, g) {
                        self.set(a, g, bg);
                    }
                }
            }
        }
    }

    fn compact(&self) -> CosetTable {
        let mut to_idx = vec![0; self.nr_rows() + 1];
        let mut i = 0;
        for k in 1..=self.nr_rows() {
            if self.canon(k) == k {
                to_idx[k] = i;
                i += 1;
            }
        }
        let to_idx = to_idx;

        let mut result = vec![];
        for k in 1..=self.nr_rows() {
            if self.canon(k) == k {
                let mut row = HashMap::new();
                for g in self.all_gens() {
                    row.insert(
                        g,
                        to_idx[self.canon(*self.get(k, g).unwrap())]
                    );
                }
                result.push(row);
            }
        }

        result
    }
}


impl fmt::Display for DynamicCosetTable {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for c in 1..=self.nr_rows() {
            for g in self.all_gens() {
                write!(f, "{} -> {}, ", g, self.get(c, g).unwrap())?;
            }
            write!(f, "({})\n", self.part.find(&c))?;
        }
        Ok(())
    }
}


fn scan(table: &DynamicCosetTable, w: &FreeWord, start: usize, limit: usize)
    -> (usize, usize)
{
    let mut row = start;

    for index in 0..limit {
        if let Some(&next) = table.get(row, w[index]) {
            row = next;
        } else {
            return (row, index);
        }
    }
    (row, limit)
}


fn scan_both_ways(table: &DynamicCosetTable, w: &FreeWord, start: usize)
    -> (usize, usize, usize, isize)
{
    let n = w.len();
    let (head, i) = scan(table, w, start, n);
    let (tail, j) = scan(table, &w.inverse(), start, n - i);
    (head, tail, n - i - j, if i < n { w[i] } else { 0 })
}


fn scan_and_connect(table: &mut DynamicCosetTable, w: &FreeWord, start: usize) {
    let (head, tail, gap, c) = scan_both_ways(table, w, start);

    if gap == 1 {
        table.join(head, tail, c);
    } else if gap == 0 && head != tail {
        table.merge(head, tail);
    }
}


pub fn coset_table(
    nr_gens: usize, relators: &Vec<Relator>, subgroup_gens: &Vec<FreeWord>
) -> CosetTable
{
    let mut rels = BTreeSet::new();
    for rel in relators {
        rels.extend(rel.permutations());
    }
    let rels = rels;

    let mut table = DynamicCosetTable::new(nr_gens);
    let mut i = 0;

    while i < table.nr_rows() {
        i += 1;

        for g in table.all_gens() {
            if i != table.canon(i) {
                break;
            }
            if table.get(i, g).is_none() {
                let n = table.nr_rows() + 1;
                assert!(n < 100_000, "Reached coset table limit of 100_000");

                table.join(i, n, g);
                for w in &rels {
                    if w[0] == g {
                        scan_and_connect(&mut table, w, i);
                    }
                }

                for w in subgroup_gens {
                    let c = table.canon(1);
                    scan_and_connect(&mut table, w, c);
                }
            }
        }
    }

    table.compact()
}


pub fn coset_representative(table: &CosetTable)
    -> Vec<FreeWord>
{
    let mut queue = VecDeque::from([0 as usize]);
    let mut result = vec![None; table.len()];
    result[0] = Some(FreeWord::empty());

    while let Some(i) = queue.pop_front() {
        let w = result[i].clone().unwrap();

        for (&g, &k) in table[i].iter() {
            if result[k].is_none() {
                result[k] = Some(&w * FreeWord::from([g]));
                queue.push_back(k);
            }
        }
    }

    result.iter().flatten().cloned().collect()
}


fn induced_table<T, F>(nr_gens: usize, img: F, start: &T) -> CosetTable
    where
        T: Clone + Eq + std::hash::Hash,
        F: Fn(&T, isize) -> T
{
    let mut table = DynamicCosetTable::new(nr_gens);
    let mut o2n = HashMap::from([(start.clone(), 1)]);
    let mut n2o = HashMap::from([(1, start.clone())]);
    let mut i = 0;

    while i < table.nr_rows() {
        i += 1;
        for g in table.all_gens() {
            let k = img(&n2o[&i], g);
            let n = *o2n.entry(k.clone()).or_insert(table.nr_rows() + 1);
            n2o.insert(n, k);
            table.join(i, n, g);
        }
    }

    table.compact()
}


#[cfg(test)]
mod coset_tests {
    use super::*;

    fn make_table(nr_gens: usize, rels: &[&[isize]], subgens: &[&[isize]])
        -> CosetTable
    {
        coset_table(
            nr_gens,
            &rels.iter().map(|r| Relator::from(*r)).collect(),
            &subgens.iter().map(|g| FreeWord::from(*g)).collect(),
        )
    }

    #[test]
    fn test_coset_table() {
        assert_eq!(
            make_table(2, &[&[1, 1], &[2, 2], &[1, 2, 1, 2]], &[]),
            vec![
                HashMap::from([(1, 1), (2, 2), (-1, 1), (-2, 2)]),
                HashMap::from([(1, 0), (2, 3), (-1, 0), (-2, 3)]),
                HashMap::from([(1, 3), (2, 0), (-1, 3), (-2, 0)]),
                HashMap::from([(1, 2), (2, 1), (-1, 2), (-2, 1)]),
            ]
        );
        assert_eq!(
            make_table(2, &[&[1, 1], &[2, 2], &[1, 2, 1, 2]], &[&[1]]),
            vec![
                HashMap::from([(1, 0), (2, 1), (-1, 0), (-2, 1)]),
                HashMap::from([(1, 1), (2, 0), (-1, 1), (-2, 0)]),
            ]
        );
        assert_eq!(
            make_table(2, &[&[1, 1], &[2, 2], &[1, 2, 1, 2]], &[&[2]]),
            vec![
                HashMap::from([(1, 1), (2, 0), (-1, 1), (-2, 0)]),
                HashMap::from([(1, 0), (2, 1), (-1, 0), (-2, 1)]),
            ]
        );
        assert_eq!(
            make_table(2, &[&[1, 1], &[2, 2], &[1, 2, 1, 2]], &[&[1], &[2]]),
            vec![
                HashMap::from([(1, 0), (2, 0), (-1, 0), (-2, 0)]),
            ]
        );
        assert_eq!(
            make_table(
                3,
                &[
                    &[1, 1], &[2, 2], &[3, 3],
                    &[1, 2, 1, 2, 1, 2], &[1, 3, 1, 3], &[3, 2, 3, 2, 3, 2]
                ],
                &[&[1, 2]]
            ).len(),
            8
        );
    }

    #[test]
    fn test_coset_representatives() {
        assert_eq!(
            coset_representative(
                &make_table(2, &[&[1, 1], &[2, 2], &[1, 2, 1, 2]], &[])
            ).len(),
            4
        );
        assert_eq!(
            coset_representative(
                &make_table(2, &[&[1, 1], &[2, 2], &[1, 2, 1, 2]], &[&[1]])
            ).len(),
            2
        );
        assert_eq!(
            coset_representative(&make_table(
                2, &[&[1, 1], &[2, 2], &[1, 2, 1, 2]], &[&[1], &[2]]
            )).len(),
            1
        );
        assert_eq!(
            coset_representative(&make_table(
                3,
                &[
                    &[1, 1], &[2, 2], &[3, 3],
                    &[1, 2, 1, 2, 1, 2], &[1, 3, 1, 3], &[3, 2, 3, 2, 3, 2]
                ],
                &[&[1, 2]]
            )).len(),
            8
        );
    }
}
