use core::fmt;
use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};

use crate::util::backtrack::{BackTrackIterator, BackTracking};
use crate::util::partitions::IntPartition;

use super::free_words::{relator_permutations, FreeWord};


#[derive(Clone)]
pub struct CosetTable {
    nr_gens: usize,
    table: Vec<Vec<isize>>,
    part: IntPartition
}


impl CosetTable {
    pub fn new(nr_gens: usize) -> Self {
        Self {
            nr_gens,
            table: vec![vec![-1; nr_gens * 2 + 1]],
            part: IntPartition::new(),
        }
    }

    pub fn nr_gens(&self) -> usize {
        self.nr_gens
    }

    pub fn all_gens(&self) -> Vec<isize> {
        let tmp: Vec<_> = (1..=self.nr_gens).map(|g| g as isize).collect();
        tmp.iter().cloned().chain(tmp.iter().map(|g| -g)).collect()
    }

    pub fn len(&self) -> usize {
        self.table.len()
    }

    pub fn get(&self, c: usize, g: isize) -> Option<usize> {
        if c < self.len() {
            let r = self.table[c][(g + self.nr_gens as isize) as usize];
            if r >= 0 {
                Some(self.canon(r as usize))
            } else {
                None
            }
        } else {
            None
        }
    }

    pub fn set(&mut self, c: usize, g: isize, d: usize) {
        while c >= self.len() {
            self.table.push(vec![-1; self.nr_gens * 2 + 1]);
        }
        self.table[c][(g + self.nr_gens as isize) as usize] = d as isize;
    }

    fn join(&mut self, c: usize, d: usize, g: isize) {
        self.set(c, g, d);
        self.set(d, -g, c);
    }

    fn canon(&self, c: usize) -> usize {
        self.part.find(c)
    }

    fn merge(&mut self, a: usize, b: usize) {
        let mut queue: VecDeque<(usize, usize)> = VecDeque::from([(a, b)]);

        while let Some((a, b)) = queue.pop_front() {
            let a = self.canon(a);
            let b = self.canon(b);

            if a != b {
                for g in self.all_gens() {
                    if let Some(ag) = self.get(a, g) {
                        if let Some(bg) = self.get(b, g) {
                            queue.push_back((ag, bg));
                        } else {
                            self.set(b, g, ag);
                        }
                    } else if let Some(bg) = self.get(b, g) {
                        self.set(a, g, bg);
                    }
                }
                self.part.unite(a, b);
            }
        }
    }

    fn compact(&self) -> CosetTable {
        let mut n = 0;
        let mut old_to_new = vec![0; self.len()];
        for k in 0..self.len() {
            if self.canon(k) == k {
                old_to_new[k] = n;
                n += 1;
            }
        }

        let mut result = CosetTable::new(self.nr_gens);
        for k in 0..self.len() {
            if self.canon(k) == k {
                for g in self.all_gens() {
                    if let Some(c) = self.get(k, g) {
                        result.set(old_to_new[k], g, old_to_new[c]);
                    }
                }
            }
        }

        result
    }

    fn as_map(&self) -> Vec<BTreeMap<isize, usize>> {
        let mut result = vec![];

        for k in 0..self.len() {
            let mut row = BTreeMap::new();
            for g in self.all_gens() {
                if let Some(c) = self.get(k, g) {
                    row.insert(g, c);
                }
            }
            result.push(row);
        }

        result
    }
}


impl fmt::Display for CosetTable {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for c in 0..self.len() {
            write!(f, "{}: ", c)?;
            for g in self.all_gens() {
                if let Some(d) = self.get(c, g) {
                    write!(f, "{} -> {}, ", g, d)?;
                }
            }
            writeln!(f, "({})", self.canon(c))?;
        }
        Ok(())
    }
}


fn expanded_relator_set(relators: &Vec<FreeWord>) -> BTreeSet<FreeWord> {
    let mut rels = BTreeSet::new();
    for rel in relators {
        rels.extend(relator_permutations(&rel));
    }
    rels
}


fn scan(table: &CosetTable, w: &FreeWord, start: usize, limit: usize)
    -> (usize, usize)
{
    let mut row = start;

    for index in 0..limit {
        if let Some(next) = table.get(row, w[index]) {
            row = next;
        } else {
            return (row, index);
        }
    }
    (row, limit)
}


fn scan_inverse(
    table: &CosetTable, w: &FreeWord, start: usize, limit: usize
)
    -> (usize, usize)
{
    let n = w.len();
    let mut row = start;

    for index in 0..limit {
        if let Some(next) = table.get(row, -w[n - 1 - index]) {
            row = next;
        } else {
            return (row, index);
        }
    }
    (row, limit)
}


fn scan_both_ways(table: &CosetTable, w: &FreeWord, start: usize)
    -> (usize, usize, usize, isize)
{
    let n = w.len();
    let (head, i) = scan(table, w, start, n);
    let (tail, j) = scan_inverse(table, w, start, n - i);
    (head, tail, n - i - j, if i < n { w[i] } else { w[0] })
}


fn scan_and_connect(
    table: &mut CosetTable, w: &FreeWord, start: usize
) {
    let (head, tail, gap, c) = scan_both_ways(table, w, start);

    if gap == 1 {
        table.join(head, tail, c);
    } else if gap == 0 && head != tail {
        table.merge(head, tail);
    }
}


pub fn coset_table(
    nr_gens: usize, relators: &Vec<FreeWord>, subgroup_gens: &Vec<FreeWord>
) -> CosetTable
{
    let rels = expanded_relator_set(relators);
    let mut table = CosetTable::new(nr_gens);

    for i in 0.. {
        if i >= table.len() {
            break;
        }

        for g in table.all_gens() {
            if i != table.canon(i) {
                break;
            }
            if table.get(i, g).is_none() {
                let n = table.len();
                assert!(n < 100_000, "Reached coset table limit of 100_000");

                table.join(i, n, g);
                for w in &rels {
                    if w[0] == g {
                        let c = table.canon(i);
                        scan_and_connect(&mut table, w, c);
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


pub fn coset_representative(table: &CosetTable) -> BTreeMap<usize, FreeWord> {
    let mut queue = VecDeque::from([0]);
    let mut result = BTreeMap::from([(0, FreeWord::empty())]);

    while let Some(i) = queue.pop_front() {
        let w = result[&i].clone();

        for k in 0..table.len() {
            for g in table.all_gens() {
                if !result.contains_key(&k) {
                    result.insert(k, &w * g);
                    queue.push_back(k);
                }
            }
        }
    }

    result
}


fn induced_table<T, F>(nr_gens: usize, img: F, start: &T) -> CosetTable
    where
        T: Clone + Eq + std::hash::Hash,
        F: Fn(&T, isize) -> T
{
    let mut table = CosetTable::new(nr_gens);
    let mut o2n = HashMap::from([(start.clone(), 0)]);
    let mut n2o = HashMap::from([(0, start.clone())]);

    for i in 0.. {
        if i >= table.len() {
            break;
        }
        for g in table.all_gens() {
            let k = img(&n2o[&i], g);
            let n = *o2n.entry(k.clone()).or_insert(table.len());
            n2o.insert(n, k);
            table.join(i, n, g);
        }
    }

    table.compact()
}


pub fn core_table(base: &CosetTable) -> CosetTable {
    let img = |es: &Vec<usize>, g| es.iter()
        .map(|&e| base.get(e, g).unwrap())
        .collect();
    let start = &(0..base.len()).collect();

    induced_table(base.nr_gens, img, start)
}


pub fn intersection_table(ta: &CosetTable, tb: &CosetTable) -> CosetTable {
    assert_eq!(ta.nr_gens, tb.nr_gens);

    let img = |&(a, b): &(usize, usize), g| {
        (ta.get(a, g).unwrap(), tb.get(b, g).unwrap())
    };
    let start = &(0 as usize, 0 as usize);

    induced_table(ta.nr_gens, img, start)
}


fn first_free_in_table(table: &CosetTable) -> Option<(usize, isize)> {
    for k in 0..table.len() {
        for g in table.all_gens() {
            if table.get(k, g).is_none() {
                return Some((k, g));
            }
        }
    }
    None
}


fn derived_table(
    table: &CosetTable,
    expanded_rels: &Vec<FreeWord>,
    from: usize,
    to: usize,
    g: isize
)
    -> Option<CosetTable>
{
    if table.get(from, g).is_some() || table.get(to, -g).is_some() {
        return None;
    }

    let mut result = table.clone();
    result.join(from, to, g);

    let mut q = VecDeque::from([from]);

    while let Some(row) = q.pop_front() {
        for rel in expanded_rels {
            let (head, tail, gap, c) = scan_both_ways(&result, rel, row);
            if gap == 1 {
                result.join(head, tail, c);
                q.push_back(head);
            } else if gap == 0 && head != tail {
                return None;
            }
        }
    }

    Some(result)
}


fn potential_children(
    table: &CosetTable, expanded_rels: &Vec<FreeWord>, max_rows: usize
)
    -> Vec<CosetTable>
{
    let mut result = vec![];

    if let Some((k, g)) = first_free_in_table(table) {
        let limit = max_rows.min(table.len() + 1);
        for pos in k..limit {
            if let Some(t) = derived_table(table, expanded_rels, k, pos, g) {
                result.push(t);
            }
        }
    }

    result
}


fn compare_renumbered_from(table: &CosetTable, start: usize) -> isize {
    let n = table.len();
    let mut n2o = HashMap::from([(0, start)]);
    let mut o2n = HashMap::from([(start, 0)]);

    for row in 0..table.len() {
        assert!(row < n2o.len(), "coset table is not transitive");

        for g in table.all_gens() {
            let oval = table.get(row, g).unwrap_or(n);

            let nval = if let Some(t) = table.get(n2o[&row], g) {
                if !o2n.contains_key(&t) {
                    let n = n2o.len();
                    o2n.insert(t, n);
                    n2o.insert(n, t);
                }
                *o2n.get(&t).unwrap()
            } else {
                n
            };

            let result = nval as isize - oval as isize;

            if result != 0 {
                return result;
            }
        }
    }
    0
}


fn is_canonical(table: &CosetTable) -> bool {
    for start in 1..table.len() {
        if compare_renumbered_from(table, start) < 0 {
            return false;
        }
    }
    true
}


struct CosetTableBacktracking {
    nr_gens: usize,
    expanded_relators: Vec<FreeWord>,
    max_rows: usize,
}


impl BackTracking for CosetTableBacktracking {
    type State = CosetTable;
    type Item = CosetTable;

    fn root(&self) -> Self::State {
        CosetTable::new(self.nr_gens)
    }

    fn extract(&self, state: &Self::State) -> Option<Self::Item> {
        if first_free_in_table(state).is_none() {
            Some(state.compact())
        } else {
            None
        }
    }

    fn children(&self, state: &Self::State) -> Vec<Self::State> {
        potential_children(state, &self.expanded_relators, self.max_rows)
            .iter()
            .cloned()
            .filter(is_canonical)
            .collect()
    }
}


pub struct CosetTables {
    bt: BackTrackIterator<CosetTableBacktracking>
}


impl CosetTables {
    fn new(nr_gens: usize, rels: &Vec<FreeWord>, max_rows: usize)
        -> CosetTables
    {
        let expanded_relators = expanded_relator_set(&rels)
            .iter().cloned().collect();

        CosetTables {
            bt: BackTrackIterator::new(
                CosetTableBacktracking { nr_gens, expanded_relators, max_rows }
            )
        }
    }
}


impl Iterator for CosetTables {
    type Item = CosetTable;

    fn next(&mut self) -> Option<Self::Item> {
        self.bt.next()
    }
}


pub fn coset_tables(nr_gens: usize, rels: &Vec<FreeWord>, max_rows: usize)
    -> CosetTables
{
    CosetTables::new(nr_gens, rels, max_rows)
}


#[cfg(test)]
mod coset_tests {
    use super::*;

    fn make_table(nr_gens: usize, rels: &[&[isize]], subgens: &[&[isize]])
        -> CosetTable
    {
        coset_table(
            nr_gens,
            &rels.iter().map(|&r| Vec::from(r).into()).collect(),
            &subgens.iter().map(|&g| Vec::from(g).into()).collect(),
        )
    }

    #[test]
    fn test_coset_table() {
        assert_eq!(
            make_table(2, &[&[1, 1], &[2, 2], &[1, 2, 1, 2]], &[]).as_map(),
            vec![
                BTreeMap::from([(1, 1), (2, 2), (-1, 1), (-2, 2)]),
                BTreeMap::from([(1, 0), (2, 3), (-1, 0), (-2, 3)]),
                BTreeMap::from([(1, 3), (2, 0), (-1, 3), (-2, 0)]),
                BTreeMap::from([(1, 2), (2, 1), (-1, 2), (-2, 1)]),
            ]
        );
        assert_eq!(
            make_table(2, &[&[1, 1], &[2, 2], &[1, 2, 1, 2]], &[&[1]]).as_map(),
            vec![
                BTreeMap::from([(1, 0), (2, 1), (-1, 0), (-2, 1)]),
                BTreeMap::from([(1, 1), (2, 0), (-1, 1), (-2, 0)]),
            ]
        );
        assert_eq!(
            make_table(2, &[&[1, 1], &[2, 2], &[1, 2, 1, 2]], &[&[2]]).as_map(),
            vec![
                BTreeMap::from([(1, 1), (2, 0), (-1, 1), (-2, 0)]),
                BTreeMap::from([(1, 0), (2, 1), (-1, 0), (-2, 1)]),
            ]
        );
        assert_eq!(
            make_table(2, &[&[1, 1], &[2, 2], &[1, 2, 1, 2]], &[&[1], &[2]]).as_map(),
            vec![
                BTreeMap::from([(1, 0), (2, 0), (-1, 0), (-2, 0)]),
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

    #[test]
    fn test_core_table() {
        assert_eq!(
            core_table(&make_table(
                3,
                &[
                    &[1, 1], &[2, 2], &[3, 3],
                    &[1, 2, 1, 2, 1, 2], &[1, 3, 1, 3], &[3, 2, 3, 2, 3, 2]
                ],
                &[&[1, 2]]
            )).len(),
            24
        );
    }

    #[test]
    fn test_enumerator() {
        assert_eq!(
            coset_tables(
                2,
                &vec![
                    FreeWord::from([1, 1]),
                    FreeWord::from([2, 2]),
                    FreeWord::from([1, 2, 1, 2]),
                ],
                8
            )
                .map(|t| t.len())
                .collect::<Vec<_>>(),
            vec![1, 2, 2, 2, 4]
        );
    }
}
