use std::collections::{VecDeque, HashSet, BTreeMap};

use crate::dsets::DSet;


pub struct Traversal<'a, T: DSet> {
    ds: &'a T,
    indices: Vec<usize>,
    seeds_left: VecDeque<usize>,
    seen: HashSet<(usize, Option<usize>)>,
    todo: BTreeMap<usize, VecDeque<usize>>,
}


impl<'a, T: DSet> Traversal<'a, T> {
    pub fn new(ds: &'a T, indices: &[usize], seeds: &[usize]) -> Self {
        let indices: Vec<_> = indices.into();
        let seeds_left = Vec::from(seeds).into();
        let seen = HashSet::new();
        let todo: BTreeMap<_, _> = indices.iter()
            .map(|&i| (i, VecDeque::new()))
            .collect();

        Self { ds, indices, seeds_left, seen, todo }
    }
}


impl<'a, T: DSet> Iterator for Traversal<'a, T> {
    type Item = (usize, Option<usize>, usize);

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let maybe_i = self.todo.iter()
                .find(|(_, q)| !q.is_empty())
                .and_then(|(&i, _)| Some(i));

            let maybe_d = if let Some(i) = maybe_i {
                self.todo.get_mut(&i).and_then(|q| q.pop_front())
            } else {
                self.seeds_left.pop_front()
            };

            if let Some(d) = maybe_d {
                if !self.seen.contains(&(d, maybe_i)) {
                    let di = if let Some(i) = maybe_i {
                        self.ds.op(i, d).unwrap()
                    } else {
                        d
                    };

                    for k in self.indices.iter() {
                        if let Some(q) = self.todo.get_mut(k) {
                            if *k < 2 {
                                q.push_front(di)
                            } else {
                                q.push_back(di)
                            }
                        }
                    }

                    self.seen.insert((di, None));
                    self.seen.insert((di, maybe_i));
                    self.seen.insert((d, maybe_i));

                    return Some((d, maybe_i, di))
                }
            } else {
                return None
            }
        }
    }
}


#[cfg(test)]
mod traversal_tests {
    use crate::dsyms::PartialDSym;

    use super::*;

    #[test]
    fn test_traversal() {
        let s = "<1.1:8:2 4 6 8,8 3 5 7,1 2 3 4 5 6 7 8:4,4 6 8 4>";
        let dsym: PartialDSym = s.parse().unwrap();

        assert_eq!(
            Traversal::new(&dsym, &[0, 1, 2], &[1]).collect::<Vec<_>>(),
            vec![
                (1, None, 1),
                (1, Some(0), 2),
                (2, Some(1), 3),
                (3, Some(0), 4),
                (4, Some(1), 5),
                (5, Some(0), 6),
                (6, Some(1), 7),
                (7, Some(0), 8),
                (8, Some(1), 1),
                (1, Some(2), 1),
                (2, Some(2), 2),
                (3, Some(2), 3),
                (4, Some(2), 4),
                (5, Some(2), 5),
                (6, Some(2), 6),
                (7, Some(2), 7),
                (8, Some(2), 8)            
            ]
        );

        assert_eq!(
            Traversal::new(&dsym, &[0, 2], &[1, 2, 3]).collect::<Vec<_>>(),
            vec![
                (1, None, 1),
                (1, Some(0), 2), (1, Some(2), 1),  (2, Some(2), 2),
                (3, None, 3),
                (3, Some(0), 4), (3, Some(2), 3),  (4, Some(2), 4),
            ]
        );
    }
}
