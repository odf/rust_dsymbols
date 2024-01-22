use std::cell::UnsafeCell;
use std::hash::Hash;
use std::collections::HashMap;


#[derive(Clone)]
struct PartitionImpl<T> {
    index: HashMap<T, usize>,
    elements: Vec<T>,
    rank: Vec<usize>,
    parent: Vec<usize>,
}

impl<T> PartitionImpl<T> where T: Clone + Eq + Hash {
    fn new() -> Self {
        PartitionImpl {
            index: HashMap::new(),
            elements: vec![],
            rank: vec![],
            parent: vec![],
        }
    }

    fn get_index(&mut self, a: &T) -> usize {
        if let Some(x) = self.index.get(a) {
            *x
        } else {
            let i = self.elements.len();
            self.index.insert(a.clone(), i);
            self.elements.push(a.clone());
            self.rank.push(0);
            self.parent.push(i);
            i
        }
    }

    fn root_index(&mut self, a: &T) -> usize {
        let mut x = self.get_index(a);
        let mut root = x;

        while self.parent[root] != root {
            root = self.parent[root];
        }

        while x != root {
            let t = x;
            x = self.parent[x];
            self.parent[t] = root;
        }

        root
    }

    fn find(&mut self, a: &T) -> T {
        let root = self.root_index(a);
        self.elements[root].clone()
    }

    fn unite(&mut self, a: &T, b: &T) {
        let x = self.root_index(a);
        let y = self.root_index(b);

        if x != y {
            let rx = self.rank[x];
            let ry = self.rank[y];

            if rx < ry {
                self.parent[x] = y;
            } else {
                if rx == ry {
                    self.rank[x] = rx + 1;
                }
                self.parent[y] = x;
            }
        }
    }
}


pub struct Partition<T> {
    _impl: UnsafeCell<PartitionImpl<T>>,
}


impl<T> Partition<T> where T: Clone + Eq + Hash {
    pub fn new() -> Self {
        Partition { _impl: UnsafeCell::new(PartitionImpl::new())}
    }

    pub fn find(&self, x: &T) -> T {
        unsafe { (*self._impl.get()).find(x) }
    }

    pub fn unite(&mut self, x: &T, y: &T) {
        unsafe { (*self._impl.get()).unite(x, y) };
    }

    pub fn classes(&self, elms: &[T]) -> Vec<Vec<T>> {
        let mut class_for_rep = HashMap::new();
        let mut classes = vec![];

        for e in elms {
            let rep = self.find(e);
            if let Some(cl) = class_for_rep.get(&rep) {
                let class: &mut Vec<_> = &mut classes[*cl];
                class.push(e.clone());
            } else {
                class_for_rep.insert(rep, classes.len());
                classes.push(vec![e.clone()]);
            }
        }

        classes
    }
}


impl<T> Clone for Partition<T> where T: Clone {
    fn clone(&self) -> Self {
        Self {
            _impl: UnsafeCell::new(unsafe { (*self._impl.get()).clone() })
        }
    }
}


#[test]
pub fn test_partition() {
    let p = {
        let mut p = Partition::new();
        for (a, b) in [(1, 2), (3, 4), (5, 6), (7, 8), (2, 3), (1, 6)] {
            p.unite(&a, &b);
        }
        p
    };

    let test = HashMap::from([
        (0, 0),
        (1, 1), (2, 1), (3, 1), (4, 1), (5, 1), (6, 1),
        (7, 2), (8, 2),
        (9, 3),
    ]);

    for a in 0..=9 {
        for b in 0..=9 {
            if test[&a] == test[&b] {
                assert_eq!(p.find(&a), p.find(&b));
            } else {
                assert_ne!(p.find(&a), p.find(&b));
            }
        }
    }

    let elms: Vec<_> = (0..=9).collect();
    let cl = p.classes(&elms);
    assert_eq!(cl, vec![vec![0], vec![1, 2, 3, 4, 5, 6], vec![7, 8], vec![9]]);
}
