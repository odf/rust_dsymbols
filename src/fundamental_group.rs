use std::collections::{HashMap, HashSet, VecDeque};

use crate::dsyms::*;
use crate::fpgroups::free_words::FreeWord;


type Edge = (usize, usize);
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
                    if i != j {
                        opposite.insert((d, i, j), ((d, j, i), 1));
                    }
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


fn spanning_tree<T: DSym>(ds: &T) -> Vec<Ridge> {
    let mut seen = HashSet::new();
    let mut result = Vec::new();

    // elements reversal allows for use of old JS results in testing
    for (maybe_i, d, di) in ds.traversal(0..=ds.dim(), (1..=ds.size()).rev()) {
        if !seen.contains(&di) {
            if let Some(i) = maybe_i {
                result.push((d, i, i));
            }
            seen.insert(di);
        }
    }

    result
}


pub fn inner_edges<T: DSym>(ds: &T) -> Vec<Edge> {
    let glued = Boundary::new(ds).glue_recursively(spanning_tree(ds));

    glued.iter().map(|&(d, i, _)| (d, i)).collect()
}


fn trace_word<T: DSym>(
    ds: &T,
    edge_to_word: &HashMap<Edge, FreeWord>,
    d: usize,
    i: usize,
    j: usize
)
    -> FreeWord
{
    let mut e = d;
    let mut result = FreeWord::empty();

    loop {
        result *= edge_to_word.get(&(e, i)).unwrap_or(&FreeWord::empty());
        e = ds.op(i, e).unwrap_or(e);
        result *= edge_to_word.get(&(e, j)).unwrap_or(&FreeWord::empty());
        e = ds.op(j, e).unwrap_or(e);

        if e == d {
            break;
        }
    }

    result
}


fn find_generators<T: DSym>(ds: &T)
    -> (HashMap<Edge, FreeWord>, HashMap<usize, Edge>)
{
    let mut edge_to_word = HashMap::new();
    let mut gen_to_edge = HashMap::new();

    let mut bnd = Boundary::new(ds);
    bnd.glue_recursively(spanning_tree(ds));

    for d in 1..=ds.size() {
        for i in 0..=ds.dim() {
            if (0..=ds.dim()).any(|j| bnd.opposite(d, i, j).is_some()) {
                let di = ds.op(i, d).unwrap();
                let gen = gen_to_edge.len() + 1;
                gen_to_edge.insert(gen, (d, i));

                edge_to_word.insert((d, i), FreeWord::from([gen as isize]));
                edge_to_word.insert((di, i), FreeWord::from([-(gen as isize)]));

                for (d, i, j) in bnd.glue_recursively(vec![(d, i, i)]) {
                    let w = trace_word(ds, &edge_to_word, di, j, i);
                    if w.len() > 0 {
                        edge_to_word.insert((d, i), w.inverse());
                        edge_to_word.insert((di, i), w);
                    }
                }
            }
        }
    }

    (edge_to_word, gen_to_edge)
}


#[test]
fn test_spanning_tree() {
    let tree = |s: &str|
        spanning_tree(&s.parse::<PartialDSym>().unwrap()).iter()
            .map(|&(d, i, _)| (d, i))
            .collect::<Vec<_>>();

    assert_eq!(tree("<1.1:3:1 2 3,1 3,2 3:4 8,3>"), vec![(3, 1), (2, 2)]);
    assert_eq!(tree("<1.1:2 3:2,1 2,1 2,2:6,3 2,6>"), vec![(2, 0)]);
    assert_eq!(
        tree("<1.1:8:2 4 6 8,8 3 5 7,1 2 3 4 5 6 7 8:4,4 6 8 4>"),
        vec![(8, 0), (7, 1), (6, 0), (5, 1), (4, 0), (3, 1), (2, 0)]
    );
    assert_eq!(
        tree("
            <1.1:12:2 5 7 10 11 12,1 4 6 9 7 12 11,3 5 8 6 11 12 10:4 4,6 3 3>
        "),
        vec![
            (12, 0), (9, 1), (5, 0), (3, 1), (6, 0), (10, 1), (11, 0), (5, 2),
            (2, 0), (2, 1), (4, 0),
        ]
    );
    assert_eq!(
        tree("
            <1.1:24:
            2 4 6 8 10 12 14 16 18 20 22 24,
            16 3 5 7 9 11 13 15 24 19 21 23,
            10 9 20 19 14 13 22 21 24 23 18 17:
            8 4,3 3 3 3
            >
        "),
        vec![
            (24,0), (23, 1), (22, 0), (21, 1), (20, 0), (19, 1), (18, 0),
            (24, 2), (11, 0), (12, 1), (13, 0), (14, 1), (15, 0), (16, 1),
            (1, 0), (2, 1), (3, 0), (4, 1), (5, 0), (6, 1), (7, 0), (8, 1),
            (9, 0)
        ]
    );
}


#[test]
fn test_inner_edges() {
    let inner = |s: &str| inner_edges(&s.parse::<PartialDSym>().unwrap());

    assert_eq!(inner("<1.1:3:1 2 3,1 3,2 3:4 8,3>"), vec![(3, 1), (2, 2)]);
    assert_eq!(inner("<1.1:2 3:2,1 2,1 2,2:6,3 2,6>"), vec![(2, 0)]);
    assert_eq!(
        inner("<1.1:8:2 4 6 8,8 3 5 7,1 2 3 4 5 6 7 8:4,4 6 8 4>"),
        vec![(8, 0), (7, 1), (6, 0), (5, 1), (4, 0), (3, 1), (2, 0), (8, 1)]
    );
    assert_eq!(
        inner("
            <1.1:12:2 5 7 10 11 12,1 4 6 9 7 12 11,3 5 8 6 11 12 10:4 4,6 3 3>
        "),
        vec![
            (12, 0), (9, 1), (5, 0), (3, 1), (6, 0), (10, 1), (11, 0), (5, 2),
            (2, 0), (2, 1), (4, 0), (12, 1), (3, 2)
        ]
    );
    assert_eq!(
        inner("
            <1.1:24:
            2 4 6 8 10 12 14 16 18 20 22 24,
            16 3 5 7 9 11 13 15 24 19 21 23,
            10 9 20 19 14 13 22 21 24 23 18 17:
            8 4,3 3 3 3
            >
        "),
        vec![
            (24,0), (23, 1), (22, 0), (21, 1), (20, 0), (19, 1), (18, 0),
            (24, 2), (11, 0), (12, 1), (13, 0), (14, 1), (15, 0), (16, 1),
            (1, 0), (2, 1), (3, 0), (4, 1), (5, 0), (6, 1), (7, 0), (8, 1),
            (9, 0), (24, 1), (23, 2), (11, 1)
        ]
    );
}


#[test]
fn test_find_generators() {
    let find = |s: &str| find_generators(&s.parse::<PartialDSym>().unwrap());

    assert_eq!(
        find("<1.1:3:1 2 3,1 3,2 3:4 8,3>").0,
            //(
                HashMap::from([
                    ((1, 0), FreeWord::from([-1])),
                    ((1, 1), FreeWord::from([-2])),
                    ((2, 0), FreeWord::from([-1])),
                    ((3, 0), FreeWord::from([-3])),
                    ((3, 2), FreeWord::from([-2])),
                ]),
            //    HashMap::from([(1, (1, 0)), (2, (1, 1)), (3, (3, 0))])
            //)
    );
}
