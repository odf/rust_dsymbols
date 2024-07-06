use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};

use crate::dsyms::*;
use crate::fpgroups::free_words::{relator_representative, FreeWord};


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

    fn glue_recursively(&mut self, todo: Vec<(usize, usize, Option<usize>)>)
        -> Vec<(usize, usize, Option<usize>)>
    {
        let mut todo = VecDeque::from(todo);
        let mut result = Vec::new();

        while let Some(next) = todo.pop_front() {
            let (d, i, j) = next;

            let good = if let Some(j) = j {
                let t = if self.ds.op(i, d) == Some(d) { 1 } else { 2 };
                let m = self.ds.m(i, j, d).unwrap() * t;
                self.opposite(d, i, j).is_some_and(|(_, n)| n == m)
            } else {
                true
            };

            if good {
                for (d, i, j) in self.glue(d, i) {
                    todo.push_back((d, i, Some(j)));
                }
                result.push(next);
            }
        }

        result
    }
}


fn spanning_tree<T: DSym>(ds: &T) -> Vec<(usize, usize, Option<usize>)> {
    let mut seen = HashSet::new();
    let mut result = Vec::new();

    // elements reversal allows for use of old JS results in testing
    for (maybe_i, d, di) in ds.traversal(0..=ds.dim(), (1..=ds.size()).rev()) {
        if !seen.contains(&di) {
            if let Some(i) = maybe_i {
                result.push((d, i, None));
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
    edge_to_word: &BTreeMap<Edge, FreeWord>,
    d: usize,
    i: Option<usize>,
    j: Option<usize>
)
    -> FreeWord
{
    let nil = FreeWord::empty();
    let mut result = FreeWord::empty();

    if let Some(i) = i {
        if let Some(j) = j {
            let mut e = d;
            loop {
                result *= edge_to_word.get(&(e, i)).unwrap_or(&nil);
                e = ds.op(i, e).unwrap_or(e);
                result *= edge_to_word.get(&(e, j)).unwrap_or(&nil);
                e = ds.op(j, e).unwrap_or(e);

                if e == d {
                    break;
                }
            }
        } else {
            result *= edge_to_word.get(&(d, i)).unwrap_or(&nil);
        }
    } else if let Some(j) = j {
        result *= edge_to_word.get(&(d, j)).unwrap_or(&nil);
    }

    result
}


fn find_generators<T: DSym>(ds: &T)
    -> (BTreeMap<Edge, FreeWord>, BTreeMap<usize, Edge>)
{
    let mut edge_to_word = BTreeMap::new();
    let mut gen_to_edge = BTreeMap::new();

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

                let glued = bnd.glue_recursively(vec![(d, i, None)]);

                for (e, i, j) in glued {
                    let ei = ds.op(i, e).unwrap();
                    let w = trace_word(ds, &edge_to_word, ei, j, Some(i));
                    if w.len() > 0 {
                        edge_to_word.insert((e, i), w.inverse());
                        edge_to_word.insert((ei, i), w);
                    }
                }
            }
        }
    }

    (edge_to_word, gen_to_edge)
}


pub struct FundamentalGroup {
    pub relators: Vec<FreeWord>,
    pub cones: BTreeSet<(FreeWord, usize)>,
    pub gen_to_edge: BTreeMap<usize, Edge>,
    pub edge_to_word: BTreeMap<Edge, FreeWord>,
}


pub fn fundamental_group<T: DSym>(ds: &T) -> FundamentalGroup {
    let (edge_to_word, gen_to_edge) = find_generators(ds);
    let mut cones = BTreeSet::new();
    let mut relators = BTreeSet::new();

    for i in 0..=ds.dim() {
        for j in i..=ds.dim() {
            for d in ds.orbit_reps_2d(i, j) {
                let di = ds.op(i, d).unwrap();
                let word = trace_word(ds, &edge_to_word, di, Some(j), Some(i));
                let degree = ds.v(i, j, d).unwrap();
                let rel = word.raised_to(degree as isize);

                if rel.len() > 0 {
                    relators.insert(relator_representative(&rel));
                }

                if degree > 1 {
                    cones.insert((relator_representative(&word), degree));
                }
            }
        }
    }

    let relators = relators.iter().cloned().collect();

    FundamentalGroup { relators, cones, gen_to_edge, edge_to_word }
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
fn test_fundamental_group_a() {
    let group = |s: &str| fundamental_group(&s.parse::<PartialDSym>().unwrap());

    let g = group("<1.1:1 3:1,1,1,1:4,3,4>");
    assert_eq!(
        g.gen_to_edge,
        BTreeMap::from([(1, (1, 0)), (2, (1, 1)), (3, (1, 2)), (4, (1, 3))])
    );
    assert_eq!(
        g.edge_to_word,
        BTreeMap::from([
            ((1, 0), FreeWord::from([-1])),
            ((1, 1), FreeWord::from([-2])),
            ((1, 2), FreeWord::from([-3])),
            ((1, 3), FreeWord::from([-4])),
        ])
    );
    assert_eq!(
        g.relators.iter().cloned().collect::<HashSet<_>>(),
        HashSet::from([
            FreeWord::from([1, 1]),
            FreeWord::from([1, 2]).raised_to(4),
            FreeWord::from([1, 3]).raised_to(2),
            FreeWord::from([1, 4]).raised_to(2),
            FreeWord::from([2, 2]),
            FreeWord::from([2, 3]).raised_to(3),
            FreeWord::from([2, 4]).raised_to(2),
            FreeWord::from([3, 3]),
            FreeWord::from([3, 4]).raised_to(4),
            FreeWord::from([4, 4]),
        ])
    );
    assert_eq!(
        g.cones,
        BTreeSet::from([
            (FreeWord::from([1, 2]), 4),
            (FreeWord::from([1, 3]), 2),
            (FreeWord::from([1, 4]), 2),
            (FreeWord::from([2, 3]), 3),
            (FreeWord::from([2, 4]), 2),
            (FreeWord::from([3, 4]), 4),
        ])
    );
}


#[test]
fn test_fundamental_group_b() {
    let group = |s: &str| fundamental_group(&s.parse::<PartialDSym>().unwrap());

    let g = group("<1.1:2:2,2,2:4,3>");
    assert_eq!(
        g.gen_to_edge,
        BTreeMap::from([(1, (1, 1)), (2, (1, 2))])
    );
    assert_eq!(
        g.edge_to_word,
        BTreeMap::from([
            ((1, 1), FreeWord::from([1])),
            ((1, 2), FreeWord::from([2])),
            ((2, 1), FreeWord::from([-1])),
            ((2, 2), FreeWord::from([-2])),
        ])
    );
    assert_eq!(
        g.relators.iter().cloned().collect::<HashSet<_>>(),
        HashSet::from([
            FreeWord::from([1]).raised_to(4),
            FreeWord::from([2]).raised_to(2),
            FreeWord::from([1, -2]).raised_to(3),
        ])
    );
    assert_eq!(
        g.cones,
        BTreeSet::from([
            (FreeWord::from([1]), 4),
            (FreeWord::from([2]), 2),
            (FreeWord::from([1, -2]), 3),
        ])
    );
}


#[test]
fn test_fundamental_group_c() {
    let group = |s: &str| fundamental_group(&s.parse::<PartialDSym>().unwrap());

    let g = group("<1.1:3:1 2 3,1 3,2 3:4 8,3>");
    assert_eq!(
        g.gen_to_edge,
        BTreeMap::from([(1, (1, 0)), (2, (1, 1)), (3, (3, 0))])
    );
    assert_eq!(
        g.edge_to_word,
        BTreeMap::from([
            ((1, 0), FreeWord::from([-1])),
            ((1, 1), FreeWord::from([-2])),
            ((2, 0), FreeWord::from([-1])),
            ((3, 0), FreeWord::from([-3])),
            ((3, 2), FreeWord::from([-2])),
        ])
    );
    assert_eq!(
        g.relators.iter().cloned().collect::<HashSet<_>>(),
        HashSet::from([
            FreeWord::from([1, 1]),
            FreeWord::from([2, 2]),
            FreeWord::from([3, 3]),
            FreeWord::from([1, 2, 1, 2, 1, 2, 1, 2]),
            FreeWord::from([1, 3, 1, 3, 1, 3, 1, 3]),
            FreeWord::from([2, 3, 2, 3]),
        ])
    );
    assert_eq!(
        g.cones,
        BTreeSet::from([
            (FreeWord::from([1, 2]), 4),
            (FreeWord::from([1, 3]), 4),
            (FreeWord::from([2, 3]), 2),
        ])
    );
}


#[test]
fn test_fundamental_group_d() {
    let group = |s: &str| fundamental_group(&s.parse::<PartialDSym>().unwrap());

    let g = group("
        <1.1:24:
        2 4 6 8 10 12 14 16 18 20 22 24,
        16 3 5 7 9 11 13 15 24 19 21 23,
        10 9 20 19 14 13 22 21 24 23 18 17:
        8 4,3 3 3 3
        >");
    assert_eq!(g.gen_to_edge, BTreeMap::from([(1, (1, 2)), (2, (3, 2))]));
    assert_eq!(
        g.relators.iter().cloned().collect::<HashSet<_>>(),
        HashSet::from([FreeWord::from([1, 2, -1, -2])])
    );
    assert_eq!(g.cones, BTreeSet::from([]));
}
