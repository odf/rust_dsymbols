use std::collections::{BTreeMap, BTreeSet, VecDeque};


pub fn min_edge_cut<I>(edges: I, source: usize, sink: usize)
    -> Vec<(usize, usize)>
    where I: IntoIterator<Item=(usize, usize)>
{
    let edges: BTreeSet<_> = edges.into_iter().collect();

    let neighbors = by_first(
        edges.iter().flat_map(|&(v, w)| [(v, w), (w, v)])
    );

    let mut path_edges = BTreeSet::new();

    loop {
        let (next, seen) = augment(
            &edges, &neighbors, source, sink, &path_edges
        );
        if let Some(next) = next {
            path_edges = next;
        } else {
            return edges.iter()
                .filter(|&(v, w)| seen.contains(v) && !seen.contains(w))
                .cloned()
                .collect();
        }
    }
}


pub fn min_vertex_cut<I>(edges: I, source: usize, sink: usize)
    -> Vec<usize>
    where I: IntoIterator<Item=(usize, usize)>
{
    let edges: BTreeSet<_> = edges.into_iter().collect();

    let vertices: BTreeSet<_> = edges.iter()
        .flat_map(|&(v, w)| [v, w])
        .collect();

    let offset = vertices.iter().max().unwrap_or(&0) + 1;

    let x_edges: Vec<_> = std::iter::empty()
        .chain(edges.iter().map(|&(v, w)| (v + offset, w)))
        .chain(vertices.iter().map(|&v| (v, v + offset)))
        .collect();

    min_edge_cut(x_edges, source + offset, sink).iter()
        .map(|&(v, w)| v.min(w))
        .collect()
}


fn by_first<I>(pairs: I) -> BTreeMap<usize, BTreeSet<usize>>
    where I: IntoIterator<Item=(usize, usize)>
{
    let mut result: BTreeMap<_, BTreeSet<_>> = BTreeMap::new();

    for (v, w) in pairs {
        result.entry(v)
            .and_modify(|a| { a.insert(w); })
            .or_insert(BTreeSet::from([w]));
    }

    result
}


fn augment(
    edges: &BTreeSet<(usize, usize)>,
    neighbors: &BTreeMap<usize, BTreeSet<usize>>,
    source: usize,
    sink: usize,
    path_edges: &BTreeSet<(usize, usize)>
)
    -> (Option<BTreeSet<(usize, usize)>>, BTreeSet<usize>)
{
    let mut q = VecDeque::from([source]);
    let mut seen = BTreeSet::from([source]);
    let mut back = BTreeMap::new();

    while let Some(v) = q.pop_front() {
        for &w in &neighbors[&v] {
            if !seen.contains(&w) && !path_edges.contains(&(v, w)) {
                if edges.contains(&(v, w)) || path_edges.contains(&(w, v)) {
                    back.insert(w, v);
                    seen.insert(w);
                    q.push_back(w);
                }
            }
        }
        if back.contains_key(&sink) {
            break;
        }
    }

    if back.contains_key(&sink) {
        let mut result = path_edges.clone();
        let mut w = sink;
        while w != source {
            let v = back[&w];
            if result.contains(&(w, v)) {
                result.remove(&(w, v));
            } else {
                result.insert((v, w));
            }
            w = v;
        }

        return (Some(result), seen);
    } else {
        return (None, seen);
    }
}


#[cfg(test)]
mod test {
    use super::*;

    fn example() -> Vec<(usize, usize)> {
        vec![
            (1, 2), (1, 3), (1, 4), (1, 5), (1, 6),
            (2, 7), (3, 7), (3, 9), (4, 8), (4, 9), (4, 10), (5, 9), (6, 9),
            (7, 11), (8, 11), (10, 11),
            (11, 9)
        ]
    }

    #[test]
    fn test_min_edge_cut() {
        let mut cut = min_edge_cut(example(), 1, 11);
        cut.sort();
        assert_eq!(cut, vec![(1, 4), (7, 11)]);
    }

    #[test]
    fn test_min_vertex_cut() {
        let mut cut = min_vertex_cut(example(), 1, 11);
        cut.sort();
        assert_eq!(cut, [4, 7]);
    }
}
