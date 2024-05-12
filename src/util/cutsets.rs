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
