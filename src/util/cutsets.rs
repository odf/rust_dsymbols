use std::collections::{BTreeMap, BTreeSet, VecDeque};


pub struct EdgeCut {
    pub cut_edges: Vec<(usize, usize)>,
    pub inside_vertices: Vec<usize>,
}


pub fn min_edge_cut<I>(edges: I, source: usize, sink: usize)
    -> EdgeCut
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
            let cut_edges = edges.iter()
                .filter(|&(v, w)| seen.contains(v) && !seen.contains(w))
                .cloned()
                .collect();
            let inside_vertices = seen.iter().cloned().collect();
            return EdgeCut { cut_edges, inside_vertices };
        }
    }
}


pub fn min_edge_cut_undirected<I>(edges: I, source: usize, sink: usize)
    -> EdgeCut
    where I: IntoIterator<Item=(usize, usize)>
{
    let edges: BTreeSet<_> = edges.into_iter()
        .flat_map(&|(v, w)| [(v, w), (w, v)])
        .collect();

    min_edge_cut(edges, source, sink)
}


pub struct VertexCut {
    pub cut_vertices: Vec<usize>,
    pub inside_vertices: Vec<usize>,
}


pub fn min_vertex_cut<I>(edges: I, source: usize, sink: usize)
    -> VertexCut
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

    let edge_cut = min_edge_cut(x_edges, source + offset, sink);
    let cut_vertices: Vec<_> = edge_cut.cut_edges.iter()
        .map(|&(v, w)| v.min(w))
        .collect();
    let inside_vertices = edge_cut.inside_vertices.iter()
        .filter(|&&v| v < offset && !cut_vertices.contains(&v))
        .cloned()
        .collect();

    VertexCut { cut_vertices, inside_vertices }
}


pub fn min_vertex_cut_undirected<I>(edges: I, source: usize, sink: usize)
    -> VertexCut
    where I: IntoIterator<Item=(usize, usize)>
{
    let edges: BTreeSet<_> = edges.into_iter()
        .flat_map(&|(v, w)| [(v, w), (w, v)])
        .collect();

    min_vertex_cut(edges, source, sink)
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

    fn nice_vs(vertices: Vec<usize>) -> Vec<usize>
    {
        let mut vs = vertices.clone();
        vs.sort();
        vs
    }

    fn nice_es(edges: Vec<(usize, usize)>) -> Vec<(usize, usize)>
    {
        let mut es: Vec<_> = edges.iter()
            .map(|&(v, w)| if w < v { (w, v) } else { (v, w) })
            .collect();

        es.sort();
        es
    }

    #[test]
    fn test_min_edge_cut() {
        let cut = min_edge_cut(example(), 1, 11);
        assert_eq!(nice_es(cut.cut_edges), [(1, 4), (7, 11)]);
    }

    #[test]
    fn test_min_edge_cut_undirected() {
        let cut = min_edge_cut_undirected(example(), 1, 11);
        assert_eq!(nice_es(cut.cut_edges), [(1, 4), (4, 9), (7, 11), (9, 11)]);
    }

    #[test]
    fn test_min_vertex_cut() {
        let cut = min_vertex_cut(example(), 1, 11);
        assert_eq!(nice_vs(cut.cut_vertices), [4, 7]);
    }

    #[test]
    fn test_min_vertex_cut_undirected() {
        let cut = min_vertex_cut_undirected(example(), 1, 11);
        assert_eq!(nice_vs(cut.cut_vertices), [4, 7, 9]);
    }
}
