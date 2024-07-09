use std::collections::{BTreeSet, HashMap, HashSet, VecDeque};
use std::fmt;

use crate::derived::{build_set, build_sym_using_vs, canonical};
use crate::dsets::{DSet, PartialDSet};
use crate::dsyms::{DSym, PartialDSym};
use crate::fundamental_group::inner_edges;
use crate::util::cutsets::min_vertex_cut_undirected;


struct DrawingInstructions<'a>
{
    ds: &'a PartialDSet
}


impl<'a> DrawingInstructions<'a>
{
    pub fn new(ds: &'a PartialDSet) -> Self {
        assert_eq!(canonical(&as_dsym(ds)), as_dsym(ds));
        assert!(ds.is_oriented());
        assert!(ds.is_complete());

        Self { ds }
    }
}


impl<'a> fmt::Display for DrawingInstructions<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let ds = self.ds;
        let mut seen= HashSet::new();
        let mut queue = VecDeque::from([(0, 1)]);

        let r = |d| ds.orbit([0, 1], d).len() / 2;

        while let Some((c, d)) = queue.pop_front() {
            if !seen.contains(&d) {
                write!(f, "{}-gon from {d}", r(d))?;
                if c != 0 {
                    write!(f, " connected to {c}")?;
                }
                writeln!(f, "")?;

                seen.extend(ds.orbit([0, 1], d));

                let mut e = d;
                loop {
                    queue.push_back((e, ds.op(2, e).unwrap()));
                    e = ds.walk(e, [0, 1]).unwrap();
                    if e == d {
                        break;
                    }
                }
            } else if d > c {
                writeln!(f, "connect {c} to {d}")?;
            }
        }
        writeln!(f, "")?;

        Ok(())
    }
}


#[derive(Clone)]
enum DSetOrEmpty {
    DSet(PartialDSet),
    Empty
}


fn as_dset<T: DSet>(ds: &T) -> PartialDSet {
    build_set(ds.size(), ds.dim(), |i, d| ds.op(i, d))
}


fn as_dsym<T: DSet>(ds: &T) -> PartialDSym {
    build_sym_using_vs(as_dset(ds), |_, _| Some(1))
}


fn dual(ds: &DSetOrEmpty) -> Option<DSetOrEmpty> {
    match ds {
        DSetOrEmpty::Empty => None,
        DSetOrEmpty::DSet(ds) => {
            let n = ds.dim();
            let ds_out = build_set(ds.size(), n, |i, d| ds.op(n - i, d));
            Some(DSetOrEmpty::DSet(ds_out))
        }
    }
}


fn r(ds: &PartialDSet, i: usize, j: usize, d: usize) -> usize {
    let mut e = d;
    for k in 1.. {
        e = ds.op(j, ds.op(i, e).unwrap()).unwrap();
        if e == d {
            return k;
        }
    }
    0
}


fn collapse<I>(ds: &DSetOrEmpty, remove: I, connector: usize)
    -> Option<DSetOrEmpty>
    where I: IntoIterator<Item=usize>
{
    match ds {
        DSetOrEmpty::Empty => None,
        DSetOrEmpty::DSet(ds) => {
            let remove: HashSet<_> = remove.into_iter().collect();

            if remove.is_empty() {
                None
            } else if remove.len() == ds.size() {
                Some(DSetOrEmpty::Empty)
            } else {
                let mut src2img = vec![0; ds.size() + 1];
                let mut img2src = vec![0; ds.size() + 1];
                let mut next = 1;
                for d in 1..=ds.size() {
                    if !remove.contains(&d) {
                        src2img[d] = next;
                        img2src[next] = d;
                        next += 1;
                    }
                }

                let op = |i, d| {
                    let mut e = ds.op(i, img2src[d]).unwrap();
                    if i != connector {
                        while src2img[e] == 0 {
                            e = ds.op(i, ds.op(connector, e).unwrap()).unwrap();
                        }
                    }
                    Some(src2img[e])
                };

                let ds_out = build_set(ds.size() - remove.len(), ds.dim(), op);
                Some(DSetOrEmpty::DSet(ds_out))
            }
        }
    }
}


fn reglue<I>(ds: &PartialDSet, pairs: I, index: usize) -> Option<PartialDSet>
    where I: IntoIterator<Item=(usize, usize)>
{
    let paired: HashMap<_, _> = pairs.into_iter()
        .flat_map(|(d, e)| [(d, e), (e, d)])
        .collect();

    if paired.is_empty() {
        None
    } else {
        let op = |i, d| {
            if i == index && paired.contains_key(&d) {
                Some(*paired.get(&d).unwrap())
            } else {
                ds.op(i, d)
            }
        };

        Some(build_set(ds.size(), ds.dim(), op))
    }
}


fn grow(ds: &PartialDSet, m: usize) -> PartialDSet {
    let n = ds.size();
    let op = |i, d| if d > n { Some(d) } else { ds.op(i, d) };
    build_set(n + m, ds.dim(), op)
}


fn cut_face(ds: &PartialDSet, d1: usize, d2: usize) -> PartialDSet {
    let n = ds.size();
    let mut ds = grow(ds, 8);
    let old = vec![
        d1, d2,
        ds.op(3, d2).unwrap(), ds.op(3, d1).unwrap(),
        ds.op(1, d1).unwrap(), ds.op(1, d2).unwrap(),
        ds.op(3, ds.op(1, d2).unwrap()).unwrap(),
        ds.op(3, ds.op(1, d1).unwrap()).unwrap(),
    ];
    let nu: Vec<_> = ((n + 1)..(n + 9)).collect();

    ds = reglue(
        &ds,
        [(nu[0], nu[1]), (nu[2], nu[3]), (nu[4], nu[5]), (nu[6], nu[7])],
        0
    ).unwrap();
    ds = reglue(
        &ds,
        (0..8).map(|k| (nu[k], old[k])),
        1
    ).unwrap();
    ds = reglue(
        &ds,
        [(nu[0], nu[4]), (nu[1], nu[5]), (nu[2], nu[6]), (nu[3], nu[7])],
        2
    ).unwrap();
    ds = reglue(
        &ds,
        [(nu[0], nu[3]), (nu[1], nu[2]), (nu[4], nu[7]), (nu[5], nu[6])],
        3
    ).unwrap();

    ds
}


fn cut_tile(ds: &PartialDSet, cut_chambers: &Vec<usize>) -> PartialDSet {
    let n = ds.size();
    let m = cut_chambers.len();
    assert!(m % 2 == 0);
    let opposites: Vec<_> = cut_chambers.iter()
        .map(|&d| ds.op(2, d).unwrap())
        .collect();
    let mut ds = grow(ds, 2 * m);

    let cycle_0 = |start: usize| ((0..(m / 2)))
        .map(move |i| (start + 2 * i, start + 2 * i + 1));

    let cycle_1 = |start: usize| ((0..(m / 2)))
        .map(move |i| (start + 2 * i + 1, start + (2 * i + 2) % m));

    ds = reglue(
        &ds,
        std::iter::empty().chain(cycle_0(n + 1)).chain(cycle_0(n + m + 1)),
        0
    ).unwrap();

    ds = reglue(
        &ds,
        std::iter::empty().chain(cycle_1(n + 1)).chain(cycle_1(n + m + 1)),
        1
    ).unwrap();

    ds = reglue(
        &ds,
        std::iter::empty()
            .chain((0..m).map(|i| (cut_chambers[i], n + 1 + i)))
            .chain((0..m).map(|i| (opposites[i], n + m + 1 + i))),
        2
    ).unwrap();

    ds = reglue(
        &ds,
        (0..m).map(|i| (n + 1 + i, n + m + 1 + i)),
        3
    ).unwrap();

    ds
}


fn squeeze_tile_3d(ds: &PartialDSet, d: usize, e: usize) -> PartialDSet {
    let f = ds.op(0, e).unwrap();
    let g = ds.op(0, d).unwrap();
    let f2 = ds.op(2, f).unwrap();
    let g2 = ds.op(2, g).unwrap();
    let d2 = ds.op(2, d).unwrap();
    let e2 = ds.op(2, e).unwrap();

    reglue(&ds, [(f, d), (g, e), (f2, d2), (g2, e2)], 2).unwrap()
}


fn merge_tiles(input: &DSetOrEmpty) -> Option<DSetOrEmpty> {
    match input {
        DSetOrEmpty::Empty => None,
        DSetOrEmpty::DSet(ds) => {
            let inner = inner_edges(&as_dsym(ds));
            let junk = inner.iter().cloned()
                .filter(|&(_, i)| i == 3)
                .flat_map(|(d, _)| ds.orbit([3], d));
            collapse(input, junk, 3)
        }
    }
}


fn merge_facets(input: &DSetOrEmpty) -> Option<DSetOrEmpty> {
    match input {
        DSetOrEmpty::Empty => None,
        DSetOrEmpty::DSet(ds) => {
            let reps = ds.orbit_reps([2, 3], 1..ds.size());
            let junk = reps.iter().cloned()
                .filter(|&d| r(ds, 2, 3, d) == 2)
                .flat_map(|d| ds.orbit([2, 3], d));
            collapse(input, junk, 2)
        }
    }
}


fn merge_all(ds: &DSetOrEmpty) -> Option<DSetOrEmpty> {
    let mut ds = ds.clone();

    for op in [
        merge_tiles, merge_facets, dual,
        merge_tiles, merge_facets, dual
    ] {
        if let Some(out) = op(&ds) {
            ds = out;
        }
    }

    Some(ds)
}


fn fix_local_1_vertex(input: &DSetOrEmpty) -> Option<DSetOrEmpty> {
    match input {
        DSetOrEmpty::Empty => None,
        DSetOrEmpty::DSet(ds) => {
            for c in ds.orbit_reps([1, 2], 1..ds.size()) {
                if ds.op(1, c) == ds.op(2, c) {
                    let d = ds.op(0, ds.op(1, c).unwrap()).unwrap();
                    let e = ds.op(1, ds.op(0, c).unwrap()).unwrap();
                    let f = ds.op(3, d).unwrap();
                    let g = ds.op(3, e).unwrap();

                    let d1 = ds.op(1, d).unwrap();
                    let e1 = ds.op(1, e).unwrap();
                    let f1 = ds.op(1, f).unwrap();
                    let g1 = ds.op(1, g).unwrap();

                    let tmp = reglue(
                        &ds, [(d, e1), (e, d1), (f, g1), (g, f1)], 1
                    ).unwrap();

                    let orb = tmp.orbit([0, 1, 3], c);
                    return collapse(&DSetOrEmpty::DSet(tmp), orb, 3);
                }
            }

            None
        }
    }
}


fn fix_local_2_vertex(input: &DSetOrEmpty) -> Option<DSetOrEmpty> {
    match input {
        DSetOrEmpty::Empty => None,
        DSetOrEmpty::DSet(ds) => {
            for d in ds.orbit_reps([1, 2], 1..ds.size()) {
                if r(&ds, 1, 2, d) == 2 {
                    let e = ds.op(3, ds.op(2, d).unwrap()).unwrap();
                    if
                        d == e ||
                        d == ds.op(1, ds.op(0, e).unwrap()).unwrap() ||
                        d == ds.op(0, ds.op(1, e).unwrap()).unwrap()
                    {
                        continue;
                    }

                    let mut ds = as_dset(ds);
                    let e = ds.op(2, ds.op(1, d).unwrap()).unwrap();

                    if r(&ds, 0, 1, d) > 3 {
                        ds = cut_face(
                            &ds,
                            ds.op(0, d).unwrap(),
                            ds.op(0, ds.op(1, d).unwrap()).unwrap()
                        );
                    }

                    if r(&ds, 0, 1, e) > 3 {
                        ds = cut_face(
                            &ds,
                            ds.op(0, e).unwrap(),
                            ds.op(0, ds.op(1, e).unwrap()).unwrap()
                        );
                    }

                    ds = squeeze_tile_3d(
                        &ds,
                        ds.op(1, ds.op(0, d).unwrap()).unwrap(),
                        ds.op(1, ds.op(0, e).unwrap()).unwrap(),
                    );

                    let orb = ds.orbit([0, 1, 3], d);
                    return collapse(&DSetOrEmpty::DSet(ds), orb, 3);
                }
            }

            None
        }
    }
}


fn fix_non_disk_face(input: &DSetOrEmpty) -> Option<DSetOrEmpty> {
    match input {
        DSetOrEmpty::Empty => None,
        DSetOrEmpty::DSet(ds) => {
            let mut face_rep = vec![0; ds.size() + 1];
            for d in ds.orbit_reps([0, 1], 1..=ds.size()) {
                for e in ds.orbit([0, 1], d) {
                    face_rep[e] = d;
                }
            }
            let face_rep = face_rep;

            for d in ds.orbit_reps([1, 2], 1..=ds.size()) {
                let mut e = ds.op(1, ds.op(2, d).unwrap()).unwrap();
                while e != d {
                    if face_rep[e] == face_rep[d] {
                        let f = ds.op(3, d).unwrap();
                        let g = ds.op(3, e).unwrap();

                        let d1 = ds.op(1, d).unwrap();
                        let e1 = ds.op(1, e).unwrap();
                        let f1 = ds.op(1, f).unwrap();
                        let g1 = ds.op(1, g).unwrap();

                        return Some(DSetOrEmpty::DSet(
                            reglue(
                                &ds, [(d, e1), (e, d1), (f, g1), (g, f1)], 1
                            ).unwrap()
                        ));
                    } else {
                        e = ds.op(1, ds.op(2, e).unwrap()).unwrap();
                    }
                }
            }

            None
        }
    }
}


fn split_and_glue(input: &DSetOrEmpty) -> Option<DSetOrEmpty> {
    match input {
        DSetOrEmpty::Empty => None,
        DSetOrEmpty::DSet(ds) => {
            let ds = as_dset(&canonical(&as_dsym(ds)));

            if let Some(cut) = small_tile_cut(&ds) {
                let (glue_chamber, cut_vertex_reps, inside_vertex_reps) = cut;
                let special_chamber = ds.op(3, glue_chamber).unwrap();
                let marked_vertex_reps = std::iter::empty()
                    .chain(ds.orbit([0, 1], glue_chamber))
                    .chain(cut_vertex_reps.iter().cloned())
                    .chain(inside_vertex_reps.iter().cloned())
                    .collect();
                let ordered = ordered_cut(
                    special_chamber, &marked_vertex_reps, &ds
                );
                assert_eq!(
                    ordered.len(), cut_vertex_reps.len(),
                    "split_and_glue(): got {:?} from {:?} in\n\n{}",
                    ordered, cut_vertex_reps, DrawingInstructions::new(&ds)
                );
                split_and_glue_attempt(&ds, glue_chamber, ordered)
            } else {
                None
            }
        }
    }
}


fn small_tile_cut(ds: &PartialDSet) -> Option<(usize, Vec<usize>, Vec<usize>)> {
    let (elm_to_index, reps, edges) = make_skeleton(ds);
    let source = reps.len();
    let sink = source + 1;

    let mut best = None;

    for d in ds.orbit_reps([0, 1, 3], 1..=ds.size()) {
        let d3 = ds.op(3, d).unwrap();
        let v_in: HashSet<_> = ds.orbit([0, 1], d).iter()
            .map(|&e| elm_to_index[e])
            .collect();
        let v_out: HashSet<_> = ds.orbit([0, 1], d3).iter()
            .map(|&e| elm_to_index[e])
            .collect();

        let edges = edges.iter().cloned()
            .chain(v_in.iter().map(|&v| (source, v)))
            .chain(v_out.iter().map(|&v| (v, sink)));

        let cut_raw = min_vertex_cut_undirected(edges, source, sink);
        let (n, m) = (v_in.len(), cut_raw.cut_vertices.len());

        if m < n {
            let cut_vertices: Vec<_> = cut_raw.cut_vertices.iter()
                .map(|&v| reps[v])
                .collect();
            let inside_vertices: Vec<_> = cut_raw.inside_vertices.iter()
                .filter(|&&v| v < reps.len())
                .map(|&v| reps[v])
                .collect();
            let key = (m, -(n as isize));
            if let Some((best_key, _, _, _)) = best {
                if key < best_key {
                    best = Some((key, d, cut_vertices, inside_vertices));
                }
            } else {
                best = Some((key, d, cut_vertices, inside_vertices));
            }
        }
    }

    best.and_then(|(_, d, cut_vertices, inside_vertices)|
        Some((d, cut_vertices, inside_vertices))
    )
}


fn ordered_cut(
    special_chamber: usize,
    marked_vertex_reps: &Vec<usize>,
    ds: &PartialDSet
) -> Vec<(usize, usize)>
{
    let special: HashSet<_> = ds.orbit([0, 1], special_chamber).iter()
        .cloned()
        .collect();
    let marked: HashSet<_> = marked_vertex_reps.iter()
        .flat_map(|&e| ds.orbit([1, 2], e))
        .collect();
    let start = *marked.iter()
        .find(|&&e| !marked.contains(&ds.op(0, e).unwrap()))
        .unwrap();

    let mut result = vec![];
    let mut d = start;

    while result.len() < ds.size() + 1 { // avoid looping forever in error case
        let mut e = ds.op(1, d).unwrap();
        while marked.contains(&ds.op(0, e).unwrap()) {
            e = ds.walk(e, [0, 1]).unwrap();
        }

        if special.contains(&d) {
            let mut d = d;
            while ds.op(1, d) != Some(e) {
                result.push((d, ds.walk(d, [1, 0, 1]).unwrap()));
                d = ds.walk(d, [1, 0]).unwrap();
            }
        } else if ds.op(1, d) != Some(e) {
            result.push((d, e));
        }

        d = ds.op(2, e).unwrap();

        if d == start {
            break;
        }
    }

    result
}


fn split_and_glue_attempt(
    ds: &PartialDSet, glue_chamber: usize, ordered: Vec<(usize, usize)>
) -> Option<DSetOrEmpty>
{
    eprintln!(
        "# INFO: attempting {}-split with {}-glue",
        ordered.len(), ds.orbit([0, 1], glue_chamber).len() / 2
    );

    let mut ds = as_dset(ds);
    let mut cut_chambers = vec![];

    for (d, e) in ordered {
        if ds.walk(d, [1, 0, 1]) != Some(e) {
            if ds.orbit([0, 1], d).contains(&e) {
                ds = cut_face(&ds, d, e);
            } else {
                eprintln!("# INFO: invalid split-and-glue attempt");
                return None;
            }
        }
        cut_chambers.push(ds.op(1, d).unwrap());
        cut_chambers.push(ds.op(1, e).unwrap());
    }
    eprintln!("# INFO: valid split-and-glue attempt");

    ds = cut_tile(&ds, &cut_chambers);

    let junk = ds.orbit([0, 1, 3], glue_chamber);
    collapse(&DSetOrEmpty::DSet(ds), junk, 3)
}


fn make_skeleton(ds: &PartialDSet)
    -> (Vec<usize>, Vec<usize>, Vec<(usize, usize)>)
{
    let reps = ds.orbit_reps([1, 2], 1..=ds.size());

    let mut elm_to_index = vec![0; ds.size() + 1];
    for (i, &d) in reps.iter().enumerate() {
        for e in ds.orbit([1, 2], d) {
            elm_to_index[e] = i;
        }
    }

    let edges: Vec<_> = ds.orbit_reps([0, 2], 1..=ds.size()).iter()
        .map(|&d| (elm_to_index[d], elm_to_index[ds.op(0, d).unwrap()]))
        .map(|(d, e)| (d.min(e), d.max(e)))
        .collect::<BTreeSet<_>>().iter().cloned()
        .collect();

    (elm_to_index, reps, edges)
}


pub fn simplify<T: DSet>(ds: &T) -> Option<PartialDSym> {
    // TODO add assertions to ensure input is legal

    let mut ds = DSetOrEmpty::DSet(as_dset(ds));
    ds = merge_all(&ds).or(Some(ds)).unwrap();

    loop {
        let mut changed = false;
        for op in [
            fix_local_1_vertex,
            fix_local_2_vertex,
            fix_non_disk_face,
            split_and_glue,
        ] {
            if let Some(out) = op(&ds) {
                ds = merge_all(&out).or(Some(out)).unwrap();
                changed = true;
                break;
            }
        }
        if !changed {
            return match ds {
                DSetOrEmpty::Empty => None,
                DSetOrEmpty::DSet(ds) => Some(as_dsym(&ds))
            }
        }
    }
}


#[cfg(test)]
mod test {
    use crate::covers::finite_universal_cover;
    use crate::delaney3d::pseudo_toroidal_cover;
    use crate::derived::{canonical, minimal_image};
    use crate::fpgroups::invariants::abelian_invariants;
    use crate::fundamental_group::fundamental_group;

    use super::*;

    fn dsym(s: &str) -> PartialDSym {
        s.parse::<PartialDSym>().unwrap()
    }


    #[test]
    fn test_make_skeleton() {
        let ds = finite_universal_cover(&dsym("<1.1:1:1,1,1:3,3>"));
        let (elm_to_index, reps, edges) = make_skeleton(&as_dset(&ds));
        assert_eq!(elm_to_index.len(), 25);
        assert_eq!(reps.len(), 4);
        assert_eq!(edges.len(), 6);
        assert_eq!(edges, vec![(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]);
    }

    #[test]
    fn test_simplify() {
        let good = HashSet::from([
            dsym("<1.1:1 3:1,1,1,1:4,3,4>").to_string(),
            dsym(
                "<1.1:8 3:2 4 6 8,6 3 5 7 8,2 7 8 5 6,4 3 6 8:3 4,5 3,3 4>"
            ).to_string()
        ]);

        let test = |s: &str| {
            let ds = s.parse::<PartialDSym>().unwrap();
            let cov = pseudo_toroidal_cover(&ds).unwrap();
            let out = simplify(&cov).unwrap();

            assert_eq!(out.orbit_reps([0, 1, 2], 1..out.size()).len(), 1);
            assert_eq!(out.orbit_reps([1, 2, 3], 1..out.size()).len(), 1);

            let reps = out.orbit_reps([2, 3], 1..out.size());
            assert!(reps.iter().all(|&d| out.r(2, 3, d) != Some(2)));
            let reps = out.orbit_reps([0, 1], 1..out.size());
            assert!(reps.iter().all(|&d| out.r(0, 1, d) != Some(2)));
            let reps = out.orbit_reps([1, 2], 1..out.size());
            assert!(reps.iter().all(|&d| out.r(1, 2, d) != Some(2)));

            let fg = fundamental_group(&out);
            let nr_gens = fg.gen_to_edge.len();
            let inv = abelian_invariants(nr_gens, fg.relators.clone());
            assert_eq!(inv, vec![0, 0, 0]);

            let tmp = canonical(&minimal_image(&out)).to_string();
            assert!(good.contains(&tmp));
        };

        test("<1.4:1 3:1,1,1,1:4,3,4>");
        test("<2.1:2 3:1 2,1 2,1 2,2:3 3,3 4,4>");
        test("<513.5:2 3:2,1 2,1 2,2:4,2 4,6>");
        test("<513.8:2 3:2,1 2,1 2,2:6,2 3,6>");
        test("<3.3:3 3:1 2 3,1 2 3,1 3,2 3:3 3 4,4 4,3>");
        test("<167.3:3 3:1 2 3,1 3,2 3,1 2 3:3 4,3,4 6>");
        test("<184.4:3 3:1 2 3,1 3,2 3,1 3:4 6,3,3>");
        test("<23.14:4 3:1 2 3 4,1 2 4,1 3 4,2 3 4:3 3 8,4 3,3 4>");
        test("<71.3:4 3:1 2 3 4,1 2 4,1 3 4,2 4:3 3 6,3 3,4>");
        test("<514.7:4 3:2 4,1 2 3 4,1 2 3 4,3 4:4 4,2 4 4 3,4 4>");
        test("<553.3:4 3:2 4,1 2 3 4,3 4,2 4:4 6,2 6,4>");
        test("<45.2:5 3:1 2 3 5,1 2 4 5,1 3 4 5,2 3 4 5:3 3 3,3 3 3,6 4 4>");
        test("<45.7:5 3:1 2 3 5,1 2 4 5,1 3 4 5,2 3 4 5:3 3 3,4 3 3,6 3 3>");
        test("<45.12:5 3:1 2 3 5,1 2 4 5,1 3 4 5,2 3 4 5:3 3 6,4 3 3,3 4 4>");
        test("<54.2:5 3:1 2 3 5,1 2 4 5,1 3 5,2 3 4 5:3 3 3,3 4,3 6>");
        test("<54.4:5 3:1 2 3 5,1 2 4 5,1 3 5,2 3 4 5:3 3 3,4 4,3 4>");
        test("<222.77:5 3:1 2 4 5,1 3 5,2 3 4 5,1 5 4:4 12,3 2,3 4>");
    }
}
