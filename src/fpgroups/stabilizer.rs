use std::collections::{HashMap, HashSet, VecDeque};

use crate::fpgroups::free_words::relator_representative;

use super::free_words::{relator_permutations, FreeWord};


fn relators_by_start_gen(rels: &Vec<FreeWord>)
    -> HashMap<isize, Vec<FreeWord>>
{
    let mut result = HashMap::new();

    for rel in rels {
        for w in relator_permutations(&rel) {
            result.entry(w[0])
                .and_modify(|v: &mut Vec<_>| v.push(w.clone()))
                .or_insert(vec![w]);
        }
    }

    result
}


fn trace_word<F>(
    point: usize,
    w: &FreeWord,
    edge_to_word: &HashMap<(usize, isize), FreeWord>,
    action: F
)
    -> FreeWord
    where F: Fn(usize, isize) -> usize
{
    let mut p = point;
    let mut result = FreeWord::empty();

    for &g in w.iter() {
        result *= &edge_to_word[&(p, g)];
        p = action(p, g);
    }

    result
}


fn close_relations_in_place<F>(
    edge_to_word: &mut HashMap<(usize, isize), FreeWord>,
    start_edge: (usize, isize),
    wd: &FreeWord,
    rels_by_gen: &HashMap<isize, Vec<FreeWord>>,
    action: F,
)
    where F: Fn(usize, isize) -> usize
{
    let (p, g) = start_edge;
    let mut queue = VecDeque::from([(p, g, wd.clone())]);

    while let Some((point, gen, w)) = queue.pop_front() {
        edge_to_word.insert((action(point, gen), -gen), w.inverse());
        edge_to_word.insert((point, gen), w);

        for r in rels_by_gen[&gen].iter() {
            let mut cuts = vec![];
            let mut x = point;

            for i in 0..r.len() {
                let h = r[i];
                if !edge_to_word.contains_key(&(x, h)) {
                    let w = r.rotated(i as isize + 1) * -h;
                    cuts.push((x, h, w));
                }
                x = action(x, h);
            }

            if cuts.len() == 1 {
                let (p, g, w) = cuts[0].clone();
                let w = trace_word(p, &w, &edge_to_word, &action);
                queue.push_back((p, g, w));
            }
        }
    }
}


fn spanning_tree<F>(base_point: usize, nr_gens: usize, action: F)
    -> Vec<(usize, isize)>
    where F: Fn(usize, isize) -> usize
{
    let mut edges = vec![];
    let mut queue = VecDeque::from([base_point]);
    let mut seen = HashSet::from([base_point]);

    while let Some(point) = queue.pop_front() {
        for i in 1..=nr_gens as isize {
            for gen in [i, -i] {
                let p = action(point, gen);
                if !seen.contains(&p) {
                    queue.push_back(p);
                    seen.insert(p);
                    edges.push((point, gen));
                }
            }
        }
    }
    edges
}


pub fn stabilizer<F>(
    base_point: usize,
    nr_points: usize,
    nr_gens: usize,
    rels: &Vec<FreeWord>,
    action: F,
)
    -> (Vec<FreeWord> , Vec<FreeWord>)
    where F: Fn(usize, isize) -> usize
{
    let rels_by_gen = relators_by_start_gen(rels);
    let tree = spanning_tree(base_point, nr_gens, &action);

    let mut point_to_word = HashMap::from([(base_point, FreeWord::empty())]);
    let mut edge_to_word = HashMap::new();

    for edge in tree {
        close_relations_in_place(
            &mut edge_to_word, edge, &FreeWord::empty(), &rels_by_gen, &action
        );
        let (pt, gen) = edge;
        point_to_word.insert(action(pt, gen), &point_to_word[&pt] * gen);
    }

    let mut generators = vec![];

    for px in 0..nr_points {
        let wx = &point_to_word[&px];

        for i in 1..=nr_gens as isize {
            for g in [i, -i] {
                let edge = (px, g);
                if edge_to_word.get(&edge).is_none() {
                    let wy = &point_to_word[&action(px, g)];
                    generators.push(wx * g * wy.inverse());

                    let w = FreeWord::from([generators.len() as isize]);
                    close_relations_in_place(
                        &mut edge_to_word, edge, &w, &rels_by_gen, &action
                    )
                }
            }
        }
    }

    let mut subrels = vec![];
    let mut seen = HashSet::new();

    for p in 0..nr_points {
        for r in rels {
            let w = relator_representative(
                &trace_word(p, &r, &edge_to_word, &action)
            );
            if w.len() > 0 && !seen.contains(&w) {
                seen.insert(w.clone());
                subrels.push(w);
            }
        }
    }
    subrels.sort();
    subrels.reverse();

    (generators, subrels)
}
