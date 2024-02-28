use std::collections::{HashMap, VecDeque};

use super::free_words::{relator_permutations, FreeWord};


fn relators_by_start_gen(rels: Vec<FreeWord>)
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
    w: FreeWord,
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
    wd: FreeWord,
    rels_by_gen: &HashMap<isize, Vec<FreeWord>>,
    action: F,
)
    where F: Fn(usize, isize) -> usize
{
    let (p, g) = start_edge;
    let mut queue = VecDeque::from([(p, g, wd)]);

    while let Some((point, gen, w)) = queue.pop_front() {
        edge_to_word.insert((action(point, gen), -gen), w.inverse());
        edge_to_word.insert((point, gen), w);

        for r in rels_by_gen[&gen].iter() {
            let mut cuts = vec![];
            let mut x = point;

            for i in 0..r.len() {
                let h = r[i];
                if !edge_to_word.contains_key(&(x, h)) {
                    let w = FreeWord::from([-h]) * r.rotated(i as isize);
                    cuts.push((x, h, w));
                }
                x = action(x, h);
            }

            if cuts.len() == 1 {
                let (p, g, w) = cuts[0].clone();
                let w = trace_word(p, w, &edge_to_word, &action);
                queue.push_back((p, g, w));
            }
        }
    }
}
