use crate::backtrack::BackTrackIterator;
use crate::backtrack::BackTracking;
use crate::dsets::*;


struct DSetGenState {
    dset: PartialDSet,
    is_remap_start: Vec<bool>,
}


struct DSetBackTracking {
    dim: usize,
    max_size: usize,
}

impl BackTracking for DSetBackTracking {
    type State = DSetGenState;
    type Item = SimpleDSet;

    fn root(&self) -> Self::State {
        DSetGenState {
            dset: PartialDSet::new(1, self.dim),
            is_remap_start: vec![false; self.max_size + 1],
        }
    }

    fn extract(&self, state: &Self::State) -> Option<Self::Item> {
        if state.dset.is_complete() {
            Some(SimpleDSet::from(&state.dset))
        } else {
            None
        }
    }

    fn children(&self, state: &Self::State) -> Vec<Self::State> {
        let mut result = vec![];

        let ds = &state.dset;

        if let Some((i, d)) = first_undefined(ds) {
            for e in d..=self.max_size.min(ds.size() + 1) {
                if ds.get(i, e) == None {
                    let mut dset = ds.clone();
                    let mut is_remap_start = state.is_remap_start.clone();

                    if e > dset.size() {
                        dset.grow(1);
                        is_remap_start[e] = true;
                    }

                    dset.set(i, d, e);

                    let (head, tail, gap, k) = scan_02_orbit(&dset, d);

                    if gap == 1 {
                        dset.set(k, head, tail);
                    } else if gap == 0 && head != tail {
                        continue;
                    }

                    if check_canonicity(&dset, &mut is_remap_start) {
                        result.push(Self::State { dset, is_remap_start });
                    }
                }
            }
        }

        result
    }
}


pub struct DSets(BackTrackIterator<DSetBackTracking>);

impl DSets {
    pub fn new(dim: usize, max_size: usize) -> DSets {
        DSets(BackTrackIterator::new(DSetBackTracking { dim, max_size }))
    }
}

impl Iterator for DSets {
    type Item = SimpleDSet;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next()
    }
}


fn first_undefined(ds: &PartialDSet) -> Option<(usize, usize)> {
    for d in 1..=ds.size() {
        for i in 0..=ds.dim() {
            if ds.get(i, d) == None {
                return Some((i, d));
            }
        }
    }
    None
}


fn scan_02_orbit(ds: &PartialDSet, d: usize) -> (usize, usize, usize, usize) {
    let (head, i) = scan(ds, [0, 2, 0, 2], d, 4);
    let (tail, j) = scan(ds, [2, 0, 2, 0], d, 4 - i);

    (head, tail, 4 - i - j, 2 * (i % 2))
}


fn scan(
    ds: &PartialDSet, w: [usize; 4], d: usize, limit: usize
) -> (usize, usize)
{
    let mut e = d;

    for k in 0..limit {
        if let Some(en) = ds.get(w[k], e) {
            e = en;
        } else {
            return (e, k);
        }
    }

    (e, limit)
}


fn check_canonicity(
    ds: &PartialDSet, is_remap_start: &mut [bool]
) -> bool
{
    let mut n2o = vec![];
    let mut o2n = vec![];

    n2o.resize(ds.size() + 1, 0);
    o2n.resize(ds.size() + 1, 0);

    for d in 1..=ds.size() {
        if is_remap_start[d] {
            let diff = compare_renumbered_from(ds, d, &mut n2o, &mut o2n);

            if diff < 0 {
                return false;
            } else if diff > 0 {
                is_remap_start[d] = false;
            }
        }
    }

    true
}


fn compare_renumbered_from(
    ds: &PartialDSet, d0: usize, n2o: &mut [usize], o2n: &mut [usize]
) -> i64
{
    n2o.fill(0);
    o2n.fill(0);

    n2o[1] = d0;
    o2n[d0] = 1;

    let mut next = 2;

    for d in 1..=ds.size() {
        for i in 0..=ds.dim() {
            let ei = ds.get(i, n2o[d]).unwrap_or(0);

            if ei == 0 {
                return 0;
            } else {
                if o2n[ei] == 0 {
                    o2n[ei] = next;
                    n2o[next] = ei;
                    next += 1;
                }

                let di = ds.get(i, d).unwrap_or(0);

                if di == 0 {
                    return 0;
                } else if o2n[ei] != di {
                    return (o2n[ei] as i64) - (di as i64);
                }
            }
        }
    }

    0
}
