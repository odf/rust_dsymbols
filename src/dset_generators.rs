use crate::backtrack::BackTrackIterator;
use crate::backtrack::BackTracking;
use crate::dsets::*;


struct DSetGenState {
    dset: PartialDSet,
    is_remap_start: Vec<bool>,
    next_i_d: Option<(usize, usize)>,
}


struct DSetBackTracking {
    dim: usize,
    max_size: usize,
}

impl BackTracking for DSetBackTracking {
    type State = DSetGenState;
    type Item = PartialDSet;

    fn root(&self) -> Self::State {
        DSetGenState {
            dset: PartialDSet::new(1, self.dim),
            is_remap_start: vec![false; self.max_size + 1],
            next_i_d: Some((0, 1)),
        }
    }

    fn extract(&self, state: &Self::State) -> Option<Self::Item> {
        if state.next_i_d.is_none() {
            Some(state.dset.clone())
        } else {
            None
        }
    }

    fn children(&self, state: &Self::State) -> Vec<Self::State> {
        let mut new2old = vec![0; self.max_size + 1];
        let mut old2new = vec![0; self.max_size + 1];

        let mut result = vec![];

        if let Some((i, d)) = state.next_i_d {
            let ds = &state.dset;
            let max_e = (ds.size() + 1).min(self.max_size);

            for e in d..=max_e {
                if ds.get(i, e) == None {
                    let mut dset = ds.clone();
                    let mut is_remap_start = state.is_remap_start.clone();

                    if e > dset.size() {
                        dset.grow(1);
                        is_remap_start[e] = true;
                    }

                    dset.set(i, d, e);

                    if !check_and_apply_implications(&mut dset, i, d) {
                        continue;
                    }

                    if check_canonicity(
                        &dset, &mut is_remap_start, &mut new2old, &mut old2new
                    ) {
                        let next_i_d = next_undefined(&dset, i, d);
                        result.push(
                            Self::State { dset, is_remap_start, next_i_d }
                        );
                    }
                }
            }
        }

        result
    }
}


pub struct DSets {
    bt: BackTrackIterator<DSetBackTracking>,
    counter: usize,
}

impl DSets {
    pub fn new(dim: usize, max_size: usize) -> DSets {
        DSets {
            bt: BackTrackIterator::new(
                DSetBackTracking { dim, max_size }
            ),
            counter: 0,
        }
    }
}

impl Iterator for DSets {
    type Item = SimpleDSet;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(ds) = self.bt.next() {
            self.counter += 1;
            Some(SimpleDSet::from_partial_unchecked(ds, self.counter))
        } else {
            None
        }
    }
}


fn next_undefined(ds: &PartialDSet, i0: usize, d0: usize)
    -> Option<(usize, usize)>
{
    for i in i0..=ds.dim() {
        if ds.get(i, d0) == None {
            return Some((i, d0));
        }
    }

    for d in (d0 + 1)..=ds.size() {
        for i in 0..=ds.dim() {
            if ds.get(i, d) == None {
                return Some((i, d));
            }
        }
    }
    None
}


fn check_and_apply_implications(
    dset: &mut PartialDSet, i: usize, d: usize
) -> bool {
    // TODO only works for 2d symbols right now

    let (head, tail, gap, k) = scan_02_orbit(dset, d);

    if i == 1 {
        true
    } else if gap == 0 && head != tail {
        false
    } else {
        if gap == 1 {
            dset.set(k, head, tail);
        }
        true
    }
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
    ds: &PartialDSet,
    is_remap_start: &mut [bool],
    n2o: &mut [usize],
    o2n: &mut [usize]
) -> bool
{
    for d in 1..=ds.size() {
        if is_remap_start[d] {
            let diff = compare_renumbered_from(ds, d, n2o, o2n);

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
    ds: &PartialDSet, d0: usize, new2old: &mut [usize], old2new: &mut [usize]
) -> i64
{
    new2old.fill(0);
    old2new.fill(0);

    new2old[1] = d0;
    old2new[d0] = 1;

    let mut next = 2;

    for d in 1..=ds.size() {
        for i in 0..=ds.dim() {
            let ei = ds.get(i, new2old[d]).unwrap_or(0);

            if ei == 0 {
                return 0;
            } else {
                if old2new[ei] == 0 {
                    old2new[ei] = next;
                    new2old[next] = ei;
                    next += 1;
                }

                let di = ds.get(i, d).unwrap_or(0);

                if di == 0 {
                    return 0;
                } else if old2new[ei] != di {
                    return (old2new[ei] as i64) - (di as i64);
                }
            }
        }
    }

    0
}
