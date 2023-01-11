use super::backtrack::BackTrackIterator;
use super::backtrack::BackTracking;
use super::dsets::*;


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

                    let (head, tail, gap, k) = scan02Orbit(&dset, d);

                    if gap == 1 {
                        dset.set(k, head, tail);
                    } else if gap == 0 && head != tail {
                        continue;
                    }

                    if check_canonicity(&dset, &is_remap_start) {
                        result.push(Self::State { dset, is_remap_start });
                    }
                }
            }
        }

        result
    }
}


fn first_undefined(ds: &PartialDSet) -> Option<(usize, usize)> {
    for i in 0..=ds.dim() {
        for d in 1..=ds.size() {
            if ds.get(i, d) == None {
                return Some((i, d));
            }
        }
    }
    None
}


fn scan02Orbit(ds: &PartialDSet, d: usize) -> (usize, usize, usize, usize) {
    todo!()
}


fn check_canonicity(ds: &PartialDSet, is_remap_start: &Vec<bool>) -> bool {
    todo!()
}
