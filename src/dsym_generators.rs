use num_rational::Rational64;

use crate::backtrack::BackTrackIterator;
use crate::backtrack::BackTracking;
use crate::dsets::*;
use crate::dsyms::*;


struct DSymGenState {
    vs: Vec<usize>,
    curv: Rational64,
    next: usize,
}


struct DSymBackTracking {
    dset: SimpleDSet,
    orbit_index: Vec<Vec<usize>>,
    orbit_rs: Vec<usize>,
    orbit_vmins: Vec<usize>,
    orbit_is_chain: Vec<bool>,
}

impl DSymBackTracking {
    fn new(dset: &SimpleDSet) -> DSymBackTracking {
        let (orbit_rs, orbit_index) = collect_orbits(&dset);

        let mut orbit_vmins = vec![0; orbit_rs.len()];
        for i in 0..orbit_vmins.len() {
            orbit_vmins[i] = match orbit_rs[i] {
                1 => 3,
                2 => 2,
                _ => 1,
            }
        }

        let mut orbit_is_chain = vec![false; orbit_rs.len()];
        for i in 0..dset.dim() {
            for d in 1..=dset.size() {
                if dset.get(i, d) == Some(d) || dset.get(i + 1, d) == Some(d) {
                    orbit_is_chain[orbit_index[i][d]] = true;
                }
            }
        }

        DSymBackTracking {
            dset: dset.clone(),
            orbit_index,
            orbit_rs,
            orbit_vmins,
            orbit_is_chain,
        }
    }
}

impl BackTracking for DSymBackTracking {
    type State = DSymGenState;
    type Item = PartialDSym;

    fn root(&self) -> Self::State {
        todo!()
    }

    fn extract(&self, state: &Self::State) -> Option<Self::Item> {
        todo!()
    }

    fn children(&self, state: &Self::State) -> Vec<Self::State> {
        todo!()
    }
}


pub struct DSyms {
    bt: BackTrackIterator<DSymBackTracking>,
    counter: usize,
}

impl DSyms {
    pub fn new(dset: &SimpleDSet) -> DSyms {
        DSyms {
            bt: BackTrackIterator::new(DSymBackTracking::new(dset)),
            counter: 0,
        }
    }
}

impl Iterator for DSyms {
    type Item = SimpleDSym;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(ds) = self.bt.next() {
            self.counter += 1;
            Some(SimpleDSym::from_partial(ds, self.counter))
        } else {
            None
        }
    }
}
