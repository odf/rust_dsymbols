use num_rational::Rational64;
use num_traits::Signed;

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
    orbit_maps: Option<Vec<Vec<usize>>>,
    base_curvature: Rational64,
    canonicity_test_needed: bool,
}

impl DSymBackTracking {
    fn new(dset: &SimpleDSet) -> DSymBackTracking {
        let (orbit_rs, orbit_is_chain, orbit_index) = collect_orbits(&dset);
        let orbit_vmins = compute_vmins(&orbit_rs);

        let mut curv_sixfold = -3 * dset.size() as i64;
        for i in 0..orbit_vmins.len() {
            let k = if orbit_is_chain[i] { 1 } else { 2 };
            curv_sixfold += 6 * k / orbit_vmins[i] as i64;
        }

        let base_curvature = Rational64::new(curv_sixfold, 6);
        let canonicity_test_needed = !base_curvature.is_negative();

        let orbit_maps = if canonicity_test_needed {
            Some(orbit_maps(&dset, orbit_vmins.len(), &orbit_index))
        } else {
            None
        };

        DSymBackTracking {
            dset: dset.clone(),
            orbit_index,
            orbit_rs,
            orbit_vmins,
            orbit_is_chain,
            orbit_maps,
            base_curvature,
            canonicity_test_needed,
        }
    }

    fn orbit_count(&self) -> usize {
        self.orbit_vmins.len()
    }

    fn is_minimally_hyperbolic(&self, vs: &[usize], curv: Rational64) -> bool {
        if !curv.is_negative() {
            false
        } else {
            for i in 0..self.orbit_count() {
                if vs[i] > self.orbit_vmins[i] {
                    let v = vs[i] as i64;
                    let k = if self.orbit_is_chain[i] { 1 } else { 2 };
                    let c = curv
                        - Rational64::new(k, v)
                        + Rational64::new(k, v - 1);
                    if c.is_negative() {
                        return false;
                    }
                }
            }
            true
        }
    }

    fn is_canonical(&self, vs: &[usize]) -> bool {
        if self.canonicity_test_needed {
            let maps = self.orbit_maps.as_ref().unwrap();

            for m in maps {
                let ws: Vec<_> = (0..vs.len()).map(|i| vs[m[i]]).collect();
                if &ws[..] > vs {
                    return false;
                }
            }
        }

        true
    }

    fn is_good(&self, vs: &[usize], curv: Rational64) -> bool {
        if !curv.is_positive() {
            true
        } else {
            let ds = &self.dset;
            let mut cones = vec![];
            let mut corners = vec![];

            for d in ds.orbit_reps_2d(0, 2) {
                let d0 = ds.get(0, d).unwrap();
                let d2 = ds.get(2, d).unwrap();

                if d0 == d && d2 == d {
                    corners.push(2)
                } else if d0 != d && d2 == d0 {
                    cones.push(2)
                }
            }

            for i in 0..self.orbit_count() {
                if vs[i] > 1 {
                    if self.orbit_is_chain[i] {
                        corners.push(vs[i]);
                    } else {
                        cones.push(vs[i]);
                    }
                }
            }

            cones.sort();
            cones.reverse();
            corners.sort();
            corners.reverse();

            let front = cones.iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join("");
            let middle = String::from(
                if self.dset.is_loopless() { "" } else { "*" }
            );
            let back = corners.iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join("");
            let cross = String::from(
                if self.dset.is_weakly_oriented() { "" } else { "x" }
            );
            let key = [ front, middle, back, cross ].join("");

            [
                "", "*", "x",
                "532", "432", "332",
                "422", "322", "222",
                "44", "33", "22",
                "*532", "*432", "*332", "3*2",
                "*422", "*322", "*222", "2*4", "2*3", "2*2",
                "*44", "*33", "*22", "4*", "3*", "2*", "4x", "3x", "2x"
            ].contains(&&key[..])
        }
    }

    fn make_dsym(&self, vs: &[usize]) -> PartialDSym {
        PartialDSym::from(
            self.dset.clone(),
            self.orbit_index.clone(),
            self.orbit_rs.clone(),
            vs.iter().cloned().collect()
        )
    }
}


fn compute_vmins(orbit_rs: &[usize]) -> Vec<usize> {
    let mut orbit_vmins = vec![0; orbit_rs.len()];

    for i in 0..orbit_vmins.len() {
        orbit_vmins[i] = match orbit_rs[i] {
            1 => 3,
            2 => 2,
            _ => 1,
        }
    }

    orbit_vmins
}


fn orbit_maps(
    dset: &SimpleDSet, orbit_count: usize, orbit_index: &Vec<Vec<usize>>
) -> Vec<Vec<usize>> {
    let mut orbit_maps =vec![];

    for map in dset.automorphisms() {
        let mut m = vec![0; orbit_count];
        for i in 0..dset.dim() {
            for d in 1..=dset.size() {
                m[orbit_index[i][d]] = orbit_index[i][map[d]];
            }
        }
        orbit_maps.push(m);
    }

    orbit_maps
}


impl BackTracking for DSymBackTracking {
    type State = DSymGenState;
    type Item = PartialDSym;

    fn root(&self) -> Self::State {
        let vs = self.orbit_vmins.clone();
        let curv = self.base_curvature;
        let next = if curv.is_negative() { self.orbit_count() } else { 0 };

        Self::State { vs, curv, next }
    }

    fn extract(&self, state: &Self::State) -> Option<Self::Item> {
        if state.next >= self.orbit_count()
            && self.is_good(&state.vs, state.curv)
            && self.is_canonical(&state.vs)
        {
            Some(self.make_dsym(&state.vs))
        } else {
            None
        }
    }

    fn children(&self, state: &Self::State) -> Vec<Self::State> {
        let mut result = vec![];

        if state.next < self.orbit_count() {
            let n = state.next;
            let vmin = state.vs[n];

            for v in vmin..=7 {
                let mut vs = state.vs.clone();
                vs[n] = v;

                let k = if self.orbit_is_chain[n] { 1 } else { 2 };
                let curv = state.curv
                    - Rational64::new(k, vmin as i64)
                    + Rational64::new(k, v as i64);

                if curv.is_negative() {
                    if self.is_minimally_hyperbolic(&vs, curv) {
                        let next = self.orbit_count();
                        result.push(Self::State { vs, curv, next });
                    }
                    break;
                } else {
                    let next = state.next + 1;
                    result.push(Self::State { vs, curv, next });
                }
            }
        }

        result
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
