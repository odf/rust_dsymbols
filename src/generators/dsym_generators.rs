use std::collections::VecDeque;
use std::str::FromStr;

use crate::util::backtrack::BackTrackIterator;
use crate::util::backtrack::BackTracking;
use crate::dsets::*;
use crate::dsyms::*;


#[derive(Clone, Copy)]
pub enum Geometries {
    Spherical,
    Euclidean,
    Hyperbolic,
    All,
}

impl Geometries {
    fn min_curvature(&self) -> i64 {
        match self {
            Geometries::Spherical => 1,
            Geometries::Euclidean => 0,
            Geometries::Hyperbolic => i64::MIN,
            Geometries::All => i64::MIN,
        }
    }

    fn max_curvature(&self) -> i64 {
        match self {
            Geometries::Spherical => 4 * CURV_FAC,
            Geometries::Euclidean => 0,
            Geometries::Hyperbolic => -1,
            Geometries::All => 4 * CURV_FAC,
        }
    }
}

impl FromStr for Geometries {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match &s.to_lowercase()[..] {
            "spherical" | "s" => Ok(Geometries::Spherical),
            "euclidean" | "e" => Ok(Geometries::Euclidean),
            "hyperbolic" | "h" => Ok(Geometries::Hyperbolic),
            "all" | "a" => Ok(Geometries::All),
            _ => Err("unknown geometry '".to_string() + s + "'"),
        }
    }
}


struct DSymGenState {
    vs: Vec<usize>,
    curv: i64,
    next: usize,
}


const CURV_FAC: i64 = 420;


struct DSymBackTracking {
    dset: SimpleDSet,
    orbit_index: Vec<Vec<usize>>,
    orbit_rs: Vec<usize>,
    orbit_vmins: Vec<usize>,
    orbit_is_chain: Vec<bool>,
    orbit_maps: Option<Vec<Vec<usize>>>,
    base_curvature: i64,
    min_curvature: i64,
    max_curvature: i64,
}

impl DSymBackTracking {
    fn new(dset: &SimpleDSet, geoms: Geometries) -> DSymBackTracking {
        let (orbit_rs, orbit_is_chain, orbit_index) = collect_orbits(&dset);
        let orbit_vmins = compute_vmins(&orbit_rs);

        let mut base_curvature = -CURV_FAC / 2 * dset.size() as i64;
        for i in 0..orbit_vmins.len() {
            let k = if orbit_is_chain[i] { 1 } else { 2 };
            base_curvature += CURV_FAC * k / orbit_vmins[i] as i64;
        }

        let min_curvature = geoms.min_curvature().max(if base_curvature < 0 {
            base_curvature
        } else {
            -CURV_FAC // implied by minimal hyperbolicity
        });
        let max_curvature = geoms.max_curvature();

        let orbit_maps = if base_curvature >= 0 {
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
            min_curvature,
            max_curvature,
        }
    }

    fn orbit_count(&self) -> usize {
        self.orbit_vmins.len()
    }

    fn is_minimally_hyperbolic(&self, vs: &[usize], curv: i64) -> bool {
        if curv < 0 {
            for i in 0..self.orbit_count() {
                if vs[i] > self.orbit_vmins[i] {
                    let v = vs[i] as i64;
                    let k = if self.orbit_is_chain[i] { 1 } else { 2 };
                    let c = curv - k * CURV_FAC / v + k * CURV_FAC / (v - 1);
                    if c < 0 {
                        return false;
                    }
                }
            }
            true
        } else {
            false
        }
    }

    fn is_canonical(&self, vs: &[usize]) -> bool {
        for m in self.orbit_maps.as_ref().unwrap() {
            let ws: Vec<_> = (0..vs.len()).map(|i| vs[m[i]]).collect();
            if &ws[..] > vs {
                return false;
            }
        }

        true
    }

    fn is_good(&self, vs: &[usize], curv: i64) -> bool {
        if curv <= 0 {
            true
        } else {
            let key = self.orbifold_symbol(vs);

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

    fn orbifold_symbol(&self, vs: &[usize]) -> String {
        let ds = &self.dset;
        let mut cones = vec![];
        let mut corners = vec![];

        for d in ds.orbit_reps_2d(0, 2) {
            let d0 = ds.op_unchecked(0, d);
            let d2 = ds.op_unchecked(2, d);

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

        let front = degree_list_as_string(cones);
        let middle = String::from(
            if self.dset.is_loopless() { "" } else { "*" }
        );
        let back = degree_list_as_string(corners);
        let cross = String::from(
            if self.is_weakly_oriented() { "" } else { "x" }
        );

        [ front, middle, back, cross ].join("")
    }

    fn is_weakly_oriented(&self) -> bool {
        let mut sgn = vec![0; self.dset.size() + 1];
        let mut queue = VecDeque::new();

        sgn[1] = 1;
        queue.push_back(1);

        while let Some(d) = queue.pop_front() {
            for i in 0..=self.dset.dim() {
                if let Some(di) = self.dset.op(i, d) {
                    if sgn[di] == 0 {
                        sgn[di] = -sgn[d];
                        queue.push_back(di);
                    } else if di != d && sgn[di] != -sgn[d] {
                        return false;
                    }
                }
            }
        }

        true
    }
}


fn degree_list_as_string(vs: Vec<usize>) -> String {
    vs.iter().map(|v| v.to_string()).collect::<Vec<_>>().join("")
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
        let next = if curv < 0 { self.orbit_count() } else { 0 };

        Self::State { vs, curv, next }
    }

    fn extract(&self, state: &Self::State) -> Option<Self::Item> {
        let good =
            state.curv >= self.min_curvature &&
            state.curv <= self.max_curvature && (
                self.base_curvature < 0 || (
                    state.next >= self.orbit_count() &&
                    self.is_good(&state.vs, state.curv) &&
                    self.is_canonical(&state.vs)
                )
            );

        if good {
            Some(PartialDSym::from_fields(
                self.dset.clone(),
                self.orbit_index.clone(),
                self.orbit_rs.clone(),
                state.vs.iter().cloned().collect()
            ))
        } else {
            None
        }
    }

    fn children(&self, state: &Self::State) -> Vec<Self::State> {
        let n = state.next;

        if n >= self.orbit_count() || self.base_curvature < 0 {
            Vec::new()
        } else {
            let vmin = state.vs[n] as i64;
            let mut result = Vec::new();

            for v in vmin..=7 {
                let mut vs = state.vs.clone();
                vs[n] = v as usize;

                let k = if self.orbit_is_chain[n] { 1 } else { 2 };
                let curv = state.curv - k * CURV_FAC / vmin + k * CURV_FAC / v;

                if curv >= self.min_curvature {
                    if curv < 0 {
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
}


pub struct DSyms {
    bt: BackTrackIterator<DSymBackTracking>,
    counter: usize,
}

impl DSyms {
    pub fn new(dset: &SimpleDSet, geoms: Geometries) -> DSyms {
        DSyms {
            bt: BackTrackIterator::new(DSymBackTracking::new(dset, geoms)),
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
