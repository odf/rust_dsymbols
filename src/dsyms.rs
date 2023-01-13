use crate::dsets::*;

pub trait DSym : DSet {
    fn r(&self, i: usize, j: usize, d: usize) -> usize;
    fn v(&self, i: usize, j: usize, d: usize) -> usize;

    fn m(&self, i: usize, j: usize, d: usize) -> usize {
        self.r(i, j, d) * self.v(i, j, d)
    }
}


pub fn collect_orbits(ds: &SimpleDSet) -> (Vec<usize>, Vec<Vec<usize>>) {
    let mut orbit_rs = vec![];
    let mut orbit_index = vec![vec![0; ds.size() + 1]; ds.dim()];
    let mut seen = vec![false; ds.size() + 1];

    for i in 0..ds.dim() {
        seen.fill(false);

        for d in 1..=ds.size() {
            if !seen[d] {
                let orbit_nr = orbit_rs.len();

                orbit_index[i][d] = orbit_nr;
                seen[d] = true;

                let mut e = d;
                let mut k = i;
                let mut steps = 0;

                loop {
                    e = ds.get(k, e).unwrap();
                    k = i + (i + 1) - k;
                    steps += 1;

                    orbit_index[i][e] = orbit_nr;
                    seen[e] = true;

                    if e == d && k == i {
                        break;
                    }
                }

                orbit_rs.push(steps / 2);
            }
        }
    }

    (orbit_rs, orbit_index)
}
