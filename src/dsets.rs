pub trait AbstractDelaneySet {
    fn set_count(&self) -> usize { 1 }
    fn symbol_count(&self) -> usize { 1 }
    fn m(&self, _i: usize, _j: usize, _d: usize) -> usize { 0 }
}


struct DelaneySet {
    size: usize,
    dim: usize,
    op: Vec<usize>,
}

impl DelaneySet {
    fn new(size: usize, dim: usize) -> DelaneySet {
        let mut op = Vec::new();
        op.resize(size * (dim + 1), 0);
        DelaneySet { size, dim, op }
    }

    fn idx(&self, i: usize, d: usize) -> usize {
        i * self.size + d - 1
    }

    fn get(&self, i: usize, d: usize) -> Option<usize> {
        if i > self.dim || d < 1 || d > self.size {
            None
        } else {
            let di = *self.op.get(self.idx(i, d)).unwrap_or(&0);
            if di == 0 { None } else { Some(di) }
        }
    }

    fn set(&mut self, i: usize, d: usize, e: usize) {
        if i > self.dim || d < 1 || d > self.size || e < 1 || e > self.size {
            panic!("parameter out of range");
        }
        let kd = self.idx(i, d);
        let ke = self.idx(i, e);
        self.op[kd] = e;
        self.op[ke] = d;
    }
}

impl AbstractDelaneySet for DelaneySet {}
