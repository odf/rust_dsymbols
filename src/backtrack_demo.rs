use super::backtrack::BackTrackIterator;
use super::backtrack::BackTracking;


struct PartitionState {
    xs: Vec<u32>,
    left: u32,
    top: u32,
}

struct PartitionBackTracking {
    n: u32,
}

impl BackTracking for PartitionBackTracking {
    type State = PartitionState;
    type Item = Vec<u32>;

    fn root(&self) -> Self::State {
        PartitionState { xs: vec![], left: self.n, top: 1 }
    }

    fn extract(&self, state: &Self::State) -> Option<Self::Item> {
        if state.left == 0 { Some(state.xs.clone()) } else { None }
    }

    fn children(&self, state: &Self::State) -> Vec<Self::State> {
        (state.top..=state.left).map(|i| {
            PartitionState {
                xs: [state.xs.clone(), vec![i]].concat(),
                left: state.left - i,
                top: state.top.max(i),
            }
        }).collect()
    }
}


pub struct IntPartitions(BackTrackIterator<PartitionBackTracking>);

impl IntPartitions {
    pub fn new(n: u32) -> IntPartitions {
        IntPartitions(BackTrackIterator::new(PartitionBackTracking { n }))
    }
}

impl Iterator for IntPartitions {
    type Item = Vec<u32>;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next()
    }
}
