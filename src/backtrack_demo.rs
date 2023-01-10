use super::backtrack::BackTrackIterator;
use super::backtrack::BackTracker;


struct PartitionState {
    xs: Vec<u32>,
    left: u32,
    top: u32,
}

struct Partitions {
    n: u32,
}

impl BackTracker for Partitions {
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

pub fn print_integer_partitions(n: u32) {
    for xs in BackTrackIterator::new(Partitions { n }) {
        println!("{:?}", xs);
    }
}
