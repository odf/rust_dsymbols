// Generic backtracker

pub trait BackTracker {
    type State;
    type Item;

    fn root(&self) -> Self::State;
    fn extract(&self, state: &Self::State) -> Option<Self::Item>;
    fn children(&self, state: &Self::State) -> Vec<Self::State>;
}

pub struct BackTrackIterator<T: BackTracker> {
    bt: T,
    stack: Vec<Vec<T::State>>,
}

impl<T: BackTracker> BackTrackIterator<T> {
    pub fn new(bt: T) -> BackTrackIterator<T> {
        BackTrackIterator { stack: vec![vec![T::root(&bt)]], bt }
    }
}

impl<T: BackTracker> Iterator for BackTrackIterator<T> {
    type Item = T::Item;

    fn next(&mut self) -> Option<Self::Item> {
        while self.stack.len() > 0 {
            let list = self.stack.last().unwrap();
            let current = list.last().unwrap();
            let value = T::extract(&self.bt, current);
            let todo = T::children(&self.bt, current);

            if todo.len() > 0 {
                let mut todo = todo;
                todo.reverse();
                self.stack.push(todo);
            } else {
                while self.stack.len() > 0 &&
                      self.stack.last().unwrap().len() < 2
                {
                    self.stack.pop();
                }

                if self.stack.len() > 0 {
                    let mut list = self.stack.pop().unwrap();
                    list.pop();
                    self.stack.push(list);
                }
            }

            if let Some(value) = value {
                return Some(value);
            }
        }

        None
    }
}
