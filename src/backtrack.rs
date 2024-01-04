pub trait BackTracking {
    type State;
    type Item;

    fn root(&self) -> Self::State;
    fn extract(&self, state: &Self::State) -> Option<Self::Item>;
    fn children(&self, state: &Self::State) -> Vec<Self::State>;
}

pub struct BackTrackIterator<T: BackTracking> {
    bt: T,
    stack: Vec<Vec<T::State>>,
}

impl<T: BackTracking> BackTrackIterator<T> {
    pub fn new(bt: T) -> BackTrackIterator<T> {
        BackTrackIterator { stack: vec![vec![T::root(&bt)]], bt }
    }
}

impl<T: BackTracking> Iterator for BackTrackIterator<T> {
    type Item = T::Item;

    fn next(&mut self) -> Option<Self::Item> {
        while self.stack.len() > 0 {
            let current = self.stack.last().unwrap().last().unwrap();
            let value = T::extract(&self.bt, current);
            let todo = T::children(&self.bt, current);

            if todo.len() > 0 {
                self.stack.push(todo.into_iter().rev().collect());
            } else {
                while self.stack.len() > 0 &&
                      self.stack.last().unwrap().len() < 2
                {
                    self.stack.pop();
                }

                if let Some(todo) = self.stack.last_mut() {
                    todo.pop();
                }
            }

            if value.is_some() {
                return value;
            }
        }

        None
    }
}
