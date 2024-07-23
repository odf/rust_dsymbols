use std::collections::VecDeque;

use std::collections::HashSet;
use std::fmt;

use crate::derived::as_dsym;
use crate::derived::canonical;
use crate::dsets::DSet;
use crate::dsets::PartialDSet;

pub(crate) struct DrawingInstructions<'a>
{
    pub(crate) ds: &'a PartialDSet
}

impl<'a> DrawingInstructions<'a>
{
    pub fn new(ds: &'a PartialDSet) -> Self {
        assert_eq!(canonical(&as_dsym(ds)), as_dsym(ds));
        assert!(ds.is_oriented());
        assert!(ds.is_complete());

        Self { ds }
    }
}

impl<'a> fmt::Display for DrawingInstructions<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let ds = self.ds;
        let mut seen= HashSet::new();
        let mut queue = VecDeque::from([(0, 1)]);

        let r = |d| ds.r(0, 1, d).unwrap();

        while let Some((c, d)) = queue.pop_front() {
            if !seen.contains(&d) {
                write!(f, "{}-gon from {d}", r(d))?;
                if c != 0 {
                    write!(f, " connected to {c}")?;
                }
                writeln!(f, "")?;

                seen.extend(ds.orbit([0, 1], d));

                let mut e = d;
                loop {
                    queue.push_back((e, ds.op(2, e).unwrap()));
                    e = ds.walk(e, [0, 1]).unwrap();
                    if e == d {
                        break;
                    }
                }
            } else if d > c {
                writeln!(f, "connect {c} to {d}")?;
            }
        }

        for d in ds.orbit_reps([0, 1, 3], 1..ds.size()) {
            writeln!(f, "3-connect {d} to {}", ds.op(3, d).unwrap())?;
        }
        writeln!(f, "")?;

        Ok(())
    }
}
