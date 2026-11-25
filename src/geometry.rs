use std::cmp::Ordering;

/// Cartesian coordinates. Predicates are exact over the entire i32 range.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Point {
    pub x: i32,
    pub y: i32,
}

impl Point {
    pub const fn new(x: i32, y: i32) -> Self {
        Self { x, y }
    }

    /// Order after the orientation-preserving shear y' = y + epsilon*x.
    pub(crate) fn height(self) -> (i32, i32) {
        (self.y, self.x)
    }
}

pub(crate) fn turn(a: Point, b: Point, c: Point) -> i128 {
    let dx = i128::from(b.x) - i128::from(a.x);
    let dy = i128::from(b.y) - i128::from(a.y);
    dx * (i128::from(c.y) - i128::from(a.y)) - dy * (i128::from(c.x) - i128::from(a.x))
}

#[derive(Clone, Copy)]
pub(crate) struct Ring<'a>(pub &'a [Point]);
impl Ring<'_> {
    pub fn next(self, v: usize) -> usize {
        (v + 1) % self.0.len()
    }
    pub fn prev(self, v: usize) -> usize {
        (v + self.0.len() - 1) % self.0.len()
    }
    pub fn rising(self, e: usize) -> bool {
        self.0[e].height() < self.0[self.next(e)].height()
    }
    pub fn ends(self, e: usize) -> (Point, Point) {
        let (a, b) = (self.0[e], self.0[self.next(e)]);
        if self.rising(e) { (a, b) } else { (b, a) }
    }

    /// Order of two disjoint edge interiors on a common transverse scanline.
    /// At the later lower endpoint the other edge is already present. The
    /// order there persists throughout the overlap by polygon simplicity.
    pub fn compare_edges(self, e: usize, f: usize) -> Ordering {
        if e == f {
            return Ordering::Equal;
        }
        let (a, b) = self.ends(e);
        let (c, d) = self.ends(f);
        let order = match a.height().cmp(&c.height()) {
            Ordering::Less => turn(a, b, c),
            Ordering::Greater => -turn(c, d, a),
            Ordering::Equal => turn(a, b, d),
        };
        order.cmp(&0)
    }

    pub fn edge_vs_vertex(self, e: usize, v: usize) -> Ordering {
        let (a, b) = self.ends(e);
        turn(a, b, self.0[v]).cmp(&0)
    }
}
