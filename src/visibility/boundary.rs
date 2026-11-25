//! Coordinates on the clockwise double boundary of an open polygonal chain.
//! The clockwise tour uses the left side when traversing the chain forward.
use crate::geometry::{Ring, turn};
use std::cmp::Ordering;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub(super) enum Feature {
    Vertex,
    Edge(usize),
}
impl Feature {
    pub fn cmp_x(self, other: Self, ring: Ring, level: usize) -> Ordering {
        match (self, other) {
            (Self::Vertex, Self::Vertex) => Ordering::Equal,
            (Self::Edge(e), Self::Vertex) => ring.edge_vs_vertex(e, level),
            (Self::Vertex, Self::Edge(e)) => ring.edge_vs_vertex(e, level).reverse(),
            (Self::Edge(e), Self::Edge(f)) => ring.compare_edges(e, f),
        }
    }
}
/// The beam travels right from 'from', wrapping through infinity if needed.
/// Every height used by the construction names an original input vertex.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub(super) enum Chord {
    Beam {
        level: usize,
        from: Feature,
        to: Feature,
    },
    Null(usize),
}
/// An interior offset measures height from the start of the oriented step.
/// At a corner offset zero, copies zero and one bound the turn in tour order.
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub(super) struct Site {
    pub step: usize,
    offset: u64,
    copy: u8,
}
impl Site {
    pub fn at_step(self, step: usize) -> Self {
        Self { step, ..self }
    }
    pub fn in_edge(self) -> bool {
        self.offset != 0
    }
    pub fn last_copy(self) -> Self {
        Self { copy: 1, ..self }
    }
    pub fn corner(step: usize, copy: u8) -> Self {
        Self {
            step,
            offset: 0,
            copy,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub(super) struct Chain {
    pub first: usize,
    pub last: usize,
}

fn height_code(ring: Ring, v: usize) -> u64 {
    let (y, x) = ring.0[v].height();
    let biased = |a: i32| (a as u32) ^ (1 << 31);
    (u64::from(biased(y)) << 32) | u64::from(biased(x))
}
impl Chain {
    pub fn period(self) -> usize {
        2 * (self.last - self.first)
    }
    pub fn vertex(self, step: usize) -> usize {
        self.first + step.min(self.period() - step)
    }
    pub fn edge(self, step: usize) -> usize {
        self.vertex(step).min(self.vertex(step + 1))
    }
    pub fn steps(self, v: usize) -> [usize; 2] {
        let f = v - self.first;
        [f, (self.period() - f) % self.period()]
    }
    // The outside of a right turn is on the forward (left) strand.
    fn outer_step(self, ring: Ring, v: usize) -> usize {
        let [forward, reverse] = self.steps(v);
        if forward == reverse || turn(ring.0[v - 1], ring.0[v], ring.0[v + 1]) < 0 {
            forward
        } else {
            reverse
        }
    }
    pub fn extremum(self, ring: Ring, v: usize) -> bool {
        self.first < v && v < self.last && ring.rising(v - 1) != ring.rising(v)
    }

    pub fn site(self, ring: Ring, level: usize, f: Feature, right: bool) -> Site {
        if let Feature::Edge(e) = f {
            let forward = ring.rising(e) != right;
            let step = if forward {
                e - self.first
            } else {
                self.period() - 1 - (e - self.first)
            };
            return Site {
                step,
                offset: height_code(ring, self.vertex(step)).abs_diff(height_code(ring, level)),
                copy: 0,
            };
        }
        let [fwd, rev] = self.steps(level);
        if fwd != rev && !self.extremum(ring, level) {
            return Site::corner(
                if ring.rising(level) != right {
                    fwd
                } else {
                    rev
                },
                0,
            );
        }
        let step = self.outer_step(ring, level);
        let leaving_rises = ring.0[self.vertex(step + 1)].height() > ring.0[level].height();
        Site::corner(step, u8::from(leaving_rises != right))
    }
    pub fn sites(self, ring: Ring, c: Chord) -> [Site; 2] {
        match c {
            Chord::Null(v) => {
                let [fwd, rev] = self.steps(v);
                let inner = if self.outer_step(ring, v) == fwd {
                    rev
                } else {
                    fwd
                };
                [Site::corner(inner, 0), Site::corner(inner, 1)]
            }
            Chord::Beam { level, from, to } => [
                self.site(ring, level, from, true),
                self.site(ring, level, to, false),
            ],
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct Arc {
    pub start: Site,
    pub end: Site,
}
impl Arc {
    /// Return the unrolled step interval containing nonzero-length pieces.
    pub fn steps(self, chain: Chain) -> std::ops::Range<usize> {
        let end = self.end.step
            + usize::from(self.end.offset != 0)
            + if self.end <= self.start {
                chain.period()
            } else {
                0
            };
        self.start.step..end
    }
    pub fn contains(self, p: Site) -> bool {
        if self.start < self.end {
            self.start < p && p < self.end
        } else {
            self.start < p || p < self.end
        }
    }
    /// Classify the open interval immediately following p.
    pub fn after(self, p: Site) -> bool {
        self.start == p || self.contains(p)
    }
    pub fn covers(self, p: Site) -> bool {
        self.start == p || self.end == p || self.contains(p)
    }
    pub fn full() -> Self {
        let at = Site::corner(0, 0);
        Self { start: at, end: at }
    }
}
