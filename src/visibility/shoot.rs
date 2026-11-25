//! A horizontal shot is specified by one of the two ports of a feature.
//! Comparing intersected edges uses orientation, never rational intersection
//! coordinates. Behind the origin, the order continues through infinity.
use super::boundary::{Arc, Chain, Chord, Feature, Site};
use crate::geometry::Ring;
use std::cmp::Ordering;

#[derive(Clone, Copy)]
pub(super) struct Shot {
    pub level: usize,
    pub origin: Feature,
    pub right: bool,
    pub target: Feature,
}
impl Shot {
    pub fn new(level: usize, origin: Feature, right: bool) -> Self {
        Self {
            level,
            origin,
            right,
            target: origin,
        }
    }
    fn along(&self, ring: Ring, a: Feature, b: Feature) -> Ordering {
        let ord = a.cmp_x(b, ring, self.level);
        if self.right { ord } else { ord.reverse() }
    }
    pub fn consider(&mut self, ring: Ring, f: Feature) {
        if f == self.origin || f == self.target {
            return;
        }
        let direct = self.along(ring, self.origin, f) == Ordering::Less;
        let old_direct = self.along(ring, self.origin, self.target) == Ordering::Less;
        if self.target == self.origin
            || (direct && !old_direct)
            || (direct == old_direct && self.along(ring, f, self.target) == Ordering::Less)
        {
            self.target = f
        }
    }
    pub fn edge(&mut self, ring: Ring, e: usize) {
        let (lo, hi) = ring.ends(e);
        let h = ring.0[self.level].height();
        if lo.height() < h && h < hi.height() {
            self.consider(ring, Feature::Edge(e))
        }
    }
    pub fn chord(self) -> Chord {
        let (from, to) = if self.right {
            (self.origin, self.target)
        } else {
            (self.target, self.origin)
        };
        Chord::Beam {
            level: self.level,
            from,
            to,
        }
    }
}

pub(super) fn scan_arc(ring: Ring, chain: Chain, arc: Arc, shot: &mut Shot) {
    for k in arc.steps(chain) {
        let step = k % chain.period();
        let edge = chain.edge(step);
        let (lo, hi) = ring.ends(edge);
        let h = ring.0[shot.level].height();
        if !(lo.height() < h && h < hi.height()) {
            continue;
        }
        let f = Feature::Edge(edge);
        let right = (step < chain.period() / 2) != ring.rising(edge);
        let site = chain.site(ring, shot.level, f, right);
        if arc.covers(site) {
            shot.consider(ring, f)
        }
    }
    if (chain.first..=chain.last).contains(&shot.level) {
        for step in chain.steps(shot.level) {
            if arc.covers(Site::corner(step, 0)) || arc.covers(Site::corner(step, 1)) {
                shot.consider(ring, Feature::Vertex)
            }
        }
    }
}
