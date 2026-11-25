//! Grade table and the arc cutter of §4.1. Intervals are aligned relative to
//! the start of the input curve. A short final interval replaces padding.
use super::boundary::{Arc, Chain, Feature};
use super::decomposition::Node;
use super::hierarchy::{base, merge};
use super::oracle::Oracle;
use super::shoot::{Shot, scan_arc};
use super::submap::Submap;
use crate::geometry::Ring;
use std::cell::OnceCell;

pub(super) const BASE: usize = 4;
pub(super) fn grade(m: usize) -> u32 {
    m.next_power_of_two().ilog2()
}
pub(super) fn granularity(m: usize) -> usize {
    1 << grade(m).div_ceil(5)
}

pub(super) struct Piece {
    pub map: Submap,
    pub tree: OnceCell<Node>,
    highest: usize,
    oracle: OnceCell<Oracle>,
}
pub(super) struct Catalog {
    origin: usize,
    levels: Vec<Vec<Piece>>,
}
impl Catalog {
    pub fn build(ring: Ring, whole: Chain) -> Self {
        let size = whole.last - whole.first;
        let mut table = Self {
            origin: whole.first,
            levels: Vec::new(),
        };
        for g in BASE.ilog2()..=grade(size).max(BASE.ilog2()) {
            let span = 1 << g;
            let mut level = Vec::new();
            for first in (whole.first..whole.last).step_by(span) {
                let chain = Chain {
                    first,
                    last: (first + span).min(whole.last),
                };
                let len = chain.last - first;
                if span > BASE && len <= span / 2 {
                    break;
                }
                let (map, highest) = if span == BASE {
                    (
                        base(ring, chain),
                        (first..=chain.last)
                            .max_by_key(|&v| ring.0[v].height())
                            .unwrap(),
                    )
                } else {
                    let middle = first + span / 2;
                    let left = table.get(Chain {
                        first,
                        last: middle,
                    });
                    let right = table.get(Chain {
                        first: middle,
                        last: chain.last,
                    });
                    let high = if ring.0[left.highest].height() > ring.0[right.highest].height() {
                        left.highest
                    } else {
                        right.highest
                    };
                    (
                        merge(ring, &table, &left.map, &right.map, granularity(len)),
                        high,
                    )
                };
                // The retained maps support both phases. Build a map's
                // query structures only when a cutter first requests it.
                level.push(Piece {
                    tree: OnceCell::new(),
                    map,
                    highest,
                    oracle: OnceCell::new(),
                });
            }
            table.levels.push(level);
        }
        table
    }
    pub fn get(&self, chain: Chain) -> &Piece {
        let g = grade(chain.last - chain.first).max(BASE.ilog2());
        let p = &self.levels[(g - BASE.ilog2()) as usize][(chain.first - self.origin) >> g];
        debug_assert_eq!(
            (p.map.chain.first, p.map.chain.last),
            (chain.first, chain.last)
        );
        p
    }
    /// Partition an integral interval into O(log length) aligned pieces.
    pub fn chunks(&self, mut first: usize, last: usize, mut visit: impl FnMut(Chain)) {
        while first < last {
            let g = (last - first)
                .ilog2()
                .min((first - self.origin).trailing_zeros());
            let end = first + (1 << g);
            visit(Chain { first, last: end });
            first = end
        }
    }
    /// Vertex-to-vertex pieces, split at cap turns and dyadic boundaries.
    /// A vertex between two clipped edges can be a singleton piece.
    pub fn cut_arc(&self, chain: Chain, arc: Arc, mut visit: impl FnMut(Chain, bool)) {
        let steps = arc.steps(chain);
        let mut begin = steps.start + usize::from(arc.start.in_edge());
        let end = steps.end.saturating_sub(usize::from(arc.end.in_edge()));
        let half = chain.period() / 2;
        if begin == end {
            let v = chain.vertex(begin % chain.period());
            visit(Chain { first: v, last: v }, true)
        }
        while begin < end {
            let stop = end.min((begin / half + 1) * half);
            let a = chain.vertex(begin % chain.period());
            let b = chain.vertex(stop % chain.period());
            self.chunks(a.min(b), a.max(b), |piece| visit(piece, a < b));
            begin = stop;
        }
    }
    pub fn shoot_arc(&self, ring: Ring, chain: Chain, arc: Arc, shot: &mut Shot) {
        let steps = arc.steps(chain);
        if steps.len() <= BASE {
            scan_arc(ring, chain, arc, shot);
            return;
        }
        let begin = steps.start;
        let end = steps.end;
        // Only the two extreme steps can be clipped in the middle of an edge.
        for (k, partial) in [(begin, arc.start.in_edge()), (end - 1, arc.end.in_edge())] {
            if partial {
                let step = k % chain.period();
                let mut clipped = *shot;
                clipped.target = shot.origin;
                clipped.edge(ring, chain.edge(step));
                if clipped.target != shot.origin {
                    let r = (step < chain.period() / 2) != ring.rising(chain.edge(step));
                    if arc.covers(chain.site(ring, shot.level, clipped.target, r)) {
                        shot.consider(ring, clipped.target)
                    }
                }
            }
        }
        self.cut_arc(chain, arc, |piece, _| {
            if piece.last - piece.first <= BASE {
                if (piece.first..=piece.last).contains(&shot.level) {
                    shot.consider(ring, Feature::Vertex)
                }
                for e in piece.first..piece.last {
                    shot.edge(ring, e)
                }
            } else {
                let p = self.get(piece);
                p.oracle
                    .get_or_init(|| Oracle::new(ring, &p.map, p.highest))
                    .query(ring, &p.map, shot)
            }
        });
    }
}
