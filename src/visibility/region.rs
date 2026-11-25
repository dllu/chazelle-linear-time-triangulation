//! Local visibility inside a region. To find a conformality chord, search
//! the separator chords of each aligned piece of a candidate arc (Lemma 3.2).
use super::boundary::{Arc, Chain, Chord, Feature, Site};
use super::catalog::{BASE, Catalog, Piece};
use super::decomposition::Node;
use super::shoot::Shot;
use crate::geometry::Ring;

pub(super) struct Region<'a> {
    pub ring: Ring<'a>,
    pub table: &'a Catalog,
    pub chain: Chain,
    pub arcs: &'a [Arc],
}
impl Region<'_> {
    pub fn locate(&self, p: Site) -> Option<usize> {
        self.arcs
            .iter()
            .position(|a| a.start == a.end || a.contains(p))
    }
    pub fn beam(&self, level: usize, origin: Feature, right: bool) -> Chord {
        let mut shot = Shot::new(level, origin, right);
        for &arc in self.arcs {
            self.table.shoot_arc(self.ring, self.chain, arc, &mut shot)
        }
        shot.chord()
    }
    /// All ports at a vertex whose chords fall in this region.
    pub fn at_vertex(&self, v: usize, mut emit: impl FnMut(Chord)) {
        if self.chain.extremum(self.ring, v) {
            let c = Chord::Null(v);
            if self
                .chain
                .sites(self.ring, c)
                .iter()
                .all(|&p| self.locate(p).is_some())
            {
                emit(c)
            }
        }
        for dir in [false, true] {
            if self
                .locate(self.chain.site(self.ring, v, Feature::Vertex, dir))
                .is_some()
            {
                emit(self.beam(v, Feature::Vertex, dir))
            }
        }
    }
    fn vertex_witness(&self, v: usize, a: usize, b: usize) -> Option<Chord> {
        let mut found = None;
        self.at_vertex(v, |c| {
            let [p, q] = self.chain.sites(self.ring, c);
            let owners = [self.locate(p), self.locate(q)];
            if owners.contains(&Some(a)) && owners.contains(&Some(b)) {
                found = Some(c)
            }
        });
        found
    }

    fn piece_witness(&self, piece: &Piece, forward: bool, a: usize, b: usize) -> Option<Chord> {
        let map = &piece.map;
        let local = map.chain;
        let len = local.last - local.first;
        let start = if forward { 0 } else { len };
        let strand = Arc {
            start: Site::corner(start, 1),
            end: Site::corner(start + len, 0),
        };
        let parent_start = if forward {
            local.first - self.chain.first
        } else {
            self.chain.period() - (local.last - self.chain.first)
        };
        let lift = |p: Site| p.at_step((parent_start + (p.step - start)) % self.chain.period());
        let mut node = piece.tree.get_or_init(|| Node::build(map));
        loop {
            match node {
                Node::Leaf(region) => {
                    for arc in map.arcs(&map.cells[*region]) {
                        let steps = arc.steps(local);
                        for k in steps.start..=steps.end {
                            let v = local.vertex(k % local.period());
                            if let Some(c) = self.vertex_witness(v, a, b) {
                                return Some(c);
                            }
                        }
                    }
                    return None;
                }
                Node::Branch {
                    end,
                    inside,
                    outside,
                } => {
                    let chord = map.chords[map.chord_at[*end]];
                    let pair = local.sites(self.ring, chord);
                    let mut bounds = vec![strand.start, strand.end];
                    let mut barriers = Vec::new();
                    for (i, p) in pair.into_iter().enumerate() {
                        if !strand.contains(p) {
                            continue;
                        }
                        bounds.push(p);
                        let actual = match chord {
                            Chord::Null(_) => chord,
                            Chord::Beam { level, from, to } => {
                                self.beam(level, if i == 0 { from } else { to }, i == 0)
                            }
                        };
                        let [u, v] = self.chain.sites(self.ring, actual);
                        let source = lift(p);
                        debug_assert!(u == source || v == source);
                        let target = if u == source { v } else { u };
                        if self.locate(target) == Some(b) {
                            return Some(actual);
                        }
                        barriers.push(Arc { start: u, end: v });
                    }
                    // The visible segment from each separator port splits the
                    // parent's boundary circle. Keep intervals on the target
                    // side of every such cut. Lemma 2.4 puts all survivors
                    // in a single component of the piece's separator chord.
                    bounds.sort_unstable();
                    let component = Arc {
                        start: map.sites[*end],
                        end: map.sites[map.mate[*end]],
                    };
                    let target = self.arcs[b].start;
                    let mut branch = None;
                    for &p in &bounds[..bounds.len() - 1] {
                        if barriers
                            .iter()
                            .all(|cut| cut.after(lift(p)) == cut.after(target))
                        {
                            let side = component.after(p);
                            debug_assert!(
                                branch.is_none_or(|prev| prev == side),
                                "shielding failed"
                            );
                            branch = Some(side)
                        }
                    }
                    node = if branch? { inside } else { outside };
                }
            }
        }
    }

    fn pair_witness(&self, a: usize, b: usize) -> Option<Chord> {
        let arc = self.arcs[a];
        let mut found = None;
        self.table.cut_arc(self.chain, arc, |chain, forward| {
            if found.is_some() {
                return;
            }
            if chain.last - chain.first <= BASE {
                found = (chain.first..=chain.last).find_map(|w| self.vertex_witness(w, a, b));
            } else {
                found = [chain.first, chain.last]
                    .into_iter()
                    .find_map(|w| self.vertex_witness(w, a, b))
                    .or_else(|| self.piece_witness(self.table.get(chain), forward, a, b))
            }
        });
        found
    }
    pub fn split_witness(&self) -> Chord {
        let n = self.arcs.len();
        for a in 0..n {
            for offset in 2..n - 1 {
                if let Some(c) = self.pair_witness(a, (a + offset) % n) {
                    return c;
                }
            }
        }
        unreachable!("Lemma 3.3 requires visibility between nonconsecutive arcs")
    }
}

#[cfg(test)]
mod tests {
    use super::super::{refine, submap::Submap};
    use super::*;
    use crate::Point;
    #[test]
    fn shielding_in_sparse_regions() {
        let mut state = 334u32;
        let mut random = || {
            state = state.wrapping_mul(1664525).wrapping_add(1013904223);
            state
        };
        for trial in 0..64 {
            let p: Vec<_> = (0..129)
                .map(|v| {
                    if trial % 2 == 0 {
                        Point::new(v, (random() % 150) as i32)
                    } else {
                        let a = std::f64::consts::TAU * v as f64 / 129.0;
                        let r = 10000.0 + f64::from(random() % 90000);
                        Point::new((r * a.cos()).round() as i32, (r * a.sin()).round() as i32)
                    }
                })
                .collect();
            let chain = Chain {
                first: 0,
                last: p.len() - 1,
            };
            let ring = Ring(&p);
            let table = Catalog::build(ring, chain);
            let map = &table.get(chain).map;
            let mut chords = refine::full_map(ring, &table, map);
            chords.retain(|_| random() % 21 == 0);
            let map = Submap::new(ring, chain, chords);
            for cell in &map.cells {
                if cell.len() < 4 {
                    continue;
                }
                let arcs: Vec<_> = map.arcs(cell).collect();
                let region = Region {
                    ring,
                    table: &table,
                    chain,
                    arcs: &arcs,
                };
                for a in 0..arcs.len() {
                    for offset in 2..arcs.len() - 1 {
                        let b = (a + offset) % arcs.len();
                        // Lemma 3.2 promises a hit whenever an original vertex
                        // of A sees B. It may also return a Steiner endpoint.
                        let expected = (0..p.len()).any(|v| {
                            [false, true].iter().any(|&dir| {
                                let port = chain.site(ring, v, Feature::Vertex, dir);
                                if !arcs[a].contains(port) {
                                    return false;
                                }
                                let mut shot = Shot::new(v, Feature::Vertex, dir);
                                for e in 0..chain.last {
                                    shot.edge(ring, e)
                                }
                                arcs[b].contains(chain.site(ring, v, shot.target, !dir))
                            })
                        });
                        let witness = region.pair_witness(a, b);
                        assert!(!expected || witness.is_some(), "{trial} {a} {b}");
                        if let Some(Chord::Beam { level, from, to }) = witness {
                            let mut shot = Shot::new(level, from, true);
                            for e in 0..chain.last {
                                shot.edge(ring, e)
                            }
                            shot.consider(ring, Feature::Vertex);
                            assert_eq!(shot.target, to)
                        }
                    }
                }
            }
        }
    }
}
