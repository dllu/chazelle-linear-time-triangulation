//! Chazelle's global horizontal ray shooter (Lemma 3.6). A separator hit
//! identifies the only remaining component that can contain a closer hit.
use super::boundary::{Chord, Feature};
use super::mesh::Mesh;
use super::separator::Division;
use super::shoot::{Shot, scan_arc};
use super::submap::Submap;
use crate::geometry::Ring;
use std::cmp::Ordering;

pub(super) struct Oracle {
    division: Division,
    /// (height, region below, region above), ordered on the meridian at infinity.
    meridian: Vec<((i32, i32), usize, usize)>,
    pole: usize,
}
impl Oracle {
    pub fn new(ring: Ring, map: &Submap, highest: usize) -> Self {
        let division = Division::new(&Mesh::new(ring, map));
        let mut meridian = Vec::new();
        for (i, &site) in map.sites.iter().enumerate() {
            let Chord::Beam { level, from, to } = map.chords[map.chord_at[i]] else {
                continue;
            };
            if from.cmp_x(to, ring, level) != Ordering::Less
                && site == map.chain.site(ring, level, from, true)
            {
                meridian.push((ring.0[level].height(), map.owner[i], map.owner[map.mate[i]]))
            }
        }
        meridian.sort_unstable();
        let pole = map.region_after(map.chain.site(ring, highest, Feature::Vertex, true));
        Self {
            division,
            meridian,
            pole,
        }
    }
    fn scan(ring: Ring, map: &Submap, region: usize, shot: &mut Shot) {
        for arc in map.arcs(&map.cells[region]) {
            scan_arc(ring, map.chain, arc, shot)
        }
    }
    pub fn query(&self, ring: Ring, map: &Submap, result: &mut Shot) {
        let mut shot = Shot::new(result.level, result.origin, result.right);
        let c = map.chain;
        let contains = match shot.origin {
            Feature::Vertex => (c.first..=c.last).contains(&shot.level),
            Feature::Edge(e) => (c.first..c.last).contains(&e),
        };
        if contains {
            let site = c.site(ring, shot.level, shot.origin, shot.right);
            let cell = map.region_after(site);
            Self::scan(ring, map, cell, &mut shot);
            result.consider(ring, shot.target);
            return;
        }
        for &r in &self.division.separator {
            Self::scan(ring, map, r, &mut shot)
        }
        let probe = |region: usize, shot: &mut Shot| {
            if let Some(g) = self.division.group_of[region] {
                for &r in &self.division.groups[g] {
                    Self::scan(ring, map, r, shot)
                }
            }
        };
        if shot.target == shot.origin {
            let h = ring.0[shot.level].height();
            let i = self.meridian.partition_point(|&(level, _, _)| level < h);
            let region = if i > 0 {
                self.meridian[i - 1].2
            } else {
                self.meridian.first().map_or(self.pole, |entry| entry.1)
            };
            probe(region, &mut shot);
        } else {
            let site = c.site(ring, shot.level, shot.target, !shot.right);
            probe(map.region_after(site), &mut shot);
        }
        result.consider(ring, shot.target)
    }
}

#[cfg(test)]
mod tests {
    use super::super::{boundary::Chain, catalog::Catalog};
    use super::*;
    use crate::Point;

    #[test]
    fn meridian_and_foreign_edge_queries_across_components() {
        let mut state = 3089u32;
        let p: Vec<_> = (0..4100)
            .map(|x| {
                state = state.wrapping_mul(1664525).wrapping_add(1013904223);
                Point::new(x, (state % 100000) as i32)
            })
            .collect();
        let ring = Ring(&p);
        let chain = Chain {
            first: 1,
            last: p.len() - 2,
        };
        let table = Catalog::build(ring, chain);
        let map = &table.get(chain).map;
        let highest = (chain.first..=chain.last)
            .max_by_key(|&v| p[v].height())
            .unwrap();
        let shooter = Oracle::new(ring, map, highest);
        assert!(shooter.division.groups.len() > 1);
        assert!(!shooter.division.separator.is_empty());
        for e in [0, chain.last] {
            let (lo, hi) = ring.ends(e);
            for v in (chain.first..=chain.last).step_by(17) {
                if !(lo.height() < p[v].height() && p[v].height() < hi.height()) {
                    continue;
                }
                for right in [false, true] {
                    let mut query = Shot::new(v, Feature::Edge(e), right);
                    let mut expected = query;
                    expected.consider(ring, Feature::Vertex);
                    for w in chain.first..chain.last {
                        expected.edge(ring, w)
                    }
                    shooter.query(ring, map, &mut query);
                    assert_eq!(query.target, expected.target)
                }
            }
        }
    }
}
