//! Horizontal visibility in the coordinates induced by the symbolic shear.
mod boundary;
mod catalog;
mod decomposition;
mod fusion;
mod hierarchy;
mod mesh;
mod oracle;
mod refine;
mod region;
mod separator;
mod shoot;
mod submap;

use crate::geometry::Ring;
use boundary::{Chain, Chord, Feature};
use std::cmp::Ordering;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct Sight {
    /// Nearest finite hit on each side, excluding incident edges.
    pub hit: [Option<usize>; 2],
}
impl Sight {
    pub fn consider(&mut self, ring: Ring, v: usize, e: usize) {
        let (lo, hi) = ring.ends(e);
        let level = ring.0[v].height();
        if !(lo.height() < level && level < hi.height()) {
            return;
        }
        let left = ring.edge_vs_vertex(e, v) == Ordering::Less;
        let nearest = &mut self.hit[usize::from(!left)];
        let closer = if left {
            Ordering::Greater
        } else {
            Ordering::Less
        };
        if nearest.is_none_or(|f| ring.compare_edges(e, f) == closer) {
            *nearest = Some(e)
        }
    }
}
pub(crate) fn polygon(ring: Ring) -> Vec<Sight> {
    let n = ring.0.len();
    let chain = Chain {
        first: 0,
        last: n - 1,
    };
    let table = catalog::Catalog::build(ring, chain);
    let submap = &table.get(chain).map;
    let mut sights = vec![Sight::default(); n];
    for chord in refine::full_map(ring, &table, submap) {
        if let Chord::Beam { level, from, to } = chord {
            match (from, to) {
                (Feature::Vertex, Feature::Edge(e)) | (Feature::Edge(e), Feature::Vertex) => {
                    sights[level].consider(ring, level, e)
                }
                _ => {}
            }
        }
    }
    for (v, s) in sights.iter_mut().enumerate() {
        s.consider(ring, v, n - 1)
    }
    sights
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Point;
    use std::collections::BTreeSet;
    fn check_chain(points: &[Point], first: usize, last: usize) {
        let ring = Ring(points);
        let chain = Chain { first, last };
        let table = catalog::Catalog::build(ring, chain);
        let map = &table.get(chain).map;
        let gamma = 1 << ((last - first).next_power_of_two().ilog2().div_ceil(5));
        for cell in &map.cells {
            assert!(cell.len() <= 4);
            for arc in map.arcs(cell) {
                assert!(arc.steps(chain).len() <= gamma)
            }
        }
        let highest = (first..=last).max_by_key(|&i| ring.0[i].height()).unwrap();
        let shooter = oracle::Oracle::new(ring, map, highest);
        for v in 0..points.len() {
            for right in [false, true] {
                let mut actual = shoot::Shot::new(v, Feature::Vertex, right);
                let mut wanted = actual;
                for e in first..last {
                    wanted.edge(ring, e)
                }
                shooter.query(ring, map, &mut actual);
                assert_eq!(
                    actual.target, wanted.target,
                    "{chain:?} vertex={v} right={right}"
                );
            }
        }
        // Augmented fusion endpoints can start in the middle of a foreign
        // edge at a height belonging to the queried chain.
        for e in (0..points.len() - 1).filter(|e| !(first..last).contains(e)) {
            let (lo, hi) = ring.ends(e);
            for v in (first..=last).step_by(3) {
                if !(lo.height() < points[v].height() && points[v].height() < hi.height()) {
                    continue;
                }
                for right in [false, true] {
                    let mut actual = shoot::Shot::new(v, Feature::Edge(e), right);
                    let mut wanted = actual;
                    wanted.consider(ring, Feature::Vertex);
                    for edge in first..last {
                        wanted.edge(ring, edge)
                    }
                    shooter.query(ring, map, &mut actual);
                    assert_eq!(
                        actual.target, wanted.target,
                        "foreign edge {e}, {v}, {right}, {chain:?}"
                    );
                }
            }
        }
        let full: BTreeSet<_> = refine::full_map(ring, &table, map).into_iter().collect();
        assert_complete(ring, chain, full);
        if (last - first).is_multiple_of(11) {
            // Force several refinement rounds on small, exhaustively checked
            // chains by starting with a much coarser conformal map.
            let coarse = map.clone().coarsen(ring, chain.period());
            assert_complete(
                ring,
                chain,
                refine::full_map(ring, &table, &coarse)
                    .into_iter()
                    .collect(),
            );
        }
    }
    fn assert_complete(ring: Ring, chain: Chain, full: BTreeSet<Chord>) {
        let Chain { first, last } = chain;
        for v in first..=last {
            for right in [false, true] {
                let mut ray = shoot::Shot::new(v, Feature::Vertex, right);
                for e in first..last {
                    ray.edge(ring, e)
                }
                assert!(
                    full.contains(&ray.chord()),
                    "missing {:?} {chain:?}",
                    ray.chord()
                );
            }
            if chain.extremum(ring, v) {
                assert!(full.contains(&Chord::Null(v)))
            }
        }
        for c in full {
            if let Chord::Beam { level, from, to } = c {
                let mut ray = shoot::Shot::new(level, from, true);
                ray.consider(ring, Feature::Vertex);
                for e in first..last {
                    ray.edge(ring, e)
                }
                assert_eq!(ray.target, to, "spurious chord {c:?} in {chain:?}")
            }
        }
    }
    #[test]
    fn open_chains_match_all_pairs_visibility() {
        let mut state = 0x12345678u32;
        for n in 4..144 {
            let mut p: Vec<_> = (0..n)
                .map(|i| {
                    state = state.wrapping_mul(1664525).wrapping_add(1013904223);
                    let a = std::f64::consts::TAU * i as f64 / n as f64;
                    let r = 1000.0 + f64::from(state % 1000000);
                    Point::new((a.cos() * r).round() as i32, (a.sin() * r).round() as i32)
                })
                .collect();
            p.rotate_left(state as usize % n);
            check_chain(&p, 0, n - 1);
            check_chain(&p, n / 8, n - 2);
        }
    }
    #[test]
    fn folded_chains() {
        let mut p = Vec::new();
        for i in 0..90 {
            p.push(Point::new(i, if i % 2 == 0 { -15 } else { 15 }))
        }
        for _ in 0..2 {
            check_chain(&p, 0, p.len() - 1);
            p.reverse()
        }
    }
    #[test]
    fn contractions_can_unlock_a_hub_chord() {
        let points: Vec<_> = [2, 7, -1, 5, 0, 8, 3, 9, -2]
            .into_iter()
            .enumerate()
            .map(|(i, y)| Point::new(i as i32, y))
            .collect();
        let ring = Ring(&points);
        let chain = Chain {
            first: 0,
            last: points.len() - 1,
        };
        let table = catalog::Catalog::build(ring, chain);
        let full = submap::Submap::new(
            ring,
            chain,
            refine::full_map(ring, &table, &table.get(chain).map),
        );
        let example = (0..full.sites.len())
            .find(|&i| {
                let a = &full.cells[full.owner[i]];
                let b = &full.cells[full.owner[full.mate[i]]];
                a.len() == 3
                    && b.len() >= 3
                    && a.iter()
                        .any(|&e| full.cells[full.owner[full.mate[e]]].len() == 1)
            })
            .expect("adjacent hubs with a leaf");
        assert_eq!(
            full.chords[full.chord_at[example]],
            Chord::Beam {
                level: 5,
                from: Feature::Edge(7),
                to: Feature::Vertex
            }
        );
        assert!(full.coarsen(ring, chain.period()).chords.is_empty());
    }
    #[test]
    fn maximality_changes_with_granularity() {
        let points = [Point::new(1, 2), Point::new(3, 4)];
        let ring = Ring(&points);
        let chain = Chain { first: 0, last: 1 };
        let map = hierarchy::base(ring, chain);
        assert_eq!(map.chords.len(), 2);
        for i in 0..map.chords.len() {
            let mut retained = map.chords.clone();
            retained.remove(i);
            let relaxed = submap::Submap::new(ring, chain, retained);
            assert!(
                relaxed
                    .cells
                    .iter()
                    .flat_map(|c| relaxed.arcs(c))
                    .any(|a| a.steps(chain).len() > 1)
            )
        }
        assert!(map.coarsen(ring, 2).chords.is_empty());
    }
}
