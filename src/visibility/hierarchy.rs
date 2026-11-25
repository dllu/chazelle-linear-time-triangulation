//! Canonical submaps (Chazelle §§2–4). The two endpoints of an old chord
//! are independent sources during fusion, including when both lie inside
//! edges. A join also introduces its own ports.
use super::boundary::{Chain, Chord, Feature};
use super::catalog::{Catalog, granularity};
use super::shoot::Shot;
use super::submap::Submap;
use crate::geometry::Ring;

pub(super) fn base(ring: Ring, chain: Chain) -> Submap {
    let mut out = Vec::new();
    for v in chain.first..=chain.last {
        for right in [false, true] {
            let mut shot = Shot::new(v, Feature::Vertex, right);
            for e in chain.first..chain.last {
                shot.edge(ring, e)
            }
            out.push(shot.chord());
        }
        if chain.extremum(ring, v) {
            out.push(Chord::Null(v))
        }
    }
    Submap::new(ring, chain, out).coarsen(ring, granularity(chain.last - chain.first))
}

pub(super) fn merge(
    ring: Ring,
    table: &Catalog,
    left: &Submap,
    right: &Submap,
    gamma: usize,
) -> Submap {
    let chain = Chain {
        first: left.chain.first,
        last: right.chain.last,
    };
    let join = left.chain.last;
    let mut out = Vec::new();
    for (source, other) in [(left, right), (right, left)] {
        super::fusion::pass(ring, table, source, other, join, &mut out);
    }
    for dir in [false, true] {
        let mut shot = Shot::new(join, Feature::Vertex, dir);
        for map in [left, right] {
            let cell = map.region_after(map.chain.site(ring, join, Feature::Vertex, dir));
            for arc in map.arcs(&map.cells[cell]) {
                table.shoot_arc(ring, map.chain, arc, &mut shot)
            }
        }
        out.push(shot.chord());
    }
    if chain.extremum(ring, join) {
        out.push(Chord::Null(join))
    }
    conform(ring, table, Submap::new(ring, chain, out)).coarsen(ring, gamma)
}
fn conform(ring: Ring, table: &Catalog, mut map: Submap) -> Submap {
    loop {
        let mut splits = Vec::new();
        for cell in &map.cells {
            if cell.len() > 4 {
                let arcs: Vec<_> = map.arcs(cell).collect();
                let region = super::region::Region {
                    ring,
                    table,
                    chain: map.chain,
                    arcs: &arcs,
                };
                splits.push(region.split_witness());
            }
        }
        if splits.is_empty() {
            return map;
        }
        map.chords.extend(splits);
        map = Submap::new(ring, map.chain, map.chords)
    }
}
