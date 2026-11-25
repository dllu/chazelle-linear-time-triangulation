//! Region-local down phase (§4.2). Reuse the cached submaps of the region's
//! arc pieces to discover finer exit chords and restore degree four.
//! Each child is described directly by arcs of the original double boundary.
use super::boundary::{Arc, Chord, Site};
use super::catalog::{Catalog, grade};
use super::region::Region;
use super::submap::Submap;
use crate::geometry::Ring;

fn between(arcs: &[Arc], i: usize, p: Site, j: usize, q: Site) -> Vec<Arc> {
    let origin = arcs[i].start;
    let key = |s: Site| (s < origin, s);
    if i == j && key(p) < key(q) {
        return vec![Arc { start: p, end: q }];
    }
    let mut path = vec![Arc {
        start: p,
        end: arcs[i].end,
    }];
    let mut k = (i + 1) % arcs.len();
    while k != j {
        path.push(arcs[k]);
        k = (k + 1) % arcs.len()
    }
    path.push(Arc {
        start: arcs[j].start,
        end: q,
    });
    path
}
fn split(arcs: &[Arc], p: Site, q: Site) -> Option<[Vec<Arc>; 2]> {
    if arcs.len() == 1 && arcs[0].start == arcs[0].end {
        return Some([
            vec![Arc { start: p, end: q }],
            vec![Arc { start: q, end: p }],
        ]);
    }
    let i = arcs.iter().position(|a| a.contains(p))?;
    let j = arcs.iter().position(|a| a.contains(q))?;
    Some([between(arcs, i, p, j, q), between(arcs, j, q, i, p)])
}

/// Fuse the region-facing ports of finer piece maps. Retaining the old exits
/// lets the region oracle resolve each port directly, even at a clipped end.
fn finer_chords(region: &Region, gamma: usize) -> Vec<Chord> {
    let mut cuts = Vec::new();
    for (i, &arc) in region.arcs.iter().enumerate() {
        if arc.steps(region.chain).len() <= gamma {
            continue;
        }
        region.table.cut_arc(region.chain, arc, |piece, forward| {
            for v in [piece.first, piece.last] {
                region.at_vertex(v, |c| cuts.push(c))
            }
            let len = piece.last - piece.first;
            if len <= gamma {
                return;
            }
            let reset = region
                .table
                .get(piece)
                .map
                .clone()
                .coarsen(region.ring, gamma);
            let start = if forward { 0 } else { len };
            let strand = Arc {
                start: Site::corner(start, 1),
                end: Site::corner(start + len, 0),
            };
            for c in reset.chords {
                let sites = piece.sites(region.ring, c);
                for (port, site) in sites.into_iter().enumerate() {
                    if !strand.contains(site) {
                        continue;
                    }
                    match c {
                        Chord::Null(v) => region.at_vertex(v, |c| cuts.push(c)),
                        Chord::Beam { level, from, to } => {
                            let right = port == 0;
                            let origin = if right { from } else { to };
                            let parent = region.chain.site(region.ring, level, origin, right);
                            if region.locate(parent) == Some(i) {
                                cuts.push(region.beam(level, origin, right))
                            }
                        }
                    }
                }
            }
        })
    }
    cuts.sort_unstable();
    cuts.dedup();
    cuts
}

pub(super) fn full_map(ring: Ring, table: &Catalog, map: &Submap) -> Vec<Chord> {
    let mut result = map.chords.clone();
    let chain = map.chain;
    let mut pending: Vec<_> = map
        .cells
        .iter()
        .map(|c| map.arcs(c).collect::<Vec<_>>())
        .collect();
    while let Some(arcs) = pending.pop() {
        let region = Region {
            ring,
            table,
            chain,
            arcs: &arcs,
        };
        let weight = arcs.iter().map(|a| a.steps(chain).len()).max().unwrap();
        if weight <= 4 {
            for arc in arcs.iter() {
                let steps = arc.steps(chain);
                for k in steps.start..=steps.end {
                    region.at_vertex(chain.vertex(k % chain.period()), |c| result.push(c));
                }
            }
            continue;
        }
        let exponent = grade(weight);
        // Four factors of (weight / next_weight) in the restoration cost leave
        // slack below the 1/15 saving from the shooting and cutting oracles.
        let next_weight = 1usize << (exponent - exponent.div_ceil(128));
        let cuts = finer_chords(&region, next_weight);
        let mut cells = vec![arcs];
        for c in cuts {
            let [p, q] = chain.sites(ring, c);
            let (i, children) = cells
                .iter()
                .enumerate()
                .find_map(|(i, arcs)| split(arcs, p, q).map(|s| (i, s)))
                .expect("a new visible pair must split a single region");
            cells.swap_remove(i);
            cells.extend(children);
            result.push(c)
        }
        while let Some(cell) = cells.pop() {
            if cell.len() <= 4 {
                debug_assert!(cell.iter().all(|a| a.steps(chain).len() <= next_weight));
                pending.push(cell);
                continue;
            }
            let region = Region {
                ring,
                table,
                chain,
                arcs: &cell,
            };
            let c = region.split_witness();
            let [p, q] = chain.sites(ring, c);
            cells.extend(split(&cell, p, q).unwrap());
            result.push(c);
        }
    }
    result
}
