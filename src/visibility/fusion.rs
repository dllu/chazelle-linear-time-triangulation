//! Section 3.1. Walk the sites of one operand in double-boundary order,
//! carrying a region of the other operand. Visibility of a region exit on
//! the current arc is the certificate for changing the carried region.
use super::boundary::{Arc, Chord, Site};
use super::catalog::Catalog;
use super::shoot::Shot;
use super::submap::Submap;
use crate::geometry::{Ring, turn};

fn in_region(
    ring: Ring,
    table: &Catalog,
    map: &Submap,
    region: usize,
    mut shot: Shot,
) -> Option<Shot> {
    let mut arcs = map.arcs(&map.cells[region]);
    for arc in arcs.clone() {
        table.shoot_arc(ring, map.chain, arc, &mut shot)
    }
    if shot.origin == shot.target {
        return None;
    }
    let site = map.chain.site(ring, shot.level, shot.target, !shot.right);
    arcs.any(|a| a.covers(site)).then_some(shot)
}

/// The three local sectors of a chain endpoint are separated by its two
/// horizontal ports. A missing port just identifies the neighboring sectors.
fn starting_region(ring: Ring, other: &Submap, join: usize, source_neighbor: usize) -> usize {
    let c = other.chain;
    let step = join - c.first;
    let next = c.vertex(step + 1);
    let before = Site::corner(step, 0);
    let after = Site::corner(step, 1);
    let rises = |v: usize| ring.0[join].height() < ring.0[v].height();
    if rises(next) != rises(source_neighbor) {
        return other.region_after(before);
    }
    if turn(ring.0[join], ring.0[next], ring.0[source_neighbor]) > 0 {
        other.region_after(after)
    } else if other.sites.is_empty() {
        0
    } else {
        let i = other.sites.partition_point(|&s| s < before);
        other.owner[(i + other.sites.len() - 1) % other.sites.len()]
    }
}

fn site_shot(ring: Ring, map: &Submap, i: usize) -> Option<Shot> {
    let Chord::Beam { level, from, to } = map.chords[map.chord_at[i]] else {
        return None;
    };
    let right = map.sites[i] == map.chain.site(ring, level, from, true);
    let mut shot = Shot::new(level, if right { from } else { to }, right);
    shot.target = if right { to } else { from };
    Some(shot)
}

struct Crossing {
    site: Site,
    region: usize,
    chord: Chord,
}
fn last_crossing(
    table: &Catalog,
    ring: Ring,
    source: &Submap,
    arc: Arc,
    current: usize,
    other: &Submap,
) -> Option<Crossing> {
    let mut last = None;
    for &exit in &other.cells[current] {
        for i in [exit, other.mate[exit]] {
            let Some(limit) = site_shot(ring, other, i) else {
                continue;
            };
            let mut shot = Shot::new(limit.level, limit.origin, limit.right);
            table.shoot_arc(ring, source.chain, arc, &mut shot);
            if shot.target == shot.origin {
                continue;
            }
            let site = source
                .chain
                .site(ring, shot.level, shot.target, !shot.right);
            if !arc.contains(site) && site != arc.end {
                continue;
            }
            let hit = shot.target;
            shot.consider(ring, limit.target);
            if shot.target != hit {
                continue;
            }
            let mut reverse = Shot::new(shot.level, hit, !shot.right);
            reverse.consider(ring, limit.origin);
            for local in source.arcs(&source.cells[source.region_after(site)]) {
                table.shoot_arc(ring, source.chain, local, &mut reverse)
            }
            if reverse.target != limit.origin {
                continue;
            }
            let c = source.chain;
            let rises =
                ring.0[c.vertex(site.step + 1)].height() > ring.0[c.vertex(site.step)].height();
            let right_port = if limit.right { i } else { other.mate[i] };
            let region = other.owner[if rises {
                other.mate[right_port]
            } else {
                right_port
            }];
            let distance = |p: Site| (p <= arc.start, p);
            if last
                .as_ref()
                .is_none_or(|old: &Crossing| distance(site) > distance(old.site))
            {
                last = Some(Crossing {
                    site,
                    region,
                    chord: shot.chord(),
                });
            }
        }
    }
    last
}

pub(super) fn pass(
    ring: Ring,
    table: &Catalog,
    source: &Submap,
    other: &Submap,
    join: usize,
    out: &mut Vec<Chord>,
) {
    let step = join - source.chain.first;
    let mut pos = Site::corner(step, 1);
    let neighbor = source.chain.vertex(step + 1);
    let mut region = starting_region(ring, other, join, neighbor);
    let start = source.sites.partition_point(|&s| s <= pos);
    for offset in 0..source.sites.len() {
        let i = (start + offset) % source.sites.len();
        let end = source.sites[i];
        let retained = site_shot(ring, source, i);
        if end == Site::corner(step, 0) || end == Site::corner(step, 1) {
            // The merge inserts both join rays explicitly.
            continue;
        }
        loop {
            if let Some(mut limit) = retained {
                let start = Shot::new(limit.level, limit.origin, limit.right);
                if let Some(hit) = in_region(ring, table, other, region, start) {
                    limit.consider(ring, hit.target);
                    if limit.target == hit.target {
                        out.push(limit.chord());
                        if limit.level == join {
                            let site = other.chain.site(ring, join, hit.target, !hit.right);
                            if let Ok(j) = other.sites.binary_search(&site) {
                                let rises = ring.0[source.chain.vertex(end.step + 1)].height()
                                    > ring.0[source.chain.vertex(end.step)].height();
                                let right_port = if hit.right { other.mate[j] } else { j };
                                region = other.owner[if rises {
                                    other.mate[right_port]
                                } else {
                                    right_port
                                }];
                            }
                        }
                        break;
                    }
                }
            }
            if let Some(cross) =
                last_crossing(table, ring, source, Arc { start: pos, end }, region, other)
            {
                pos = cross.site;
                region = cross.region;
                if pos == end {
                    out.push(cross.chord);
                    break;
                }
                continue;
            }
            if let Some(limit) = retained {
                out.push(limit.chord())
            } else {
                out.push(source.chords[source.chord_at[i]])
            }
            break;
        }
        pos = end;
    }
}
