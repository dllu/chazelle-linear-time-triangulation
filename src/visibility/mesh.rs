//! A triangulated sphere for the separator in §3.4. Crossbars through the
//! thin chain join paired sides between successive projected chord endpoints.
//! Thus two regions incident along a polygon edge are linked by a path of
//! length three, even if the side opposite a long arc has many chord ends.
use super::boundary::{Chord, Feature, Site};
use super::submap::Submap;
use crate::geometry::Ring;

pub(super) struct Mesh {
    /// Each oriented triangle has three vertex indices and three neighbors.
    /// The neighbor in slot i is across the edge i -> (i+1) mod 3.
    pub triangles: Vec<[usize; 3]>,
    pub neighbors: Vec<[usize; 3]>,
    /// Exterior region represented by each original region center.
    pub centers: Vec<usize>,
    /// Regions to remove if a mesh vertex belongs to a separator.
    pub touches: Vec<Vec<usize>>,
}

#[derive(Clone, Copy)]
struct Event {
    site: Site,
    shift: i8,
    pair: usize,
    outside: bool,
    after: usize,
    before: usize,
}

impl Mesh {
    pub fn new(ring: Ring, map: &Submap) -> Self {
        let chain = map.chain;
        // Project the endpoints to the forward strand and introduce a pair of
        // auxiliary sites in the interval immediately following each projection.
        let pair_at = |level: usize, f: Feature| -> (Site, Site) {
            match f {
                Feature::Vertex => (
                    Site::corner(level - chain.first, 0),
                    Site::corner(chain.period() - (level - chain.first), 0),
                ),
                Feature::Edge(e) => {
                    let r = ring.rising(e);
                    (
                        chain.site(ring, level, f, !r),
                        chain.site(ring, level, f, r),
                    )
                }
            }
        };
        let mut knots = vec![pair_at(chain.first, Feature::Vertex)];
        for &c in &map.chords {
            match c {
                Chord::Null(v) => knots.push(pair_at(v, Feature::Vertex)),
                Chord::Beam { level, from, to } => {
                    knots.push(pair_at(level, from));
                    knots.push(pair_at(level, to));
                }
            }
        }
        knots.sort_unstable();
        knots.dedup();
        knots.retain(|&(a, _)| a < Site::corner(chain.last - chain.first, 0));
        let n = map.sites.len();
        let mut events: Vec<_> = map
            .sites
            .iter()
            .enumerate()
            .map(|(i, &site)| Event {
                site,
                shift: 0,
                pair: map.chord_at[i],
                outside: true,
                after: map.owner[i],
                before: map.owner[(i + n - 1) % n],
            })
            .collect();
        for (j, (forward, reverse)) in knots.into_iter().enumerate() {
            for (site, shift) in [(forward, 1), (reverse, -1)] {
                let owner = if shift > 0 {
                    map.region_after(site.last_copy())
                } else {
                    let i = map.sites.partition_point(|&p| p < site);
                    if n == 0 {
                        0
                    } else {
                        map.owner[(i + n - 1) % n]
                    }
                };
                events.push(Event {
                    site: if shift > 0 { site.last_copy() } else { site },
                    shift,
                    pair: map.chords.len() + j,
                    outside: false,
                    after: owner,
                    before: owner,
                });
            }
        }
        events.sort_unstable_by_key(|e| (e.site, e.shift));
        let m = events.len();
        let mut mates = vec![0; m];
        let mut at = vec![None; m / 2];
        for (i, e) in events.iter().enumerate() {
            if let Some(j) = at[e.pair].replace(i) {
                mates[i] = j;
                mates[j] = i
            }
        }
        // Dart slots at vertex i are forward boundary, backward boundary,
        // and chord. Their rotations depend only on the side of the chord.
        let twin = |d: usize| match d % 3 {
            0 => 3 * ((d / 3 + 1) % m) + 1,
            1 => 3 * ((d / 3 + m - 1) % m),
            _ => 3 * mates[d / 3] + 2,
        };
        let next = |d: usize| {
            let t = twin(d);
            let rotation = if events[t / 3].outside {
                [1, 2, 0]
            } else {
                [2, 0, 1]
            };
            3 * (t / 3) + rotation[t % 3]
        };
        let mut centers = vec![usize::MAX; map.cells.len()];
        let mut touches: Vec<_> = events
            .iter()
            .map(|e| {
                if e.after == e.before {
                    vec![e.after]
                } else {
                    vec![e.after, e.before]
                }
            })
            .collect();
        let mut triangles = vec![[0; 3]; 3 * m];
        let mut neighbors = triangles.clone();
        let mut used = vec![false; 3 * m];
        let mut face = Vec::new();
        for start in 0..3 * m {
            if used[start] {
                continue;
            }
            let center = touches.len();
            face.clear();
            let mut d = start;
            while !used[d] {
                used[d] = true;
                face.push(d);
                d = next(d)
            }
            // The outer faces always contain a forward boundary dart.
            let regions = if let Some(&dart) = face.iter().find(|&&d| d % 3 == 0) {
                let region = events[dart / 3].after;
                debug_assert_eq!(centers[region], usize::MAX);
                centers[region] = center;
                vec![region]
            } else {
                let mut incident: Vec<_> = face
                    .iter()
                    .flat_map(|&d| touches[d / 3].iter().copied())
                    .collect();
                incident.sort_unstable();
                incident.dedup();
                incident
            };
            touches.push(regions);
            for (k, &dart) in face.iter().enumerate() {
                triangles[dart] = [dart / 3, twin(dart) / 3, center];
                neighbors[dart] = [
                    twin(dart),
                    face[(k + 1) % face.len()],
                    face[(k + face.len() - 1) % face.len()],
                ];
            }
        }
        debug_assert!(centers.iter().all(|&c| c != usize::MAX));
        Self {
            triangles,
            neighbors,
            centers,
            touches,
        }
    }
}
