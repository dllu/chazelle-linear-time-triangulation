//! A submap is a noninterleaving matching on the boundary circle (Lemma 2.2).
//! An arc is owned by the region on its exterior, indexed by its first site.
use super::boundary::{Arc, Chain, Chord, Site};
use crate::geometry::Ring;
use std::collections::VecDeque;

#[derive(Clone)]
pub(super) struct Submap {
    pub chain: Chain,
    pub chords: Vec<Chord>,
    pub sites: Vec<Site>,
    pub cells: Vec<Vec<usize>>,
    pub mate: Vec<usize>,
    pub owner: Vec<usize>,
    pub chord_at: Vec<usize>,
}
impl Submap {
    pub fn new(ring: Ring, chain: Chain, mut chords: Vec<Chord>) -> Self {
        chords.sort_unstable();
        chords.dedup();
        let mut ends: Vec<_> = chords
            .iter()
            .enumerate()
            .flat_map(|(c, &chord)| chain.sites(ring, chord).map(|site| (site, c)))
            .collect();
        ends.sort_unstable();
        debug_assert!(ends.windows(2).all(|p| p[0].0 != p[1].0));
        let (sites, chord_at): (Vec<_>, Vec<_>) = ends.into_iter().unzip();
        let n = sites.len();
        let mut first = vec![None; chords.len()];
        let mut mate = vec![0; n];
        for (i, &c) in chord_at.iter().enumerate() {
            if let Some(j) = first[c].replace(i) {
                mate[i] = j;
                mate[j] = i
            }
        }
        let mut seen = vec![false; n];
        let mut cells = Vec::new();
        for start in 0..n {
            if seen[start] {
                continue;
            }
            let mut i = start;
            let mut cell = Vec::new();
            while !seen[i] {
                seen[i] = true;
                cell.push(i);
                i = mate[(i + 1) % n];
            }
            cells.push(cell);
        }
        if n == 0 {
            cells.push(Vec::new())
        }
        debug_assert_eq!(
            cells.len(),
            chords.len() + 1,
            "crossing chords on the double boundary"
        );
        let mut owner = vec![0; n];
        for (id, cell) in cells.iter().enumerate() {
            for &i in cell {
                owner[i] = id
            }
        }
        Self {
            chain,
            owner,
            chords,
            sites,
            cells,
            mate,
            chord_at,
        }
    }
    pub fn region_after(&self, p: Site) -> usize {
        if self.sites.is_empty() {
            return 0;
        }
        let end = self.sites.partition_point(|&s| s <= p);
        self.owner[if end == 0 {
            self.sites.len() - 1
        } else {
            end - 1
        }]
    }
    pub fn arcs(&self, cell: &[usize]) -> impl Iterator<Item = Arc> + Clone {
        (0..cell.len().max(1)).map(|k| {
            if cell.is_empty() {
                return Arc::full();
            }
            let i = cell[k];
            Arc {
                start: self.sites[i],
                end: self.sites[(i + 1) % self.sites.len()],
            }
        })
    }

    /// Section 3.3. Merge the two incident regions when every resulting arc
    /// satisfies the weight bound. Requeue surviving exits after degree drops.
    pub fn coarsen(mut self, ring: Ring, gamma: usize) -> Self {
        let n = self.sites.len();
        if n == 0 {
            return self;
        }
        let mut owner = self.owner;
        let mut next: Vec<_> = (0..n).map(|i| (i + 1) % n).collect();
        let mut prev: Vec<_> = (0..n).map(|i| (i + n - 1) % n).collect();
        let mut removed = vec![false; self.chords.len()];
        let mut queue: VecDeque<_> = (0..n).filter(|&i| i < self.mate[i]).collect();
        while let Some(i) = queue.pop_front() {
            let c = self.chord_at[i];
            if removed[c] {
                continue;
            }
            let j = self.mate[i];
            let (a, b) = (owner[i], owner[j]);
            if self.cells[a].len().min(self.cells[b].len()) >= 3 {
                continue;
            }
            let exits = self.cells[a]
                .iter()
                .chain(&self.cells[b])
                .copied()
                .filter(|&k| k != i && k != j);
            let degree = self.cells[a].len() + self.cells[b].len() - 2;
            if degree == 0 && self.chain.period() > gamma {
                continue;
            }
            if exits.clone().any(|k| {
                let mut end = next[k];
                while end == i || end == j {
                    end = next[end]
                }
                Arc {
                    start: self.sites[k],
                    end: self.sites[end],
                }
                .steps(self.chain)
                .len()
                    > gamma
            }) {
                continue;
            }
            removed[c] = true;
            for k in [i, j] {
                next[prev[k]] = next[k];
                prev[next[k]] = prev[k]
            }
            let mut merged = std::mem::take(&mut self.cells[a]);
            merged.append(&mut self.cells[b]);
            merged.retain(|&k| k != i && k != j);
            for &k in &merged {
                owner[k] = a;
                queue.push_back(k)
            }
            self.cells[a] = merged;
        }
        let chords = self
            .chords
            .into_iter()
            .zip(removed)
            .filter_map(|(c, gone)| (!gone).then_some(c))
            .collect();
        Self::new(ring, self.chain, chords)
    }
}
