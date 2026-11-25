//! The binary edge-centroid decomposition in §2.3. A branch records the
//! boundary interval of one component after cutting its separator chord.
use super::submap::Submap;

pub(super) enum Node {
    Leaf(usize),
    Branch {
        end: usize,
        inside: Box<Node>,
        outside: Box<Node>,
    },
}
impl Node {
    pub fn build(map: &Submap) -> Self {
        let mut cut = vec![false; map.chords.len()];
        let mut parent = vec![0; map.cells.len()];
        let mut size = vec![0; map.cells.len()];
        Self::component(map, 0, &mut cut, &mut parent, &mut size)
    }
    fn component(
        map: &Submap,
        root: usize,
        cut: &mut [bool],
        parent: &mut [usize],
        size: &mut [usize],
    ) -> Self {
        let mut order = vec![root];
        parent[root] = root;
        let mut at = 0;
        while at < order.len() {
            let v = order[at];
            size[v] = 1;
            for &e in &map.cells[v] {
                let w = map.owner[map.mate[e]];
                if w != parent[v] && !cut[map.chord_at[e]] {
                    parent[w] = v;
                    order.push(w)
                }
            }
            at += 1;
        }
        if order.len() == 1 {
            return Self::Leaf(root);
        }
        let mut best = (usize::MAX, 0);
        for &v in order[1..].iter().rev() {
            size[parent[v]] += size[v];
            best = best.min((size[v].max(order.len() - size[v]), v));
        }
        let v = best.1;
        let e = *map.cells[v]
            .iter()
            .find(|&&e| map.owner[map.mate[e]] == parent[v])
            .unwrap();
        let end = e.min(map.mate[e]);
        cut[map.chord_at[end]] = true;
        Self::Branch {
            end,
            inside: Box::new(Self::component(map, map.owner[end], cut, parent, size)),
            outside: Box::new(Self::component(
                map,
                map.owner[map.mate[end]],
                cut,
                parent,
                size,
            )),
        }
    }
}
