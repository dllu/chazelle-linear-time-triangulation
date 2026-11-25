//! Comparison sweep used as an independent O(n log n) reference. The AVL
//! tree holds just the segments crossing the current open horizontal slab.
use crate::geometry::Ring;
use crate::visibility::Sight;
use std::cmp::Ordering;

type Link = Option<Box<Node>>;
struct Node {
    edge: usize,
    child: [Link; 2],
    height: i32,
}
fn height(p: &Link) -> i32 {
    p.as_ref().map_or(0, |v| v.height)
}
impl Node {
    fn fix(&mut self) {
        self.height = 1 + height(&self.child[0]).max(height(&self.child[1]));
    }
    fn rotate(mut self: Box<Self>, side: usize) -> Box<Self> {
        let mut top = self.child[side].take().unwrap();
        self.child[side] = top.child[1 - side].take();
        self.fix();
        top.child[1 - side] = Some(self);
        top.fix();
        top
    }
    fn balance(mut self: Box<Self>) -> Link {
        self.fix();
        for side in 0..2 {
            if height(&self.child[side]) > height(&self.child[1 - side]) + 1 {
                let child = self.child[side].as_mut().unwrap();
                if height(&child.child[1 - side]) > height(&child.child[side]) {
                    self.child[side] = Some(self.child[side].take().unwrap().rotate(1 - side));
                }
                return Some(self.rotate(side));
            }
        }
        Some(self)
    }
}
fn change(node: Link, edge: usize, insert: bool, ring: Ring) -> Link {
    let Some(mut n) = node else {
        assert!(insert);
        return Some(Box::new(Node {
            edge,
            child: [None, None],
            height: 1,
        }));
    };
    let ord = ring.compare_edges(edge, n.edge);
    if ord == Ordering::Equal {
        assert!(!insert);
        if n.child[0].is_none() {
            return n.child[1].take();
        }
        if n.child[1].is_none() {
            return n.child[0].take();
        }
        let mut successor = n.child[1].as_ref().unwrap();
        while let Some(left) = successor.child[0].as_ref() {
            successor = left
        }
        n.edge = successor.edge;
        n.child[1] = change(n.child[1].take(), n.edge, false, ring);
    } else {
        let side = usize::from(ord == Ordering::Greater);
        n.child[side] = change(n.child[side].take(), edge, insert, ring);
    }
    n.balance()
}
pub(crate) fn visibility(ring: Ring) -> Vec<Sight> {
    let mut order: Vec<_> = (0..ring.0.len()).collect();
    order.sort_unstable_by_key(|&v| ring.0[v].height());
    let mut active = None;
    let mut sights = vec![Sight::default(); order.len()];
    for v in order {
        let incident = [
            (ring.prev(v), ring.rising(ring.prev(v))),
            (v, !ring.rising(v)),
        ];
        for (e, ending) in incident {
            if ending {
                active = change(active, e, false, ring);
            }
        }
        let mut p = active.as_ref();
        while let Some(node) = p {
            let left = ring.edge_vs_vertex(node.edge, v) == Ordering::Less;
            sights[v].hit[usize::from(!left)] = Some(node.edge);
            p = node.child[usize::from(left)].as_ref();
        }
        for (e, ending) in incident {
            if !ending {
                active = change(active, e, true, ring);
            }
        }
    }
    sights
}
