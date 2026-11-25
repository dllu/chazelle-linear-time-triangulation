//! A separator division for Lemma 3.6. In a triangulated sphere, the dual
//! edges complementary to a breadth-first spanning tree form a tree. Its
//! face bags (the three root paths) form a tree decomposition. Slicing the
//! bags into BFS layers bounds the work of extracting each separator.
use super::mesh::Mesh;
use std::collections::VecDeque;

struct RootedTree {
    up: Vec<Vec<usize>>,
    depth: Vec<usize>,
    tin: Vec<usize>,
    end: Vec<usize>,
}
impl RootedTree {
    fn new(parent: Vec<usize>, depth: Vec<usize>, root: usize) -> Self {
        let n = parent.len();
        let mut children = vec![Vec::new(); n];
        for (v, &p) in parent.iter().enumerate() {
            if v != root {
                children[p].push(v)
            }
        }
        let mut tin = vec![0; n];
        let mut end = vec![0; n];
        let mut time = 0;
        let mut stack = vec![(root, false)];
        while let Some((v, finished)) = stack.pop() {
            if finished {
                end[v] = time;
                continue;
            }
            tin[v] = time;
            time += 1;
            stack.push((v, true));
            stack.extend(children[v].iter().map(|&w| (w, false)));
        }
        let mut up = vec![parent];
        while (1usize << up.len()) < n {
            let p = up.last().unwrap();
            up.push(p.iter().map(|&v| p[v]).collect());
        }
        Self {
            up,
            depth,
            tin,
            end,
        }
    }
    fn ancestor(&self, a: usize, b: usize) -> bool {
        self.tin[a] <= self.tin[b] && self.tin[b] < self.end[a]
    }
    fn lca(&self, mut a: usize, b: usize) -> usize {
        if self.ancestor(a, b) {
            return a;
        }
        for p in self.up.iter().rev() {
            if !self.ancestor(p[a], b) {
                a = p[a]
            }
        }
        self.up[0][a]
    }
    fn mark_path(&self, mut v: usize, lo: usize, hi: usize, removed: &mut [bool]) {
        for p in self.up.iter().rev() {
            if self.depth[p[v]] >= hi {
                v = p[v]
            }
        }
        if self.depth[v] >= hi {
            v = self.up[0][v]
        }
        while self.depth[v] >= lo {
            removed[v] = true;
            if self.depth[v] == 0 {
                break;
            }
            v = self.up[0][v]
        }
    }
}

/// Pair each primal tree edge with the dual edge forbidden in the cotree.
fn spanning_trees(mesh: &Mesh, adj: &[Vec<(usize, usize)>]) -> (RootedTree, RootedTree) {
    let n = adj.len();
    let root = mesh.centers[0];
    let mut parent = vec![usize::MAX; n];
    let mut depth = vec![0; n];
    let mut blocked = vec![[false; 3]; mesh.triangles.len()];
    parent[root] = root;
    let mut bfs = VecDeque::from([root]);
    while let Some(v) = bfs.pop_front() {
        for &(w, d) in &adj[v] {
            if parent[w] != usize::MAX {
                continue;
            }
            bfs.push_back(w);
            parent[w] = v;
            depth[w] = depth[v] + 1;
            let (t, j) = (d / 3, d % 3);
            let neighbor = mesh.neighbors[t][j];
            let reverse = (0..3)
                .find(|&i| mesh.neighbors[neighbor][i] == t && mesh.triangles[neighbor][i] == w)
                .unwrap();
            blocked[t][j] = true;
            blocked[neighbor][reverse] = true;
        }
    }
    let primal = RootedTree::new(parent, depth, root);
    let f = mesh.triangles.len();
    let mut parent = vec![usize::MAX; f];
    let mut depth = vec![0; f];
    parent[0] = 0;
    bfs.push_back(0);
    while let Some(t) = bfs.pop_front() {
        for (j, &w) in mesh.neighbors[t].iter().enumerate() {
            if !blocked[t][j] && parent[w] == usize::MAX {
                parent[w] = t;
                depth[w] = depth[t] + 1;
                bfs.push_back(w)
            }
        }
    }
    (primal, RootedTree::new(parent, depth, 0))
}

fn band_separator(
    mesh: &Mesh,
    primal: &RootedTree,
    dual: &RootedTree,
    face: &[usize],
    k: usize,
) -> Vec<bool> {
    // k = ceil(cuberoot(n)), up to a constant factor. Delete the lightest
    // residue class of layers, leaving bands of fewer than k layers.
    let mut sizes = vec![0; k];
    for &d in &primal.depth {
        sizes[d % k] += 1
    }
    let residue = (0..k).min_by_key(|&r| sizes[r]).unwrap();
    let mut removed: Vec<_> = primal.depth.iter().map(|d| d % k == residue).collect();
    let mut bands = vec![Vec::new(); primal.depth.iter().max().unwrap() / k + 2];
    for &v in &mesh.centers {
        let d = primal.depth[v];
        if !removed[v] {
            bands[(d + k - 1 - residue) / k].push(face[v])
        }
    }
    let mut weight = vec![0; mesh.triangles.len()];
    for (band, mut terminals) in bands.into_iter().enumerate() {
        if terminals.len() <= k * k {
            continue;
        }
        for &t in &terminals {
            weight[t] += 1
        }
        terminals.sort_unstable_by_key(|&t| dual.tin[t]);
        terminals.dedup();
        let joins: Vec<_> = terminals.windows(2).map(|w| dual.lca(w[0], w[1])).collect();
        terminals.extend(joins);
        terminals.sort_unstable_by_key(|&t| dual.tin[t]);
        terminals.dedup();
        // The minimal cotree spanning the terminals, with unary paths
        // suppressed. Removing a heavy bag accounts for k² distinct terminals.
        let mut stack = Vec::new();
        let mut edges = Vec::new();
        for &t in &terminals {
            while stack.last().is_some_and(|&p| !dual.ancestor(p, t)) {
                stack.pop();
            }
            edges.push((t, stack.last().copied()));
            stack.push(t);
        }
        let hi = band * k + residue;
        let lo = (hi + 1).saturating_sub(k);
        for (t, p) in edges.into_iter().rev() {
            let w = std::mem::take(&mut weight[t]);
            if w > k * k {
                for &v in &mesh.triangles[t] {
                    primal.mark_path(v, lo, hi, &mut removed)
                }
            } else if let Some(p) = p {
                weight[p] += w
            }
        }
    }
    removed
}

pub(super) struct Division {
    pub separator: Vec<usize>,
    pub groups: Vec<Vec<usize>>,
    pub group_of: Vec<Option<usize>>,
}
impl Division {
    pub fn new(mesh: &Mesh) -> Self {
        let n = mesh.touches.len();
        let mut adj = vec![Vec::new(); n];
        let mut face = vec![0; n];
        for (t, triangle) in mesh.triangles.iter().enumerate() {
            for j in 0..3 {
                let v = triangle[j];
                adj[v].push((triangle[(j + 1) % 3], 3 * t + j));
                face[v] = t;
            }
        }
        let k = 1usize << n.ilog2().div_ceil(3);
        let (primal, dual) = spanning_trees(mesh, &adj);
        let removed = band_separator(mesh, &primal, &dual, &face, k);
        let mut bfs = VecDeque::new();
        let mut separated = vec![false; mesh.centers.len()];
        for (v, &gone) in removed.iter().enumerate() {
            if gone {
                for &r in &mesh.touches[v] {
                    separated[r] = true
                }
            }
        }
        let separator = (0..separated.len()).filter(|&r| separated[r]).collect();
        let mut component = vec![None; n];
        let mut groups = Vec::new();
        let mut group_of = vec![None; separated.len()];
        for (r, &v) in mesh.centers.iter().enumerate() {
            if separated[r] {
                continue;
            }
            let id = *component[v].get_or_insert(groups.len());
            if id == groups.len() {
                groups.push(Vec::new());
                bfs.push_back(v);
                while let Some(u) = bfs.pop_front() {
                    for &(w, _) in &adj[u] {
                        if !removed[w] && component[w].is_none() {
                            component[w] = Some(id);
                            bfs.push_back(w)
                        }
                    }
                }
            }
            groups[id].push(r);
            group_of[r] = Some(id);
        }
        debug_assert!(groups.iter().all(|g| g.len() <= k * k));
        Self {
            separator,
            groups,
            group_of,
        }
    }
}
