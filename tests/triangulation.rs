use chazelle::{Point, Triangle, triangulate, triangulate_sweep};
use std::collections::BTreeMap;

fn cross(a: Point, b: Point, c: Point) -> i128 {
    (i128::from(b.x) - i128::from(a.x)) * (i128::from(c.y) - i128::from(a.y))
        - (i128::from(b.y) - i128::from(a.y)) * (i128::from(c.x) - i128::from(a.x))
}
fn verify(p: &[Point], tris: &[Triangle]) {
    let n = p.len();
    assert_eq!(tris.len(), n - 2);
    let mut coverage = BTreeMap::new();
    let mut area = 0;
    for &[a, b, c] in tris {
        let signed = cross(p[a], p[b], p[c]);
        assert!(signed > 0);
        area += signed;
        for (x, y) in [(a, b), (b, c), (c, a)] {
            *coverage.entry((x, y)).or_insert(0usize) += 1;
        }
    }
    let original: i128 = (0..n).map(|i| cross(p[0], p[i], p[(i + 1) % n])).sum();
    assert_eq!(area, original.abs());
    for i in 0..n {
        let j = (i + 1) % n;
        let edge = if original > 0 { (i, j) } else { (j, i) };
        assert_eq!(coverage.remove(&edge), Some(1));
    }
    for (&(a, b), &v) in &coverage {
        assert_eq!(v, 1);
        assert_eq!(coverage.get(&(b, a)), Some(&1));
    }
    // Even a diagonal passing through a nonincident boundary vertex is invalid.
    for &(a, b) in coverage.keys() {
        for i in 0..n {
            let j = (i + 1) % n;
            if [a, b].contains(&i) || [a, b].contains(&j) {
                continue;
            }
            let s = cross(p[a], p[b], p[i]);
            let t = cross(p[a], p[b], p[j]);
            if s == 0 {
                assert!(
                    p[i].x < p[a].x.min(p[b].x)
                        || p[i].x > p[a].x.max(p[b].x)
                        || p[i].y < p[a].y.min(p[b].y)
                        || p[i].y > p[a].y.max(p[b].y)
                );
            }
            let u = cross(p[i], p[j], p[a]);
            let v = cross(p[i], p[j], p[b]);
            assert!(!(s.signum() * t.signum() < 0 && u.signum() * v.signum() < 0));
        }
    }
}
fn check(p: &[Point]) {
    for f in [triangulate, triangulate_sweep] {
        verify(p, &f(p).unwrap_or_else(|e| panic!("{p:?}: {e}")));
    }
}
#[test]
fn degeneracies_and_paper_counterexample() {
    for coords in [
        vec![(0, 0), (1, 0), (2, 0), (3, 0), (3, 2), (1, 2), (0, 2)],
        vec![(0, 8), (6, 6), (4, 4), (1, 2), (0, 0)],
        vec![
            (0, 0),
            (8, 0),
            (8, 5),
            (6, 5),
            (6, 1),
            (5, 1),
            (5, 4),
            (4, 4),
            (4, 1),
            (3, 1),
            (3, 5),
            (0, 5),
        ],
        vec![
            (i32::MIN, i32::MIN),
            (i32::MAX, i32::MIN),
            (i32::MAX, i32::MAX),
            (i32::MIN, i32::MAX),
        ],
    ] {
        let mut p: Vec<_> = coords.into_iter().map(|(x, y)| Point::new(x, y)).collect();
        for _ in 0..2 {
            for _ in 0..p.len() {
                p.rotate_left(1);
                check(&p)
            }
            p.reverse();
        }
    }
}
#[test]
fn generated_stars() {
    let mut state = 321u32;
    for n in 3..90 {
        for _ in 0..3 {
            let p: Vec<_> = (0..n)
                .map(|i| {
                    state = state.wrapping_mul(1664525).wrapping_add(1013904223);
                    let r = 1000. + (state % 30000) as f64;
                    let a = std::f64::consts::TAU * (i as f64) / (n as f64);
                    Point::new((a.cos() * r).round() as i32, (a.sin() * r).round() as i32)
                })
                .collect();
            check(&p);
        }
    }
}

#[test]
fn bent_spiral_corridors() {
    // A thin open corridor following successive sides of shrinking squares.
    // Nonneighboring turns are separated by at least six lattice units.
    for turns in 2..32 {
        let (mut x, mut y) = (0i32, 0i32);
        let mut path = vec![(x, y)];
        for j in 0..turns {
            let (dx, dy) = [(1, 0), (0, 1), (-1, 0), (0, -1)][j % 4];
            let distance = 8 * (turns - j) as i32;
            x += distance * dx;
            y += distance * dy;
            path.push((x, y))
        }
        let normals: Vec<_> = path
            .windows(2)
            .map(|w| ((w[0].1 - w[1].1).signum(), (w[1].0 - w[0].0).signum()))
            .collect();
        let mut p = Vec::new();
        for side in [1, -1] {
            for i in 0..=turns {
                let j = if side == 1 { i } else { turns - i };
                let (mut dx, mut dy) = (0, 0);
                if j > 0 {
                    dx += normals[j - 1].0;
                    dy += normals[j - 1].1
                }
                if j < turns {
                    dx += normals[j].0;
                    dy += normals[j].1
                }
                p.push(Point::new(path[j].0 + side * dx, path[j].1 + side * dy));
            }
        }
        for _ in 0..3 {
            p.rotate_left(turns / 2);
            check(&p)
        }
    }
}

#[test]
fn large_comb_matches_sweep() {
    let n = 4099;
    let mut p: Vec<_> = (0..n)
        .map(|i| Point::new(i, ((i % 2) * 100000) + (137 * i) % 73))
        .collect();
    p.extend([Point::new(n - 1, -1000), Point::new(0, -1000)]);
    for i in 0..3 {
        p.rotate_left(123);
        let a = triangulate(&p).unwrap();
        let b = triangulate_sweep(&p).unwrap();
        assert_eq!(a, b, "rotation {i}");
    }
}

#[test]
fn untangled_tours() {
    let mut state = 17u32;
    let mut sample = || {
        state = state.wrapping_mul(1664525).wrapping_add(1013904223);
        (state >> 3) as i32
    };
    for n in 5..82 {
        let mut p: Vec<_> = (0..n).map(|_| Point::new(sample(), sample())).collect();
        // Removing a proper crossing by reversing the intervening path
        // strictly reduces tour length, so this produces a simple boundary.
        loop {
            let mut changed = false;
            for i in 0..n - 2 {
                for j in i + 2..n {
                    let (a, b, c, d) = (p[i], p[i + 1], p[j], p[(j + 1) % n]);
                    if cross(a, b, c).signum() * cross(a, b, d).signum() < 0
                        && cross(c, d, a).signum() * cross(c, d, b).signum() < 0
                    {
                        p[i + 1..=j].reverse();
                        changed = true;
                    }
                }
            }
            if !changed {
                break;
            }
        }
        check(&p);
    }
}

#[test]
fn input_errors() {
    use chazelle::Error;
    let a = Point::new(1, 2);
    let b = Point::new(2, 3);
    for f in [triangulate, triangulate_sweep] {
        assert_eq!(f(&[a, b]), Err(Error::TooFewVertices));
        assert_eq!(f(&[a, b, a]), Err(Error::DuplicateNeighbor));
        assert_eq!(f(&[a, b, Point::new(-1, 0)]), Err(Error::ZeroArea));
    }
}
