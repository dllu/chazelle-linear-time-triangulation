//! Deterministic end-to-end timings. Input generation and verification are
//! outside the timed region. Run with cargo run --release --example scaling.
use chazelle::{Point, triangulate, triangulate_sweep};
use std::hint::black_box;
use std::time::Instant;

const MAX_POWER: u32 = 20;

// Link primitive lattice directions in angular order. Quarter turns close
// the polygon, and choosing short steps keeps the coordinates within i32.
fn convex_polygon(n: usize) -> Vec<Point> {
    fn gcd(mut a: i32, mut b: i32) -> i32 {
        while b != 0 {
            (a, b) = (b, a % b);
        }
        a
    }

    assert!(n >= 8 && n.is_multiple_of(8));
    let mut quarter = vec![(1, 0), (1, 1)];
    'directions: for x in 2.. {
        for y in 1..x {
            if quarter.len() == n / 4 {
                break 'directions;
            }
            if gcd(x, y) == 1 {
                quarter.extend([(x, y), (y, x)]);
            }
        }
    }
    quarter.sort_unstable_by(|&(x, y), &(u, v)| (y * u).cmp(&(v * x)));
    let mut vertices = Vec::with_capacity(n);
    let mut point = Point::new(0, 0);
    for _ in 0..4 {
        for step in &mut quarter {
            vertices.push(point);
            let (dx, dy) = *step;
            point.x += dx;
            point.y += dy;
            *step = (-dy, dx);
        }
    }
    vertices
}

fn polygon(n: usize, family: &str) -> Vec<Point> {
    if family == "comb" {
        let teeth = (n - 2) as i32;
        let mut vertices: Vec<_> = (0..teeth)
            .map(|i| Point::new(i, (i % 2) * 100000 + ((137 * i) % 73)))
            .collect();
        vertices.extend([Point::new(teeth - 1, -1000), Point::new(0, -1000)]);
        return vertices;
    }
    if family == "convex" {
        return convex_polygon(n);
    }
    // Varying radii create many simultaneous active edges in the sweep.
    let mut state = 41599u32;
    (0..n)
        .map(|i| {
            let angle = std::f64::consts::TAU * i as f64 / n as f64;
            state = state.wrapping_mul(1664525).wrapping_add(1013904223);
            let r = 4e8 + (state % 500000000) as f64;
            Point::new(
                (r * angle.cos()).round() as i32,
                (r * angle.sin()).round() as i32,
            )
        })
        .collect()
}
fn measure(mut run: impl FnMut()) -> u128 {
    run();
    let mut samples = [0; 5];
    for nanos in &mut samples {
        let start = Instant::now();
        run();
        *nanos = start.elapsed().as_nanos()
    }
    samples.sort_unstable();
    samples[2]
}
fn main() {
    println!("family,vertices,submap_ns,sweep_ns");
    for family in ["convex", "star", "comb"] {
        for k in 7..=MAX_POWER {
            let n = 1 << k;
            let p = polygon(n, family);
            assert_eq!(triangulate(&p).unwrap(), triangulate_sweep(&p).unwrap());
            let submap = measure(|| {
                black_box(triangulate(black_box(&p)).unwrap());
            });
            let sweep = measure(|| {
                black_box(triangulate_sweep(black_box(&p)).unwrap());
            });
            println!("{family},{n},{submap},{sweep}");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn largest_convex_fixture_has_a_strict_left_turn_at_every_vertex() {
        let n = 1 << MAX_POWER;
        let p = convex_polygon(n);
        assert_eq!(p.len(), n);
        for i in 0..n {
            let (a, b, c) = (p[i], p[(i + 1) % n], p[(i + 2) % n]);
            let dx = i128::from(b.x) - i128::from(a.x);
            let dy = i128::from(b.y) - i128::from(a.y);
            assert!(
                dx * (i128::from(c.y) - i128::from(b.y)) > dy * (i128::from(c.x) - i128::from(b.x))
            );
        }
    }
}
