//! Deterministic end-to-end timings. Input generation and verification are
//! outside the timed region. Run with cargo run --release --example scaling.
use chazelle::{Point, triangulate, triangulate_sweep};
use std::hint::black_box;
use std::time::Instant;

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
        let half = (n / 2) as i32;
        let bottom = |i: i32| {
            let x = i - half / 2;
            x * x - (half / 2) * (half / 2)
        };
        return (0..half)
            .map(|i| Point::new(i, bottom(i)))
            .chain((0..half).rev().map(|i| Point::new(i, 1 - bottom(i))))
            .collect();
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
        for k in 7..=17 {
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
