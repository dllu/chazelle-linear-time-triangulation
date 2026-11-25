# Chazelle triangulation

A Rust implementation of Chazelle’s deterministic coarse-to-fine visibility algorithm for simple polygons. It accepts either winding, horizontal edges and collinear boundary vertices. Coordinates are signed 32-bit integers. The result contains `n − 2` counterclockwise triangles indexed into the original boundary.

The implementation uses canonical submaps, two directional fusion walks, separator-based horizontal ray shooting, and Fournier–Montuno triangulation. The separator construction and the scheduling of down-phase cuts differ from the published procedures. The bounds and counterexamples are in [DEVIATIONS.md](DEVIATIONS.md). [Chazelle’s paper](https://www.cs.princeton.edu/~chazelle/pubs/polygon-triang.pdf)

## Run

Pass a simple boundary to `chazelle::triangulate(&[Point])`, omitting the repeated closing vertex. The command-line program reads one `x y` pair per line and writes zero-based triangle indices.

```sh
printf '0 0\n4 0\n4 3\n2 1\n0 3\n' | cargo run --quiet --release
cargo test --release
```

Simplicity is a caller precondition. The library checks the vertex count, adjacent duplicates, and signed area. It takes one boundary without holes.

## Scaling

![Median microseconds per vertex for Chazelle submaps and a comparison sweep on convex, star and comb polygons](bench/scaling.svg)

Both paths include triangle extraction. Each panel has its own vertical scale. From 16,384 to 131,072 vertices on the star input, time per vertex grew from 15.6 to 19.2 µs for submaps and from 0.46 to 0.61 µs for the O(n log n) AVL sweep. The sweep was faster throughout this experiment. Timings illustrate finite-size scaling. The worst-case O(n) bound follows from the operation counts in [DEVIATIONS.md](DEVIATIONS.md).

Reproduce the [raw medians](bench/scaling.csv) with `cargo run --release --example scaling > bench/scaling.csv`. Run `python3 bench/plot.py` with Matplotlib installed to redraw the chart. The recorded run used Rust 1.100 on one Cortex-X925 core, one warmup, and five timed runs per point.
