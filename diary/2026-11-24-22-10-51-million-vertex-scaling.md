# Scaling through 1,048,576 vertices

Replaced the convex benchmark’s parabolic boundary with a polygon built from primitive integer edge directions. The original coordinates would overflow i32 at the larger sizes. A regression test checks that the million-vertex boundary makes a strict left turn at every vertex.

Reran the release benchmark on core 5 for all 14 powers of two from 128 to 1,048,576 on each of the three families. The submap and sweep triangulations agreed at all 42 inputs. Each chart point uses the median of five runs after warmup. At the largest star, the medians were 20.68 and 0.91 µs per vertex, respectively.
