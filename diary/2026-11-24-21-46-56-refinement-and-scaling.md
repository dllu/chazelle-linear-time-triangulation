# Refinement and scaling checks

The down phase now cuts a long region arc into aligned pieces from the up-phase table. It resets the larger piece maps to the next granularity and resolves their region-facing ports by local shooting. Tests force several refinement rounds on small inputs so the entire resulting map can be compared with exhaustive visibility.

The final scaling run covers convex, star and comb polygons at eleven powers of two. At 131,072 vertices the star took 2.52 s through submaps and 79.8 ms through the AVL sweep. The README chart uses microseconds per vertex and separate scales. Regenerating it from the CSV reproduced the SVG byte for byte.

The debug and release suites exercise the corrected contraction rule, centroid shielding, component selection at the meridian, and triangle geometry. A separate deterministic stress run checked identical output against the sweep on 4,000 radial polygons and narrow corridors. The standard binary reads a boundary from stdin with the command in the README.
