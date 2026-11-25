# From open-chain visibility to triangles

Started with the paper’s spherical open-chain map. Distinguishing the two copies at a turning vertex made it possible to encode a null chord, a wrap chord, and a finite ray hit with the same boundary traversal. The independent AVL sweep and the all-pairs map checks were useful at different scales. The all-pairs checks also see the exterior chords that triangle counts miss.

During fusion a star-shaped input exposed a missed region transition. A target exit landed exactly at the end of the current source arc on the join level. Recording that event on the inclusive end, then choosing the region entered by the next step, fixed the map.

Implemented the canonical grade table, centroid search, planar separator ray shooter and region-local down phase. DEVIATIONS.md records the separator and refinement arguments. The sparse-region tests exercise the shielding lemma on much larger arcs than appear in the smallest end-to-end tests. Spiral corridors, equal-height combs and untangled tours exercise the finished triangulator. The release benchmark compares identical triangle output on two families through 131,072 vertices.
