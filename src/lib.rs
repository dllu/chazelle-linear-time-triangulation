//! Triangulation of simple, hole-free polygons with exact integral predicates.
//!
//! ```
//! use chazelle::{Point, triangulate};
//!
//! let boundary = [Point::new(0, 0), Point::new(4, 0), Point::new(0, 3)];
//! let faces = triangulate(&boundary)?;
//! assert_eq!(faces.len(), 1);
//! # Ok::<(), chazelle::Error>(())
//! ```
mod geometry;
mod sweep;
mod triangles;
mod visibility;

pub use geometry::Point;
use geometry::{Ring, turn};
/// Indices of three input vertices, in counterclockwise order.
pub type Triangle = [usize; 3];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
/// An input check or a subdivision invariant failed.
pub enum Error {
    TooFewVertices,
    DuplicateNeighbor,
    ZeroArea,
    InvalidSubdivision,
}
impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::TooFewVertices => "a polygon needs at least three vertices",
            Self::DuplicateNeighbor => "neighboring vertices coincide",
            Self::ZeroArea => "the polygon has zero signed area",
            Self::InvalidSubdivision => "inconsistent visibility subdivision",
        })
    }
}
impl std::error::Error for Error {}

/// Triangulate a simple boundary, in either winding, with no repeated closing
/// vertex. Output triples index the input and are counterclockwise. Simplicity
/// is a precondition. Time and additional space are O(n).
pub fn triangulate(polygon: &[Point]) -> Result<Vec<Triangle>, Error> {
    run(polygon, visibility::polygon)
}
/// Independent O(n log n) sweep for comparison, using the same output stage.
pub fn triangulate_sweep(polygon: &[Point]) -> Result<Vec<Triangle>, Error> {
    run(polygon, sweep::visibility)
}
fn run(p: &[Point], map: fn(Ring) -> Vec<visibility::Sight>) -> Result<Vec<Triangle>, Error> {
    if p.len() < 3 {
        return Err(Error::TooFewVertices);
    }
    let ring = Ring(p);
    let mut area = 0;
    for i in 0..p.len() {
        let next = ring.next(i);
        if p[i] == p[next] {
            return Err(Error::DuplicateNeighbor);
        }
        area += turn(p[0], p[i], p[next]);
    }
    if area == 0 {
        return Err(Error::ZeroArea);
    }
    let reversed;
    let ring = if area < 0 {
        reversed = p.iter().copied().rev().collect::<Vec<_>>();
        Ring(&reversed)
    } else {
        ring
    };
    let mut triangles = triangles::extract(ring, &map(ring))?;
    if area < 0 {
        for t in &mut triangles {
            for v in t {
                *v = p.len() - 1 - *v
            }
        }
    }
    Ok(triangles)
}
