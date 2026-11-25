//! Read one integer x y pair per line and emit zero-based triangle indices.
use chazelle::{Point, triangulate};
use std::io::{self, BufRead, Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut vertices = Vec::new();
    for line in io::stdin().lock().lines() {
        let line = line?;
        let mut fields = line.split_whitespace();
        let Some(x) = fields.next() else { continue };
        let y = fields.next().ok_or("expected x y")?;
        if fields.next().is_some() {
            return Err("expected exactly two coordinates".into());
        }
        vertices.push(Point::new(x.parse()?, y.parse()?));
    }
    let mut out = io::BufWriter::new(io::stdout().lock());
    for [a, b, c] in triangulate(&vertices)? {
        writeln!(out, "{a} {b} {c}")?
    }
    out.flush()?;
    Ok(())
}
