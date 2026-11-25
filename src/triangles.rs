//! Fournier–Montuno Algorithms 2 and 3. All cuts have their endpoints
//! on the original boundary circle, so twin indices suffice to walk faces.
use crate::geometry::{Ring, turn};
use crate::visibility::Sight;
use crate::{Error, Triangle};

#[derive(Clone, Copy)]
struct Event {
    vertex: Option<usize>,
    chord: Option<usize>,
}

/// Visit the faces of a disk slit by a noncrossing matching. A face follows
/// the boundary to the next event, crosses its chord, then continues forward.
fn faces(
    events: &[Event],
    cuts: usize,
    mut visit: impl FnMut(&[usize]) -> Result<(), Error>,
) -> Result<(), Error> {
    let mut first = vec![None; cuts];
    let mut twin: Vec<_> = (0..events.len()).collect();
    for (i, event) in events.iter().enumerate() {
        if let Some(c) = event.chord
            && let Some(j) = first[c].replace(i)
        {
            twin[i] = j;
            twin[j] = i;
        }
    }
    let mut seen = vec![false; events.len()];
    let mut vertices = Vec::new();
    for start in 0..events.len() {
        if seen[start] {
            continue;
        }
        vertices.clear();
        let mut i = start;
        while !seen[i] {
            seen[i] = true;
            if let Some(v) = events[i].vertex
                && vertices.last() != Some(&v)
            {
                vertices.push(v)
            }
            i = (i + 1) % events.len();
            // The two sides of a chord can carry different vertex labels
            // when one end is a Steiner point on a polygon edge.
            if let Some(v) = events[i].vertex
                && vertices.last() != Some(&v)
            {
                vertices.push(v)
            }
            i = twin[i];
        }
        if vertices.len() > 1 && vertices.first() == vertices.last() {
            vertices.pop();
        }
        visit(&vertices)?;
    }
    Ok(())
}

pub(crate) fn extract(ring: Ring, sights: &[Sight]) -> Result<Vec<Triangle>, Error> {
    let n = sights.len();
    let mut from = vec![Vec::new(); n];
    let mut before = vec![Vec::new(); n];
    let mut after = vec![Vec::new(); n];
    let mut cuts = 0;
    for v in (0..n).rev() {
        let a = ring.prev(v);
        let incoming = usize::from(!ring.rising(a));
        let outgoing = usize::from(!ring.rising(v));
        let extremum = incoming != outgoing;
        if extremum && turn(ring.0[a], ring.0[v], ring.0[ring.next(v)]) > 0 {
            continue;
        }
        let directions = [incoming, outgoing];
        for &side in &directions[..1 + usize::from(extremum)] {
            let edge = sights[v].hit[side].ok_or(Error::InvalidSubdivision)?;
            from[v].push(cuts);
            // On the target edge the sources occur in reverse cyclic order:
            // interleaving two pairs of boundary endpoints would cross their rays.
            (if v < edge { &mut before } else { &mut after })[edge].push(cuts);
            cuts += 1;
        }
    }
    let mut events = Vec::with_capacity(n + 2 * cuts);
    for v in 0..n {
        if from[v].is_empty() {
            events.push(Event {
                vertex: Some(v),
                chord: None,
            })
        }
        events.extend(from[v].iter().map(|&c| Event {
            vertex: Some(v),
            chord: Some(c),
        }));
        events.extend(before[v].iter().chain(&after[v]).map(|&c| Event {
            vertex: None,
            chord: Some(c),
        }));
    }
    let mut diagonals = vec![Vec::new(); n];
    let mut count = 0;
    faces(&events, cuts, |corners| {
        let &[a, b] = corners else {
            return Err(Error::InvalidSubdivision);
        };
        if ring.next(a) != b && ring.next(b) != a {
            diagonals[a].push((b, count));
            diagonals[b].push((a, count));
            count += 1;
        }
        Ok(())
    })?;
    events.clear();
    for (v, incident) in diagonals.iter_mut().enumerate() {
        // At most four trapezoids terminate at any one vertex.
        incident.sort_unstable_by_key(|&(w, _)| (v + n - w) % n);
        if incident.is_empty() {
            events.push(Event {
                vertex: Some(v),
                chord: None,
            });
        }
        events.extend(incident.iter().map(|&(_, c)| Event {
            vertex: Some(v),
            chord: Some(c),
        }));
    }
    let mut triangles = Vec::with_capacity(n - 2);
    let mut stack = Vec::new();
    faces(&events, count, |piece| {
        let k = piece.len();
        if k < 3 {
            return Err(Error::InvalidSubdivision);
        }
        let low = (0..k).min_by_key(|&i| ring.0[piece[i]].height()).unwrap();
        let high = (0..k).max_by_key(|&i| ring.0[piece[i]].height()).unwrap();
        let start = if (low + 1) % k == high {
            high
        } else if (high + 1) % k == low {
            low
        } else {
            return Err(Error::InvalidSubdivision);
        };
        stack.clear();
        for offset in 0..k {
            let c = piece[(start + offset) % k];
            while stack.len() >= 2 {
                let b = *stack.last().unwrap();
                let a = stack[stack.len() - 2];
                if turn(ring.0[a], ring.0[b], ring.0[c]) <= 0 {
                    break;
                }
                triangles.push([a, b, c]);
                stack.pop();
            }
            stack.push(c);
        }
        if stack.len() != 2 {
            return Err(Error::InvalidSubdivision);
        }
        Ok(())
    })?;
    if triangles.len() != n - 2 {
        return Err(Error::InvalidSubdivision);
    }
    Ok(triangles)
}
