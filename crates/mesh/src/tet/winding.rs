//! How many times the tetrahedra round each edge wind about it (#685). Positive tetrahedra
//! round an interior edge fill exactly one turn; round a boundary edge they leave a gap. A star
//! of positive tetrahedra can instead wind twice round an edge while every face is held at most
//! twice and its volume equals the winding-weighted volume its boundary encloses, which none of
//! the gate's other checks sees.
//!
//! The count is an integer, from exact predicates alone. Round the edge `ab`, each tetrahedron
//! `abcd` (its corners ordered so it is positive) turns from the half-plane through `c` to the
//! one through `d`, by less than half a turn. Those steps chain into the edge's link: cycles
//! round an interior edge, and one path from one boundary face to the other round a boundary
//! edge (the gate has already refused an edge in more or fewer than two boundary faces). The
//! winding number of a closed chain is the number of its steps that cross a fixed half-plane
//! through the edge, each step counted on its half-open interval `(c, d]` (the crossing count of
//! the point-in-polygon winding number, Hormann & Agathos, *The point in polygon problem for
//! arbitrary polygons*, CGTA 20(3):131–144, 2001). Measured from the half-plane through a vertex
//! of the link (for a boundary edge, where its path starts), an interior edge's steps cross it
//! once and a boundary edge's never; a star wound twice crosses it twice, and a boundary fan
//! that closes onto itself crosses it once. Whether a step crosses is decided by the signs of
//! Shewchuk's exact `orient3d` ([`crate::predicates`]), so there is no tolerance.
//!
//! The work is two predicates per edge of each tetrahedron and one sort of the six edge
//! incidences per tetrahedron. The certificate is local: it does not prove that unrelated parts
//! of a mesh do not intersect.

use fieldcraft_model::progress;

use super::sorted2;
use crate::predicates::{Point3, Sign, snap_to_domain, tet_volume_sign};
use crate::{Error, Result};

/// Each edge of a tetrahedron with its other two corners, the four an even permutation of the
/// corners: a positive tetrahedron is positive in this order too.
const EDGES: [[usize; 4]; 6] = [[0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2], [1, 2, 0, 3], [1, 3, 2, 0], [2, 3, 0, 1]];

/// One tetrahedron round an edge: the edge (sorted) and the corners it turns from and to.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Step {
    edge: (u32, u32),
    from: u32,
    to: u32,
}

/// An error unless the tetrahedra wind exactly once round each interior edge and less than once
/// round each boundary edge. Only called once the rest of the gate passed: every tetrahedron
/// positive, no face held by more than two, every boundary edge in exactly two boundary faces.
pub(super) fn check(points: &[[f64; 3]], tets: &[[u32; 4]]) -> Result<()> {
    let mut steps = Vec::with_capacity(tets.len().saturating_mul(6));
    for t in tets {
        for [i, j, k, l] in EDGES {
            let [a, b, c, d] = [i, j, k, l].map(|n| t.get(n).copied().unwrap_or_default());
            // Swapping the edge's ends reverses the turn.
            let (from, to) = if a < b { (c, d) } else { (d, c) };
            steps.push(Step { edge: sorted2(a, b), from, to });
        }
    }
    steps.sort_unstable();
    let point = |v: u32| -> Result<Point3> {
        let p = points.get(v as usize).ok_or_else(|| Error::Failed("gate: vertex id".into()))?;
        Point3::try_from(p.map(snap_to_domain))
    };
    for (n, round) in steps.chunk_by(|x, y| x.edge == y.edge).enumerate() {
        if n % 4096 == 0 {
            progress::check()?;
        }
        let Some(first) = round.first() else { continue };
        let (a, b) = first.edge;
        // A path's start is a corner turned from and never to. More than one, or a corner turned
        // from or to twice, is not one chain (only a mesh the gate refused has that).
        let mut froms: Vec<u32> = round.iter().map(|s| s.from).collect();
        let mut tos: Vec<u32> = round.iter().map(|s| s.to).collect();
        froms.sort_unstable();
        tos.sort_unstable();
        let unique = |v: &[u32]| v.windows(2).all(|w| w.first() != w.get(1));
        let starts: Vec<u32> = froms.iter().copied().filter(|f| tos.binary_search(f).is_err()).collect();
        let chained = unique(&froms) && unique(&tos) && starts.len() <= 1;
        let start = starts.first().copied();
        let boundary = start.is_some();
        let crossings = if chained {
            let (pa, pb) = (point(a)?, point(b)?);
            let r = point(start.unwrap_or(first.from))?;
            let mut crossings = 0usize;
            for s in round {
                let (c, d) = (point(s.from)?, point(s.to)?);
                let side = |x: &Point3, y: &Point3| tet_volume_sign(&pa, &pb, x, y);
                let crosses = match (side(&r, &d), side(&c, &r)) {
                    // `d` on the plane through the edge and `r`: crossing when on `r`'s side of
                    // the edge (`c` then sees both on its positive side), else half a turn away.
                    (Sign::Zero, from_c) => from_c == Sign::Positive,
                    // `c` on that plane and `d` not: the step starts on `r`'s half-plane (the
                    // interval is open there) or half a turn from it, and turns less than half a
                    // turn.
                    (_, Sign::Zero) => false,
                    (to_d, from_c) => to_d == Sign::Positive && from_c == Sign::Positive,
                };
                crossings += usize::from(crosses);
            }
            Some(crossings)
        } else {
            None
        };
        if crossings != Some(usize::from(!boundary)) {
            let kind = if boundary { "boundary" } else { "interior" };
            let expected = if boundary { "less than one turn" } else { "exactly one turn" };
            let found = match crossings {
                Some(_) if boundary => "wind a whole turn or more".to_string(),
                Some(n) => format!("wind {n} turns"),
                None => "do not chain into one fan".to_string(),
            };
            return Err(Error::failed(crate::message![
                format!("mesh validity gate: the {} tetrahedra round {kind} edge {a}–{b} near ", round.len()),
                crate::message::point(points.get(a as usize).copied().unwrap_or_default()),
                format!(" {found} where they must fill {expected}: the mesh overlaps itself there")
            ]));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
