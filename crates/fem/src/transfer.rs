//! Field transfer between non-matching meshes: the value of a nodal field of one mesh (the
//! source, such as a moulding or thermal analysis) at the nodes of another (the target, the
//! structural mesh), by consistent interpolation.
//!
//! Each target point is located in a source solid element, and the element's own shape
//! functions are evaluated there: u(p) = Σ Nᵢ(ξ(p)) uᵢ, with ξ(p) the natural coordinates of the
//! point. This is the consistent (pointwise) interpolation of A. de Boer, A. H. van Zuijlen,
//! H. Bijl, *Review of coupling methods for non-matching meshes*, Comput. Methods Appl. Mech.
//! Engrg. 196(8):1515–1525, 2007, §2.1: it reproduces every field the source's shape functions
//! reproduce, so a uniform field and a field linear in x, y, z transfer exactly (every
//! isoparametric element contains the complete linear polynomials), whatever the two meshes.
//!
//! The natural coordinates are found by Newton's method on x(ξ) = p, the inverse of the
//! isoparametric map x(ξ) = Σ Nᵢ(ξ) xᵢ (Zienkiewicz, Taylor & Zhu 2013, §6.5), started from the
//! reference element's centroid. A point inside no element (a node of a curved surface just
//! outside the faceted source, say) takes the value at the closest point of the closest
//! element found within a tolerance: Newton's ξ clamped to the reference element, which for the
//! small gaps between two meshes of one part is the nearest point of the element's surface.
//! Candidates come from a bucket grid of the element boxes ([`crate::contact::Grid`]), so each
//! point tries a handful of elements.

use fieldcraft_model::ElementKind;

use crate::contact::{Aabb, Grid, bounds};
use crate::dense::inverse3;
use crate::element::{ElementType, Family, element_type};
use crate::{Error, Result};

/// Newton iterations of the inverse map.
const MAX_NEWTON: usize = 30;

/// Newton stops when no natural coordinate moves by more than this.
const NEWTON_TOLERANCE: f64 = 1e-13;

/// How far outside the reference element (in natural coordinates) a point still counts as
/// inside: round-off on a face shared by two elements.
const INSIDE: f64 = 1e-10;

/// A gap below this fraction of the element's size is round-off, not a point outside.
const ROUND_OFF_GAP: f64 = 1e-9; // unit-factors: relative to the element's size

/// Most grid cells of a locator.
const MAX_CELLS: usize = 1 << 22;

/// A source mesh's solid elements, ready to locate points in.
#[derive(Debug, Clone)]
pub struct Locator {
    elements: Vec<(&'static ElementType, Vec<[f64; 3]>)>,
    grid: Grid,
    mean_size: f64,
}

/// Where a point lies in the source mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct Hit {
    /// The element, by its index in the list the locator was made from.
    pub element: usize,
    /// The element's shape functions at the point (at the closest point of the element when the
    /// point is outside it): the weights of its nodal values, summing to one.
    pub weights: Vec<f64>,
    /// How far the point is from the element: zero inside it.
    pub distance: f64,
}

impl Locator {
    /// A locator over solid elements given by kind and node coordinates (in the model's node
    /// order). Only the solid kinds the kernels have are accepted (tet4, tet10, wedge6, pyramid5,
    /// hex8, hex20); a node count that does not match the kind, or a coordinate that is not
    /// finite, is refused.
    pub fn new(elements: Vec<(ElementKind, Vec<[f64; 3]>)>) -> Result<Locator> {
        let mut out = Vec::with_capacity(elements.len());
        let mut boxes: Vec<Aabb> = Vec::with_capacity(elements.len());
        let mut size = 0.0;
        for (i, (kind, xyz)) in elements.into_iter().enumerate() {
            let ty = element_type(kind)?;
            if !matches!(ty.family, Family::Tetrahedron | Family::Wedge | Family::Pyramid | Family::Hexahedron) {
                return Err(Error::NotSupported(format!("locating points in {} elements (solid elements only)", kind.name())));
            }
            if xyz.len() != ty.nodes {
                return Err(Error::Invalid(format!("source element {i}: {} nodes for a {}", xyz.len(), kind.name())));
            }
            if xyz.iter().flatten().any(|c| !c.is_finite()) {
                return Err(Error::Invalid(format!("source element {i}: a coordinate is not finite")));
            }
            let b = bounds(&xyz, 0.0);
            size += (0..3).map(|k| b[1][k] - b[0][k]).fold(0.0, f64::max);
            boxes.push(b);
            out.push((ty, xyz));
        }
        let mean_size = if out.is_empty() { 0.0 } else { size / out.len() as f64 };
        Ok(Locator { grid: Grid::new(&boxes, MAX_CELLS.min(8 * boxes.len().max(1))), elements: out, mean_size })
    }

    /// The number of elements.
    pub fn len(&self) -> usize {
        self.elements.len()
    }

    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }

    /// The mean size of the elements (the largest side of their boxes).
    pub fn mean_size(&self) -> f64 {
        self.mean_size
    }

    /// The element holding `p`, the first in element order when it is on a face they share;
    /// else the closest element within `tolerance` (a length, zero or more); `None` when there
    /// is none or `p` is not finite.
    pub fn locate(&self, p: [f64; 3], tolerance: f64) -> Option<Hit> {
        if p.iter().any(|c| !c.is_finite()) {
            return None;
        }
        let tolerance = if tolerance.is_finite() { tolerance.max(0.0) } else { 0.0 };
        let mut candidates = Vec::new();
        self.grid.candidates_in(&bounds(&[p], tolerance), &mut candidates);
        let mut best: Option<Hit> = None;
        for index in candidates {
            let i = index as usize;
            let Some((ty, xyz)) = self.elements.get(i) else { continue };
            let Some((xi, distance)) = natural(ty, xyz, p) else { continue };
            if distance > tolerance || best.as_ref().is_some_and(|b| b.distance <= distance) {
                continue;
            }
            let weights = (ty.shape)(xi).values().to_vec();
            best = Some(Hit { element: i, weights, distance });
            if distance == 0.0 {
                break;
            }
        }
        best
    }
}

/// The natural coordinates of `p` in an element (clamped to the reference element) and the
/// distance from `p` to the point they map to (zero when `p` is inside); `None` when the map
/// cannot be inverted there (a degenerate element).
pub fn natural(ty: &ElementType, xyz: &[[f64; 3]], p: [f64; 3]) -> Option<([f64; 3], f64)> {
    let position = |xi: [f64; 3]| -> ([f64; 3], crate::shape::Shape) {
        let s = (ty.shape)(xi);
        let mut x = [0.0; 3];
        for (n, c) in s.values().iter().zip(xyz) {
            for k in 0..3 {
                x[k] += n * c[k];
            }
        }
        (x, s)
    };
    let mut xi = centroid(ty.family);
    for _ in 0..MAX_NEWTON {
        let (x, s) = position(xi);
        // J[a][b] = ∂x_b/∂ξ_a; its inverse inv[b][a] = ∂ξ_a/∂x_b.
        let mut j = [[0.0; 3]; 3];
        for (d, c) in s.derivatives().iter().zip(xyz) {
            for a in 0..3 {
                for b in 0..3 {
                    j[a][b] += d[a] * c[b];
                }
            }
        }
        let (inv, _) = inverse3(&j)?;
        let r = [0, 1, 2].map(|k| p[k] - x[k]);
        let step = [0, 1, 2].map(|a| (0..3).map(|b| inv[b][a] * r[b]).sum::<f64>());
        if step.iter().any(|v| !v.is_finite()) {
            return None;
        }
        // A point far outside a curved element can send Newton away: keep ξ near the element.
        xi = [0, 1, 2].map(|a| (xi[a] + step[a]).clamp(-3.0, 3.0));
        if step.iter().all(|v| v.abs() < NEWTON_TOLERANCE) {
            break;
        }
    }
    let clamped = project(ty.family, xi);
    let outside = (0..3).map(|k| (clamped[k] - xi[k]).abs()).fold(0.0, f64::max);
    let (x, _) = position(clamped);
    let gap = (0..3).map(|k| (x[k] - p[k]).powi(2)).sum::<f64>().sqrt();
    if outside > INSIDE {
        return Some((clamped, gap));
    }
    // Inside, unless Newton did not converge (then the clamped point says how far).
    let scale = bounds(xyz, 0.0);
    let size = (0..3).map(|k| scale[1][k] - scale[0][k]).fold(0.0, f64::max);
    Some((clamped, if gap <= ROUND_OFF_GAP * size { 0.0 } else { gap }))
}

/// The centroid of a family's reference element.
fn centroid(family: Family) -> [f64; 3] {
    match family {
        Family::Tetrahedron => [0.25; 3],
        Family::Wedge => [1.0 / 3.0, 1.0 / 3.0, 0.0],
        Family::Pyramid => [0.0, 0.0, 0.25],
        Family::Hexahedron => [0.0; 3],
        Family::Triangle => [1.0 / 3.0, 1.0 / 3.0, 0.0],
        Family::Quadrilateral => [0.0; 3],
    }
}

/// `xi` moved into a family's reference element: the closest point for the tetrahedron, the
/// wedge and the hexahedron; for the pyramid, ζ clamped and then ξ, η to that height's square.
fn project(family: Family, xi: [f64; 3]) -> [f64; 3] {
    let [a, b, c] = xi;
    match family {
        Family::Hexahedron | Family::Quadrilateral => xi.map(|v| v.clamp(-1.0, 1.0)),
        Family::Tetrahedron => simplex(&xi),
        Family::Wedge => {
            let t = simplex(&[a, b]);
            [t[0], t[1], c.clamp(-1.0, 1.0)]
        }
        Family::Triangle => {
            let t = simplex(&[a, b]);
            [t[0], t[1], 0.0]
        }
        Family::Pyramid => {
            let z = c.clamp(0.0, 1.0);
            let h = 1.0 - z;
            [a.clamp(-h, h), b.clamp(-h, h), z]
        }
    }
}

/// The closest point of the corner simplex {x ≥ 0, Σ x ≤ 1} (two or three coordinates) to `x`:
/// the positive part when its sum is at most one, else the projection onto the face Σ x = 1
/// (by the sorted threshold of Held, Wolfe & Crowder 1974).
fn simplex<const N: usize>(x: &[f64; N]) -> [f64; N] {
    let positive = x.map(|v| v.max(0.0));
    if positive.iter().sum::<f64>() <= 1.0 {
        return positive;
    }
    let mut sorted = *x;
    sorted.sort_by(|a, b| b.total_cmp(a));
    let mut sum = 0.0;
    let mut theta = 0.0;
    for (i, v) in sorted.iter().enumerate() {
        sum += v;
        let t = (sum - 1.0) / (i + 1) as f64;
        if v - t > 0.0 {
            theta = t;
        }
    }
    x.map(|v| (v - theta).max(0.0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::ELEMENT_TYPES;

    /// A structured mesh of the box [0, a] × [0, b] × [0, c] in `n` cells per side, of `kind`
    /// (hex8 or hex20 cells, each hex split into six tet4 or tet10, or into two wedge6), with
    /// the interior nodes moved by `wobble` so no face is planar.
    fn box_mesh(kind: ElementKind, n: [usize; 3], size: [f64; 3], wobble: f64) -> Vec<(ElementKind, Vec<[f64; 3]>)> {
        let at = |i: f64, j: f64, k: f64| -> [f64; 3] {
            let u = [i / n[0] as f64, j / n[1] as f64, k / n[2] as f64];
            let inner = u.iter().all(|v| *v > 1e-9 && *v < 1.0 - 1e-9);
            let w = if inner { wobble * (7.0 * u[0] + 3.0 * u[1]).sin() * (5.0 * u[2]).cos() } else { 0.0 };
            [size[0] * u[0] + w, size[1] * u[1] - 0.5 * w, size[2] * u[2] + 0.3 * w]
        };
        let mut out = Vec::new();
        for i in 0..n[0] {
            for j in 0..n[1] {
                for k in 0..n[2] {
                    let (i, j, k) = (i as f64, j as f64, k as f64);
                    let c = [
                        at(i, j, k),
                        at(i + 1.0, j, k),
                        at(i + 1.0, j + 1.0, k),
                        at(i, j + 1.0, k),
                        at(i, j, k + 1.0),
                        at(i + 1.0, j, k + 1.0),
                        at(i + 1.0, j + 1.0, k + 1.0),
                        at(i, j + 1.0, k + 1.0),
                    ];
                    // Straight-sided mid-side nodes (the mean of the ends).
                    let mid = |a: usize, b: usize| [0, 1, 2].map(|q| 0.5 * (c[a][q] + c[b][q]));
                    match kind {
                        ElementKind::Hex8 => out.push((kind, c.to_vec())),
                        ElementKind::Hex20 => {
                            let edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)];
                            let mut x = c.to_vec();
                            x.extend(edges.iter().map(|(a, b)| mid(*a, *b)));
                            out.push((kind, x));
                        }
                        ElementKind::Wedge6 => {
                            out.push((kind, vec![c[0], c[1], c[2], c[4], c[5], c[6]]));
                            out.push((kind, vec![c[0], c[2], c[3], c[4], c[6], c[7]]));
                        }
                        ElementKind::Tet4 | ElementKind::Tet10 => {
                            // Six tetrahedra around the diagonal 0–6, positively oriented.
                            for (a, b) in [(1, 2), (2, 3), (3, 7), (7, 4), (4, 5), (5, 1)] {
                                let corners = [c[0], c[a], c[b], c[6]];
                                let mut x = corners.to_vec();
                                if kind == ElementKind::Tet10 {
                                    for (p, q) in [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)] {
                                        x.push([0, 1, 2].map(|r| 0.5 * (corners[p][r] + corners[q][r])));
                                    }
                                }
                                out.push((kind, x));
                            }
                        }
                        _ => {}
                    }
                }
            }
        }
        out
    }

    fn linear(p: [f64; 3]) -> f64 {
        3.0 - 2.0 * p[0] + 0.7 * p[1] + 5.5 * p[2]
    }

    #[test]
    fn the_box_meshes_are_valid() {
        for kind in [ElementKind::Tet4, ElementKind::Tet10, ElementKind::Wedge6, ElementKind::Hex8, ElementKind::Hex20] {
            for (k, xyz) in box_mesh(kind, [3, 2, 2], [1.0, 0.8, 0.6], 0.03) {
                let v = crate::element::Solid::new(k, &xyz, crate::element::Formulation::Full).unwrap().volume().unwrap();
                assert!(v > 0.0, "{kind:?}: {v}");
            }
        }
    }

    /// Consistent interpolation is exact for uniform and linear fields between two different
    /// meshes of one box, of every solid kind, with distorted interiors.
    #[test]
    fn uniform_and_linear_fields_transfer_exactly_between_non_matching_meshes() {
        let target = box_mesh(ElementKind::Hex20, [5, 3, 4], [1.0, 0.8, 0.6], 0.02);
        let points: Vec<[f64; 3]> = target.iter().flat_map(|(_, x)| x.clone()).collect();
        for kind in [ElementKind::Tet4, ElementKind::Tet10, ElementKind::Wedge6, ElementKind::Hex8, ElementKind::Hex20] {
            let source = box_mesh(kind, [3, 2, 2], [1.0, 0.8, 0.6], 0.03);
            let locator = Locator::new(source.clone()).unwrap();
            for p in &points {
                let hit = locator.locate(*p, 1e-9).unwrap_or_else(|| panic!("{kind:?}: {p:?} not found"));
                assert_eq!(hit.distance, 0.0, "{kind:?}: {p:?}");
                let xyz = &source[hit.element].1;
                let uniform: f64 = hit.weights.iter().map(|w| w * 42.0).sum();
                let value: f64 = hit.weights.iter().zip(xyz).map(|(w, x)| w * linear(*x)).sum();
                assert!((uniform - 42.0).abs() < 1e-12, "{kind:?}: {uniform}");
                assert!((value - linear(*p)).abs() < 1e-11, "{kind:?} at {p:?}: {value} vs {}", linear(*p));
            }
        }
    }

    #[test]
    fn a_point_just_outside_takes_the_closest_element_within_the_tolerance() {
        let locator = Locator::new(box_mesh(ElementKind::Hex8, [2, 2, 2], [1.0, 1.0, 1.0], 0.0)).unwrap();
        assert!((locator.mean_size() - 0.5).abs() < 1e-12);
        let p = [0.3, 0.4, 1.0 + 1e-3];
        assert!(locator.locate(p, 1e-4).is_none());
        let hit = locator.locate(p, 1e-2).unwrap();
        assert!((hit.distance - 1e-3).abs() < 1e-12, "{}", hit.distance);
        // The weights are those of the nearest surface point: a linear field there.
        let source = box_mesh(ElementKind::Hex8, [2, 2, 2], [1.0, 1.0, 1.0], 0.0);
        let value: f64 = hit.weights.iter().zip(&source[hit.element].1).map(|(w, x)| w * linear(*x)).sum();
        assert!((value - linear([0.3, 0.4, 1.0])).abs() < 1e-12);
        assert!(locator.locate([5.0, 5.0, 5.0], 0.1).is_none());
        assert!(locator.locate([f64::NAN, 0.0, 0.0], 1.0).is_none());
        assert!(locator.locate([0.5, 0.5, 0.5], f64::INFINITY).is_some());
    }

    #[test]
    fn a_pyramid_reproduces_a_linear_field() {
        let ty = ELEMENT_TYPES.iter().find(|t| t.kind == ElementKind::Pyramid5).unwrap();
        let xyz = vec![[0.0, 0.0, 0.0], [2.0, 0.1, 0.0], [2.1, 2.0, 0.1], [0.0, 1.9, 0.0], [1.0, 1.0, 1.5]];
        for p in [[1.0, 1.0, 0.5], [0.3, 0.3, 0.05], [1.0, 1.0, 1.4]] {
            let (xi, d) = natural(ty, &xyz, p).unwrap();
            assert_eq!(d, 0.0);
            let v: f64 = (ty.shape)(xi).values().iter().zip(&xyz).map(|(w, x)| w * linear(*x)).sum();
            assert!((v - linear(p)).abs() < 1e-10, "{p:?}: {v}");
        }
    }

    #[test]
    fn bad_sources_are_refused() {
        assert!(Locator::new(vec![(ElementKind::Hex8, vec![[0.0; 3]; 4])]).is_err());
        assert!(Locator::new(vec![(ElementKind::Quad4, vec![[0.0; 3]; 4])]).is_err());
        assert!(Locator::new(vec![(ElementKind::Tet4, vec![[f64::NAN; 3]; 4])]).is_err());
        let empty = Locator::new(Vec::new()).unwrap();
        assert!(empty.is_empty() && empty.locate([0.0; 3], 1.0).is_none());
        // A collapsed element is skipped, not a panic.
        let flat = Locator::new(vec![(ElementKind::Tet4, vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]])]).unwrap();
        assert!(flat.locate([0.2, 0.2, 0.0], 1.0).is_none());
    }

    #[test]
    fn the_simplex_projection_is_the_closest_point() {
        assert_eq!(simplex(&[0.2, 0.3, 0.1]), [0.2, 0.3, 0.1]);
        assert_eq!(simplex(&[-0.2, 0.3, 0.1]), [0.0, 0.3, 0.1]);
        let p = simplex(&[1.0, 1.0, 1.0]);
        assert!(p.iter().all(|v| (v - 1.0 / 3.0).abs() < 1e-15));
        assert_eq!(simplex(&[2.0, -1.0]), [1.0, 0.0]);
    }
}
