//! `mesh.revolve`: a 2D mesh of a profile in the XY plane (x the radius, Y the axis, as
//! axisymmetric sections take it) swept about the axis into a solid mesh.
//!
//! The profile's point (r, y) at angle θ becomes (r cos θ, r sin θ, y): the axis becomes Z (up)
//! and the profile's half-plane at θ = `start` the first ring of nodes. Each quad4 becomes a hex8
//! per sector and each tri3 a wedge6, the profile's node order reversed so the swept elements are
//! positively oriented. The 2D elements and their nodes are replaced (new labels after the
//! largest ones). Geometry tags carry over: an element filling face F of body B fills body B, an
//! element edge on edge E of body B gives the swept faces F named E of body B, and when the
//! sweep is not a whole turn its first and last sides are faces `start` and `end` of the body,
//! so sets named by geometry resolve on the solid. A whole turn joins its last ring to the first.

use std::collections::BTreeMap;

use fieldcraft_model::{ElementKind, Entity, GeometryRef};
use serde_json::{Value, json};

use super::mesh::optional_count;
use super::optional_number;
use crate::{CommandKind, CommandSpec, EngineError, Result, Session};

pub const REVOLVE: CommandSpec = CommandSpec {
    id: "mesh.revolve",
    label: "Revolve 2D Mesh",
    params: r#"{"angle": 360 (degrees, up to a whole turn), "sectors": 36 (elements around, 1..=3600), "start": 0 (degrees: where the profile's half-plane lies, from +X towards +Y)} — sweeps the 2D mesh (quad4 and tri3 in the XY plane, x ≥ 0 the radius, Y the axis) about the axis into hex8 and wedge6 with the axis along Z: (r, y) at angle θ becomes (r cos θ, r sin θ, y). Replaces the 2D elements and their nodes; geometry tags carry over (a face's elements fill its body, an edge's element edges sweep to faces of that name, and a partial turn's first and last sides are faces `start` and `end`), so sets named by geometry resolve on the solid"#,
    form: Some(super::forms::REVOLVE),
    kind: CommandKind::Model,
    run: revolve,
};

/// Most sectors of a sweep.
const MAX_SECTORS: u64 = 3600;

/// The hex8 face (1-based) swept by each quad4 edge, and the wedge6 face swept by each tri3
/// edge, with the profile's nodes reversed (see [`swept`]); faces 1 and 2 are the first and
/// last sides of a sector.
const QUAD_EDGE_FACE: [u8; 4] = [6, 5, 4, 3];
const TRI_EDGE_FACE: [u8; 3] = [5, 4, 3];

/// A swept element's nodes: the profile's nodes reversed on ring `k`, then on ring `k + 1`.
fn swept(profile: &[u32], ring: impl Fn(u32, usize) -> u32, k: usize) -> Vec<u32> {
    let reversed: Vec<u32> = profile.first().copied().into_iter().chain(profile.iter().skip(1).rev().copied()).collect();
    reversed.iter().map(|n| ring(*n, k)).chain(reversed.iter().map(|n| ring(*n, k + 1))).collect()
}

fn revolve(s: &mut Session, p: &Value) -> Result<Value> {
    let angle = optional_number(p, "angle", 360.0)?;
    if !(angle > 0.0 && angle <= 360.0) {
        return Err(EngineError::Params(format!("`angle` must be in (0, 360] degrees, got {angle}")));
    }
    let start = optional_number(p, "start", 0.0)?;
    if !start.is_finite() {
        return Err(EngineError::Params("`start` must be a number of degrees".into()));
    }
    let sectors = optional_count(p, "sectors", 36, 1..=MAX_SECTORS)?;
    let whole = angle >= 360.0;
    if whole && sectors < 3 {
        return Err(EngineError::Params("a whole turn takes at least 3 sectors".into()));
    }
    // The profile: every element a linear 2D element in the XY plane, off the axis.
    let profile: Vec<(u32, ElementKind, Vec<u32>)> = s.model.elements().iter().map(|(l, e)| (*l, e.kind, e.nodes.clone())).collect();
    if profile.is_empty() {
        return Err(EngineError::Params("there is no 2D mesh to revolve".into()));
    }
    if let Some((l, k, _)) = profile.iter().find(|(_, k, _)| !matches!(k, ElementKind::Quad4 | ElementKind::Tri3)) {
        return Err(EngineError::Params(format!(
            "not supported yet: revolving {} elements (element {l}): mesh.revolve takes a mesh of quad4 and tri3 elements",
            k.name()
        )));
    }
    let used: Vec<u32> = {
        let mut v: Vec<u32> = profile.iter().flat_map(|(_, _, n)| n.iter().copied()).collect();
        v.sort_unstable();
        v.dedup();
        v
    };
    let mut points = Vec::with_capacity(used.len());
    for n in &used {
        let xyz = s.model.node(*n).map(|n| n.xyz).ok_or_else(|| EngineError::Params(format!("node {n} is missing")))?;
        if xyz[2] != 0.0 {
            return Err(EngineError::Params(format!("node {n} is off the XY plane (z = {})", xyz[2])));
        }
        if !(xyz[0] > 0.0) {
            return Err(EngineError::Params(format!(
                "not supported yet: revolving a profile that reaches the axis (node {n} at x = {}): its elements would collapse; keep the profile at x > 0",
                xyz[0]
            )));
        }
        points.push(xyz);
    }
    let rings = if whole { sectors } else { sectors + 1 };
    let count = |a: usize, b: usize| a.checked_mul(b).ok_or_else(|| EngineError::Params("the revolved mesh is too large".into()));
    let (node_count, element_count) = (count(used.len(), rings)?, count(profile.len(), sectors)?);
    if element_count > fieldcraft_model::MAX_ELEMENTS {
        return Err(EngineError::Params(format!("the revolved mesh would have {element_count} elements, more than {}", fieldcraft_model::MAX_ELEMENTS)));
    }
    // The tags to carry over, read before the profile is removed.
    let tags = s.model.classification().clone();
    let entity = |tag: u32| tags.entity(tag).cloned();
    let (first_node, first_element) = super::fresh_labels(&s.model, node_count, element_count)?;
    let index: BTreeMap<u32, usize> = used.iter().enumerate().map(|(i, n)| (*n, i)).collect();
    let ring = |n: u32, k: usize| -> u32 {
        let i = index.get(&n).copied().unwrap_or(0);
        first_node + ((k % rings) * used.len() + i) as u32
    };
    let (start_rad, step) = (start.to_radians(), angle.to_radians() / sectors as f64);
    for k in 0..rings {
        let (sin, cos) = (start_rad + step * k as f64).sin_cos();
        for (i, x) in points.iter().enumerate() {
            s.model.set_node(first_node + (k * used.len() + i) as u32, [x[0] * cos, x[0] * sin, x[1]])?;
        }
    }
    let mut label = first_element;
    let mut new_tags: Vec<(u32, u8, GeometryRef)> = Vec::new();
    let mut regions: Vec<(u32, GeometryRef)> = Vec::new();
    for k in 0..sectors {
        for (old, kind, nodes) in &profile {
            let (solid, edge_face) = match kind {
                ElementKind::Quad4 => (ElementKind::Hex8, &QUAD_EDGE_FACE[..]),
                _ => (ElementKind::Wedge6, &TRI_EDGE_FACE[..]),
            };
            s.model.set_element(label, solid, swept(nodes, ring, k))?;
            if let Some(t) = tags.elements().get(old) {
                if let Some(region) = t.region.and_then(entity) {
                    let body = GeometryRef::body(region.body.clone())?;
                    regions.push((label, body.clone()));
                    if !whole && k == 0 {
                        new_tags.push((label, 1, GeometryRef::face(body.body.clone(), "start".to_string())?));
                    }
                    if !whole && k + 1 == sectors {
                        new_tags.push((label, 2, GeometryRef::face(body.body.clone(), "end".to_string())?));
                    }
                }
                for (edge, tag) in &t.edges {
                    let (Some(e), Some(face)) = (entity(*tag), edge_face.get(usize::from(*edge).saturating_sub(1))) else { continue };
                    if let Entity::Edge(name) = &e.entity {
                        new_tags.push((label, *face, GeometryRef::face(e.body.clone(), name.clone())?));
                    }
                }
            }
            label += 1;
        }
    }
    // The profile goes; its tags with it.
    for (old, _, _) in &profile {
        s.model.remove_element(*old)?;
    }
    for n in &used {
        s.model.remove_node(*n)?;
    }
    for (element, body) in regions {
        let tag = s.model.geometry_tag(&body)?;
        s.model.tag_element(element, tag)?;
    }
    for (element, face, reference) in new_tags {
        let tag = s.model.geometry_tag(&reference)?;
        s.model.tag_element_face(element, face, tag)?;
    }
    Ok(json!({
        "nodes": node_count,
        "elements": element_count,
        "first_node": first_node,
        "first_element": first_element,
        "sectors": sectors,
        "angle": angle,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A ring of rectangular section (1 ≤ r ≤ 2, 0 ≤ y ≤ 1) as two quads and two triangles,
    /// revolved: its volume is π (2² − 1²) · 1 = 3π, every element valid, the tags carried.
    #[test]
    fn a_revolved_profile_is_a_valid_solid_of_the_right_volume() {
        for (angle, sectors) in [(360.0f64, 48usize), (90.0, 12)] {
            let mut s = Session::default();
            s.execute("mesh.planar", &json!({"size": 0.5, "polygon": [[1, 0], [2, 0], [2, 1], [1, 1]], "name": "Ring"})).unwrap();
            let profile = s.model.elements().len();
            let r = s.execute("mesh.revolve", &json!({"angle": angle, "sectors": sectors})).unwrap();
            assert_eq!(r["elements"].as_u64(), Some((profile * sectors) as u64));
            let q = s.execute("mesh.quality", &json!({})).unwrap();
            assert_eq!(q["invalid"]["count"].as_u64(), Some(0), "{q}");
            let volume: f64 = q["volume"].as_f64().unwrap();
            // Chords: a polygon of `sectors` sides over the arc.
            let chord = (angle.to_radians() / sectors as f64).sin() * sectors as f64 / angle.to_radians();
            let exact = 3.0 * std::f64::consts::PI * angle / 360.0;
            assert!((volume - exact * chord).abs() < 1e-9 * exact, "{volume} {}", exact * chord);
            // The bottom edge (e1, y = 0) sweeps to face e1: a set on it holds the z = 0 faces.
            s.execute("model.set.create", &json!({"name": "bottom", "geometry": [{"body": "Ring", "face": "e1"}]})).unwrap();
            let set = s.execute("model.set.resolve", &json!({"name": "bottom"})).unwrap();
            assert_eq!(set["faces"]["count"].as_u64(), Some(2 * sectors as u64), "two element edges along e1: {set}");
            if angle < 360.0 {
                s.execute("model.set.create", &json!({"name": "cut", "geometry": [{"body": "Ring", "face": "start"}]})).unwrap();
                let set = s.execute("model.set.resolve", &json!({"name": "cut"})).unwrap();
                assert_eq!(set["faces"]["count"].as_u64(), Some(profile as u64), "{set}");
            }
        }
    }

    #[test]
    fn what_mesh_revolve_refuses() {
        let mut s = Session::default();
        assert!(s.execute("mesh.revolve", &json!({})).is_err(), "no mesh");
        s.execute("mesh.planar", &json!({"size": 0.5, "polygon": [[0, 0], [2, 0], [2, 1], [0, 1]]})).unwrap();
        let e = s.execute("mesh.revolve", &json!({})).unwrap_err().to_string();
        assert!(e.contains("axis"), "{e}");
        let mut s = Session::default();
        s.execute("mesh.planar", &json!({"size": 0.5, "kind": "quad8", "polygon": [[1, 0], [2, 0], [2, 1], [1, 1]]})).unwrap();
        assert!(s.execute("mesh.revolve", &json!({})).unwrap_err().to_string().contains("quad8"));
        for bad in [json!({"angle": 0}), json!({"angle": 400}), json!({"sectors": 0}), json!({"angle": 360, "sectors": 2}), json!({"start": "x"})] {
            assert!(s.execute("mesh.revolve", &bad).is_err(), "{bad}");
        }
    }
}
