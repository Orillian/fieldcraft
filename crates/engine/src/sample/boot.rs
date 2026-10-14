//! The boot seal: a conical rubber bellows of our own design, its six rounded convolutions
//! squeezed axially by a rigid collar until they fold onto each other.
//!
//! **Profile.** The wall's mid-surface in the (r, y) half-plane, y the axis: a base cylinder
//! (r = [`BASE_RADIUS`]), then six V-shaped convolutions zig-zagging out to a crest and back in to
//! a root, the crests and roots tapering linearly towards the top, and the top collar, a cylinder
//! at the last root's radius. Every corner is rounded with radius [`BEND`]. The wall is the
//! mid-surface offset by ± [`WALL`]/2 along its normal: a structured quad mesh, [`LAYERS`]
//! elements through the wall and about [`STEP`] mm along it.

use serde_json::{Value, json};

/// Mid-surface radius of the base cylinder (mm).
pub const BASE_RADIUS: f64 = 27.0;
/// Wall thickness (mm).
pub const WALL: f64 = 1.8;
/// Radius of the rounded corners of the mid-surface (mm).
pub const BEND: f64 = 1.6;
/// Convolutions, and their pitch along the axis (mm).
pub const FOLDS: usize = 6;
const PITCH: f64 = 9.5;
/// Height of the straight base, and the boot's height (mm).
const BASE_LENGTH: f64 = 4.0;
pub const HEIGHT: f64 = 70.0;
/// Crest and root radii of the first and last convolutions (mm).
const CREST: [f64; 2] = [38.0, 26.0];
const ROOT: [f64; 2] = [28.0, 17.0];
/// Elements through the wall and the target element length along it (mm).
pub const LAYERS: usize = 3;
pub const STEP: f64 = 0.6;

/// A piece of the mid-surface: a straight line, or an arc about `center` of `radius` from angle
/// `start` through `sweep` (radians, counter-clockwise positive).
#[derive(Debug, Clone, Copy)]
enum Piece {
    Line { from: [f64; 2], to: [f64; 2] },
    Arc { center: [f64; 2], radius: f64, start: f64, sweep: f64 },
}

impl Piece {
    fn length(&self) -> f64 {
        match *self {
            Piece::Line { from, to } => (to[0] - from[0]).hypot(to[1] - from[1]),
            Piece::Arc { radius, sweep, .. } => radius * sweep.abs(),
        }
    }

    /// The point and unit tangent at fraction `s` of the piece.
    fn at(&self, s: f64) -> ([f64; 2], [f64; 2]) {
        match *self {
            Piece::Line { from, to } => {
                let l = self.length();
                ([from[0] + s * (to[0] - from[0]), from[1] + s * (to[1] - from[1])], [(to[0] - from[0]) / l, (to[1] - from[1]) / l])
            }
            Piece::Arc { center, radius, start, sweep } => {
                let a = start + s * sweep;
                let (sin, cos) = a.sin_cos();
                let turn = sweep.signum();
                ([center[0] + radius * cos, center[1] + radius * sin], [-turn * sin, turn * cos])
            }
        }
    }

}

/// The crest or root radius of convolution `k` (linear from the first to the last).
fn taper(ends: [f64; 2], k: usize) -> f64 {
    let f = k as f64 / (FOLDS - 1) as f64;
    ends[0] + f * (ends[1] - ends[0])
}

/// The corners of the mid-surface's polyline, base to collar top.
fn corners() -> Vec<[f64; 2]> {
    let mut out = vec![[BASE_RADIUS, 0.0], [BASE_RADIUS, BASE_LENGTH]];
    for k in 0..FOLDS {
        let y = BASE_LENGTH + k as f64 * PITCH;
        out.push([taper(CREST, k), y + 0.5 * PITCH]);
        out.push([taper(ROOT, k), y + PITCH]);
    }
    let collar = out.last().map_or(0.0, |c| c[0]);
    out.push([collar, HEIGHT]);
    out
}

/// The mid-surface, base to collar: the polyline of [`corners`] with each inner corner rounded
/// by an arc of radius [`BEND`] tangent to both sides (`None` when a side is too short for its
/// arcs).
fn mid_surface() -> Option<Vec<Piece>> {
    let c = corners();
    let unit = |a: [f64; 2], b: [f64; 2]| {
        let l = (b[0] - a[0]).hypot(b[1] - a[1]);
        [(b[0] - a[0]) / l, (b[1] - a[1]) / l]
    };
    let mut out = Vec::new();
    let mut from = *c.first()?;
    for w in c.windows(3) {
        let (a, v, b) = (w[0], w[1], w[2]);
        let (d1, d2) = (unit(a, v), unit(v, b));
        let cross = d1[0] * d2[1] - d1[1] * d2[0];
        let turn = cross.atan2(d1[0] * d2[0] + d1[1] * d2[1]);
        let tangent = BEND * (0.5 * turn.abs()).tan();
        let enter = [v[0] - tangent * d1[0], v[1] - tangent * d1[1]];
        // The side must hold both of its arcs' tangent lengths.
        if tangent > (v[0] - from[0]).hypot(v[1] - from[1]) {
            return None;
        }
        out.push(Piece::Line { from, to: enter });
        // The centre on the inside of the turn: left of d1 for a left (counter-clockwise) turn.
        let side = turn.signum();
        let center = [enter[0] - side * BEND * d1[1], enter[1] + side * BEND * d1[0]];
        let start = (enter[1] - center[1]).atan2(enter[0] - center[0]);
        out.push(Piece::Arc { center, radius: BEND, start, sweep: turn });
        from = [v[0] + tangent * d2[0], v[1] + tangent * d2[1]];
    }
    out.push(Piece::Line { from, to: *c.last()? });
    Some(out)
}

/// The mesh of the wall in the (r, y) half-plane: nodes row by row along the mid-surface
/// (`LAYERS + 1` per row, from the outer face in), and the number of rows.
pub fn wall_nodes() -> (Vec<[f64; 2]>, usize) {
    let mut rows: Vec<([f64; 2], [f64; 2])> = Vec::new();
    for piece in mid_surface().unwrap_or_default() {
        let n = (piece.length() / STEP).round().max(1.0) as usize;
        let first = if rows.is_empty() { 0 } else { 1 };
        for i in first..=n {
            rows.push(piece.at(i as f64 / n as f64));
        }
    }
    let mut nodes = Vec::with_capacity(rows.len() * (LAYERS + 1));
    for (p, t) in &rows {
        // The left normal of the tangent: through the wall from the outer face in, so the quads
        // run counter-clockwise.
        let n = [-t[1], t[0]];
        for j in 0..=LAYERS {
            let d = WALL * (j as f64 / LAYERS as f64 - 0.5);
            nodes.push([p[0] + d * n[0], p[1] + d * n[1]]);
        }
    }
    (nodes, rows.len())
}

/// The top of the collar (mm): the boot's height.
pub fn height() -> f64 {
    HEIGHT
}

/// The wall's quad4 mesh in the (r, y) = (x, y) plane, tagged as face `wall` of body `Boot`
/// with its edges `base`, `collar` (the ends), `outer` and `inner` (the faces of the wall).
fn profile(steps: &mut Vec<Value>) {
    let (nodes, rows) = wall_nodes();
    let per = LAYERS + 1;
    let label = |i: usize, j: usize| (i * per + j + 1) as u64;
    for (k, p) in nodes.iter().enumerate() {
        steps.push(json!({"command": "mesh.node.add", "params": {"label": k + 1, "xyz": [p[0], p[1], 0.0]}}));
    }
    let mut element = 0u64;
    let mut edges: [Vec<Value>; 4] = Default::default();
    for i in 0..rows.saturating_sub(1) {
        for j in 0..LAYERS {
            element += 1;
            let quad = [label(i, j), label(i + 1, j), label(i + 1, j + 1), label(i, j + 1)];
            steps.push(json!({"command": "mesh.element.add", "params": {"label": element, "kind": "quad4", "nodes": quad}}));
            // Edge 1 runs along the outer face (j = 0), 3 along the inner one, 4 and 2 across
            // the wall at the base and the collar.
            let mut on = |k: usize, edge: u8| edges[k].push(json!({"element": element, "edge": edge}));
            if i == 0 {
                on(0, 4);
            }
            if i + 2 == rows {
                on(1, 2);
            }
            if j == 0 {
                on(2, 1);
            }
            if j + 1 == LAYERS {
                on(3, 3);
            }
        }
    }
    steps.push(json!({"command": "mesh.tag", "params": {"geometry": {"body": "Boot", "face": "wall"}, "elements": (1..=element).collect::<Vec<_>>()}}));
    for (name, list) in ["base", "collar", "outer", "inner"].iter().zip(edges) {
        steps.push(json!({"command": "mesh.tag", "params": {"geometry": {"body": "Boot", "edge": name}, "edges": list}}));
    }
}

/// The rubber, its section, general contact on every exterior face or edge, the clamped base
/// and the collar pushed down by `travel` along the axis (`axis`: 1 for Y in 2D, 2 for Z in 3D)
/// and held across it (a rigid collar).
fn analysis(steps: &mut Vec<Value>, section: Value, entity: &str, axis: usize, travel: f64, extra: &[&str]) {
    // The elements: the profile's face, or the solid body.
    let region = if axis == 1 { json!({"body": "Boot", "face": "wall"}) } else { json!({"body": "Boot"}) };
    let mut push = [json!(0), json!(0), json!(0), Value::Null, Value::Null, Value::Null];
    if axis == 1 {
        push[2] = Value::Null;
    }
    push[axis] = json!(-travel);
    let mut held = vec!["Base", "Push"];
    held.extend_from_slice(extra);
    steps.extend([
        json!({"command": "model.set.create", "params": {"name": "boot", "geometry": [region]}}),
        json!({"command": "model.set.create", "params": {"name": "base", "geometry": [{"body": "Boot", entity: "base"}]}}),
        json!({"command": "model.set.create", "params": {"name": "collar", "geometry": [{"body": "Boot", entity: "collar"}]}}),
        json!({"command": "model.material.create", "params": {"name": "Rubber", "elastic": {"type": "isotropic", "young": RUBBER[0], "poisson": RUBBER[1]}, "density": 1.1e-9,
            "hyperelastic": {"type": "mooney_rivlin", "c01": RUBBER[2]}}}),
        json!({"command": "model.section.create", "params": {"name": "Boot", "set": "boot", "material": "Rubber", "kind": section, "formulation": {"hybrid": true}}}),
        json!({"command": "interaction.contact.general", "params": {"name": "Folds"}}),
        json!({"command": "model.boundary_condition.create", "params": {"name": "Base", "type": "fixed", "target": {"set": "base"}}}),
        json!({"command": "model.boundary_condition.create", "params": {"name": "Push", "type": "displacement", "target": {"set": "collar"}, "values": push}}),
        json!({"command": "model.step.create", "params": {"name": "Squeeze", "procedure": {"type": "nonlinear_static", "initial_increment": 0.02, "max_increment": 0.05, "max_increments": 1000, "line_search": true, "stabilization": STABILIZATION}, "boundary_conditions": held}}),
    ]);
}

/// The rubber: Young's modulus (MPa), Poisson's ratio and Mooney–Rivlin C01 (MPa), so
/// μ = E / (2 (1 + ν)) = 1.0 MPa, C10 = μ/2 − C01 = 0.4 MPa: a soft generic rubber, nearly
/// incompressible (κ/μ ≈ 830).
pub const RUBBER: [f64; 3] = [3.0, 0.499, 0.1];

/// The squeeze's stabilization (the fraction of its work the damping may dissipate).
pub const STABILIZATION: f64 = 2e-3;

/// The axisymmetric boot (hybrid quad4) squeezed by `travel` mm.
pub fn axisymmetric(travel: f64) -> Value {
    let mut steps = vec![json!({"command": "model.units", "params": {"system": "mm-t-s-N-MPa"}})];
    profile(&mut steps);
    analysis(&mut steps, json!({"type": "axisymmetric"}), "edge", 1, travel, &[]);
    Value::Array(steps)
}

/// The boot in 3D: the axisymmetric wall mesh revolved about the axis (`mesh.revolve`, Z up)
/// into hybrid hex8, `sectors` around the circumference, or around the half x ≥ 0 (from −90° to
/// 90°) with the cut held on the symmetry plane x = 0.
pub fn solid(travel: f64, sectors: usize, half: bool) -> Value {
    let mut steps = vec![json!({"command": "model.units", "params": {"system": "mm-t-s-N-MPa"}})];
    profile(&mut steps);
    let (angle, start) = if half { (180, -90) } else { (360, 0) };
    steps.push(json!({"command": "mesh.revolve", "params": {"angle": angle, "sectors": sectors, "start": start}}));
    let mut extra = Vec::new();
    if half {
        steps.extend([
            json!({"command": "model.set.create", "params": {"name": "cut", "geometry": [{"body": "Boot", "face": "start"}, {"body": "Boot", "face": "end"}]}}),
            json!({"command": "model.boundary_condition.create", "params": {"name": "Symmetry", "type": "symmetry", "target": {"set": "cut"}, "plane": "x"}}),
        ]);
        extra.push("Symmetry");
    }
    analysis(&mut steps, json!({"type": "solid"}), "face", 2, travel, &extra);
    Value::Array(steps)
}
