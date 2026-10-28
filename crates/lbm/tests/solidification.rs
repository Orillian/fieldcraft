//! Melting and freezing on the temperature lattice (the total-enthalpy scheme): the two-phase
//! Stefan problem against Neumann's closed form, the latent heat's balance, and a channel half
//! frozen by a cold wall, whose liquid layer carries the Poiseuille flow.

mod common;

use common::*;
use fieldcraft_lbm::Flow;
use fieldcraft_model::{BoxFace, FaceCondition, FluidThermal, Material, Model, PhaseChange, Step, ThermalCondition};

const PERIODIC: FaceCondition = FaceCondition::Periodic {};

/// [`hot_model`] with a phase change: latent heat `latent` (`L / c`, the specific heat is 1)
/// between `solidus` and `liquidus` (°C).
fn freezing_model(nu: f64, alpha: f64, phase: [f64; 3], domain: fieldcraft_model::FluidDomain, step: fieldcraft_model::Procedure) -> Model {
    let [latent, solidus, liquidus] = phase;
    let mut m = Model::default();
    let material = Material {
        phase_change: Some(PhaseChange { latent_heat: latent, solidus: kelvin(solidus), liquidus: kelvin(liquidus) }),
        ..hot_fluid(1.0, nu, alpha)
    };
    m.add_material("Fluid", material).unwrap();
    m.add_fluid_domain("Flow", domain).unwrap();
    m.add_step(Step::new("Run", step), None).unwrap();
    m
}

/// A column one cell wide along y, periodic in x and z.
fn column(n: usize, length: f64) -> fieldcraft_model::FluidDomain {
    let dx = length / n as f64;
    let mut d = domain([0.0; 3], [dx, length, dx], dx);
    d.faces.x_min = PERIODIC;
    d.faces.x_max = PERIODIC;
    d.faces.z_min = PERIODIC;
    d.faces.z_max = PERIODIC;
    d
}

/// The error function by its Maclaurin series (|x| ≤ 3 here, to round-off).
fn erf(x: f64) -> f64 {
    let (mut term, mut sum) = (x, x);
    for n in 1..200 {
        term *= -x * x / n as f64;
        let next = term / (2 * n + 1) as f64;
        sum += next;
        if next.abs() < 1e-17 * sum.abs() {
            break;
        }
    }
    2.0 / std::f64::consts::PI.sqrt() * sum
}

/// Neumann's λ of the two-phase Stefan problem with equal properties in both phases: a liquid at
/// `t_i` above the melting temperature `t_m` freezes from a wall held at `t_w` below it, the front
/// at `s = 2 λ √(α t)`, where `St_s e^{−λ²} / erf λ − St_l e^{−λ²} / erfc λ = λ √π` with
/// `St_s = c (t_m − t_w) / L` and `St_l = c (t_i − t_m) / L` (H. S. Carslaw, J. C. Jaeger,
/// *Conduction of Heat in Solids*, 2nd ed., Oxford 1959, §11.2).
fn neumann(latent: f64, t_w: f64, t_m: f64, t_i: f64) -> f64 {
    let (st_s, st_l) = ((t_m - t_w) / latent, (t_i - t_m) / latent);
    let residual = |l: f64| {
        let e = (-l * l).exp();
        st_s * e / erf(l) - st_l * e / (1.0 - erf(l)) - l * std::f64::consts::PI.sqrt()
    };
    let (mut lo, mut hi) = (1e-6, 3.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if residual(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// The total enthalpy on the lattice (the sum of every scalar population).
fn enthalpy(flow: &Flow, cells: usize) -> f64 {
    (0..cells).map(|c| (0..7).map(|q| flow.thermal_population(q, c)).sum::<f64>()).sum()
}

/// The two-phase Stefan problem: a liquid at +0.5 above its melting point 0 freezes from a wall
/// held at −1, with L / c = 2 (Stefan numbers 0.5 in the solid, 0.25 in the liquid). The frozen
/// thickness `Σ (1 − f_l) Δx` follows Neumann's `2 λ √(α t)`, the temperatures his profiles, and
/// the heat the cold wall takes is the fall in the lattice's total enthalpy to round-off.
#[test]
fn a_freezing_front_follows_neumanns_solution() {
    let (latent, t_w, t_m, t_i, alpha) = (2.0, -1.0, 0.0, 0.5, 0.05);
    let lambda = neumann(latent, t_w, t_m, t_i);
    let mut last = 0.0;
    for n in [64usize, 128] {
        // a column twice as long as the depth the cold reaches, so its far end stays at t_i
        let length = 2.0;
        let dx = 1.0 / n as f64;
        let cells = (length / dx).round() as usize;
        let mut d = column(cells, length);
        let mut t = FluidThermal { initial_temperature: kelvin(t_i), ..Default::default() };
        t.faces.y_min = ThermalCondition::Temperature { value: kelvin(t_w) };
        d.thermal = Some(t);
        // τ_T = 1: Δt = Δx² / (8 α); the frozen layer grows to 0.3
        let dt = dx * dx / (8.0 * alpha);
        let end = (0.3 / (2.0 * lambda)).powi(2) / alpha;
        let m = freezing_model(0.1, alpha, [latent, t_m, t_m], d, step(end, dt, None));
        let (setup, mut flow) = prepare(&m, &[]);
        assert!(flow.phase_change());
        let e0 = enthalpy(&flow, cells);
        let mut taken = 0.0;
        let mut worst: f64 = 0.0;
        for k in 1..=setup.steps {
            flow.step().unwrap();
            taken += flow.heat()[BoxFace::YMin.index()];
            if k % (setup.steps / 4).max(1) == 0 {
                let time = k as f64 * dt;
                let s = 2.0 * lambda * (alpha * time).sqrt();
                let frozen: f64 = (0..cells).map(|c| 1.0 - flow.liquid_fraction(c).unwrap()).sum::<f64>() * dx;
                println!("n {n}: t {time:.4} front {frozen:.5} vs Neumann {s:.5} ({:+.2} %)", 100.0 * (frozen / s - 1.0));
                worst = worst.max((frozen / s - 1.0).abs());
            }
        }
        let time = setup.steps as f64 * dt;
        let root = 2.0 * (alpha * time).sqrt();
        let s = lambda * root;
        let mut profile: f64 = 0.0;
        for c in 0..cells {
            let y = (c as f64 + 0.5) * dx;
            let exact =
                if y < s { t_w + (t_m - t_w) * erf(y / root) / erf(lambda) } else { t_i - (t_i - t_m) * (1.0 - erf(y / root)) / (1.0 - erf(lambda)) };
            profile = profile.max((flow.temperature(c).unwrap() - exact).abs());
        }
        let balance = (e0 - enthalpy(&flow, cells) - taken).abs() / taken;
        println!("n {n}: λ {lambda:.6}, front within {:.2} %, temperatures within {profile:.2e}, heat balance {balance:.1e}", 100.0 * worst);
        assert!(worst < [0.012, 0.004][usize::from(n == 128)], "{n}: front {worst}");
        // measured: 0.6 % and 0.11 % of the range
        assert!(profile < [0.008, 0.0015][usize::from(n == 128)] * (t_i - t_w), "{n}: temperature {profile}");
        assert!(balance < 1e-10, "{n}: the latent heat is accounted for: {balance}");
        if n == 128 {
            assert!(worst < last, "the front converges: {worst} after {last}");
        }
        last = worst;
    }
}

/// A channel between a wall held at −1 and one held at +1, with a melting point of −0.5: the
/// steady temperature is linear and the fluid freezes on the cold quarter of the gap. A body
/// force along x drives the liquid three quarters as plane Poiseuille flow between the front and
/// the hot wall, `Q = a (H − y*)³ / (12 ν)`, and the frozen layer holds still.
#[test]
fn a_frozen_layer_stops_the_flow() {
    let (h, nu, alpha, a, latent): (f64, f64, f64, f64, f64) = (1.0, 0.1, 0.05, 1e-3, 1.0);
    let front = 0.25 * h;
    let expected = a * (h - front).powi(3) / (12.0 * nu);
    let mut last = f64::INFINITY;
    for n in [16usize, 32] {
        let dx = h / n as f64;
        let mut d = column(n, h);
        d.acceleration = [a, 0.0, 0.0];
        let mut t = FluidThermal::default();
        t.faces.y_min = ThermalCondition::Temperature { value: kelvin(-1.0) };
        t.faces.y_max = ThermalCondition::Temperature { value: kelvin(1.0) };
        d.thermal = Some(t);
        let dt = dx * dx / (8.0 * alpha);
        let m = freezing_model(nu, alpha, [latent, -0.5, -0.5], d, step(4.0 * (h - front).powi(2) / nu, dt, None));
        let (setup, mut flow) = prepare(&m, &[]);
        // start from the steady temperature; the flow develops from rest
        flow.set_temperature(|c| -1.0 + 2.0 * (c as f64 + 0.5) * dx / h);
        run_all(&setup, &mut flow);
        let u = |c: usize| flow.moments(c).unwrap().1[0] * dx / dt;
        let q: f64 = (0..n).map(u).sum::<f64>() * dx;
        let still = (0..n).filter(|&c| flow.liquid_fraction(c).unwrap() == 0.0).map(|c| u(c).abs()).fold(0.0, f64::max);
        let frozen = (0..n).filter(|&c| flow.liquid_fraction(c).unwrap() < 1.0).count() as f64 * dx;
        let fastest = (0..n).map(|c| u(c).abs()).fold(0.0, f64::max);
        let error = (q / expected - 1.0).abs();
        println!(
            "n {n}: frozen {frozen:.4} (front {front}), flow {q:.6e} vs {expected:.6e} ({:+.2} %), frozen cells move at {:.1e} of the fastest",
            100.0 * (q / expected - 1.0),
            still / fastest
        );
        assert!((frozen - front).abs() < 0.5 * dx, "{n}: the front at {frozen}");
        assert!(still < 5e-7 * fastest, "{n}: frozen cells move at {still}");
        assert!(error < [0.095, 0.048][usize::from(n == 32)], "{n}: flow {q} vs {expected}");
        assert!(error < last, "the flow converges: {error} after {last}");
        last = error;
    }
}

/// Without a phase change nothing changes: the populations sum to the temperature and no
/// liquid fraction is kept.
#[test]
fn without_a_phase_change_there_is_no_liquid_fraction() {
    let mut d = column(4, 1.0);
    d.thermal = Some(FluidThermal::default());
    let m = hot_model(1.0, 0.1, 0.05, d, step(1.0, 0.01, None));
    let (_, flow) = prepare(&m, &[]);
    assert!(!flow.phase_change());
    assert_eq!(flow.liquid_fraction(0), None);
}

/// A phase change without heat in the fluid is refused by name, not ignored: nothing would ever
/// freeze.
#[test]
fn a_phase_change_without_heat_is_refused() {
    let m = freezing_model(0.1, 0.05, [1.0, 0.0, 0.0], column(4, 1.0), step(1.0, 0.01, None));
    let e = fieldcraft_lbm::Setup::new(&m, None, &[], &fieldcraft_lbm::Limits::default()).err().map(|e| e.to_string());
    assert!(e.as_deref().is_some_and(|e| e.contains("phase change") && e.contains("thermal") && e.contains("`Flow`")), "{e:?}");
}
