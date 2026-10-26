//! Shear-dependent (generalised Newtonian) viscosity in the flow collision: power law,
//! Carreau–Yasuda, Cross and Cross-WLF ([`fieldcraft_model::Fluid`]).
//!
//! # Local shear rate
//! The lattice Boltzmann method carries the strain rate in each cell's populations: the
//! Chapman–Enskog expansion ties the non-equilibrium momentum flux (with Guo's forcing,
//! [`crate::turbulence::flux`]) to it, `P ≈ −2 ρ c_s² τ S` (lattice units, Δx = Δt = 1; Krüger
//! et al. 2017, §4.1 and §6.3). The shear rate `γ̇ = √(2 S:S)` is then
//!
//! `γ̇ = |P| / (√2 ρ c_s² τ)`, with `|P| = √(P:P)`,
//!
//! which depends on the relaxation time τ the cell relaxes with, while τ depends on the viscosity
//! at γ̇: `τ = ½ + ν(γ̇) / c_s²`. The product `σ = τ γ̇ = 3 |P| / (√2 ρ)` is known (a scaled shear
//! stress), so each cell solves the one equation
//!
//! `F(τ) = τ − ½ − 3 ν(σ / τ) = 0`
//!
//! for its own τ, within the stable range `[τ_min, τ_max]` (J. Boyd, J. Buick, S. Green,
//! *A second-order accurate lattice Boltzmann non-Newtonian flow model*, J. Phys. A 39:14241,
//! 2006, take γ̇ from the same moments; S. Gabbanelli, G. Drazer, J. Koplik, *Lattice Boltzmann
//! method for non-Newtonian (power-law) fluids*, PRE 72:046312, 2005, bound τ the same way). For
//! every model here the shear stress `η(γ̇) γ̇` grows with γ̇ (n > 0), so `τ γ̇` does too and
//! the root is unique: F is negative below it and positive above, at every σ. A clamped τ is
//! counted ([`Clamp`]) and reported by the solve: a shear-thinning fluid at rest (the middle of
//! a channel) reaches τ_max, a shear-thickening one τ_min, where the stress is small and the
//! viscosity matters little.
//!
//! The root is found by regula falsi with the Illinois modification (M. Dowell, P. Jarratt,
//! *A modified regula falsi method for computing the root of an equation*, BIT 11:168–174, 1971)
//! on `[τ_min, τ_max]`, bracketed from the start and bounded in iterations. For a constant
//! viscosity F is linear and the first step lands on the Newtonian τ.
//!
//! The two relaxation times keep their Λ (τ⁻ follows τ⁺ as with the eddy viscosity), and the
//! recursive regularised collision takes the cell's rate in place of the fluid's
//! ([`crate::regularized::collide`]).
//!
//! # Temperature
//! Cross-WLF reads the cell's temperature from the scalar lattice (the previous step's, which the
//! collision of the momentum precedes); the other models ignore it.
//!
//! No code of any lattice Boltzmann package was read; this is written from the papers.

use fieldcraft_model::Fluid;

/// The largest relaxation time τ⁺ a cell takes: the viscosity is capped there. Beyond it the
/// scheme stays stable but its time accuracy degrades; a shear-thinning fluid at rest reaches it.
pub const MAX_TAU: f64 = 5.0;

/// Most iterations of a cell's relaxation time solve (the regula falsi converges in a handful;
/// a cell that has not stops at its bracket's best point).
const MAX_ITERATIONS: usize = 64;

/// A cell's relaxation time is converged when a step moves it by less than this fraction.
const TOLERANCE: f64 = 1e-13;

/// Where a cell's relaxation time ended up.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Clamp {
    /// Within the range: the model's viscosity.
    None,
    /// At τ_min: the viscosity is floored (a shear-thinning fluid at high shear, a shear-thickening
    /// one at rest).
    Low,
    /// At τ_max: the viscosity is capped (a shear-thinning fluid at rest, a shear-thickening one at
    /// high shear, a frozen melt).
    High,
}

/// How many cells' relaxation times were clamped, `[low, high]`: in the last collision and over
/// the run, with the number of collisions.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Clamps {
    pub last: [u64; 2],
    pub total: [u64; 2],
    pub collisions: u64,
}

impl Clamps {
    /// Count one collision's clamped cells.
    pub fn add(&mut self, last: [u64; 2]) {
        self.last = last;
        self.total = [self.total[0].saturating_add(last[0]), self.total[1].saturating_add(last[1])];
        self.collisions = self.collisions.saturating_add(1);
    }

    /// The report of a run of `fluid_cells` cells, when a relaxation time was clamped, for the
    /// solve's warnings; `tau` the range.
    pub fn warning(&self, fluid_cells: usize, tau: [f64; 2]) -> Option<String> {
        if self.total == [0, 0] {
            return None;
        }
        let updates = (fluid_cells as f64) * (self.collisions as f64);
        let share = |n: u64| if updates > 0.0 { 100.0 * n as f64 / updates } else { 0.0 };
        Some(format!(
            "the shear-dependent viscosity left the lattice's stable range of relaxation times [{}, {}] and was held at its ends: capped (τ = {}) in {:.2}% and floored (τ = {}) in {:.2}% of the fluid cell updates, {} and {} of the {fluid_cells} fluid cells in the last step; a shear-thinning fluid at rest (a channel's middle) reaches the cap by nature, a shear-thickening one the floor, where the stress and so the viscosity's effect are small; a longer time step raises every τ, a shorter one lowers it",
            tau[0],
            tau[1],
            tau[1],
            share(self.total[1]),
            tau[0],
            share(self.total[0]),
            self.last[1],
            self.last[0],
        ))
    }
}

/// A shear-dependent fluid in lattice units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rheology {
    /// The model, in model (SI) units.
    pub fluid: Fluid,
    /// A lattice shear rate in model units: `1 / Δt`.
    pub rate_scale: f64,
    /// A dynamic viscosity (model units) as a lattice kinematic one: `Δt / (ρ₀ Δx²)`.
    pub viscosity_scale: f64,
    /// The range of τ⁺, `[τ_min, τ_max]`.
    pub tau: [f64; 2],
    /// The temperature, in kelvin, of the lattice's 0 ([`crate::ThermalSetup::zero`]), or the
    /// temperature of a flow without heat (which only Cross-WLF would read, and it needs heat).
    pub zero: f64,
}

impl Rheology {
    /// The unbounded τ⁺ of the model at lattice shear rate `rate` and kelvin `kelvin` (infinite
    /// where the viscosity is).
    #[inline(always)]
    fn tau_at(&self, rate: f64, kelvin: f64) -> f64 {
        let tau = 0.5 + 3.0 * self.fluid.viscosity(rate * self.rate_scale, kelvin) * self.viscosity_scale;
        if tau.is_nan() { f64::INFINITY } else { tau }
    }

    /// The unbounded τ⁺ of the fluid's largest viscosity at `kelvin` over every shear rate: at
    /// rest or at an unbounded rate, since every model here is monotone in γ̇ (infinite where it
    /// has no bound: a power law with n ≠ 1, a shear-thickening Carreau–Yasuda fluid).
    pub fn largest_tau(&self, kelvin: f64) -> f64 {
        [0.0, f64::MAX].into_iter().map(|rate| self.tau_at(rate, kelvin)).fold(0.0, f64::max)
    }

    /// The relaxation time τ⁺ of a cell with non-equilibrium flux `p` (as
    /// [`crate::turbulence::flux`] gives it), density `rho` and lattice temperature
    /// `temperature`, and whether it was clamped.
    #[inline]
    pub fn relax(&self, p: &[f64; 6], rho: f64, temperature: f64) -> (f64, Clamp) {
        let [lo, hi] = self.tau;
        let kelvin = temperature + self.zero;
        let pp = crate::turbulence::flux_norm2(p);
        // σ = τ γ̇, from P = −2 ρ c_s² τ S with c_s² = 1/3.
        let sigma = 3.0 * pp.sqrt() / (std::f64::consts::SQRT_2 * rho);
        if !(sigma.is_finite() && sigma >= 0.0) {
            // A broken cell (the run stops on its density); relax it as slowly as allowed.
            return (hi, Clamp::None);
        }
        let f = |tau: f64| tau - self.tau_at(sigma / tau, kelvin);
        let (mut a, mut fa) = (lo, f(lo));
        if fa >= 0.0 {
            return (lo, if fa > 0.0 { Clamp::Low } else { Clamp::None });
        }
        let (mut b, mut fb) = (hi, f(hi));
        if !(fb > 0.0) {
            return (hi, if fb < 0.0 || fb.is_nan() { Clamp::High } else { Clamp::None });
        }
        // Illinois: halve the end that stayed twice in a row, so both ends move.
        let mut last = 0i8;
        let mut c = a;
        for _ in 0..MAX_ITERATIONS {
            let mut next = (a * fb - b * fa) / (fb - fa);
            if !(next > a && next < b) {
                next = 0.5 * (a + b);
            }
            let moved = (next - c).abs();
            c = next;
            let fc = f(c);
            if fc == 0.0 || moved <= TOLERANCE * c {
                break;
            }
            if fc < 0.0 {
                (a, fa) = (c, fc);
                if last == -1 {
                    fb *= 0.5;
                }
                last = -1;
            } else {
                (b, fb) = (c, fc);
                if last == 1 {
                    fa *= 0.5;
                }
                last = 1;
            }
        }
        (c, Clamp::None)
    }
}

/// The shear rate (model units) at which `fluid` at `kelvin` carries the shear stress `stress`:
/// the root of `η(γ̇) γ̇ = stress`, which grows with γ̇. Bisection on log γ̇ after bracketing by
/// decades; `None` when there is none in 10⁻³⁰⁰ to 10³⁰⁰.
pub fn shear_rate_at_stress(fluid: &Fluid, stress: f64, kelvin: f64) -> Option<f64> {
    const DECADES: i32 = 300;
    const BISECTIONS: usize = 200;
    if !(stress > 0.0 && stress.is_finite()) {
        return None;
    }
    let excess = |log_rate: f64| {
        let rate = 10f64.powf(log_rate);
        let s = fluid.viscosity(rate, kelvin) * rate;
        if s.is_nan() { f64::INFINITY } else { s - stress }
    };
    // Bracket [lo, hi] in log₁₀ γ̇, one decade at a time from 1/s.
    let (mut lo, mut hi) = (0.0f64, 0.0f64);
    if excess(0.0) > 0.0 {
        while excess(lo) > 0.0 {
            lo -= 1.0;
            if lo < -f64::from(DECADES) {
                return None;
            }
        }
        hi = lo + 1.0;
    } else {
        while !(excess(hi) > 0.0) {
            hi += 1.0;
            if hi > f64::from(DECADES) {
                return None;
            }
        }
        lo = hi - 1.0;
    }
    for _ in 0..BISECTIONS {
        let mid = 0.5 * (lo + hi);
        if excess(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
        if hi - lo < 1e-14 {
            break;
        }
    }
    Some(10f64.powf(0.5 * (lo + hi)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::CS2;

    fn rheology(fluid: Fluid) -> Rheology {
        Rheology { fluid, rate_scale: 1.0, viscosity_scale: 1.0, tau: [0.505, MAX_TAU], zero: 0.0 }
    }

    /// The flux of simple shear at rate γ̇ relaxed with τ: `P_xy = −2 ρ c_s² τ (γ̇/2)`.
    fn shear(rate: f64, tau: f64, rho: f64) -> [f64; 6] {
        [0.0, 0.0, 0.0, -rho * CS2 * tau * rate, 0.0, 0.0]
    }

    /// In simple shear the cell's τ is the model's at the shear rate it implies.
    #[test]
    fn simple_shear_gets_the_models_viscosity() {
        let fluids = [
            Fluid::PowerLaw { consistency: 0.1, index: 0.5, reference_shear_rate: None },
            Fluid::PowerLaw { consistency: 0.1, index: 1.5, reference_shear_rate: None },
            Fluid::CarreauYasuda { zero_shear_viscosity: 0.5, infinite_shear_viscosity: 0.01, time_constant: 100.0, a: 2.0, index: 0.4 },
            Fluid::Cross { zero_shear_viscosity: 0.4, time_constant: 50.0, index: 0.3 },
        ];
        for fluid in fluids {
            let r = rheology(fluid);
            for (rate, rho) in [(1e-3, 1.0), (1e-2, 1.02), (3e-2, 0.97)] {
                let tau = 0.5 + 3.0 * fluid.viscosity(rate, 0.0);
                if !(0.505..=MAX_TAU).contains(&tau) {
                    continue;
                }
                let (got, clamp) = r.relax(&shear(rate, tau, rho), rho, 0.0);
                assert_eq!(clamp, Clamp::None, "{fluid:?} at {rate}");
                assert!((got - tau).abs() < 1e-11 * tau, "{fluid:?} at {rate}: {got} vs {tau}");
            }
        }
    }

    /// A constant viscosity gives the Newtonian τ to round-off.
    #[test]
    fn a_newtonian_limit_relaxes_with_the_newtonian_time() {
        let tau = 0.8;
        let nu = (tau - 0.5) * CS2;
        for fluid in [
            Fluid::PowerLaw { consistency: nu, index: 1.0, reference_shear_rate: None },
            Fluid::CarreauYasuda { zero_shear_viscosity: nu, infinite_shear_viscosity: nu, time_constant: 3.0, a: 2.0, index: 0.3 },
        ] {
            for rate in [0.0, 1e-4, 0.05] {
                let (got, clamp) = rheology(fluid).relax(&shear(rate, tau, 1.0), 1.0, 0.0);
                assert_eq!(clamp, Clamp::None);
                assert!((got - tau).abs() < 1e-14, "{fluid:?} at {rate}: {got}");
            }
        }
    }

    #[test]
    fn rest_and_extreme_shear_clamp_and_say_so() {
        let thinning = rheology(Fluid::PowerLaw { consistency: 0.1, index: 0.5, reference_shear_rate: None });
        assert_eq!(thinning.relax(&[0.0; 6], 1.0, 0.0), (MAX_TAU, Clamp::High), "infinite viscosity at rest");
        assert_eq!(thinning.relax(&shear(1e7, 0.6, 1.0), 1.0, 0.0), (0.505, Clamp::Low));
        let thickening = rheology(Fluid::PowerLaw { consistency: 0.1, index: 1.5, reference_shear_rate: None });
        assert_eq!(thickening.relax(&[0.0; 6], 1.0, 0.0), (0.505, Clamp::Low), "no viscosity at rest");
        // A frozen melt (below D2 − A2) is as viscous as allowed.
        let melt =
            Rheology { zero: 300.0, ..rheology(Fluid::CrossWlf { d1: 1e12, d2: 373.0, a1: 28.0, a2: 51.6, critical_stress: 3e4, index: 0.3 }) };
        assert_eq!(melt.relax(&shear(1e-3, 1.0, 1.0), 1.0, 0.0), (MAX_TAU, Clamp::High));
        // Hostile fluxes keep the run going without a number that is not one.
        for p in [[f64::NAN; 6], [f64::INFINITY; 6]] {
            let (tau, _) = thinning.relax(&p, 1.0, 0.0);
            assert!(tau.is_finite());
        }
        assert!(thinning.relax(&shear(1e-2, 1.0, 1.0), 0.0, 0.0).0.is_finite());
    }

    #[test]
    fn the_stress_inverts_to_its_shear_rate() {
        let fluid = Fluid::PowerLaw { consistency: 2.0, index: 0.5, reference_shear_rate: None };
        // K γ̇ⁿ = 8 at γ̇ = 16.
        let rate = shear_rate_at_stress(&fluid, 8.0, 0.0).unwrap();
        assert!((rate - 16.0).abs() < 1e-9, "{rate}");
        let newtonian = Fluid::Newtonian { viscosity: 1e-3 };
        assert!((shear_rate_at_stress(&newtonian, 2.0, 0.0).unwrap() - 2000.0).abs() < 1e-8);
        assert_eq!(shear_rate_at_stress(&fluid, 0.0, 0.0), None);
        assert_eq!(shear_rate_at_stress(&fluid, f64::NAN, 0.0), None);
    }
}
