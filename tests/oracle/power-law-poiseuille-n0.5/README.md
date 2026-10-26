# Plane Poiseuille flow of a power-law fluid, n = 0.5

A shear-thinning power-law fluid, η = K γ̇^(n−1), between two plates H = 1 m apart, driven by
a uniform body acceleration G = 0.4 m/s² along x, with K = 0.2 Pa·sⁿ (at γ̇_ref = 1/s) and ρ = 1 kg/m³. The
closed form (Bird, Armstrong & Hassager 1987, §4.2; from the stress ρ G (h − y) and the power
law, h = H/2):

u(y) = n/(n+1) (ρG/K)^(1/n) (h^((n+1)/n) − |y − h|^((n+1)/n)),

peak 1/6 m/s; the wall shear rate (ρ G h / K)^(1/n) is 1/s, where the viscosity is 0.2 Pa·s.

**2D in a 3D solver:** the domain is one cell long and one cell thick, periodic along x and z
(the lattice tests' convention), with walls at y = 0 and y = H. The solver takes each cell's
viscosity at the cell's own shear rate from the non-equilibrium momentum flux
(`fieldcraft_lbm::rheology`).

**Run:** 32 cells across, a time step of Δx² (τ = 1.1 at the walls), TRT, for 48 s (about ten
viscous times at the wall viscosity); 49 152 steps, under 2 s in release.
The middle of a shear-thinning channel is at rest, where the power law's viscosity is
infinite: the 4 middle cells end the run at the largest τ the solver allows (5), and the solve
says so in its warnings; the stress there is small, so the cap moves the profile little.

**Measures:** `solve.probe` at three cell centres (the probe reads the cell a point lies in):
beside the mid-plane, beside the quarter of the gap and the second cell from the lower plate;
each reference is the closed form at that centre.

**Grid refinement** (the same case at 16 and 64 cells across, the same three cells counted from
the plate and the mid-plane; ours, with its error from the closed form at its centre):

| Quantity | y (32 cells) | closed form | 32 cells (this case) | 16 cells | 64 cells |
|---|---|---|---|---|---|
| u_centre | 0.484375 | 0.166662 | 0.16669 (+2.5e-5) | 0.16616 (−4.6e-4, y = 0.46875) | 0.16682 (+1.5e-4, y = 0.4921875) |
| u_quarter | 0.234375 | 0.141678 | 0.14147 (−2.1e-4) | 0.13617 (−8.3e-4, y = 0.21875) | 0.14377 (−5.1e-5, y = 0.2421875) |
| u_wall | 0.046875 | 0.042618 | 0.042381 (−2.4e-4) | 0.076355 (−9.2e-4, y = 0.09375) | 0.022296 (−6.0e-5, y = 0.0234375) |

The largest error is 1.4e-3 of the peak at 32 cells, 5.5e-3 at 16 and 9e-4 at 64; the
tolerance, 2e-3 of the peak, bounds the 32-cell error.
