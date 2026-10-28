# Two-phase Stefan problem (Neumann's solution)

A liquid at +0.5 °C, above its melting point of 0 °C, fills a column whose bottom wall (y = 0) is
suddenly held at −1 °C. The liquid freezes from the wall and the front moves as s = 2 λ √(α t),
with λ the root of Neumann's transcendental equation for equal properties in both phases,
St_s e^{−λ²} / erf λ − St_l e^{−λ²} / erfc λ = λ √π, where St_s = c (T_m − T_w) / L = 0.5 and
St_l = c (T_i − T_m) / L = 0.25 (Carslaw and Jaeger, §11.2). The temperatures in the solid and the
liquid are erf profiles through the wall, the front and the far field.

The model is dimensionless in SI: ρ = 1, c = 1, k = 0.05 (α = 0.05), L = 2, ν = 0.1 (the liquid
never moves: no force acts on it). λ = 0.3834098626712249 (bisection to round-off).

**1D in a 3D solver:** the column is one cell wide and periodic across x and z, 2 long so that its
far end stays at the initial temperature over the run (the cold reaches about 1).

**Measures:** the solve's `frozen` is the mean solid fraction over the fluid cells, so the frozen
thickness Σ (1 − f_l) Δx is twice it. The temperatures are read with `solve.probe` at two cell
centres, y = 0.1484375 (in the solid) and y = 0.5078125 (in the liquid).

**Run:** 64 cells per unit (128 along the column), Δt = Δx² / (8 α) (thermal relaxation time 1),
5000 steps to t = 3.0517578125, when the front is near 0.3.

**Grid refinement** (Δt scaled with Δx², the same end time). `solve.probe` reads the cell that
contains the point, so on the other grids each value is compared with Neumann's solution at that
cell's centre:

| Quantity | 32 cells per unit | 64 (this case) | 128 |
|---|---|---|---|
| front, error | +0.94 % | +0.27 % | +0.19 % |
| T in the solid, error | −2.3·10⁻³ | +3.6·10⁻⁴ | < 10⁻⁵ |
| T in the liquid, error | +3.2·10⁻³ | +2.2·10⁻³ | +1.1·10⁻³ |

The lattice's own test (`crates/lbm/tests/solidification.rs`) checks the same problem along the
whole profile and the latent heat's balance.
