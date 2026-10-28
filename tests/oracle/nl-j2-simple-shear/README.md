# Non-coaxial logarithmic plasticity in simple shear

In SI units, a one-metre hex8 cube follows F = I + γ e₁ ⊗ e₂ from γ = 0 to 8. The bottom is fixed in x and y,
the top moves in x while y is held, and z = 0 is a symmetry plane. The upper z face is
traction-free. E = 2·10⁷ Pa, ν = 0.3, Y = 250 Pa and H = 0 approach rigid-perfect plasticity.
`kinematic` selects the additive Lagrangian logarithmic formulation, even at H = 0.

The independent reference is M. Itskov, *On the applicability of generalized strain measures
in large strain plasticity*, in *Trends and Applications of Mathematics to Mechanics*,
pp. 101–113 (2005), [doi:10.1007/88-470-0354-7_10](https://doi.org/10.1007/88-470-0354-7_10).
Equation (60), with f₀(Γ) = ln Γ from (45), gives the r = 0 curve of Figure 1:

```text
s = 4 + γ²,  e = asinh(γ/2)
τ₁₂/(Y/√3) = 2 √(s + 4e²)/s.
```

This response decreases smoothly with large shear, unlike multiplicative perfect plasticity.
The rotating principal axes distinguish the formulations, which agree in a coaxial test.
The top-face x reaction is the shear stress since its reference area and J are both one.
The four samples are evaluated from the published equation, not from solver results or a
plot digitisation. The 0.005% reference tolerance allows Δγ = 1/512 integration error and
Y/E = 1.25·10⁻⁵ elasticity; material-point tests verify step and rigidity refinement, also with
nonzero Prager hardening. Once checked, the fixture keeps the usual 10⁻⁹ relative band.
