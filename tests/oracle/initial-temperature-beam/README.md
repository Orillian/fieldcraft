# A beam ejected with a temperature linear through its depth (#160)

A steel beam 200 × 20 × 10 of 8 × 2 × 2 hex20 elements (E = 200000, ν = 0.3, α = 1.2·10⁻⁵).
A heat transfer step `Mould` holds its bottom face at 100 and its top face at 140, so the
temperature is T₀ = 100 + g z with g = 4 per length. `model.initial_field.map` takes that
result as the initial temperature `Ejected`. The static step `Cool` then holds the beam
statically determinately at its x = 0 end (corner a fixed, b at (0, 20, 0) held in x and z, c at
(0, 0, 10) held in x) and brings it to a uniform 20 by a temperature load.

The thermal strain α (20 − T₀(z)) is isotropic and linear in z, so it is compatible
(Timoshenko and Goodier §150): the beam is left stress-free, bent to the curvature α g, and its
free end rises by α g L²/2 = 0.96 while the bottom fibre shortens by α · 80 · L = 0.192. Both
fields are polynomials the hex20 elements represent exactly, so the tolerances are rounding-level.
