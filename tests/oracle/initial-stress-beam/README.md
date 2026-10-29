# A residual stress linear through a beam's depth, released (#160)

The beam, material and mould of `initial-temperature-beam`. The static step `Hold` cools it
from the mapped initial temperature T₀ = 100 + 4 z to a uniform 20 with both ends held axially
(symmetry planes at x = 0 and x = 200), the y = 0 face on a symmetry plane and corner a held in
z. Every fibre keeps its length and contracts freely sideways, so the held beam carries
σ_xx = −E α (20 − T₀(z)) = 192 + 9.6 z and nothing else.

`model.initial_field.map` takes that stress as the initial stress `Residual`; the initial
temperature is deleted, and the static step `Release` holds the beam statically determinately
(as in `initial-temperature-beam`) with no loads. The free beam's equilibrium with σ = C ε + σ₀
is the strain ε = −C⁻¹σ₀, linear in z and so compatible (Timoshenko and Goodier §150): the
residual stress is released entirely, leaving no stress, and the beam takes the shape the
eigenstrain gave in `initial-temperature-beam` (tip deflection 0.96, bottom-fibre shortening
0.192). In the terms of `crates/solver/tests/initial.rs` this is the linear residual stress
s (2z/H − 1) about the mean 240, with s = 48 and curvature 2 s/(E H).
