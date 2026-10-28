# Finite-strain kinematic hardening cycle

A unit hex8 with E = 200000, ν = 0.3, yield stress 250 and Prager modulus
H = 2000 is stretched to an axial logarithmic strain of 0.2 and returned to
zero axial displacement. Its sides contract freely. The prescribed peak
stretch is exp(0.2), and increments end at the reversal.

For `E_t = E H/(E + H)`, the peak Kirchhoff stress is
`τ_a = σ_y + E_t (ε_a − σ_y/E)` and the residual is
`τ_r = τ_a − 2σ_y − E_t (ε_a − 2σ_y/E)`.
The reported Cauchy stress is `σ_r = τ_r/exp(τ_r/(3K))`.
The accumulated plastic strain is the sum of the forward and reverse increments:
`ε̄_p = (E ε_a − σ_y)/(E + H) + (E ε_a − 2σ_y)/(E + H)`.
`expected.json` records both references and their tolerances.
