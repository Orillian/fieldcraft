# Hertz line contact on tet10 faces: cylinder on a flat (implicit mortar contact, curved faces)

The model of [nl-hertz-cylinder](../nl-hertz-cylinder/) meshed with tet10: the same steel cylinder of
radius R = 10 mm on a steel block (E = 210 GPa, ν = 0.3; mm, N, MPa), a quarter model one element
thick with every node held in z (plane strain), loaded by the same pressure for Hertz's half-width
b = 0.5 mm and peak pressure p₀ = 2884.6 MPa. The blocks have half the divisions of the hex8 case
in x and y (10 × 5 on the indenter, 12 × 7 on the foundation, the same grading), and `mesh.block`
puts the indenter's mid-side nodes on the circular arc, so the secondary faces are curved and the
facets of the mortar segmentation are not the faces. One nonlinear static step without nonlinear
geometry: small sliding, mortar contact of the round face (secondary) on the flat (main) by an
augmented Lagrangian.

Checked against the closed form: the peak pressure p₀, the force across the contact, the half-width
(the last secondary node inside b, at x = 0.42 = 0.85 b, is closed; the next one, at x = 0.51 = 1.03
b, is open by 6.1·10⁻⁴ mm, which no quantity checks: Hertz's separation there, 1.1·10⁻⁴ mm, is too
small for these faces) and, beyond the contact, zero pressure and Hertz's separation of the surfaces
at x = 0.72 = 1.44 b.
