# Hertz line contact: cylinder on a flat (implicit mortar contact)

A steel cylinder of radius R = 10 mm pressed on a steel block (E = 210 GPa, ν = 0.3; mm, N, MPa)
in plane strain: a quarter model of the cylinder's lower part (a 2.5 mm wide block whose bottom
is the circular arc, graded towards the contact) on a 5 × 5 mm foundation, symmetric about
x = 0, one element thick with every node held in z. A pressure on the indenter's top carries
the load P = 2265.57 N/mm (per length of the whole cylinder) chosen so that Hertz's half-width is
b = 0.5 mm and the peak pressure p₀ = 2884.6 MPa. Soft springs (0.01 N/mm per node) hold the
indenter's top vertically until the contact closes. One nonlinear static step without nonlinear
geometry: small sliding, mortar contact of the round face (secondary) on the flat (main) by an
augmented Lagrangian.

`normal_force` sums nodal contact forces projected onto their local overlapping-pair normals.
In this small-sliding case those normals follow the **reference cylinder**, while the vector
`force` is the resultant to compare with vertical equilibrium. The two quantities differ on a
curved contact; the grounded hold springs also carry part of the applied load. Hertz's
reference load is a vertical resultant; it is retained with its original
0.2% tolerance. A continuum estimate explains the size of the difference: with
`p(x) = p₀ sqrt(1 - x²/b²)`, the pressure-weighted mean of `x²` is `b²/4`, so projecting a
radius-`R` curved normal onto the vertical reduces its component by about `b²/(8R²)`.
Here that is 0.03125%. This is an estimate of the geometric distinction, not a replacement
reference or a change to the solver's equilibrium equations.

The reporting correction excludes uncontacted faces from the nodal normal and uses the same
reference-pair normals as the small-sliding mortar constraint. The former report used all
secondary faces' **deformed** area vectors. Its slightly smaller scalar sum was closer to the
vertical load by coincidence; assembly, displacements, stresses and the vector resultant are
unchanged. The refrozen values record this deliberate reporting definition, with the existing
Hertz values and tolerance unchanged.
