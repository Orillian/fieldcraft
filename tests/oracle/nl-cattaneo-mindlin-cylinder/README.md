# Cattaneo–Mindlin partial slip of a cylinder on a flat (implicit mortar contact with friction, plane strain)

The Hertz line contact of `nl-hertz-cylinder-plane-strain` (steel cylinder of radius R = 10 mm on
a steel block, E = 210 GPa, ν = 0.3; mm, N, MPa; plane strain, unit thickness, quad4), now with
Coulomb friction μ = 0.3 and a tangential force Q below μP. The load per length P = 2265.57 N/mm
gives Hertz's half-width a = 0.5 mm and p₀ = 2884.6 MPa. Cattaneo (1938) and Mindlin (1949):
the contact sticks over |x| < c = a √(1 − Q/(μP)) and slips at its edges, where q = μp; in the
stick zone q = μp₀ (√(1 − x²/a²) − (c/a) √(1 − x²/c²)) (Johnson 1985, §7.2).

The whole contact is modelled (no symmetry, since the shear is not symmetric): the cylinder's
lower part is three blocks 2.5 mm to each side and 3 mm tall with the circular arc at the bottom,
on a 10 × 5 mm foundation of three blocks. The middle blocks (|x| ≤ 0.75) are uniform, h = 0.025
= a/20, and the side blocks grade from h outwards: 5000 quad4, 5252 nodes. The round edge is the
secondary of a frictional pair on the flat; soft springs hold the indenter's top until the contact
closes.

Two nonlinear static steps without nonlinear geometry, `Shear50` (Q = 0.5 μP) and `Shear75`
(Q = 0.75 μP), each with per-load paths (`load_paths`): the normal pressure on the top ramps to
P over the first half of the step and holds, while the shear traction on the top grows from 0 to
Q over the second half. The shear on the top, 3 mm above the contact, would also turn the
indenter by Q·H; a couple of vertical tractions on the outer parts of the top edge (the side
blocks') carries that moment, so the contact sees the pressure of P and the force Q alone.

| h (middle) | Q/μP | slipping / closed nodes | 1 − c/a (0.2929, 0.5) | q(0)/p(0) (0.0879, 0.15) | p₀ (2884.6) |
|---|---|---|---|---|---|
| 0.05 | 0.5 | 5 / 19 | 0.263 | 0.0907 | 2887.4 |
| 0.025 | 0.5 | 12 / 39 | 0.308 | 0.0906 | 2887.9 |
| 0.0125 | 0.5 | 24 / 79 | 0.304 | 0.0904 | 2888.1 |
| 0.05 | 0.75 | 10 / 19 | 0.526 | 0.1551 | 2887.7 |
| 0.025 | 0.75 | 19 / 39 | 0.487 | 0.1544 | 2888.2 |
| 0.0125 | 0.75 | 40 / 79 | 0.506 | 0.1543 | 2888.4 |

The table was measured with the earlier node normal (its faces' area vectors summed). Since #731
the node normal comes from its overlaps, and both the friction law's tangent plane and the report
use it (see below). At h = 0.025 that moves q(0)/p(0) to 0.0901 (Q = 0.5 μP) and 0.1534
(Q = 0.75 μP), and p₀ to 2888.6 at Q = 0.75 μP. The slipping counts are unchanged. The other
meshes were not rerun.

The case uses h = 0.025 (about 25 s for both steps in a release build). The share of slipping
nodes is within one node of 1 − c/a and settles on it with the mesh; the shear at the centre stays 2–3% above μ(1 − c/a) on every mesh: the bodies are a few
half-widths deep rather than half-planes, and the stick zone sits slightly towards the leading
edge (at h = 0.025 and Q = 0.5 μP it spans −0.30 ≤ x ≤ 0.35 against |x| < 0.354). The ratio of
the tangential to the normal force across the contact is not among the quantities: the soft
springs hold the indenter, so equilibrium makes the contact carry the applied Q whatever its
stick and slip zones are, and Q/P would only restate the applied ratio.

Before #731, in the slip zones the probe's nodal shear over μ times its nodal pressure was
0.92–0.95 at the trailing edge and 1.07–1.09 at the leading edge on every mesh (mean 1), a
departure of about x/(2μR). The cause was the frame of the report. Without nonlinear geometry
the contact is small sliding: the normal force and the slip operator P = I − n n use each
secondary face's reference normal, tilted by x/R on the round edge. The probe split the nodal
force with the node's current normal, built from the deformed faces, which tilts by about x/(2R)
on an interface of two bodies of one material. That share of the pressure, p x/(2R), landed in
the reported shear, positive at one edge and negative at the other. Split with the reference
normal (same probe output, h = 0.025), q/(μp) at the slipping nodes was 0.9997–1.0004.

Since #731 a node's normal is its overlaps' normals weighted by its support: in small sliding,
the reference pairs' normals. The report and the friction law's projection of the committed
traction both use it, so the probe's per-node shear and pressure are split in the frame the law
enforces the cone in.
