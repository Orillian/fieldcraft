# Large-deflection cantilever (elastica), incompatible-mode hex8

The cantilever of nl-elastica-cantilever (L = 10, b = 1, h = 0.5, L/h = 20; E = 10⁴, ν = 0, an
elastic material in a nonlinear step: St Venant–Kirchhoff, exact for large rotations with these
small strains) meshed with 40 × 1 × 1 incompatible-mode hex8 elements: one element through the
depth, where a fully integrated hex8 locks in shear (−9.3 % on the deflection at P L²/(E I) = 1).
Clamped at x = 0, a dead vertical tip force P spread over the end face, four steps at
P L²/(E I) = 1, 2, 5 and 10 (tip rotations 0.46 to 1.43 rad). The tip's neutral axis is the mean
of the bottom and top tip nodes (10, 0, 0) and (10, 0, 0.5), compared with the elastica
(Bisshopp–Drucker).

Refinement (relative error of δ_v / δ_u against the elastica, measured on this recipe with the
divisions changed):

| mesh        | P L²/EI = 1     | 2               | 5               | 10              |
|-------------|-----------------|-----------------|-----------------|-----------------|
| 20 × 1 × 1  | 0.096 / 0.215 % | 0.141 / 0.290 % | 0.224 / 0.374 % | 0.348 / 0.420 % |
| 40 × 1 × 1  | 0.144 / 0.296 % | 0.191 / 0.374 % | 0.278 / 0.459 % | 0.413 / 0.509 % |
| 80 × 1 × 1  | 0.156 / 0.317 % | 0.204 / 0.396 % | 0.293 / 0.482 % | 0.432 / 0.535 % |
| 40 × 1 × 2  | 0.142 / 0.293 % | 0.187 / 0.365 % | 0.271 / 0.445 % | 0.404 / 0.491 % |

The error grows slightly as the mesh is refined along the length and barely moves through the
depth: the solid converges to a continuum answer just above the elastica (the hex20 case
nl-elastica-cantilever gives 0.179 % at P L²/(E I) = 1), the difference being shear deformation
(0.6 (h/L)² ≈ 0.15 % of the deflection) and the clamp's local restraint, both absent from the
elastica. The coarse mesh's slight stiffness offsets part of that difference.
