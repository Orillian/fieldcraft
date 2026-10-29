# Nonlinear patch test, incompatible-mode hex8

The mesh, boundary and material of nl-patch-hex8 with incompatible-mode hex8 elements: a unit
cube of 3 × 3 × 3 elements whose non-corner nodes are moved off the lattice, so no element is a
parallelepiped. Every boundary node is prescribed u = (F − I) X with
F = [[1.3, 0.2, -0.1], [-0.15, 0.9, 0.25], [0.05, -0.2, 1.1]] (stretch, shear and rotation,
J = 1.389000); the eight inner nodes are free. The neo-Hookean material (E = 100, ν = 0.3) must
take the homogeneous deformation with every enhanced mode at zero: the inner node 22 moves by
(F − I) X and the Cauchy stress is uniform. The mode gradients are mapped with the Jacobian at
the element centre and scaled by det J₀/det J (Taylor, Beresford & Wilson 1976), which keeps
their integral over a distorted element zero, the condition for passing.
