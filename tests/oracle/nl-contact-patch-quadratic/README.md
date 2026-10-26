# Contact patch test with quadratic faces (implicit mortar contact, non-matching meshes)

The case of [nl-contact-patch](../nl-contact-patch/) with hex20 elements: two unit cubes stacked
along z, meshed independently, 3 × 3 × 2 hex20 below (graded in x and y) and 2 × 2 × 2 above, so
the nodes of the interface at z = 1, corners and mid-side nodes, do not match. A frictionless
contact pair joins the lower block's top (main) and the upper block's bottom (secondary). The
faces x = 0 and y = 0 of both blocks are symmetry planes, the base z = 0 is held in z, and the top
z = 2 is pushed down by 0.02 (E = 1000, ν = 0.3) in two nonlinear static steps: `Small` without
nonlinear geometry (small sliding) and `Finite` with it (finite sliding). The exact state is
uniaxial stress, uniform in both blocks, with a uniform contact pressure equal to it, at the
quad8 faces' corners as at their mid-side nodes.
