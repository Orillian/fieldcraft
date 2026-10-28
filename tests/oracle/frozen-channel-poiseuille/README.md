# A channel frozen on its cold quarter: Poiseuille flow in the liquid layer

A gap H = 1 between a wall held at −1 °C (y = 0) and one held at +1 °C (y = H), filled with a
substance that melts at −0.5 °C (L / c = 1). At steady state the temperature is linear across the
gap whatever the latent heat, so the substance is frozen on the cold quarter, y < y* = H / 4. A body
force a = 10⁻³ along x drives the liquid between the front and the hot wall as plane Poiseuille
flow, u = a (y − y*)(H − y) / (2 ν), largest at the middle of the liquid layer, y = (H + y*) / 2 =
0.625: u_max = a (H − y*)² / (8 ν) = 7.03125·10⁻⁴. The frozen layer holds still.

The model is dimensionless in SI: ρ = 1, c = 1, k = 0.05 (α = 0.05), ν = 0.1, everything starting
liquid at 0 °C and at rest.

**1D in a 3D solver:** the channel is one cell wide and periodic across x and z.

**Measures:** the frozen share is the solve's `frozen` (the mean solid fraction; H = 1), the
velocity profile a `solve.probe` line across the gap through the cell centres, its largest value
and its height located by the parabola through the largest sample and its neighbours.

**Run:** 64 cells across, Δt = Δx² / (8 α) (thermal relaxation time 1), to t = 40 (65536 steps):
twenty times the gap's thermal time H² / (π² α) and seventy times the liquid layer's viscous
time (H − y*)² / (π² ν), so both the front and the flow are steady.

**Grid refinement** (Δt scaled with Δx², the same end time):

| Quantity | 16 cells | 32 | 64 (this case) | 128 |
|---|---|---|---|---|
| frozen share | 0.25 | 0.25 | 0.25 | 0.25 |
| u_max, error | +5.81 % | +2.92 % | +1.47 % | +0.74 % |
| height of u_max, error | −1.1·10⁻² | −5.4·10⁻³ | −2.7·10⁻³ | −1.4·10⁻³ |

The front lies on a cell face on every grid, so each cell is wholly frozen or wholly liquid. The
flow's error halves with the cell (first order): the frozen cells are held still at their centres
by a drag rather than by a wall on the front, so the effective no-slip wall sits inside the frozen
layer, widening the liquid gap by about 0.35 Δx on every grid.
