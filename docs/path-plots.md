# Path plots

After solving a solid step, open **Results > Path Plot**. Add reference-coordinate points by
typing their coordinates or enabling **Pick surface points** and clicking the visible mesh.
The editor keeps SI coordinates internally: changing display units does not round a pick.
Picking a deformed or curved face maps the actual displayed triangle back to the element's
reference geometry. It does not snap to a nearby node.

Choose a nodal field and the number of equally spaced arc-length stations, then **Sample
path**. Sampling runs in a cancellable background job. The chart retains the sampled step,
field and frame or mode independently of the current contour. Hovering it marks the station
in the viewport. Its overlay follows the viewport's reference or deformed configuration;
coordinates and distances in the query remain reference coordinates. **Export CSV** samples
the draft through the same command and writes the result. Cancellation or an error preserves
the preceding successful plot.

## Commands and units

`results.path` works headlessly, through MCP and in the desktop control channel. Programmatic
calls never open a picker. For example, on the built-in solved cantilever:

```json
{"step":"Load","field":"stress_xx","points":[["20 mm","5 mm","20 mm"],["80 mm","5 mm","20 mm"]],"samples":101,"frame":null,"mode":null,"path":"top-fibre.csv"}
```

The points define a polyline with 2–1000 vertices; repeated vertices are allowed, but its
total length must be positive. `samples` is 2–10000, default 101. Bare point coordinates use
the session's display length unit. Explicit unit strings avoid dependence on that setting.
The reply carries `distance`, `xyz`, `values`, containing `elements`, `units`, `field_title`,
`frame`, `mode` and `configuration: "reference"`. Values and coordinates are converted once
at the display boundary. Modal displacement fields retain the existing dimensionless
normalisation rather than being labelled as a physical length.

The optional CSV has exactly `distance,x,y,z,value` as its header and uses the same numeric
values as the reply, CRLF line endings and empty values at gaps. Units are supplied by the
reply and the chart, rather than extra CSV columns. Writes respect the session's file roots.

## Interpolation and limits

Sampling uses the containing solid's shape functions: tet4, tet10, wedge6, pyramid5, hex8 and
hex20, including curved quadratic geometry. Shell, beam, rod and 2D elements in the same
model are left out: a point is located in the solids only. The lowest element label wins on
shared faces.
This interpolates the existing nodal field; it does not recompute stresses at integration
points or make a discontinuous stress field continuous. Numerical contact fields are
available where their contributing nodal values are defined. Categorical `contact_status`
is refused because fractional status codes would be misleading.

Proven outside points and undefined nodal contributions give `null` and visible chart gaps;
no nearest-node extrapolation or gap bridging is used. An inverse search that cannot establish
containment or exclusion within its work limit returns an error, rather than calling the
point outside. A straight segment between two picks on a concave surface can leave the solid.
Positive local Jacobian checks do not certify that an entire element is globally injective.

Stale results are refused. Omitted/null `frame` selects the final state; stored frames start
at zero. Nonlinear stored frames give displacements, and stress and contact fields where the
step saved them (a frame that did not is refused, never given the final state's); temperature
requires the final state. Heat flux is final-only;
stored heat frames provide temperatures. Frequency results use `mode` from one (default one)
and reject frames. Other procedures reject non-null modes. Fluid cell fields are outside this
command's scope.

The index is cached for an exact result and model definition. Construction is bounded by
500000 mesh elements, 2000000 result nodes, 96 MiB index storage and 14000000 basis/enclosure
work units. A query separately permits 4000000 basis/enclosure work units, 2000000 BVH visits
and 4096 candidate elements per station. Index construction and sampling check cancellation
and report progress. Higher-order element bounds use conservative Bernstein enclosures,
including outward rounding; a nodal bounding box alone is insufficient for curved elements.
