# PD-PWB-SMF Complex Centerline Design

## Goal

Add a new PD-PWB-SMF simulation path mode for a more complex 2D PWB bending structure. The new mode should support a smooth centerline in the `x-z` plane while keeping the existing Lumerical/FDTD workflow and current `_1/_2` code available for comparison.

This design targets the structure shown by the user as a more complex version of the existing simple 90-degree bend. The path remains planar: `y = 0`.

## Current State

`PD-PWB-SMF/pwb_core.py` currently has two main path/structure variants:

- `_1`: builds only the simple bending section.
- `_2`: builds the bending section plus the vertical taper2 section.

Both variants use a hard-coded centerline shape:

```text
straight taper1 + quarter-circle bend + vertical taper2
```

The structure generation also assumes fixed index ranges for the bend. That makes it difficult to reuse the code for arbitrary smooth centerlines.

## Recommended Approach

Add a third variant based on a reusable centerline builder:

- `generate_pwb_path_3(params)`: generate a smooth planar complex centerline.
- `generate_pwb_structure_3(fdtd, params)`: build the FDTD geometry by placing short circular cylinders along the full centerline.
- `setup_fdtd_simulation_3(fdtd, params)`: set the simulation region and monitors from path bounds instead of hard-coded bend bounds.
- `test_setup_complex.py`: save the geometry without running FDTD, for fast visual inspection in Lumerical.

The existing `_1/_2` functions should remain unchanged unless a small helper can be reused without changing their behavior.

## Centerline Model

Use a cubic Bezier or piecewise Bezier centerline for the complex bend. The first implementation should stay dependency-light and use only `numpy`.

The path should still be composed conceptually as:

```text
input taper/straight section -> smooth complex bend -> output vertical section
```

The bend control points should be derived from existing physical parameters where possible:

- `L`: total horizontal extent.
- `h`: total vertical drop.
- `l1`: input straight/taper section length.
- `R`: approximate bend scale or clearance parameter.
- `curve_points`: sampling density.

New optional parameters can be added to `PWBParameters`:

- `bend_shape`: normalized horizontal control distance for the bend.
- `bend_lift`: normalized upward or outward overshoot of the bend before it drops.
- `complex_segments`: number of cylinder segments used for geometry generation.

All generated coordinates remain in SI units.

## Radius Model

Use a separate radius function along normalized path position `s`:

```text
s = 0 at path start
s = 1 at path end
```

The initial implementation should support three regions:

- Input taper: radius changes from `r1` to `r`.
- Main bend: radius stays close to `r`.
- Output taper: radius changes from `r` to `r2`.

This keeps the radius behavior explicit and easier to modify later. The radius model should be independent of the centerline model, so future work can use arbitrary imported centerline points without rewriting radius logic.

## Geometry Generation

Create a helper that converts any sampled centerline into FDTD cylinder segments:

```python
add_centerline_segments(fdtd, path, radii, material="Vancore B", name_prefix="PWB")
```

For each adjacent point pair:

1. Compute segment center.
2. Compute segment length.
3. Compute local direction.
4. Rotate the cylinder so its local axis follows the direction.
5. Use the average radius for that segment.

This generalizes the current hard-coded bending segment loop.

## Simulation Setup

`setup_fdtd_simulation_3` should compute simulation bounds from `generate_pwb_path_3(params)`:

- `x min/max` from path coordinates plus margin.
- `z min/max` from path coordinates plus margin.
- `y min/max` centered around zero with margin based on the maximum radius.

The first implementation can place:

- Source near the start of the path, still using the current mode source style.
- Output monitor near the end of the path.
- Transmission monitor over the full path bounding box.

Monitor orientation is the main risk area because the path end is vertical while the source is horizontal. The first version should prioritize geometry generation and setup inspection, then refine monitor placement after checking the saved `.fsp`.

## Files To Change

Planned files:

- `PD-PWB-SMF/pwb_core.py`
- `PD-PWB-SMF/test_setup_complex.py`

Optional later files:

- `PD-PWB-SMF/run_complex.py`
- `PD-PWB-SMF/sweep_complex.py`

## Verification

Lightweight verification:

- `python -m py_compile PD-PWB-SMF/pwb_core.py PD-PWB-SMF/test_setup_complex.py`

Geometry verification when Lumerical is available:

- `python PD-PWB-SMF/test_setup_complex.py`
- Open the saved `.fsp` and visually inspect whether the generated structure matches the intended complex 2D bend.

Full FDTD runs should wait until the geometry and monitor placement are checked.

## Open Implementation Notes

- Keep `_1/_2` behavior stable.
- Avoid adding `scipy` for the first version.
- Keep path parameters in SI units, while comments and output filenames can mention micrometers.
- Treat the existing absolute path layout as real for this project unless the user asks to refactor paths.

