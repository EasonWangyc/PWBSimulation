# PWB Complex Centerline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a planar complex-centerline PWB geometry mode for PD-PWB-SMF while preserving the existing `_1/_2` simple bend code.

**Architecture:** Keep the current `pwb_core.py` style, but add reusable pure-geometry helpers for Bezier centerlines, normalized path length, radius profiles, and FDTD cylinder segment placement. Add a small setup script that saves the complex geometry for visual inspection without running FDTD.

**Tech Stack:** Python, NumPy, Lumerical/Ansys `lumapi` for geometry generation, standard-library `unittest` for pure geometry tests.

---

## File Structure

- Modify `PD-PWB-SMF/pwb_core.py`: add complex path parameters, pure path/radius helpers, generic centerline segment placement, `_3` generation/setup/data/visualization wrappers.
- Create `PD-PWB-SMF/test_complex_geometry.py`: pure unit tests for path shape, normalized arc length, radius profile, and segment descriptor generation.
- Create `PD-PWB-SMF/test_setup_complex.py`: Lumerical entry point that builds and saves the complex geometry without running a full simulation.

## Tasks

### Task 1: Add Pure Geometry Tests

**Files:**
- Create: `PD-PWB-SMF/test_complex_geometry.py`

- [ ] **Step 1: Write failing tests**

Create tests that import `pwb_core`, generate the complex path, and validate shape, endpoints, radius profile, and segment descriptors.

- [ ] **Step 2: Run tests and verify RED**

Run: `python PD-PWB-SMF/test_complex_geometry.py`

Expected: FAIL because `generate_pwb_path_3`, `generate_radius_profile_3`, and `centerline_to_segments` do not exist yet.

### Task 2: Implement Complex Centerline Helpers

**Files:**
- Modify: `PD-PWB-SMF/pwb_core.py`

- [ ] **Step 1: Add parameters to `PWBParameters`**

Add Bezier control parameters and segment count defaults.

- [ ] **Step 2: Add pure helper functions**

Implement:

- `cubic_bezier(points, t)`
- `path_arc_length(path)`
- `normalized_path_position(path)`
- `generate_pwb_path_3(params)`
- `generate_radius_profile_3(params, path)`
- `centerline_to_segments(path, radii)`

- [ ] **Step 3: Run tests and verify GREEN**

Run: `python PD-PWB-SMF/test_complex_geometry.py`

Expected: PASS.

### Task 3: Add FDTD Structure and Setup Entry Points

**Files:**
- Modify: `PD-PWB-SMF/pwb_core.py`
- Create: `PD-PWB-SMF/test_setup_complex.py`

- [ ] **Step 1: Add generic FDTD segment placement**

Implement `add_centerline_segments(fdtd, segments, material="Vancore B", name_prefix="PWB_complex")`.

- [ ] **Step 2: Add `_3` wrappers**

Implement:

- `generate_pwb_structure_3(fdtd, params)`
- `setup_fdtd_simulation_3(fdtd, params)`
- `get_data_3(fdtd, params)`
- `visualize_and_save_results_3(fdtd, params)`

- [ ] **Step 3: Add complex setup script**

Create `PD-PWB-SMF/test_setup_complex.py` that imports `lumapi`, builds `_3`, sets up FDTD, saves a `.fsp`, and exits without `fdtd.run()`.

- [ ] **Step 4: Syntax-check implementation**

Run: `python -m py_compile PD-PWB-SMF/pwb_core.py PD-PWB-SMF/test_complex_geometry.py PD-PWB-SMF/test_setup_complex.py`

Expected: PASS.

### Task 4: Final Verification

**Files:**
- Verify only.

- [ ] **Step 1: Run pure geometry tests**

Run: `python PD-PWB-SMF/test_complex_geometry.py`

Expected: PASS.

- [ ] **Step 2: Do not run Lumerical by default**

Do not run `test_setup_complex.py` unless the user explicitly wants to launch Lumerical and save the `.fsp`.

