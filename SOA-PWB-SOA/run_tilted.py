r"""
Run SOA-PWB-SOA simulation with tilted SOA geometry.

Loads a pre-built tilted-SOA .fsp, builds the PWB on top with
extended reach (to fill the wedge gaps at tilted facets), sets
mesh order so PWB material wins at interface overlaps, then
configures and runs FDTD.

Usage
-----
.. code-block:: powershell

   & "D:\Program Files\Lumerical\python\python.exe" SOA-PWB-SOA\run_tilted.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import (MATERIAL_DB, SOA_DIR, SOA_RESULTS_DIR,
                        SOA_TILTED_FSP, add_lumerical_api_path)

# Ensure local pwb_core is found first (sim_config inserts PD-PWB-SMF at [0])
sys.path.insert(0, str(Path(__file__).resolve().parent))

add_lumerical_api_path()
import lumapi  # noqa: E402

from pwb_core import (  # noqa: E402
    SOAPWBParams,
    _set_segment_rotation,
    get_data,
    visualize_and_save_results,
)

# ---------------------------------------------------------------------------
# Tilted-case constants
# ---------------------------------------------------------------------------
PWB_OVERHANG  = 3e-6    # how far PWB extends beyond 0 and 250 μm  [m]
PWB_MESH_ORDER = 3       # PWB mesh order (higher than SOA material → wins at overlap)

# SOA tilt angles (rotation around global Z, in XY plane)
SOA_TILT_DEG    = 7.0    # physical tilt from +X axis  [degrees]
SOA_OUT_ROT     = 180.0 - SOA_TILT_DEG   # ≈ 173° — SOA output waveguide
SOA_IN_ROT      = SOA_TILT_DEG           # ≈   7° — SOA input  waveguide

SAVE_PATH = SOA_DIR / "SOA_PWB_SOA_tilted.fsp"
PICTURES_DIR = SOA_RESULTS_DIR / "Pictures"


# ---------------------------------------------------------------------------
# Extended path & radius helpers
# ---------------------------------------------------------------------------

def _pwb_path_extended(params, overhang):
    """Straight centreline along +X, extending *overhang* beyond both facets."""
    n = params.curve_points
    x0 = -overhang
    x1 = params.total_length + overhang
    x = np.linspace(x0, x1, n)
    return np.column_stack((x, np.zeros(n), np.zeros(n)))


def _radius_at_t(t, params, overhang, get_radius_fn):
    """Evaluate a radius function on the *extended* path's normalised t ∈ [0,1].

    The extended path runs from ``-overhang`` to ``total_length + overhang``.
    """
    total_ext = params.total_length + 2.0 * overhang
    x = t * total_ext - overhang                       # global X on extended path

    if x < 0.0:
        # left overhang: constant r_in / r_in_2
        return get_radius_fn(0.0, params)
    if x < params.taper1_length:
        local = x / params.taper1_length
        return get_radius_fn(local * params.taper1_length / params.total_length, params)
    if x < params.total_length - params.taper2_length:
        # PWB straight region
        mid_t = (params.taper1_length + 1e-12) / params.total_length
        return get_radius_fn(mid_t, params)
    if x < params.total_length:
        local = (x - (params.total_length - params.taper2_length)) / params.taper2_length
        frac = (params.total_length - params.taper2_length) / params.total_length
        base_t = frac + local * (params.taper2_length / params.total_length)
        return get_radius_fn(base_t, params)
    # right overhang: constant r_out / r_out_2
    return get_radius_fn(1.0 - 1e-12, params)


# ---------------------------------------------------------------------------
# PWB segment builder (extended + mesh-order aware)
# ---------------------------------------------------------------------------

def _add_path_segments_tilted(
    fdtd, path, params, overhang,
    start_idx, end_idx, count, name_prefix, mesh_order,
    radius_fn_1, radius_fn_2,
):
    """Like _add_path_segments but uses extended radius functions and sets mesh order."""
    path_length = len(path)
    span = end_idx - start_idx
    if span <= 0:
        return
    seg_len_idx = max(1, span // count)

    for i in range(count):
        idx0 = start_idx + i * seg_len_idx
        idx1 = min(start_idx + (i + 1) * seg_len_idx, end_idx)
        if idx0 >= idx1 or idx1 >= path_length:
            continue

        p0, p1 = path[idx0], path[idx1]
        centre = (p0 + p1) / 2.0
        direction = p1 - p0
        seg_len = np.linalg.norm(direction)
        if seg_len <= 0:
            continue
        direction = direction / seg_len

        t0 = idx0 / (path_length - 1)
        t1 = idx1 / (path_length - 1)
        r1 = (_radius_at_t(t0, params, overhang, radius_fn_1)
              + _radius_at_t(t1, params, overhang, radius_fn_1)) / 2.0
        r2 = (_radius_at_t(t0, params, overhang, radius_fn_2)
              + _radius_at_t(t1, params, overhang, radius_fn_2)) / 2.0

        fdtd.addcircle()
        fdtd.set("name", f"{name_prefix}_{i}")
        fdtd.set("material", params.material_pwb)
        fdtd.set("make ellipsoid", 1 if params.use_ellipsoid else 0)
        fdtd.set("x", centre[0])
        fdtd.set("y", centre[1])
        fdtd.set("z", centre[2])
        fdtd.set("radius", r1)
        if params.use_ellipsoid:
            fdtd.set("radius 2", r2)
        fdtd.set("z span", seg_len)
        fdtd.set("override mesh order from material database", 1)
        fdtd.set("mesh order", mesh_order)
        _set_segment_rotation(fdtd, direction)


# ---------------------------------------------------------------------------
# Structure creation
# ---------------------------------------------------------------------------

def create_tilted_structure(fdtd, params, overhang=PWB_OVERHANG):
    """Load tilted-SOA base .fsp and build the extended PWB on top.

    Returns the extended centreline path (N, 3).
    """
    from pwb_core import get_radius_at_position_1, get_radius_at_position_2

    path = _pwb_path_extended(params, overhang)

    fdtd.deleteall()
    fdtd.importmaterialdb(str(MATERIAL_DB))   # DB first → "ridge" resolved from DB on load
    fdtd.load(str(SOA_TILTED_FSP))

    cp = params.curve_points
    total_ext = params.total_length + 2.0 * overhang

    # Section boundaries on the extended path
    def _x_to_idx(x_global):
        """Map a global X coordinate to a path index."""
        t = (x_global + overhang) / total_ext
        return int(np.clip(np.round(t * (cp - 1)), 0, cp - 1))

    idx_overhang_L_start = 0
    idx_taper1_start     = _x_to_idx(0.0)
    idx_taper1_end       = _x_to_idx(params.taper1_length)
    idx_taper2_start     = _x_to_idx(params.total_length - params.taper2_length)
    idx_taper2_end       = _x_to_idx(params.total_length)
    idx_overhang_R_end   = cp - 1

    mo = PWB_MESH_ORDER
    rf1 = get_radius_at_position_1
    rf2 = get_radius_at_position_2

    # 1 — left overhang (constant r_in)
    _add_path_segments_tilted(fdtd, path, params, overhang,
        idx_overhang_L_start, idx_taper1_start,
        max(1, idx_taper1_start - idx_overhang_L_start),
        "PWB_ovL", mo, rf1, rf2)

    # 2 — taper-1 (r_in → r_pwb)
    _add_path_segments_tilted(fdtd, path, params, overhang,
        idx_taper1_start, idx_taper1_end,
        max(1, idx_taper1_end - idx_taper1_start),
        "PWB_taper1", mo, rf1, rf2)

    # 3 — PWB straight (constant r_pwb)
    _add_path_segments_tilted(fdtd, path, params, overhang,
        idx_taper1_end, idx_taper2_start,
        max(1, idx_taper2_start - idx_taper1_end),
        "PWB_straight", mo, rf1, rf2)

    # 4 — taper-2 (r_pwb → r_out)
    _add_path_segments_tilted(fdtd, path, params, overhang,
        idx_taper2_start, idx_taper2_end,
        max(1, idx_taper2_end - idx_taper2_start),
        "PWB_taper2", mo, rf1, rf2)

    # 5 — right overhang (constant r_out)
    _add_path_segments_tilted(fdtd, path, params, overhang,
        idx_taper2_end, idx_overhang_R_end,
        max(1, idx_overhang_R_end - idx_taper2_end),
        "PWB_ovR", mo, rf1, rf2)

    return path


# ---------------------------------------------------------------------------
# FDTD simulation setup (adapted for tilted case — wider Y margins)
# ---------------------------------------------------------------------------

def setup_fdtd_tilted(fdtd, params, path, overhang=PWB_OVERHANG):
    """Configure FDTD region, source, monitors for the tilted-SOA simulation.

    The source and SOA-input-side monitor are rotated to match the tilted
    waveguide axes.  Monitors inside the straight PWB section are NOT rotated.
    """
    import numpy as np

    max_radius = max(params.r_pwb, params.r_in, params.r_out,
                     params.r_in_2, params.r_pwb_2, params.r_out_2)
    margin_yz = max_radius + 2e-6          # Y / Z margin
    margin_y_tilt = 3e-6                 # extra Y margin for tilted SOA sections

    soa_lead_in = 15e-6                   # places source at -7.5 μm (inside SOA WG)
    x_min = -soa_lead_in - overhang
    x_max = params.total_length + overhang + params.taper2_length * 0.3

    # ---- FDTD region ----
    fdtd.setresource("FDTD", "GPU", True)
    fdtd.addfdtd()
    fdtd.set("express mode", True)
    fdtd.set("dimension", "3D")
    fdtd.set("simulation time", params.simulation_time)
    fdtd.set("mesh accuracy", params.mesh_accuracy)
    fdtd.set("mesh type", "auto non-uniform")
    fdtd.set("x min", x_min)
    fdtd.set("x max", x_max)
    fdtd.set("y min", -margin_yz - margin_y_tilt)
    fdtd.set("y max",  margin_yz + margin_y_tilt)
    fdtd.set("z min", -margin_yz)
    fdtd.set("z max",  margin_yz)
    for side in ("x min", "x max", "y min", "y max", "z min", "z max"):
        fdtd.set(f"{side} bc", "PML")

    # ---- Mode source (inside tilted SOA output waveguide) ----
    # Lumerical mode sources support theta / phi rotation angles.
    # For injection axis = "x-axis", local y' → global Z, so theta
    # rotates k in the XY plane.  Negative theta → tilt toward -Y.
    # Source plane is YZ (2D X-normal); the 7° cross-section mismatch
    # gives ~0.75% width error (1/cos 7° ≈ 1.0075), acceptable.
    source_x = -soa_lead_in * 0.5
    # Y on the SOA output waveguide axis (rotated by SOA_OUT_ROT ≈ 173°)
    source_y = source_x * np.tan(np.deg2rad(SOA_OUT_ROT))

    fdtd.addmode()
    fdtd.set("name", "source")
    fdtd.set("injection axis", "x-axis")
    fdtd.set("direction", "forward")
    fdtd.set("theta", -SOA_TILT_DEG)          # -7° — tilt from +X toward -Y in XY plane
    fdtd.set("phi", 0.0)                      # no out-of-plane rotation
    fdtd.set("x", source_x)
    fdtd.set("y", source_y)
    fdtd.set("z", 0.0)
    fdtd.set("y span", 2.0 * (margin_yz + margin_y_tilt))
    fdtd.set("z span", 2.0 * margin_yz)
    fdtd.set("mode selection", "fundamental TE mode")
    fdtd.set("wavelength start", params.wavelength)
    fdtd.set("wavelength stop", params.wavelength)

    # ---- Power monitors (all 2D X-normal) ----
    # Monitors don't support rotation in Lumerical.  For the output
    # monitor inside the tilted SOA input waveguide we offset Y to sit
    # on the waveguide axis; the 7° cross-section tilt gives ~1% area
    # error (cos 7° ≈ 0.993) which is acceptable.
    mon_specs = [
        # (name,                      x,       in SOA?)
        ("input_monitor",            5e-6,    False),
        ("pwb_in_monitor",           params.taper1_length, False),
        ("pwb_out_monitor",          params.total_length - params.taper2_length, False),
        ("output_monitor",           params.total_length + 5e-6,  True),
    ]

    mon_positions = {}   # name → x
    for name, mx, in_soa in mon_specs:
        mon_positions[name] = mx
        if in_soa:
            d = mx - params.total_length
            my = d * np.tan(np.deg2rad(SOA_IN_ROT))
        else:
            my = 0.0

        fdtd.addpower()
        fdtd.set("name", name)
        fdtd.set("monitor type", "2D X-normal")
        fdtd.set("x", mx)
        fdtd.set("y", my)
        fdtd.set("z", 0.0)
        fdtd.set("y span", 2.0 * (margin_yz + margin_y_tilt))
        fdtd.set("z span", 2.0 * margin_yz)

    # ---- Mode expansions (reference the power monitors above) ----
    for name, ref_mon in [("mode_exp_input",  "input_monitor"),
                           ("mode_exp_output", "output_monitor")]:
        mx = mon_positions[ref_mon]
        in_soa = (ref_mon == "output_monitor")
        if in_soa:
            d = mx - params.total_length
            my = d * np.tan(np.deg2rad(SOA_IN_ROT))
        else:
            my = 0.0

        fdtd.addmodeexpansion()
        fdtd.set("name", name)
        fdtd.setexpansion("input", ref_mon)
        fdtd.set("mode selection", "fundamental TE mode")
        fdtd.set("x", mx)
        fdtd.set("y", my)
        fdtd.set("z", 0.0)
        fdtd.set("y span", 2.0 * (margin_yz + margin_y_tilt))
        fdtd.set("z span", 2.0 * margin_yz)

    # ---- Y-normal field profile monitor (side view: X-Z plane) ----
    # Use addprofile (not addpower) so E fields are always recorded for
    # visualisation.  addpower on a full-domain 2D monitor may skip E to
    # save memory, storing only the integrated T.
    fdtd.addprofile()
    fdtd.set("name", "transmission_monitor")
    fdtd.set("monitor type", "2D Y-normal")
    fdtd.set("y", 0.0)
    fdtd.set("x", (x_min + x_max) * 0.5)
    fdtd.set("x span", x_max - x_min)
    fdtd.set("z", 0.0)
    fdtd.set("z span", 2.0 * margin_yz)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    params = SOAPWBParams()
    overhang = PWB_OVERHANG

    fdtd = lumapi.FDTD()
    try:
        path = create_tilted_structure(fdtd, params, overhang)
        setup_fdtd_tilted(fdtd, params, path, overhang)

        fdtd.save(str(SAVE_PATH))
        print(f"Structure saved → {SAVE_PATH}")
        print(f"  PWB path:  [{(-overhang)*1e6:.0f}, {(params.total_length+overhang)*1e6:.0f}] μm")
        print(f"  Overhang:  {overhang*1e6:.0f} μm each side")
        print(f"  Mesh order: {PWB_MESH_ORDER}")

        print("\nRunning FDTD simulation ...")
        fdtd.run()
        print("Simulation complete.  Extracting results ...")

        T_forward = visualize_and_save_results(fdtd, params)
        print(f"T_forward = {T_forward}")
    finally:
        fdtd.close()


if __name__ == "__main__":
    main()
