r"""
Sweep taper cross-section dimensions for SOA-PWB-SOA.

Tests each taper *independently* using a tiny FDTD window (~20 um)
that only covers the first 10 um of the taper.  A mode expansion
monitor at 5 um into the taper records fundamental TE-mode coupling
efficiency.

Taper-1:  x in [-10, 10] um,  source in SOA-out WG (forward),   monitor at 5 um
Taper-2:  x in [240, 260] um, source in SOA-in WG  (backward),  monitor at 245 um

Both use 2D grid sweeps (width x height).

Outputs
-------
- results/sweep_width_results.csv
- results/Pictures/sweep_width_{direction}_{r1:.2f}_{r2:.2f}um.png

Usage
-----
.. code-block:: powershell

   & "D:\Program Files\Lumerical\python\python.exe" SOA-PWB-SOA\sweep_taper_width.py
"""

import csv
import sys
from pathlib import Path

import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import (MATERIAL_DB, SOA_DIR,
                        SOA_RESULTS_DIR, add_lumerical_api_path)

sys.path.insert(0, str(Path(__file__).resolve().parent))

add_lumerical_api_path()
import lumapi  # noqa: E402

from pwb_core import (  # noqa: E402
    SOAPWBParams,
    _add_path_segments,
    generate_pwb_path,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SOA_BASE_FSP = SOA_DIR / "SOA_base_with_ar&cladding.fsp"

RESULTS_CSV  = SOA_RESULTS_DIR / "sweep_width_results.csv"
PICTURES_DIR = SOA_RESULTS_DIR / "Pictures"
TEMP_FSP     = SOA_DIR / "temp" / "sweep_width_current.fsp"

for d in (PICTURES_DIR, TEMP_FSP.parent):
    d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Geometry constants
# ---------------------------------------------------------------------------

# Total PWB length and taper lengths (fixed for width sweep)
TOTAL_LENGTH    = 250e-6
TAPER1_LENGTH   = 100e-6    # taper-1: 0 -> 100 um
TAPER2_LENGTH   = 100e-6    # taper-2: 150 -> 250 um
PWB_STRAIGHT    = TOTAL_LENGTH - TAPER1_LENGTH - TAPER2_LENGTH  # 50 um

# FDTD windows (tiny — only the taper entrance region)
TAPER1_X_WINDOW = (-10e-6, 10e-6)         # taper-1 test
TAPER2_X_WINDOW = (240e-6, 260e-6)        # taper-2 test

MONITOR_OFFSET  = 5e-6    # monitor placed 5 um into the taper from the facet
# Extra margin (um) added to the build window so PWB segments fully
# cover the FDTD region — the last segment extends half its z_span
# beyond its centre, so without a buffer the edge may be uncovered.
PWB_BUFFER = 2e-6
SOA_LEAD_IN     = 5e-6    # source placed 5 um inside SOA waveguide from facet


# ---------------------------------------------------------------------------
# Sweep configuration — 2D grids
# ---------------------------------------------------------------------------

# Each entry: (direction, param1, param2, label1, label2,
#              [values1 in um], [values2 in um])
# Every (v1, v2) combination is run.  ~7 pts/axis = ~49 runs/taper.
SWEEP_SPECS = [
    # ---- taper-1: r_in (width) x r_in_2 (height) ----
    ("taper1", "r_in", "r_in_2", "r_in", "r_in_2",
     [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4],
     [2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8]),

    # ---- taper-2: r_out (width) x r_out_2 (height) ----
    ("taper2", "r_out", "r_out_2", "r_out", "r_out_2",
     # centre 0.5 um, step 0.1 um
     [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1],
     # centre 1.5 um, step 0.2 um
     [0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7]),
]

BASELINE_OVERRIDES = {
    "mesh_accuracy": 1,
    "simulation_time": 2000e-15,   # tiny region — shorter sim time
    "r_pwb":   0.8e-6,            # PWB radius (axis-1) — fixed at 0.8 um
    "r_pwb_2": 0.8e-6,            # PWB radius (axis-2) — fixed at 0.8 um
}


# ---------------------------------------------------------------------------
# CSV helpers
# ---------------------------------------------------------------------------

def _write_csv_header():
    if not RESULTS_CSV.exists():
        with open(RESULTS_CSV, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "direction", "r1_attr", "r1_um", "r2_attr", "r2_um",
                "r_in_um", "r_in_2_um", "r_out_um", "r_out_2_um",
                "r_pwb_um", "r_pwb_2_um",
                "taper1_len_um", "taper2_len_um", "T_forward",
            ])


def _append_csv_row(direction, attr1, v1, attr2, v2, params, T_forward):
    with open(RESULTS_CSV, "a", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            direction, attr1, round(v1 * 1e6, 4), attr2, round(v2 * 1e6, 4),
            round(params.r_in   * 1e6, 4), round(params.r_in_2  * 1e6, 4),
            round(params.r_out  * 1e6, 4), round(params.r_out_2 * 1e6, 4),
            round(params.r_pwb  * 1e6, 4), round(params.r_pwb_2 * 1e6, 4),
            round(params.taper1_length * 1e6, 1),
            round(params.taper2_length * 1e6, 1),
            T_forward,
        ])


# ---------------------------------------------------------------------------
# Structure builder — only the taper entrance window
# ---------------------------------------------------------------------------

def _build_window_structure(fdtd, params, direction, x_window):
    """Build PWB segments covering *x_window* (with buffer).

    Uses the ACTUAL 250 um path and taper lengths so radius-profile
    calculations are correct.  A buffer is added so segments extend
    beyond the FDTD boundaries — otherwise the outermost segment
    may leave a gap at the PML edge.
    """
    path = generate_pwb_path(params)   # full 250 um path, N=200 points

    fdtd.deleteall()
    fdtd.load(str(SOA_BASE_FSP))
    fdtd.importmaterialdb(str(MATERIAL_DB))

    x_path = path[:, 0]
    x_min = x_window[0] - PWB_BUFFER
    x_max = x_window[1] + PWB_BUFFER

    mask = (x_path >= x_min) & (x_path <= x_max)
    inside_indices = np.where(mask)[0]
    if len(inside_indices) == 0:
        raise ValueError(f"No path points inside build window "
                         f"[{x_min*1e6:.0f}, {x_max*1e6:.0f}]")

    start_idx = inside_indices[0]
    end_idx   = inside_indices[-1]
    n_seg     = max(1, end_idx - start_idx)

    if direction == "taper1":
        _add_path_segments(fdtd, path, params,
                           start_idx, end_idx, n_seg, "PWB_taper1")
    else:
        _add_path_segments(fdtd, path, params,
                           start_idx, end_idx, n_seg, "PWB_taper2")

    return path


# ---------------------------------------------------------------------------
# FDTD setup for the tiny-window simulation
# ---------------------------------------------------------------------------

def _setup_window_fdtd(fdtd, params, r_start, r_start_2,
                       x_window, source_x, monitor_x, source_dir):
    """Configure FDTD region, source, and monitor for the truncated window."""
    max_radius = max(r_start, r_start_2 or 0, params.r_pwb, params.r_pwb_2)
    margin = max_radius + 1e-6

    x_min, x_max = x_window

    # ---- FDTD region ----
    fdtd.setresource("FDTD", "GPU", True)
    fdtd.addfdtd()
    fdtd.set("express mode", True)
    fdtd.set("dimension", "3D")
    fdtd.set("x min", x_min)
    fdtd.set("x max", x_max)
    fdtd.set("y min", -margin)
    fdtd.set("y max",  margin)
    fdtd.set("z min", -margin)
    fdtd.set("z max",  margin)
    for side in ("x min", "x max", "y min", "y max", "z min", "z max"):
        fdtd.set(f"{side} bc", "PML")
    fdtd.set("mesh accuracy", params.mesh_accuracy)
    fdtd.set("mesh type", "auto non-uniform")
    fdtd.set("simulation time", params.simulation_time)

    # ---- Mode source ----
    fdtd.addmode()
    fdtd.set("name", "source")
    fdtd.set("injection axis", "x-axis")
    fdtd.set("direction", source_dir)
    fdtd.set("x", source_x)
    fdtd.set("y", 0)
    fdtd.set("z", 0)
    fdtd.set("y span", 2 * margin)
    fdtd.set("z span", 2 * margin)
    fdtd.set("mode selection", "fundamental TE mode")
    fdtd.set("wavelength start", params.wavelength)
    fdtd.set("wavelength stop", params.wavelength)

    # ---- Power monitor at 5 um into the taper ----
    fdtd.addpower()
    fdtd.set("name", "taper_monitor")
    fdtd.set("monitor type", "2D X-normal")
    fdtd.set("x", monitor_x)
    fdtd.set("y", 0)
    fdtd.set("z", 0)
    fdtd.set("y span", 2 * margin)
    fdtd.set("z span", 2 * margin)

    # ---- Mode expansion at same position ----
    fdtd.addmodeexpansion()
    fdtd.set("name", "mode_exp")
    fdtd.setexpansion("input", "taper_monitor")
    fdtd.set("mode selection", "fundamental TE mode")
    fdtd.set("x", monitor_x)
    fdtd.set("y", 0)
    fdtd.set("z", 0)
    fdtd.set("y span", 2 * margin)
    fdtd.set("z span", 2 * margin)

    # ---- Y-normal profile monitor ----
    fdtd.addprofile()
    fdtd.set("name", "transmission_monitor")
    fdtd.set("monitor type", "2D Y-normal")
    fdtd.set("y", 0)
    fdtd.set("x min", x_min)
    fdtd.set("x max", x_max)
    fdtd.set("z", 0)
    fdtd.set("z span", 2 * margin)


# ---------------------------------------------------------------------------
# Single run
# ---------------------------------------------------------------------------

def _run_single_window(params, direction, r_start, r_start_2,
                       x_window, source_x, monitor_x, source_dir):
    """Build, run, extract T_forward for one (v1, v2) point."""
    fdtd = lumapi.FDTD()
    try:
        _build_window_structure(fdtd, params, direction, x_window)
        _setup_window_fdtd(fdtd, params,
                           r_start, r_start_2,
                           x_window, source_x, monitor_x, source_dir)

        fdtd.save(str(TEMP_FSP))
        print(f"  Checkpoint -> {TEMP_FSP}")

        print("  Running FDTD ...")
        fdtd.run()
        print("  Done.  Extracting ...")

        try:
            mode_data = fdtd.getresult("mode_exp", "expansion for input")
            # forward source → T_forward; backward → T_backward
            t_key = "T_backward" if source_dir == "backward" else "T_forward"
            T_forward = float(np.abs(mode_data[t_key]).flat[0])
        except Exception as exc:
            T_forward = None
            print(f"  WARNING: Could not extract T_forward: {exc}")

        print(f"  T_forward = {T_forward}")
        return T_forward

    except Exception as exc:
        print(f"  *** FAILED: {exc}")
        import traceback
        traceback.print_exc()
        try:
            fdtd.save(str(TEMP_FSP))
        except Exception:
            pass
        return None

    finally:
        try:
            fdtd.close()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# 2D grid sweep
# ---------------------------------------------------------------------------

def _sweep_2d_grid(direction, attr1, attr2, label1, label2,
                   values1_um, values2_um, base_params):
    """Run all (v1, v2) combinations."""
    values1_m = [v * 1e-6 for v in values1_um]
    values2_m = [v * 1e-6 for v in values2_um]
    total = len(values1_m) * len(values2_m)

    if direction == "taper1":
        x_window  = TAPER1_X_WINDOW
        facet_x   = 0.0                # SOA output facet
        source_x  = facet_x - SOA_LEAD_IN          # -5 um
        monitor_x = facet_x + MONITOR_OFFSET       #  5 um
        source_dir = "forward"
    else:
        x_window  = TAPER2_X_WINDOW
        facet_x   = TOTAL_LENGTH       # 250 um — SOA input facet
        source_x  = facet_x + SOA_LEAD_IN          # 255 um
        monitor_x = facet_x - MONITOR_OFFSET       # 245 um
        source_dir = "backward"

    count = 0
    for v1_m in values1_m:
        for v2_m in values2_m:
            count += 1
            run_params = SOAPWBParams()
            for key in vars(base_params):
                setattr(run_params, key, getattr(base_params, key))
            for attr, val in BASELINE_OVERRIDES.items():
                setattr(run_params, attr, val)

            # Fix taper lengths
            run_params.taper1_length = TAPER1_LENGTH
            run_params.taper2_length = TAPER2_LENGTH
            run_params.total_length  = TOTAL_LENGTH

            setattr(run_params, attr1, v1_m)
            setattr(run_params, attr2, v2_m)

            if direction == "taper1":
                r_start   = run_params.r_in
                r_start_2 = run_params.r_in_2
            else:
                r_start   = run_params.r_out
                r_start_2 = run_params.r_out_2

            r1_str = f"{v1_m*1e6:.2f}"
            r2_str = f"{v2_m*1e6:.2f}"

            print(f"\n{'='*60}")
            print(f"  [{direction}] [{count}/{total}]  "
                  f"{label1}={r1_str}um  {label2}={r2_str}um")
            print(f"  FDTD window: {x_window[0]*1e6:.0f} to {x_window[1]*1e6:.0f} um")
            print(f"  source @ {source_x*1e6:.0f} um ({source_dir})"
                  f"  monitor @ {monitor_x*1e6:.0f} um")
            print(f"  r_start={r_start*1e6:.2f}  r_start_2={r_start_2*1e6:.2f}")
            print(f"{'='*60}")

            T_forward = _run_single_window(
                run_params, direction, r_start, r_start_2,
                x_window, source_x, monitor_x, source_dir)
            _append_csv_row(direction, attr1, v1_m, attr2, v2_m,
                            run_params, T_forward)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    _write_csv_header()
    print(f"Results CSV -> {RESULTS_CSV}")
    print(f"Base .fsp  -> {SOA_BASE_FSP}")

    base = SOAPWBParams()
    base.taper1_length = TAPER1_LENGTH
    base.taper2_length = TAPER2_LENGTH
    base.total_length  = TOTAL_LENGTH
    for attr, val in BASELINE_OVERRIDES.items():
        setattr(base, attr, val)

    for (direction, attr1, attr2, label1, label2,
         values1_um, values2_um) in SWEEP_SPECS:
        n1, n2 = len(values1_um), len(values2_um)
        print(f"\n{'#'*60}")
        print(f"  2D GRID: [{direction}]  {label1} x {label2}"
              f"  ({n1} x {n2} = {n1*n2} runs)")
        print(f"  {label1}: {values1_um} um")
        print(f"  {label2}: {values2_um} um")
        print(f"{'#'*60}")
        _sweep_2d_grid(direction, attr1, attr2, label1, label2,
                       values1_um, values2_um, base)

    print(f"\n{'='*60}")
    print(f"  All sweeps complete!  Results -> {RESULTS_CSV}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
