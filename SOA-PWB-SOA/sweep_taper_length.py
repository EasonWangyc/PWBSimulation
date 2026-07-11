r"""
Sweep taper lengths for SOA-PWB-SOA.

Uses isolated single-taper testing (same approach as sweep_taper_width.py):
each taper is built with a 5 um PWB stub, and a mode expansion monitor
at the *end of the taper* records fundamental TE-mode coupling.

Taper-1:  SOA output WG -> Taper-1 (variable L) -> 5 um PWB stub
          FDTD scales with L, source forward, monitor at x = L, T_forward
Taper-2:  5 um PWB stub -> Taper-2 (variable L) -> SOA input WG
          FDTD scales with L, source backward at x=255, monitor at x=250-L,
          T_backward

Outputs
-------
- results/sweep_length_results.csv

Usage
-----
.. code-block:: powershell

   & "D:\Program Files\Lumerical\python\python.exe" SOA-PWB-SOA\sweep_taper_length.py
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

RESULTS_CSV  = SOA_RESULTS_DIR / "sweep_length_results.csv"
TEMP_FSP     = SOA_DIR / "temp" / "sweep_length_current.fsp"
TEMP_FSP.parent.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Fixed geometry
# ---------------------------------------------------------------------------
TOTAL_LENGTH  = 250e-6
TAPER1_LENGTH = 100e-6    # default — overwritten during sweep
TAPER2_LENGTH = 100e-6

PWB_STUB   = 5e-6         # short PWB after taper
SOA_LEAD   = 5e-6         # source placed this far inside SOA WG
EXTRA_MARGIN = 3e-6       # extra x-margin beyond structure for PML

# ---- optimal cross-section from width sweep ----
OPT_R_IN    = 0.5e-6
OPT_R_IN_2  = 2.8e-6
OPT_R_OUT   = 0.5e-6
OPT_R_OUT_2 = 1.5e-6

# FDTD window anchor points (fixed sides)
TAPER1_X_MIN = -10e-6                           # taper-1 test
TAPER2_X_MAX = 260e-6                           # taper-2 test

# Sweep ranges: 10–100 um, step 10 um
LENGTH_VALUES_UM = list(range(10, 110, 10))     # [10, 20, ..., 100]

# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------

def _write_csv_header():
    if not RESULTS_CSV.exists():
        with open(RESULTS_CSV, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "direction", "length_um",
                "r_in_um", "r_in_2_um", "r_out_um", "r_out_2_um",
                "r_pwb_um", "r_pwb_2_um",
                "T_forward",
            ])


def _append_csv_row(direction, L, params, T):
    with open(RESULTS_CSV, "a", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            direction, round(L * 1e6, 1),
            round(params.r_in   * 1e6, 4),
            round(params.r_in_2  * 1e6, 4),
            round(params.r_out  * 1e6, 4),
            round(params.r_out_2 * 1e6, 4),
            round(params.r_pwb  * 1e6, 4),
            round(params.r_pwb_2 * 1e6, 4),
            T,
        ])


# ---------------------------------------------------------------------------
# Build & setup (adapted from sweep_taper_width.py)
# ---------------------------------------------------------------------------

PWB_BUFFER = 3e-6  # extra margin so segments fully cover FDTD edges


def _build_structure(fdtd, params, direction, x_window):
    """Build PWB segments covering *x_window* (with buffer).

    Uses the full 250 um path so radius calculations are correct.
    A buffer is added so outermost segments reach beyond the FDTD
    PML boundaries.
    """
    path = generate_pwb_path(params)
    x_path = path[:, 0]

    fdtd.deleteall()
    fdtd.load(str(SOA_BASE_FSP))
    fdtd.importmaterialdb(str(MATERIAL_DB))

    x_min = x_window[0] - PWB_BUFFER
    x_max = x_window[1] + PWB_BUFFER

    mask = (x_path >= x_min) & (x_path <= x_max)
    inside = np.where(mask)[0]
    if len(inside) == 0:
        raise ValueError(f"No path points in [{x_min*1e6:.0f}, "
                         f"{x_max*1e6:.0f}]")

    start_idx = inside[0]
    end_idx   = inside[-1]
    n_seg     = max(1, end_idx - start_idx)
    prefix    = "PWB_taper1" if direction == "taper1" else "PWB_taper2"
    _add_path_segments(fdtd, path, params, start_idx, end_idx, n_seg, prefix)

    return path


def _setup_fdtd(fdtd, params, r_start, r_start_2, direction, L,
                x_window, source_x, monitor_x, source_dir):
    """Configure FDTD region, source, monitors."""
    max_radius = max(r_start, r_start_2 or 0, params.r_pwb, params.r_pwb_2)
    margin = max_radius + 1e-6
    x_min, x_max = x_window

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
    fdtd.set("mesh accuracy", 1)
    fdtd.set("mesh type", "auto non-uniform")
    fdtd.set("simulation time", 4000e-15)

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

    # ---- Power monitor at taper end ----
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

def _run_single(params, direction, L,
                x_window, source_x, monitor_x, source_dir,
                r_start, r_start_2):
    """Build, run, extract T_forward."""
    fdtd = lumapi.FDTD()
    try:
        _build_structure(fdtd, params, direction, x_window)
        _setup_fdtd(fdtd, params, r_start, r_start_2,
                    direction, L,
                    x_window, source_x, monitor_x, source_dir)

        fdtd.save(str(TEMP_FSP))
        print(f"  Checkpoint -> {TEMP_FSP}")
        print("  Running FDTD ...")
        fdtd.run()
        print("  Done.  Extracting ...")

        try:
            mode_data = fdtd.getresult("mode_exp", "expansion for input")
            t_key = "T_backward" if source_dir == "backward" else "T_forward"
            T = float(np.abs(mode_data[t_key]).flat[0])
        except Exception as exc:
            T = None
            print(f"  WARNING: Could not extract T: {exc}")

        print(f"  T = {T}")
        return T

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
# Sweep
# ---------------------------------------------------------------------------

def _sweep_length(direction, base_params):
    values_m = [v * 1e-6 for v in LENGTH_VALUES_UM]
    total = len(values_m)

    for i, L in enumerate(values_m, 1):
        run_params = SOAPWBParams()
        for key in vars(base_params):
            setattr(run_params, key, getattr(base_params, key))
        run_params.total_length = TOTAL_LENGTH
        run_params.taper1_length = TAPER1_LENGTH
        run_params.taper2_length = TAPER2_LENGTH

        if direction == "taper1":
            run_params.taper1_length = L
            run_params.r_in   = OPT_R_IN
            run_params.r_in_2 = OPT_R_IN_2
            r_start   = OPT_R_IN
            r_start_2 = OPT_R_IN_2
            x_window  = (TAPER1_X_MIN, L + PWB_STUB + EXTRA_MARGIN)
            source_x  = -SOA_LEAD
            monitor_x = L
            source_dir = "forward"
        else:
            run_params.taper2_length = L
            run_params.r_out   = OPT_R_OUT
            run_params.r_out_2 = OPT_R_OUT_2
            r_start   = OPT_R_OUT
            r_start_2 = OPT_R_OUT_2
            x_window  = (TOTAL_LENGTH - L - PWB_STUB - EXTRA_MARGIN, TAPER2_X_MAX)
            source_x  = TOTAL_LENGTH + SOA_LEAD     # 255 um
            monitor_x = TOTAL_LENGTH - L             # 250 - L
            source_dir = "backward"

        print(f"\n{'='*60}")
        print(f"  [{direction}] [{i}/{total}]  L = {L*1e6:.0f} um")
        print(f"  FDTD: [{x_window[0]*1e6:.0f}, {x_window[1]*1e6:.0f}] um")
        print(f"  source @ {source_x*1e6:.0f} ({source_dir})"
              f"  monitor @ {monitor_x*1e6:.0f}")
        print(f"{'='*60}")

        T = _run_single(run_params, direction, L,
                        x_window, source_x, monitor_x, source_dir,
                        r_start, r_start_2)
        _append_csv_row(direction, L, run_params, T)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    _write_csv_header()
    print(f"Results CSV -> {RESULTS_CSV}")
    print(f"Optimal params: r_in={OPT_R_IN*1e6:.1f}  r_in_2={OPT_R_IN_2*1e6:.1f}"
          f"  r_out={OPT_R_OUT*1e6:.1f}  r_out_2={OPT_R_OUT_2*1e6:.1f}")
    print(f"Lengths: {LENGTH_VALUES_UM} um")

    base = SOAPWBParams()
    base.r_pwb   = 0.8e-6
    base.r_pwb_2 = 0.8e-6

    for direction in ("taper1", "taper2"):
        _sweep_length(direction, base)

    print(f"\n{'='*60}")
    print(f"  All sweeps complete!  Results -> {RESULTS_CSV}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
