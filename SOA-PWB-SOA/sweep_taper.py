r"""
Parameter sweep for SOA-PWB-SOA taper geometry optimisation.

Sweeps taper-1 / taper-2 axis radii and lengths independently (sequential
1-D sweeps — one parameter varied at a time, others held at baseline).

Outputs
-------
- ``results/sweep_taper_results.csv``  — T_forward for every run
- ``results/Pictures/sweep_{param}_{val:.2f}um.png`` — field-distribution image per run

Usage
-----
.. code-block:: powershell

   & "D:\Program Files\Lumerical\python\python.exe" SOA-PWB-SOA\sweep_taper.py
"""

import csv
import sys
from pathlib import Path

import numpy as np

# -- project root -----------------------------------------------------------
sys.path.append(str(Path(__file__).resolve().parents[1]))
from sim_config import SOA_DIR, SOA_RESULTS_DIR, add_lumerical_api_path

sys.path.insert(0, str(Path(__file__).resolve().parent))

add_lumerical_api_path()
import lumapi  # noqa: E402

from pwb_core import (  # noqa: E402
    SOAPWBParams,
    get_data,
    visualize_and_save_results,
)

# Import the tilted-case builders from the tilted run script.
# (run_tilted is a script, but its function defs are importable.)
import run_tilted as rt  # noqa: E402

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
RESULTS_CSV = SOA_RESULTS_DIR / "sweep_taper_results.csv"
PICTURES_DIR = SOA_RESULTS_DIR / "Pictures"
TEMP_DIR = SOA_DIR / "temp"               # checkpoint .fsp files before each run

for d in (PICTURES_DIR, TEMP_DIR):
    d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Sweep configuration
# ---------------------------------------------------------------------------

# Each entry is (attribute_name, display_label, [value1, value2, ...]).
# Values are in SI (metres).  Edit the lists to adjust sweep ranges.
SWEEP_PARAMS = [
    ("r_in",            "r_in",            [1.0e-6, 1.25e-6, 1.5e-6, 2.0e-6]),
    ("r_in_2",          "r_in_2",          [2.0e-6, 3.0e-6, 4.0e-6]),
    ("taper1_length",   "taper1_length",   [15e-6, 30e-6, 50e-6, 80e-6]),
    ("r_out",           "r_out",           [0.8e-6, 1.0e-6, 1.25e-6, 1.5e-6]),
    ("r_out_2",         "r_out_2",         [1.5e-6, 2.0e-6, 3.0e-6]),
    ("taper2_length",   "taper2_length",   [15e-6, 30e-6, 50e-6, 80e-6]),
]

# Baseline overrides (applied to every run regardless of sweep).
# Leave empty to match run_tilted.py defaults.  Add keys here to override
# specific SOAPWBParams attributes (e.g. "mesh_accuracy": 3).
BASELINE_OVERRIDES = {}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fmt_um(val):
    """Format a SI-length as a compact µm string."""
    return f"{val * 1e6:.2f}um"


def _param_label(attr, val):
    return f"{attr}={_fmt_um(val)}"


def _sweep_image_path(attr, val):
    return PICTURES_DIR / f"sweep_{attr}_{_fmt_um(val)}.png"


def _write_csv_header():
    """Write CSV header if the file doesn't exist yet."""
    if not RESULTS_CSV.exists():
        with open(RESULTS_CSV, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow([
                "param", "value_m", "value_um",
                "r_in_um", "r_in_2_um", "taper1_length_um",
                "r_pwb_um", "r_pwb_2_um",
                "r_out_um", "r_out_2_um", "taper2_length_um",
                "T_forward",
            ])


def _append_csv_row(param_name, val, params, T_forward):
    with open(RESULTS_CSV, "a", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            param_name,
            val,
            val * 1e6,
            params.r_in * 1e6,
            params.r_in_2 * 1e6,
            params.taper1_length * 1e6,
            params.r_pwb * 1e6,
            params.r_pwb_2 * 1e6,
            params.r_out * 1e6,
            params.r_out_2 * 1e6,
            params.taper2_length * 1e6,
            T_forward,
        ])


def _apply_overrides(params):
    """Apply baseline overrides to a params instance."""
    for attr, val in BASELINE_OVERRIDES.items():
        setattr(params, attr, val)


# Single temp file — overwritten each iteration (no need to keep one per run)
TEMP_FSP = SOA_DIR / "temp" / "sweep_current.fsp"


def _run_single(params, tag, img_filename):
    """Build, run, extract — identical flow to ``run_tilted.main()``.

    Each call opens a **fresh** ``lumapi.FDTD()``.
    """
    overhang = rt.PWB_OVERHANG
    img_path = PICTURES_DIR / img_filename

    fdtd = lumapi.FDTD()
    try:
        path = rt.create_tilted_structure(fdtd, params, overhang)
        rt.setup_fdtd_tilted(fdtd, params, path, overhang)

        fdtd.save(str(TEMP_FSP))
        print(f"  Checkpoint saved → {TEMP_FSP}")

        print("  Running FDTD ...")
        fdtd.run()
        print("  Done.  Extracting results ...")

        T_forward = visualize_and_save_results(
            fdtd, params, output_path=str(img_path), show_figure=False,
        )
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


def _sweep_param_attr(attr, values, base_params):
    """Sweep a single parameter through *values*, keeping everything else at baseline."""
    overhang = rt.PWB_OVERHANG

    for val in values:
        # Fresh params copy from baseline
        run_params = SOAPWBParams()
        for key in vars(base_params):
            setattr(run_params, key, getattr(base_params, key))
        _apply_overrides(run_params)
        setattr(run_params, attr, val)

        label = _param_label(attr, val)
        tag = f"sweep_{attr}_{_fmt_um(val)}"
        img_filename = f"sweep_{attr}_{_fmt_um(val)}.png"

        print(f"\n{'='*60}")
        print(f"  Sweep: {label}")
        print(f"  Taper1: r_in={run_params.r_in*1e6:.2f}  r_in_2={run_params.r_in_2*1e6:.2f}  L1={run_params.taper1_length*1e6:.0f} um")
        print(f"  Taper2: r_out={run_params.r_out*1e6:.2f}  r_out_2={run_params.r_out_2*1e6:.2f}  L2={run_params.taper2_length*1e6:.0f} um")
        print(f"  Mesh accuracy: {run_params.mesh_accuracy}  |  Sim time: {run_params.simulation_time*1e15:.0f} fs")
        print(f"{'='*60}")

        T_forward = _run_single(run_params, tag, img_filename)
        _append_csv_row(attr, val, run_params, T_forward)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    _write_csv_header()
    print(f"Results CSV → {RESULTS_CSV}")

    # --- baseline run ---
    print(f"\n{'#'*60}")
    print(f"  BASELINE RUN")
    print(f"{'#'*60}")
    baseline = SOAPWBParams()
    _apply_overrides(baseline)
    print(f"  Taper1: r_in={baseline.r_in*1e6:.2f}  r_in_2={baseline.r_in_2*1e6:.2f}  L1={baseline.taper1_length*1e6:.0f} um")
    print(f"  Taper2: r_out={baseline.r_out*1e6:.2f}  r_out_2={baseline.r_out_2*1e6:.2f}  L2={baseline.taper2_length*1e6:.0f} um")
    print(f"  Mesh accuracy: {baseline.mesh_accuracy}  |  Sim time: {baseline.simulation_time*1e15:.0f} fs")

    T_base = _run_single(baseline, "baseline", "sweep_baseline.png")
    _append_csv_row("baseline", 0.0, baseline, T_base)

    # --- parameter sweeps ---
    for attr, label, values in SWEEP_PARAMS:
        print(f"\n{'#'*60}")
        print(f"  SWEEP: {label}  →  {[_fmt_um(v) for v in values]}")
        print(f"{'#'*60}")
        _sweep_param_attr(attr, values, baseline)

    print(f"\n{'='*60}")
    print(f"  All sweeps complete!  Results → {RESULTS_CSV}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
