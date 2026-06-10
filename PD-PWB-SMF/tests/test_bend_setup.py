"""
test_bend_setup.py — Build bend-only structure and save .fsp for inspection.

Does NOT run FDTD. Open the saved .fsp in Lumerical to verify:
  - Straight input waveguide (10 μm) before the Bezier bend
  - Mode source in the middle of the input straight section
  - Bezier bend arch + drop
  - Straight output waveguide (5 μm) after the bend
  - Transmission monitor on the output straight section

Usage:
    D:/Program Files/Lumerical/python/python.exe PD-PWB-SMF/tests/test_bend_setup.py
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from sim_config import PD_DIR, add_lumerical_api_path

add_lumerical_api_path()
import lumapi

from pwb_core import (
    PWBParameters,
    generate_pwb_structure_bend,
    setup_fdtd_bend,
)

SAVE_PATH = PD_DIR / "temp_bend_only.fsp"

params = PWBParameters()
params.L = 250e-6
params.h = 100e-6
params.r = 0.8e-6
params.r2 = 7e-6
params.l1 = 100e-6
params.bend_lift = 0.12
params.arch_position = 0.45
params.bend_shape = 0.5
params.drop_shape = 0.7
params.complex_segments = 300
params.mesh_accuracy = 1
params.simulation_time = 2000e-15

print("Building bend-only structure...")
fdtd = lumapi.FDTD()

generate_pwb_structure_bend(fdtd, params)
print("  Structure built.")

setup_fdtd_bend(fdtd, params)
print("  FDTD region, source, and monitor configured.")

fdtd.save(str(SAVE_PATH))
print(f"\nSaved: {SAVE_PATH}")
print("Open this file in Lumerical GUI to inspect geometry before running sweeps.")
fdtd.close()
