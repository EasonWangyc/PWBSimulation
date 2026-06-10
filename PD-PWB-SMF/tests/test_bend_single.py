"""
test_bend_single.py — Run ONE bend-only simulation to verify the monitor fix.

If this produces valid T_total data, the full sweep is ready to go.

Usage:
    D:/Program Files/Lumerical/python/python.exe PD-PWB-SMF/tests/test_bend_single.py
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from sim_config import PD_DIR, PD_RESULTS_DIR, add_lumerical_api_path

add_lumerical_api_path()
import lumapi
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from pwb_core import (
    PWBParameters,
    generate_pwb_structure_bend,
    setup_fdtd_bend,
    get_data_bend,
    get_field_bend,
)

SAVE_PATH = PD_DIR / "temp_bend_single.fsp"


def compute_loss(T_total, eps=1e-12):
    return -10.0 * np.log10(np.clip(np.abs(-1.0 * T_total), eps, None))


params = PWBParameters()
params.L = 250e-6
params.h = 100e-6
params.r = 0.8e-6
params.r2 = 5e-6
params.l1 = 100e-6
params.bend_lift = 0.12
params.arch_position = 0.45
params.bend_shape = 0.5
params.drop_shape = 0.7
params.complex_segments = 100        # enough for smooth Bezier, fewer mesh triggers
params.mesh_accuracy = 1            # coarsest = smallest memory footprint
params.simulation_time = 3000e-15

print("Single bend-only test")
print(f"  bend_lift={params.bend_lift}  arch_position={params.arch_position}")
print(f"  mesh_accuracy={params.mesh_accuracy}  complex_segments={params.complex_segments}")
print(f"  sim_time={params.simulation_time*1e15:.0f}fs")
print("-" * 50)

fdtd = lumapi.FDTD()
try:
    generate_pwb_structure_bend(fdtd, params)
    print("Structure built.")

    setup_fdtd_bend(fdtd, params)
    print("FDTD setup complete (source + power monitor + field monitor).")

    fdtd.save(str(SAVE_PATH))
    print(f"Saved: {SAVE_PATH}")

    print("Running simulation...")
    fdtd.run()
    print("Done.")

    # --- Transmission ---
    T_total = get_data_bend(fdtd, params)
    loss_db = compute_loss(T_total)
    print(f"\n  T_total = {T_total:.6f}")
    print(f"  Loss    = {loss_db:.4f} dB")

    # --- Field profile plot ---
    Ex, Ey, Ez, x, z = get_field_bend(fdtd, params)
    intensity = np.abs(Ex)**2 + np.abs(Ey)**2 + np.abs(Ez)**2
    X, Z = np.meshgrid(x, z, indexing='ij')
    vmax = np.nanmax(intensity)

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.pcolormesh(X, Z, intensity, cmap='inferno',
                       shading='auto', norm=LogNorm(vmin=vmax*1e-6, vmax=vmax))
    fig.colorbar(im, ax=ax, label='|E|² (a.u.)')
    ax.set_xlabel('X (μm)')
    ax.set_ylabel('Z (μm)')
    ax.set_title(f'Bend-only _3  |  bend_lift={params.bend_lift}  '
                 f'arch_pos={params.arch_position}\n'
                 f'T_total={T_total:.4f}  loss={loss_db:.3f} dB  '
                 f'accuracy={params.mesh_accuracy}  segs={params.complex_segments}')
    ax.set_aspect('equal')
    plt.tight_layout()

    field_path = PD_RESULTS_DIR / "bend_shape_scan" / "single_test_field.png"
    field_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(field_path), dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Field plot saved: {field_path}")

except Exception as e:
    import traceback
    traceback.print_exc()
    print(f"ERROR: {e}")
finally:
    try:
        fdtd.close()
    except Exception:
        pass
