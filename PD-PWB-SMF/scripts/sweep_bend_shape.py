"""
sweep_bend_shape.py — 2D scan of bend_lift x arch_position for bend-only _3 PWB.

Bend-only: no SMF.fsp, no taper1. Source injects directly at bend entry (x=l1).
Fixed params: r=0.8um, L=250um, h=100um, l1=100um, bend_shape=0.5, drop_shape=0.7.

Usage:
    D:/Program Files/Lumerical/python/python.exe PD-PWB-SMF/scripts/sweep_bend_shape.py
"""

import os
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from sim_config import PD_RESULTS_DIR, add_lumerical_api_path

add_lumerical_api_path()
import lumapi

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pwb_core import (
    PWBParameters,
    generate_pwb_structure_bend,
    setup_fdtd_bend,
    get_data_bend,
)


def compute_loss(T_total, eps=1e-12):
    return -10.0 * np.log10(np.clip(np.abs(-1.0 * T_total), eps, None))


def run_sweep():
    RESULTS_DIR = PD_RESULTS_DIR / "bend_shape_scan"
    os.makedirs(RESULTS_DIR, exist_ok=True)

    base = PWBParameters()
    base.L = 250e-6
    base.h = 100e-6
    base.r = 0.8e-6
    base.r2 = 7e-6
    base.l1 = 100e-6
    base.bend_shape = 0.5
    base.drop_shape = 0.7
    base.complex_segments = 100    # enough for smooth Bezier; higher = more mesh triggers → more RAM
    base.mesh_accuracy = 1
    base.simulation_time = 2000e-15

    bend_lift_vals = [0.0, 0.06, 0.12, 0.18, 0.24, 0.30]
    arch_position_vals = [0.25, 0.40, 0.55, 0.70]

    total = len(bend_lift_vals) * len(arch_position_vals)
    all_results = []

    print(f"Bend-only sweep: {len(bend_lift_vals)} x {len(arch_position_vals)} = {total} simulations")
    print(f"bend_lift: {bend_lift_vals}")
    print(f"arch_position: {arch_position_vals}")
    print(f"Fixed: r={base.r*1e6:.1f}um L={base.L*1e6:.0f}um h={base.h*1e6:.0f}um")
    print("-" * 60)

    idx = 0
    for lift in bend_lift_vals:
        for arch_pos in arch_position_vals:
            idx += 1
            params = PWBParameters()
            for attr in ['L', 'h', 'r', 'r2', 'l1', 'bend_shape', 'drop_shape',
                         'complex_segments', 'mesh_accuracy', 'simulation_time']:
                setattr(params, attr, getattr(base, attr))
            params.bend_lift = lift
            params.arch_position = arch_pos

            label = f"lift={lift:.2f}_arch={arch_pos:.2f}"
            print(f"[{idx}/{total}] {label}", end=" ", flush=True)

            fdtd = lumapi.FDTD()
            try:
                generate_pwb_structure_bend(fdtd, params)
                setup_fdtd_bend(fdtd, params)
                fdtd.run()

                T_total = get_data_bend(fdtd, params)
                loss_db = compute_loss(T_total)
                all_results.append({
                    'bend_lift': lift, 'arch_position': arch_pos,
                    'T_total': T_total, 'loss_dB': loss_db,
                })
                print(f"T_total={T_total:.6f}  loss={loss_db:.4f} dB")

            except Exception as e:
                print(f"ERROR: {e}")
                all_results.append({
                    'bend_lift': lift, 'arch_position': arch_pos,
                    'T_total': float('nan'), 'loss_dB': float('nan'),
                })
            finally:
                try:
                    fdtd.close()
                except Exception:
                    pass

    # Save results
    results_file = RESULTS_DIR / "sweep_results.txt"
    with open(results_file, 'w', encoding='utf-8') as f:
        f.write("bend_lift\tarch_position\tT_total\tloss_dB\n")
        for r in all_results:
            f.write(f"{r['bend_lift']:.3f}\t{r['arch_position']:.3f}\t"
                    f"{r['T_total']:.6f}\t{r['loss_dB']:.4f}\n")
    print(f"\nResults saved: {results_file}")

    # Heatmap
    lift_u = np.array(bend_lift_vals)
    arch_u = np.array(arch_position_vals)
    loss_matrix = np.full((len(lift_u), len(arch_u)), np.nan)
    for res in all_results:
        if np.isnan(res['loss_dB']):
            continue
        i = np.argwhere(lift_u == res['bend_lift'])[0][0]
        j = np.argwhere(arch_u == res['arch_position'])[0][0]
        loss_matrix[i, j] = res['loss_dB']

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(loss_matrix, origin='lower', aspect='auto', cmap='viridis_r')
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('loss (dB)', fontsize=12)

    ax.set_xticks(np.arange(len(arch_u)))
    ax.set_xticklabels([f"{v:.2f}" for v in arch_u], rotation=45)
    ax.set_yticks(np.arange(len(lift_u)))
    ax.set_yticklabels([f"{v:.2f}" for v in lift_u])
    ax.set_xlabel('arch_position', fontsize=12)
    ax.set_ylabel('bend_lift', fontsize=12)
    ax.set_title(f'Bend-only loss: bend_lift x arch_position\n'
                 f'r={base.r*1e6:.1f}um L={base.L*1e6:.0f}um h={base.h*1e6:.0f}um '
                 f'bs={base.bend_shape} ds={base.drop_shape}', fontsize=11)
    plt.tight_layout()

    heatmap_path = RESULTS_DIR / "loss_heatmap.png"
    plt.savefig(str(heatmap_path), dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Heatmap saved: {heatmap_path}")

    # Top results
    valid = [r for r in all_results if not np.isnan(r['loss_dB'])]
    valid.sort(key=lambda x: x['loss_dB'])
    print("\nTop 5:")
    for r in valid[:5]:
        print(f"  lift={r['bend_lift']:.2f} arch={r['arch_position']:.2f}  "
              f"loss={r['loss_dB']:.4f} dB")

    return all_results


if __name__ == "__main__":
    run_sweep()
