"""
baseline_2.py — Run _2 (simple 90-degree arc) bend-only simulations
with SMOOTH radius taper (r → r2 over second half of path).

COMPARISON PURPOSE:
  - OLD _2 (results/baseline_2/):  constant r in bend, sharp linear taper in taper2
  - NEW _2 (results/baseline_2_smooth/): smooth raised-cosine taper like _3 bend-only

Same (L, h, r, r2, l1) as the old baseline for direct A/B comparison.

Usage:
    D:/Program Files/Lumerical/python/python.exe PD-PWB-SMF/scripts/baseline_2.py
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
from matplotlib.colors import LogNorm

from pwb_core import (
    PWBParameters,
    generate_pwb_bend_path_2,
    generate_pwb_structure_bend,
    setup_fdtd_bend,
    get_data_bend,
    get_field_bend,
)

RESULTS_DIR = PD_RESULTS_DIR / "baseline_2_smooth"
SAVE_PATH = RESULTS_DIR / "temp_baseline.fsp"
os.makedirs(RESULTS_DIR, exist_ok=True)


def compute_loss(T_total, eps=1e-12):
    return -10.0 * np.log10(np.clip(np.abs(-1.0 * T_total), eps, None))


def run_baseline():
    base = PWBParameters()
    base.L = 250e-6
    base.h = 100e-6
    base.r = 0.8e-6
    base.r2 = 7e-6             # match old baseline for direct comparison
    base.l1 = 100e-6
    base.complex_segments = 100
    base.mesh_accuracy = 1
    base.simulation_time = 3000e-15

    # R must be < min(L - l1, h) = min(150, 100) = 100 um
    R_values = [20e-6, 40e-6, 60e-6, 80e-6, 95e-6]
    all_results = []

    print("_2 Baseline (smooth taper): 90-degree arc + raised-cosine radius profile")
    print(f"L={base.L*1e6:.0f}um  h={base.h*1e6:.0f}um  "
          f"r={base.r*1e6:.1f}um  r2={base.r2*1e6:.1f}um  l1={base.l1*1e6:.0f}um")
    print(f"R values: {[f'{v*1e6:.0f}' for v in R_values]} um")
    print(f"mesh_accuracy={base.mesh_accuracy}  segments={base.complex_segments}")
    print("-" * 60)

    for i, R in enumerate(R_values):
        params = PWBParameters()
        for attr in ['L', 'h', 'r', 'r2', 'l1', 'complex_segments',
                     'mesh_accuracy', 'simulation_time']:
            setattr(params, attr, getattr(base, attr))
        params.R = R

        label = f"R={R*1e6:.0f}um"
        print(f"\n[{i+1}/{len(R_values)}] {label}")

        fdtd = lumapi.FDTD()
        try:
            generate_pwb_structure_bend(fdtd, params, path_func=generate_pwb_bend_path_2)
            setup_fdtd_bend(fdtd, params, path_func=generate_pwb_bend_path_2)
            fdtd.save(str(SAVE_PATH))
            fdtd.run()

            T_total = get_data_bend(fdtd, params)
            loss_db = compute_loss(T_total)
            all_results.append({
                'R_um': R * 1e6, 'T_total': T_total, 'loss_dB': loss_db,
            })
            print(f"  T_total = {T_total:.6f}  loss = {loss_db:.4f} dB")

            # field plot from field_monitor (X-Z plane)
            Ex, Ey, Ez, x_f, z_f = get_field_bend(fdtd, params)
            intensity = np.abs(Ex)**2 + np.abs(Ey)**2 + np.abs(Ez)**2
            vmax = np.nanmax(intensity)
            X, Z = np.meshgrid(x_f, z_f, indexing='ij')

            fig, ax = plt.subplots(figsize=(8, 6))
            im = ax.pcolormesh(X, Z, intensity, cmap='inferno',
                              shading='auto', norm=LogNorm(vmin=vmax*1e-6, vmax=vmax))
            fig.colorbar(im, ax=ax, label='|E|² (a.u.)')
            ax.set_xlabel('X (μm)')
            ax.set_ylabel('Z (μm)')
            ax.set_title(f'_2 smooth taper  R={R*1e6:.0f}μm\n'
                        f'T_total={T_total:.4f}  loss={loss_db:.3f} dB')
            ax.set_aspect('equal')
            plt.tight_layout()

            fig_path = RESULTS_DIR / f"field_R_{R*1e6:.0f}um.png"
            plt.savefig(str(fig_path), dpi=200, bbox_inches='tight')
            plt.close(fig)

        except Exception as e:
            import traceback
            traceback.print_exc()
            print(f"  ERROR: {e}")
            all_results.append({
                'R_um': R * 1e6, 'T_total': float('nan'), 'loss_dB': float('nan'),
            })
        finally:
            try:
                fdtd.close()
            except Exception:
                pass

    # --- Summary ---
    print("\n" + "=" * 60)
    print("_2 Baseline (smooth taper) — Summary")
    print("-" * 60)
    for r in all_results:
        print(f"  R = {r['R_um']:5.0f} μm  |  T_total = {r['T_total']:+.6f}  |  "
              f"loss = {r['loss_dB']:.4f} dB")

    # Save text
    results_file = RESULTS_DIR / "baseline_results.txt"
    with open(results_file, 'w', encoding='utf-8') as f:
        f.write("R(um)\tT_total\tloss_dB\n")
        for r in all_results:
            f.write(f"{r['R_um']:.0f}\t{r['T_total']:.6f}\t{r['loss_dB']:.4f}\n")
    print(f"\nResults saved: {results_file}")

    # Loss vs R plot
    valid = [r for r in all_results if not np.isnan(r['T_total'])]
    if valid:
        fig, ax = plt.subplots(figsize=(7, 4))
        R_plot = [r['R_um'] for r in valid]
        loss_plot = [r['loss_dB'] for r in valid]
        ax.plot(R_plot, loss_plot, 'o-', color='coral', markersize=8, label='_2 smooth taper')
        ax.set_xlabel('R (μm)')
        ax.set_ylabel('Loss (dB)')
        ax.set_title(f'_2 Baseline (smooth): r={base.r*1e6:.1f}μm  r2={base.r2*1e6:.1f}μm')
        ax.grid(True, alpha=0.3)
        ax.legend()
        plt.tight_layout()
        fig_path = RESULTS_DIR / "baseline_loss_vs_R.png"
        plt.savefig(str(fig_path), dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f"Plot saved: {fig_path}")

    return all_results


if __name__ == "__main__":
    run_baseline()
